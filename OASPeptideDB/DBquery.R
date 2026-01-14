# app.R
# DAT-DB style Shiny app: query partitioned Parquet folder (final_enriched) + download files
# Update (v1.3 -> v1.4):
# - FIX: DuckDB decodes %2F to "/" in partition values (Disease column). Using folder codes caused 0 rows.
# - Disease filter now uses the *actual Disease values* found in partition folders, URL-decoded (e.g., "Allergy/NoSIT").
# - Summary UI + downloaded summary TXT show FULL disease names ONLY for Disease line (mapping supports both encoded and decoded).

suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinyjs)
  library(DBI)
  library(duckdb)
  library(DT)
  library(shinyFiles)
  library(data.table)
  library(Biostrings)
})

# -------------------------
# Helpers
# -------------------------
norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
empty_preview_df <- function() data.frame()

# ---- Disease full-name mapping (partition code -> full name) ----
# (Your original mapping uses encoded forms like %2F)
disease_code_to_full <- c(
  "AChR-MG"                           = "Acetylcholine receptor-positive myasthenia gravis",
  "Allergic-Rhinitis-In-Season"       = "Allergic rhinitis in season",
  "Allergic-Rhinitis-Out-Of-Season"   = "Allergic rhinitis out of  season",
  "Allergy%2FNoSIT"                   = "Allergy NoSIT",
  "Allergy%2FSIT"                     = "Allergy SIT",
  "Asthma"                            = "Asthma",
  "CLL"                               = "Chronic lymphocytic leukemia",
  "CMV"                               = "Cytomegalovirus",
  "CMV%2FEBV"                         = "Cytomegalovirus/Epstein–Barr virus",
  "Dengue"                            = "Dengue",
  "Ebola"                             = "Ebola",
  "EBV"                               = "Epstein–Barr virus",
  "HCV"                               = "Hepatitis C virus",
  "Healthy%2Fceliac-disease"          = "Celiac",
  "HIV"                               = "HIV",
  "Light-Chain-Amyloidosis"           = "Light chain amyloidosis",
  "MS"                                = "Multiple sclerosis",
  "MuSK-MG"                           = "Muscle-specific kinase myasthenia gravis",
  "Non-Dengue-Febrile-Illness"        = "Non-dengue febrile illness",
  "Obstructive-Sleep-Apnea"           = "Obstructive sleep apnea",
  "SARS-COV-2"                        = "SARS-COV-2",
  "SLE"                               = "Systemic lupus erythematosus",
  "Tonsillitis"                       = "Tonsillitis",
  "Tonsillitis%2FObstructive-Sleep-Apnea" = "Tonsillitis/Obstructive sleep apnea"
)

norm_disease_code <- function(x) trimws(as.character(x))

# Build an "extended" mapping that also matches decoded keys (e.g., "Allergy/NoSIT")
disease_code_to_full_extended <- local({
  m <- disease_code_to_full
  decoded_names <- tryCatch(utils::URLdecode(names(m)), error = function(e) names(m))
  # Add decoded keys (do not overwrite existing keys)
  add <- m
  names(add) <- decoded_names
  for (k in names(add)) {
    if (!k %in% names(m)) m[k] <- add[[k]]
  }
  m
})

# Only for summary: return full name if known, otherwise return input as-is (best-effort)
disease_fullname_for_summary <- function(code_or_decoded) {
  x <- norm_disease_code(code_or_decoded)
  if (!nzchar(x) || x %in% c("None", "All")) return(x)
  if (x %in% names(disease_code_to_full_extended)) return(unname(disease_code_to_full_extended[x]))
  x
}

# List partitions from folder tree; return Disease already URL-decoded to match DuckDB
list_partitions_from_folders <- function(root) {
  root <- norm(root)
  if (!dir.exists(root)) return(NULL)
  
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  part_dirs <- dirs[
    grepl("Disease=", dirs) &
      grepl("BSource=", dirs) &
      grepl("BType=", dirs) &
      grepl("Isotype=", dirs)
  ]
  if (length(part_dirs) == 0) return(NULL)
  
  get_val <- function(x, key) sub(paste0(".*", key, "=?([^/]+).*"), "\\1", x)
  
  disease_raw <- norm_disease_code(get_val(part_dirs, "Disease"))
  # Decode to match what DuckDB actually returns in the Disease column
  disease_dec <- tryCatch(utils::URLdecode(disease_raw), error = function(e) disease_raw)
  
  list(
    Disease = sort(unique(disease_dec)),
    BSource = sort(unique(get_val(part_dirs, "BSource"))),
    BType   = sort(unique(get_val(part_dirs, "BType"))),
    Isotype = sort(unique(get_val(part_dirs, "Isotype"))),
    n_part_dirs = length(part_dirs)
  )
}

sql_quote <- function(x) paste0("'", gsub("'", "''", x), "'")

discover_uniqueness_col <- function(con) {
  cols <- dbGetQuery(con, "DESCRIBE v;")$column_name
  if ("n_distinct_peptide" %in% cols) return("n_distinct_peptide")
  if ("n_distinct_antibodies" %in% cols) return("n_distinct_antibodies")
  NULL
}

build_where <- function(f, col_uniqueness = NULL) {
  w <- character()
  
  if (f$Disease != "All") w <- c(w, paste0("Disease = ", sql_quote(f$Disease)))
  if (f$BSource != "All") w <- c(w, paste0("BSource = ", sql_quote(f$BSource)))
  if (f$BType   != "All") w <- c(w, paste0("BType = ", sql_quote(f$BType)))
  if (f$Isotype != "All") w <- c(w, paste0("Isotype = ", sql_quote(f$Isotype)))
  
  if (f$Disease_specific == "Yes") w <- c(w, "COALESCE(Disease_presence_peptide,0) = 1")
  if (f$Disease_specific == "No")  w <- c(w, "COALESCE(Disease_presence_peptide,0) <> 1")
  
  if (!is.null(col_uniqueness)) {
    if (f$Peptide_uniqueness == "Yes") w <- c(w, paste0("COALESCE(", col_uniqueness, ",0) = 1"))
    if (f$Peptide_uniqueness == "No")  w <- c(w, paste0("COALESCE(", col_uniqueness, ",0) <> 1"))
  }
  
  if (f$CDR3 == "Yes") w <- c(w, "COALESCE(peptide_in_cdr3,0) = 1")
  if (f$CDR3 == "No")  w <- c(w, "COALESCE(peptide_in_cdr3,0) <> 1")
  
  if (length(w) == 0) return("")
  paste("WHERE", paste(w, collapse = " AND "))
}

reorder_antibody_last <- function(df) {
  if (!("Antibody" %in% names(df))) return(df)
  df[, c(setdiff(names(df), "Antibody"), "Antibody"), drop = FALSE]
}

filters_ready <- function(f) !identical(f$Disease, "None")

fmt_int <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
  format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

export_select_list <- function(con) {
  cols <- tryCatch(dbGetQuery(con, "DESCRIBE v;")$column_name, error = function(e) character())
  needed <- c("Peptide", "v_call", "d_call", "j_call", "cdr3_aa", "filename", "Antibody")
  out <- vapply(needed, function(nm) {
    if (nm %in% cols) nm else sprintf("NULL AS %s", nm)
  }, character(1))
  paste(out, collapse = ", ")
}

write_summary_txt <- function(path_txt, summary_lines) {
  dir.create(dirname(path_txt), recursive = TRUE, showWarnings = FALSE)
  writeLines(summary_lines, con = path_txt, useBytes = TRUE)
}

safe_header_token <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "NA"
  gsub("[[:space:]]+", "_", x)
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  shinyjs::useShinyjs(),
  theme = shinytheme("cerulean"),
  
  tags$head(tags$style(HTML("
    body { font-family: Arial; }

    .custom-header {
      font-size: 20px;
      font-weight: bold;
      background: #f4f4f9;
      border-left: 6px solid #ff9800;
      padding: 10px;
      margin: 0;
      border-radius: 6px;
      box-shadow: 0 2px 4px rgba(0,0,0,0.10);
      color: #333;
    }

    .panel-container {
      border: 1px solid #ddd;
      padding: 16px;
      border-radius: 8px;
      background: #f9f9f9;
      margin-top: 20px;
    }

    .vstack {
      display: flex;
      flex-direction: column;
      gap: var(--gap, 12px);
    }
    .gap-12 { --gap: 12px; }

    .vstack > * { margin: 0 !important; }
    .vstack .btn { margin: 0 !important; }

    .panel-container.nested { margin-top: 0 !important; }

    .download-row { display: flex; gap: 16px; align-items: stretch; flex-wrap: wrap; }
    .download-col-summary { flex: 0 0 34%; min-width: 260px; }
    .download-col-table { flex: 1 1 0; min-width: 0; }

    .download-row .panel-container { margin-top: 0 !important; }

    .download-col-table .dataTables_wrapper {
      width: 100% !important;
      max-width: 100% !important;
      overflow-x: auto;
    }
    .download-col-table table.dataTable { width: 100% !important; }

    .plain-title { font-weight: 700; font-size: 18px; margin-bottom: 6px !important; }

    .download-btn-row { display:flex; justify-content:flex-start; gap: 14px; flex-wrap: wrap; }

    .summary-line { margin: 2px 0 !important; }
    .summary-indent { margin-left: 18px !important; }
  "))),
  tags$title("DAT-DB"),
  
  titlePanel(HTML("<strong>D</strong>isease-specific <strong>A</strong>ntibody <strong>T</strong>ryptic peptides <strong>D</strong>atabases for <strong>B</strong>ottom-up proteomics (<strong>DAT-DB</strong>)")),
  
  tabsetPanel(
    tabPanel(
      "Introduction",
      fluidPage(
        div(class = "panel-container",
            #div(class = "custom-header", "Overview"),
            tags$div(
              tags$p(tags$strong("Welcome to DAT-DB app!")),
              tags$p("The identification of human antibodies using bottom-up proteomics relies on database searches that match experimental peptide fragments to theoretical values derived from protein sequences in databases.
Standard protein databases, such as UniProt and NCBI-RefSeq, contain a limited number of antibody sequences compared to the human body’s capacity to generate billions and currently lack disease-specific antibody sequences.
This limitation could lead to the misidentification of antibodies in samples.
To overcome this challenge, DAT-DB, a database of disease-specific antibodies, provides researchers with antibody tryptic peptides derived from next-generation sequencing of antibody repertoires."),
              tags$p(tags$strong("Author: Xuan-Tung Trinh (xttrinh1199@gmail.com)")),
              tags$hr(),
              tags$p(tags$strong("How to use this app")),
              tags$ul(
                tags$li(tags$b("Step 1:"), " Click ", tags$b("Browse…"), " to select your peptide database folder (the partitioned Parquet root)."),
                tags$li(tags$b("Step 2:"), " Choose filters. The ", tags$b("Disease"), " filter displays the actual Disease values used by DuckDB (e.g., ", tags$code("Allergy/NoSIT"), ")."),
                tags$li(tags$b("Step 3:"), " Preview results and download as ", tags$b("CSV"), ", ", tags$b("FASTA"), " or a ", tags$b("Summary TXT"), "."),
                tags$li(tags$b("Note:"), " The on-screen Summary and downloaded Summary TXT show ", tags$b("full disease names"), " for the Disease line only.")
              ),
              tags$p("If you select Disease = 'None', previews and downloads are disabled to avoid heavy queries.")
            )
        )
      )
    ),
    
    tabPanel(
      "Download Peptides",
      fluidRow(
        column(
          3,
          div(class = "panel-container",
              div(class = "vstack gap-12",
                  div(class = "custom-header", "1. Select Data"),
                  shinyDirButton("db_folder_btn", "Browse…", "Select folder"),
                  div(class = "panel-container nested",
                      div(class = "plain-title", "Database info"),
                      uiOutput("db_info")
                  )
              ),
              
              tags$hr(),
              
              # Disease filter now matches DuckDB values (decoded)
              selectInput("Disease", "Disease", c("None")),
              selectInput("BSource", "B-cells source", "All"),
              selectInput("BType", "B-cells type", "All"),
              selectInput("Isotype", "Antibody isotype", "All"),
              
              tags$hr(),
              
              selectInput("Disease_specific", "Disease-specific", c("Yes","No"), "Yes"),
              selectInput("Peptide_uniqueness", "Peptide uniqueness", c("Yes","No"), "Yes"),
              selectInput("CDR3", "CDR3", c("All","Yes","No"), "All"),
              
              tags$hr(),
              
              selectInput("max_preview", "Preview rows", c(20,100,500,1000), 20)
          )
        ),
        column(
          9,
          div(class = "panel-container",
              div(class = "vstack gap-12",
                  div(class = "custom-header", "2. Download Data"),
                  div(class = "download-btn-row", uiOutput("download_btns")),
                  div(class = "download-row",
                      div(class = "panel-container download-col-summary",
                          div(class = "plain-title", "Summary"),
                          uiOutput("summary")
                      ),
                      div(class = "panel-container download-col-table",
                          div(class = "plain-title", "Preview table"),
                          DTOutput("preview")
                      )
                  )
              )
          )
        )
      )
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  # ---- Controlled DuckDB temp dir (and cleanup on close) ----
  duck_tmp_dir <- file.path(tempdir(), "DATDB_duckdb_tmp")
  if (dir.exists(duck_tmp_dir)) unlink(duck_tmp_dir, recursive = TRUE, force = TRUE)
  dir.create(duck_tmp_dir, showWarnings = FALSE, recursive = TRUE)
  
  con <- dbConnect(duckdb(), ":memory:")
  dbExecute(con, sprintf("PRAGMA temp_directory='%s';", gsub("\\\\", "/", duck_tmp_dir)))
  dbExecute(con, "PRAGMA threads=6;")
  dbExecute(con, "PRAGMA enable_progress_bar=false;")
  
  session$onSessionEnded(function() {
    try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    if (dir.exists(duck_tmp_dir)) unlink(duck_tmp_dir, recursive = TRUE, force = TRUE)
  })
  
  volumes <- c("Computer" = shinyFiles::getVolumes()())
  shinyDirChoose(input, "db_folder_btn", roots = volumes)
  
  db_loaded <- reactiveVal(FALSE)
  db_root   <- reactiveVal("")
  part_count_rv <- reactiveVal(NA_integer_)
  uniq_col_rv <- reactiveVal(NULL)
  has_antibody_rv <- reactiveVal(FALSE)
  
  output$db_info <- renderUI({
    tags$div(class = "summary-line", "No database loaded.")
  })
  
  observeEvent(input$db_folder_btn, {
    path_vec <- shinyFiles::parseDirPath(volumes, input$db_folder_btn)
    if (length(path_vec) == 0) return()
    path <- norm(as.character(path_vec[1]))
    if (!nzchar(path) || !dir.exists(path)) return()
    
    p <- list_partitions_from_folders(path)
    if (is.null(p)) {
      db_loaded(FALSE)
      db_root("")
      part_count_rv(NA_integer_)
      
      output$db_info <- renderUI({
        tags$div(
          tags$div(class = "summary-line", tags$b("No database was loaded.")),
          tags$div(class = "summary-line", paste0("Selected folder: ", path)),
          tags$div(class = "summary-line", "Reason: no partition folders found (Disease=/BSource=/BType=/Isotype=).")
        )
      })
      return()
    }
    
    # IMPORTANT: Disease choices are now decoded (match DuckDB output)
    updateSelectInput(session, "Disease", choices = c("None", "All", p$Disease), selected = "None")
    updateSelectInput(session, "BSource", choices = c("All", p$BSource), selected = "All")
    updateSelectInput(session, "BType",   choices = c("All", p$BType),   selected = "All")
    updateSelectInput(session, "Isotype", choices = c("All", p$Isotype), selected = "All")
    
    glob <- paste0(gsub("\\\\", "/", path), "/**/*.parquet")
    dbExecute(con, "DROP VIEW IF EXISTS v;")
    dbExecute(con, sprintf("CREATE VIEW v AS SELECT * FROM read_parquet('%s');", glob))
    
    uniq_col_rv(tryCatch(discover_uniqueness_col(con), error = function(e) NULL))
    cols <- tryCatch(dbGetQuery(con, "DESCRIBE v;")$column_name, error = function(e) character())
    has_antibody_rv("Antibody" %in% cols)
    
    db_loaded(TRUE)
    db_root(path)
    part_count_rv(p$n_part_dirs)
    
    output$db_info <- renderUI({
      tags$div(
        tags$div(class = "summary-line", "Database was loaded from selected folder:"),
        tags$div(class = "summary-indent summary-line", db_root()),
        tags$div(class = "summary-line", paste0("Detected ", fmt_int(part_count_rv()), " partition folders."))
      )
    })
  }, ignoreInit = TRUE)
  
  filters <- reactive({
    list(
      Disease = input$Disease,
      BSource = input$BSource,
      BType   = input$BType,
      Isotype = input$Isotype,
      Disease_specific   = input$Disease_specific,
      Peptide_uniqueness = input$Peptide_uniqueness,
      CDR3               = input$CDR3
    )
  })
  
  make_summary_lines <- function(f) {
    where <- build_where(f, uniq_col_rv())
    
    if (isTRUE(has_antibody_rv())) {
      cnt <- dbGetQuery(con, sprintf("SELECT COUNT(*) AS n_pep, COUNT(DISTINCT Antibody) AS n_ab FROM v %s;", where))
      n_pep <- cnt$n_pep[1]; n_ab <- cnt$n_ab[1]
    } else {
      cnt <- dbGetQuery(con, sprintf("SELECT COUNT(*) AS n_pep FROM v %s;", where))
      n_pep <- cnt$n_pep[1]; n_ab <- NA
    }
    
    fn_df <- tryCatch(
      dbGetQuery(con, sprintf("SELECT COUNT(DISTINCT filename) AS n_fn FROM v %s;", where)),
      error = function(e) data.frame(n_fn = NA)
    )
    n_fn <- fn_df$n_fn[1]
    
    export_cols <- c("Peptide","v_call","d_call","j_call","cdr3_aa","filename","Antibody")
    
    disease_full <- disease_fullname_for_summary(f$Disease)
    
    c(
      "DAT-DB export summary",
      paste0("Created_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("Database_root: ", db_root()),
      paste0("Partition_folders_detected: ", fmt_int(part_count_rv())),
      paste0("Export_columns: ", paste(export_cols, collapse = ", ")),
      paste0("Rows_exported_peptides: ", fmt_int(n_pep)),
      paste0("Distinct_antibodies: ", if (isTRUE(has_antibody_rv())) fmt_int(n_ab) else "Not available"),
      paste0("Distinct_filenames: ", fmt_int(n_fn)),
      "",
      "Filters:",
      paste0("  Disease: ", disease_full),
      paste0("  BSource: ", f$BSource),
      paste0("  BType  : ", f$BType),
      paste0("  Isotype: ", f$Isotype),
      paste0("  Disease-specific: ", f$Disease_specific),
      paste0("  Peptide uniqueness: ", f$Peptide_uniqueness),
      paste0("  CDR3: ", f$CDR3)
    )
  }
  
  output$summary <- renderUI({
    if (!isTRUE(db_loaded())) return(tags$div("No database loaded."))
    
    f <- filters()
    
    if (!filters_ready(f)) {
      return(tags$div(
        tags$div(class = "summary-line", tags$b("Current filters:")),
        tags$div(class = "summary-indent summary-line", paste0("Disease: ", f$Disease))
      ))
    }
    
    withProgress(message = "Running filter query…", value = 0, {
      incProgress(0.30, detail = "Counting")
      where <- build_where(f, uniq_col_rv())
      
      if (isTRUE(has_antibody_rv())) {
        cnt <- dbGetQuery(con, sprintf("SELECT COUNT(*) AS n_pep, COUNT(DISTINCT Antibody) AS n_ab FROM v %s;", where))
        n_pep <- cnt$n_pep[1]; n_ab <- cnt$n_ab[1]
      } else {
        cnt <- dbGetQuery(con, sprintf("SELECT COUNT(*) AS n_pep FROM v %s;", where))
        n_pep <- cnt$n_pep[1]; n_ab <- NA
      }
      incProgress(0.70, detail = "Done")
      
      disease_full <- disease_fullname_for_summary(f$Disease)
      
      tags$div(
        tags$div(class = "summary-line", tags$b("Number of peptides:"), paste0(" ", fmt_int(n_pep))),
        tags$div(class = "summary-line", tags$b("Number of antibodies:"), paste0(" ", if (isTRUE(has_antibody_rv())) fmt_int(n_ab) else "Not available")),
        tags$br(),
        tags$div(class = "summary-line", tags$b("Current filters:")),
        tags$div(class = "summary-indent summary-line", paste0("Disease: ", disease_full)),
        tags$div(class = "summary-indent summary-line", paste0("B cells source: ", f$BSource)),
        tags$div(class = "summary-indent summary-line", paste0("B cells type: ", f$BType)),
        tags$div(class = "summary-indent summary-line", paste0("Antibody isotype: ", f$Isotype)),
        tags$div(class = "summary-indent summary-line", paste0("Disease-specific: ", f$Disease_specific)),
        tags$div(class = "summary-indent summary-line", paste0("Peptide uniqueness: ", f$Peptide_uniqueness)),
        tags$div(class = "summary-indent summary-line", paste0("CDR3: ", f$CDR3))
      )
    })
  })
  
  output$preview <- renderDT({
    if (!isTRUE(db_loaded())) {
      return(datatable(empty_preview_df(), options = list(scrollX = TRUE)))
    }
    
    f <- filters()
    if (!filters_ready(f)) {
      return(datatable(empty_preview_df(), options = list(scrollX = TRUE)))
    }
    
    limit_n <- as.integer(input$max_preview)
    
    dat <- withProgress(message = "Running preview query…", value = 0, {
      incProgress(0.25, detail = "Building WHERE clause")
      where <- build_where(f, uniq_col_rv())
      
      incProgress(0.55, detail = paste0("Fetching first ", limit_n, " rows"))
      
      sel <- export_select_list(con)
      sql <- sprintf("SELECT %s FROM v %s LIMIT %d;", sel, where, limit_n)
      out <- tryCatch(dbGetQuery(con, sql), error = function(e) empty_preview_df())
      
      incProgress(0.20, detail = "Formatting")
      reorder_antibody_last(out)
    })
    
    datatable(dat, options = list(scrollX = TRUE, pageLength = min(20, nrow(dat))))
  })
  
  # --- Three download buttons (CSV | FASTA | Summary) ---
  output$download_btns <- renderUI({
    ready <- isTRUE(db_loaded()) && filters_ready(filters())
    if (!ready) {
      tagList(
        tags$button("CSV", class = "btn btn-default", disabled = "disabled"),
        tags$button("FASTA", class = "btn btn-default", disabled = "disabled"),
        tags$button("Summary", class = "btn btn-default", disabled = "disabled")
      )
    } else {
      tagList(
        downloadButton("download_csv", "CSV"),
        downloadButton("download_fasta", "FASTA"),
        downloadButton("download_summary", "Summary")
      )
    }
  })
  
  # --- CSV download ---
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("OAS_peptides_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      if (!isTRUE(db_loaded())) {
        writeLines("ERROR: No database loaded. Please Browse… to a valid folder first.", con = file)
        return()
      }
      f <- filters()
      if (!filters_ready(f)) {
        writeLines("ERROR: Please set Disease from 'None' before downloading.", con = file)
        return()
      }
      
      out_csv <- gsub("\\\\", "/", file)
      where   <- build_where(f, uniq_col_rv())
      sel     <- export_select_list(con)
      
      withProgress(message = "Exporting CSV…", value = 0, {
        incProgress(0.60, detail = "Writing CSV")
        sql_copy <- sprintf("
          COPY (
            SELECT %s FROM v
            %s
          )
          TO '%s'
          (FORMAT CSV, HEADER, DELIMITER ',');
        ", sel, where, out_csv)
        dbExecute(con, sql_copy)
        incProgress(0.40, detail = "Done")
      })
    }
  )
  
  # --- FASTA download ---
  output$download_fasta <- downloadHandler(
    filename = function() {
      paste0("OAS_peptides_", format(Sys.Date(), "%Y%m%d"), ".fasta")
    },
    content = function(file) {
      if (!isTRUE(db_loaded())) {
        writeLines("ERROR: No database loaded. Please Browse… to a valid folder first.", con = file)
        return()
      }
      f <- filters()
      if (!filters_ready(f)) {
        writeLines("ERROR: Please set Disease from 'None' before downloading.", con = file)
        return()
      }
      
      out_fa <- gsub("\\\\", "/", file)
      where  <- build_where(f, uniq_col_rv())
      
      withProgress(message = "Exporting FASTA…", value = 0, {
        incProgress(0.20, detail = "Fetching Peptide + filename")
        
        sql <- sprintf("
          SELECT
            Peptide,
            filename
          FROM v
          %s
          ORDER BY filename;
        ", where)
        
        dt <- tryCatch(as.data.table(dbGetQuery(con, sql)), error = function(e) data.table())
        if (nrow(dt) == 0) {
          writeLines("No rows matched the current filters (FASTA is empty).", con = out_fa)
          return()
        }
        
        if (!all(c("Peptide", "filename") %in% names(dt))) {
          writeLines("ERROR: Required columns (Peptide, filename) not found in the dataset.", con = out_fa)
          return()
        }
        
        dt <- dt[!is.na(Peptide) & nzchar(Peptide)]
        if (nrow(dt) == 0) {
          writeLines("No non-empty Peptide sequences matched the current filters (FASTA is empty).", con = out_fa)
          return()
        }
        
        dt[, filename := safe_header_token(filename)]
        dt[, abno := paste0("ab", seq_len(.N)), by = filename]
        dt[, hdr := paste0("OAS|", filename, "_", abno, "|")]
        
        incProgress(0.50, detail = "Writing FASTA (Biostrings)")
        aa <- AAStringSet(dt$Peptide)
        names(aa) <- dt$hdr
        writeXStringSet(aa, filepath = out_fa, format = "fasta")
        
        incProgress(0.30, detail = "Done")
      })
    }
  )
  
  # --- Summary TXT download ---
  output$download_summary <- downloadHandler(
    filename = function() {
      paste0("OAS_peptides_summary_", format(Sys.Date(), "%Y%m%d"), ".txt")
    },
    content = function(file) {
      if (!isTRUE(db_loaded())) {
        writeLines("ERROR: No database loaded. Please Browse… to a valid folder first.", con = file)
        return()
      }
      f <- filters()
      if (!filters_ready(f)) {
        writeLines("ERROR: Please set Disease from 'None' before downloading.", con = file)
        return()
      }
      
      out_txt <- gsub("\\\\", "/", file)
      
      withProgress(message = "Exporting summary…", value = 0, {
        incProgress(0.35, detail = "Computing summary")
        summary_lines <- make_summary_lines(f)
        
        incProgress(0.65, detail = "Writing TXT")
        write_summary_txt(out_txt, summary_lines)
      })
    }
  )
}

shinyApp(ui, server)
