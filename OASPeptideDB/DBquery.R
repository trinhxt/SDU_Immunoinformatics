# ==============================================================================
# DAT-DB Shiny Viewer (v4.4 - About Text Update)
# ------------------------------------------------------------------------------
# Features:
# 1. Auto-loads "db_stage3" recursively.
# 2. Implements filters: Disease, Metadata, Specificity, Uniqueness, Patients.
# 3. CDR3 Spanning Filter.
# 4. FASTA/CSV/Summary Export: Located below Live Summary.
# 5. UI: Updated About tab with project description.
# ==============================================================================

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
norm_path <- function(p) {
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

fmt_int <- function(x) {
  if (is.null(x) || is.na(x)) return("0")
  format(as.numeric(x), big.mark = ",", scientific = FALSE)
}

# --- Disease Mapping Dictionary ---
disease_code_to_full <- c(
  "AChR-MG"                           = "Acetylcholine receptor-positive myasthenia gravis",
  "Allergic-Rhinitis-In-Season"       = "Allergic rhinitis in season",
  "Allergic-Rhinitis-Out-Of-Season"   = "Allergic rhinitis out of season",
  "Allergy%2FNoSIT"                   = "Allergy NoSIT",
  "Allergy/NoSIT"                     = "Allergy NoSIT",
  "Allergy%2FSIT"                     = "Allergy SIT",
  "Allergy/SIT"                       = "Allergy SIT",
  "Asthma"                            = "Asthma",
  "CLL"                               = "Chronic lymphocytic leukemia",
  "CMV"                               = "Cytomegalovirus",
  "CMV%2FEBV"                         = "Cytomegalovirus/Epstein–Barr virus",
  "CMV/EBV"                           = "Cytomegalovirus/Epstein–Barr virus",
  "Dengue"                            = "Dengue",
  "Ebola"                             = "Ebola",
  "EBV"                               = "Epstein–Barr virus",
  "HCV"                               = "Hepatitis C virus",
  "Healthy%2Fceliac-disease"          = "Celiac",
  "Healthy/celiac-disease"            = "Celiac",
  "HIV"                               = "HIV",
  "Light-Chain-Amyloidosis"           = "Light chain amyloidosis",
  "MS"                                = "Multiple sclerosis",
  "MuSK-MG"                           = "Muscle-specific kinase myasthenia gravis",
  "Non-Dengue-Febrile-Illness"        = "Non-dengue febrile illness",
  "Obstructive-Sleep-Apnea"           = "Obstructive sleep apnea",
  "SARS-COV-2"                        = "SARS-COV-2",
  "SLE"                               = "Systemic lupus erythematosus",
  "Tonsillitis"                       = "Tonsillitis",
  "Tonsillitis%2FObstructive-Sleep-Apnea" = "Tonsillitis/Obstructive sleep apnea",
  "Tonsillitis/Obstructive-Sleep-Apnea"   = "Tonsillitis/Obstructive sleep apnea"
)

get_full_disease_name <- function(code) {
  if (is.null(code) || code %in% c("None", "All")) return(code)
  if (code %in% names(disease_code_to_full)) return(disease_code_to_full[[code]])
  decoded <- tryCatch(utils::URLdecode(code), error = function(e) code)
  if (decoded %in% names(disease_code_to_full)) return(disease_code_to_full[[decoded]])
  return(decoded)
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  useShinyjs(),
  theme = shinytheme("cerulean"),
  
  tags$head(tags$style(HTML("
    body { font-family: Arial, sans-serif; }
    .panel-container { 
      border: 1px solid #ddd; 
      padding: 15px; 
      border-radius: 5px; 
      background: #fcfcfc; 
      margin-bottom: 15px; 
      box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }
    .custom-header { 
      font-size: 16px; 
      font-weight: bold; 
      border-left: 4px solid #007bc2; 
      padding-left: 10px; 
      margin-bottom: 15px; 
      color: #333; 
    }
    .status-ok { color: green; font-weight: bold; }
    .status-err { color: red; font-weight: bold; }
    .summary-stat { 
      font-size: 13px; 
      margin-bottom: 2px; 
      color: #444;
      display: flex;
      justify-content: space-between;
      border-bottom: 1px dotted #eee;
    }
    .summary-stat strong { color: #222; margin-right: 10px; text-align: right;}
    .btn-download { width: auto; margin-right: 10px; margin-bottom: 5px; }
  "))),
  
  titlePanel(HTML("<strong>DAT-DB</strong>: Disease-Specific Antibody Peptide Database")),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # --- 1. Data Selection ---
      div(class = "panel-container",
          div(class = "custom-header", "1. Database Source"),
          uiOutput("db_status_ui"),
          hidden(shinyDirButton("db_folder_btn", "Browse Folder", 
                                "Select the 'db_stage3' folder"))
      ),
      
      # --- 2. Clinical Filters ---
      div(class = "panel-container",
          div(class = "custom-header", "2. Clinical Filters"),
          selectInput("filter_disease", "Disease", choices = c("None"), selected = "None"),
          selectInput("filter_bsource", "B-Cell Source", choices = "PBMC", selected = "PBMC"),
          selectInput("filter_btype", "B-Cell Type", choices = "All", selected = "All"),
          selectInput("filter_isotype", "Isotype", choices = "All", selected = "All")
      ),
      
      # --- 3. Advanced Logic Filters ---
      div(class = "panel-container",
          div(class = "custom-header", "3. Specificity & Quality"),
          
          # CDR3 Spanning Filter
          radioButtons("logic_cdr3", "CDR3 Spanning",
                       choices = c("All" = "all",
                                   "None (0 AA)" = "none",
                                   "Low (1-3 AA)" = "low",
                                   "High (>3 AA)" = "high"),
                       selected = "high"),
          
          # Disease Specificity
          radioButtons("logic_specificity", "Disease Specificity",
                       choices = c("All" = "all",
                                   "Disease-Specific (N_Diseases=1)" = "specific",
                                   "Shared (N_Diseases>1)" = "shared"),
                       selected = "specific"),
          
          # Uniqueness
          radioButtons("logic_uniqueness", "Peptide Uniqueness",
                       choices = c("All" = "all",
                                   "Unique (N_Antibodies=1)" = "unique",
                                   "Common (N_Antibodies>1)" = "common"),
                       selected = "all"),
          
          # Patient Count
          radioButtons("logic_patients", "Patient Frequency",
                       choices = c("All" = "all",
                                   "Singleton (N_Patients=1)" = "single",
                                   "Recurring (N_Patients>1)" = "multi"),
                       selected = "multi")
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Data Preview", 
                 br(),
                 div(class = "panel-container",
                     div(class = "custom-header", "Live Summary"),
                     uiOutput("summary_stats"),
                     hr(),
                     h4("Download Results", style="font-size:14px; font-weight:bold; color:#555; margin-bottom:10px;"),
                     uiOutput("download_ui")
                 ),
                 div(class = "panel-container",
                     DTOutput("preview_table")
                 )
        ),
        tabPanel("About",
                 br(),
                 # --- NEW INTRODUCTION TEXT ---
                 div(class = "panel-container",
                     h4("Welcome to DAT-DB app!"),
                     p("The identification of human antibodies using bottom-up proteomics relies on database searches that match experimental peptide fragments to theoretical values derived from protein sequences in databases. Standard protein databases, such as UniProt and NCBI-RefSeq, contain a limited number of antibody sequences compared to the human body’s capacity to generate billions and currently lack disease-specific antibody sequences. This limitation could lead to the misidentification of antibodies in samples. To overcome this challenge, DAT-DB, a database of disease-specific antibodies, provides researchers with antibody tryptic peptides derived from next-generation sequencing of antibody repertoires.")
                 ),
                 # -----------------------------
                 div(class = "panel-container",
                     h4("Database Status"),
                     textOutput("debug_path"),
                     hr(),
                     h4("Filter Definitions"),
                     tags$ul(
                       tags$li(strong("CDR3 Spanning:"), " Length of peptide overlap with the hypervariable CDR3 region."),
                       tags$li(strong("Disease-Specific:"), " Found ONLY in the selected disease."),
                       tags$li(strong("Unique Clone:"), " Found in exactly one antibody sequence."),
                       tags$li(strong("Recurring:"), " Found in at least 2 different patients.")
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
  
  output$debug_path <- renderText({
    paste("Current Working Directory:", getwd())
  })
  
  con <- dbConnect(duckdb())
  dbExecute(con, "PRAGMA memory_limit='4GB'")
  dbExecute(con, "PRAGMA threads=4")
  
  session$onSessionEnded(function() {
    dbDisconnect(con, shutdown = TRUE)
  })
  
  vals <- reactiveValues(
    db_loaded = FALSE,
    db_path = NULL
  )
  
  # --- Core Loader Function ---
  load_database <- function(path) {
    withProgress(message = 'Connecting to Database...', value = 0.1, {
      tryCatch({
        dbExecute(con, "DROP VIEW IF EXISTS data")
        
        incProgress(0.2, detail = "Scanning parquet files...")
        clean_path <- gsub("\\\\", "/", path)
        glob <- file.path(clean_path, "**", "*.parquet")
        
        incProgress(0.3, detail = "Creating virtual views...")
        dbExecute(con, sprintf("CREATE VIEW data AS SELECT * FROM read_parquet('%s')", glob))
        
        incProgress(0.5, detail = "Reading Disease list...")
        dis <- dbGetQuery(con, "SELECT DISTINCT Disease FROM data WHERE Disease IS NOT NULL ORDER BY Disease")
        
        if (nrow(dis) == 0) {
          showNotification("Database connected but contains no data/diseases.", type="error")
          return()
        }
        
        incProgress(0.8, detail = "Reading Metadata...")
        updateSelectInput(session, "filter_disease", choices = c("None", dis$Disease), selected = "None")
        
        iso <- dbGetQuery(con, "SELECT DISTINCT Isotype FROM data WHERE Isotype IS NOT NULL ORDER BY Isotype")
        updateSelectInput(session, "filter_isotype", choices = c("All", iso$Isotype), selected = "All")
        
        bs <- dbGetQuery(con, "SELECT DISTINCT BSource FROM data WHERE BSource IS NOT NULL ORDER BY BSource")
        updateSelectInput(session, "filter_bsource", choices = c("All", bs$BSource), selected = "PBMC")
        
        bt <- dbGetQuery(con, "SELECT DISTINCT BType FROM data WHERE BType IS NOT NULL ORDER BY BType")
        updateSelectInput(session, "filter_btype", choices = c("All", bt$BType), selected = "All")
        
        vals$db_loaded <- TRUE
        vals$db_path <- path
        shinyjs::hide("db_folder_btn") 
        showNotification("Database Loaded Successfully!", type="message")
        
      }, error = function(e) {
        showNotification(paste("Error loading DB:", e$message), type = "error")
        vals$db_loaded <- FALSE
        shinyjs::show("db_folder_btn")
      })
    })
  }
  
  observe({
    cat("Checking for DB in:", getwd(), "\n")
    candidates <- c("db_stage3", "OASpepDB", "OASpeptideDB")
    found_path <- NULL
    for (p in candidates) {
      if (dir.exists(p)) {
        files <- list.files(p, pattern = "\\.parquet$", recursive = TRUE)
        if (length(files) > 0) {
          found_path <- p
          break
        }
      }
    }
    if (!is.null(found_path)) {
      cat("Auto-loading found DB:", found_path, "\n")
      load_database(found_path)
    } else {
      cat("No DB found. Waiting for manual selection.\n")
      shinyjs::show("db_folder_btn")
    }
  })
  
  roots <- c(Computer = "/")
  if (.Platform$OS.type == "windows") roots <- shinyFiles::getVolumes()()
  shinyDirChoose(input, "db_folder_btn", roots = roots)
  
  observeEvent(input$db_folder_btn, {
    req(input$db_folder_btn)
    selected <- shinyFiles::parseDirPath(roots, input$db_folder_btn)
    if (length(selected) > 0) {
      load_database(norm_path(selected))
    }
  })
  
  output$db_status_ui <- renderUI({
    if (vals$db_loaded) {
      tagList(
        div(class="status-ok", icon("check-circle"), " Database Loaded"),
        div(style="font-size:11px; color:#666; word-break:break-all;", vals$db_path)
      )
    } else {
      div(class="status-err", icon("exclamation-circle"), " No Database Found (Select Folder)")
    }
  })
  
  # --- Query Builder ---
  get_sql_query <- function(count_only = FALSE, limit = NULL) {
    req(vals$db_loaded)
    req(input$filter_disease != "None")
    
    if (count_only) {
      sel <- "COUNT(*) as n, COUNT(DISTINCT Peptide) as n_pep"
    } else {
      sel <- "Peptide, Disease, Patient, Isotype, 
              CDR3_spanning_count, N_Patients, N_Diseases, N_Antibodies,
              cdr3_aa, Antibody, filename"
    }
    
    q <- sprintf("SELECT %s FROM data WHERE 1=1", sel)
    
    if (input$filter_disease != "All") q <- paste0(q, sprintf(" AND Disease = '%s'", input$filter_disease))
    if (input$filter_bsource != "All") q <- paste0(q, sprintf(" AND BSource = '%s'", input$filter_bsource))
    if (input$filter_btype   != "All") q <- paste0(q, sprintf(" AND BType   = '%s'", input$filter_btype))
    if (input$filter_isotype != "All") q <- paste0(q, sprintf(" AND Isotype = '%s'", input$filter_isotype))
    
    if (input$logic_cdr3 == "none") {
      q <- paste0(q, " AND CDR3_spanning_count = 0")
    } else if (input$logic_cdr3 == "low") {
      q <- paste0(q, " AND CDR3_spanning_count BETWEEN 1 AND 3")
    } else if (input$logic_cdr3 == "high") {
      q <- paste0(q, " AND CDR3_spanning_count > 3")
    }
    
    if (input$logic_specificity == "specific") q <- paste0(q, " AND N_Diseases = 1")
    if (input$logic_specificity == "shared")   q <- paste0(q, " AND N_Diseases > 1")
    
    if (input$logic_uniqueness == "unique")    q <- paste0(q, " AND N_Antibodies = 1")
    if (input$logic_uniqueness == "common")    q <- paste0(q, " AND N_Antibodies > 1")
    
    if (input$logic_patients == "single")      q <- paste0(q, " AND N_Patients = 1")
    if (input$logic_patients == "multi")       q <- paste0(q, " AND N_Patients > 1")
    
    if (!is.null(limit) && !count_only) q <- paste0(q, " LIMIT ", limit)
    
    return(q)
  }
  
  # --- Outputs ---
  output$summary_stats <- renderUI({
    if (!vals$db_loaded) return(div("Please load a database."))
    if (input$filter_disease == "None") return(div("Select a Disease to begin."))
    
    withProgress(message = 'Filtering data...', value = 0.5, {
      sql <- get_sql_query(count_only = TRUE)
      res <- tryCatch(dbGetQuery(con, sql), error = function(e) data.frame(n=0, n_pep=0))
    })
    
    map_dict <- list(
      "all" = "All", "none" = "None (0 AA)", "low" = "Low (1-3 AA)", "high" = "High (>3 AA)",
      "specific" = "Disease-Specific", "shared" = "Shared",
      "unique" = "Unique Clone", "common" = "Common",
      "single" = "Singleton", "multi" = "Recurring"
    )
    
    txt_cdr3 <- map_dict[[input$logic_cdr3]] %||% input$logic_cdr3
    txt_spec <- map_dict[[input$logic_specificity]] %||% input$logic_specificity
    txt_uniq <- map_dict[[input$logic_uniqueness]] %||% input$logic_uniqueness
    txt_pat  <- map_dict[[input$logic_patients]] %||% input$logic_patients
    full_disease <- get_full_disease_name(input$filter_disease)
    
    tagList(
      h4(strong("Distinct Peptides: "), fmt_int(res$n_pep), style="color: #007bc2; margin-top:0;"),
      div(style="font-size: 11px; color:#666; margin-bottom:8px;", paste("Total Rows Matched:", fmt_int(res$n))),
      hr(style="margin: 8px 0;"),
      div(class="summary-stat", span("Disease:"), strong(full_disease)),
      div(class="summary-stat", span("B-Source:"), strong(input$filter_bsource)),
      div(class="summary-stat", span("B-Type:"), strong(input$filter_btype)),
      div(class="summary-stat", span("Isotype:"), strong(input$filter_isotype)),
      hr(style="margin: 8px 0;"),
      div(class="summary-stat", span("CDR3 Overlap:"), strong(txt_cdr3)),
      div(class="summary-stat", span("Specificity:"), strong(txt_spec)),
      div(class="summary-stat", span("Uniqueness:"), strong(txt_uniq)),
      div(class="summary-stat", span("Patient Freq:"), strong(txt_pat))
    )
  })
  
  output$preview_table <- renderDT({
    if (!vals$db_loaded || input$filter_disease == "None") return(NULL)
    withProgress(message = 'Fetching preview...', value = 0.5, {
      sql <- get_sql_query(limit = 100)
      df <- dbGetQuery(con, sql)
    })
    datatable(df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # --- Output: Download Buttons ---
  output$download_ui <- renderUI({
    if (!vals$db_loaded || input$filter_disease == "None") return(div("Select Disease to enable downloads."))
    tagList(
      downloadButton("dl_csv", "Download CSV", class="btn-download"),
      downloadButton("dl_fasta", "Download FASTA", class="btn-download"),
      downloadButton("dl_summary", "Download Summary", class="btn-download")
    )
  })
  
  # --- 1. CSV Download ---
  output$dl_csv <- downloadHandler(
    filename = function() { paste0("OAS_Export_", Sys.Date(), ".csv") },
    content = function(file) {
      withProgress(message = "Exporting CSV...", value = 0.1, {
        sql <- get_sql_query()
        incProgress(0.2, detail = "Executing DuckDB COPY command...")
        tmp_sql <- sprintf("COPY (%s) TO '%s' (FORMAT CSV, HEADER)", sql, norm_path(file))
        dbExecute(con, tmp_sql)
        incProgress(0.7, detail = "Done.")
      })
    }
  )
  
  # --- 2. FASTA Download ---
  output$dl_fasta <- downloadHandler(
    filename = function() { paste0("OAS_Export_", Sys.Date(), ".fasta") },
    content = function(file) {
      withProgress(message = "Generating FASTA...", value = 0, {
        incProgress(0.2, detail = "Fetching data from DB...")
        sql <- get_sql_query()
        df <- dbGetQuery(con, sql)
        
        if (nrow(df) == 0) {
          writeLines("No data found.", file)
          return()
        }
        
        incProgress(0.4, detail = "Aggregating unique peptides...")
        dt <- as.data.table(df)
        dt_agg <- dt[, .(All_Filenames = paste(unique(filename), collapse = ";")), by = Peptide]
        dt_agg[, N := .I]
        dt_agg[, header := paste0("pep_", N, "|", All_Filenames, "|")]
        
        incProgress(0.3, detail = "Writing file...")
        if (nrow(dt_agg) > 0) {
          seqs <- AAStringSet(dt_agg$Peptide)
          names(seqs) <- dt_agg$header
          writeXStringSet(seqs, file)
        } else {
          writeLines("No sequences to export.", file)
        }
        incProgress(0.1, detail = "Done.")
      })
    }
  )
  
  # --- 3. Summary Text Download ---
  output$dl_summary <- downloadHandler(
    filename = function() { paste0("OAS_Summary_", Sys.Date(), ".txt") },
    content = function(file) {
      sql <- get_sql_query(count_only = TRUE)
      counts <- dbGetQuery(con, sql)
      
      map_dict <- list(
        "all" = "All", "none" = "None (0 AA)", "low" = "Low (1-3 AA)", "high" = "High (>3 AA)",
        "specific" = "Disease-Specific", "shared" = "Shared",
        "unique" = "Unique Clone", "common" = "Common",
        "single" = "Singleton", "multi" = "Recurring"
      )
      get_lab <- function(x) map_dict[[x]] %||% x
      
      lines <- c(
        "================================================================",
        " DAT-DB EXPORT SUMMARY",
        "================================================================",
        paste("Date:              ", Sys.time()),
        paste("Database Path:     ", vals$db_path),
        "",
        "----------------------------------------------------------------",
        " STATISTICS",
        "----------------------------------------------------------------",
        paste("Distinct Peptides: ", fmt_int(counts$n_pep)),
        paste("Total Rows Matched:", fmt_int(counts$n)),
        "",
        "----------------------------------------------------------------",
        " FILTERS APPLIED",
        "----------------------------------------------------------------",
        paste("Disease:           ", get_full_disease_name(input$filter_disease)),
        paste("B-Cell Source:     ", input$filter_bsource),
        paste("B-Cell Type:       ", input$filter_btype),
        paste("Isotype:           ", input$filter_isotype),
        "",
        paste("CDR3 Spanning:     ", get_lab(input$logic_cdr3)),
        paste("Disease Specificity:", get_lab(input$logic_specificity)),
        paste("Peptide Uniqueness:", get_lab(input$logic_uniqueness)),
        paste("Patient Frequency: ", get_lab(input$logic_patients)),
        "================================================================"
      )
      
      writeLines(lines, file)
    }
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui, server)