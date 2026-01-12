# app.R
# Shiny app: query a partitioned Parquet folder (final_enriched) and download CSV
# Filters:
#   - Disease, BSource, BType, Isotype (from partition folders)
#   - Disease_specific:
#       Yes  => Disease_presence_peptide = 1
#       No   => Disease_presence_peptide != 1 (includes 0/NA/other)
#       All  => no filter
#   - Peptide_uniqueness:
#       Yes  => n_distinct_peptide = 1  (or falls back to n_distinct_antibodies if present)
#       No   => != 1
#       All  => no filter
#   - CDR3:
#       Yes  => peptide_in_cdr3 = 1
#       No   => peptide_in_cdr3 = 0 (or != 1)
#       All  => no filter
#
# Notes:
#   - Uses DuckDB for fast predicate pushdown + export.
#   - Gets partition choices by scanning folder names (fast; no parquet scan).
#   - Download uses DuckDB COPY TO CSV for speed (does not pull all rows into R).

suppressPackageStartupMessages({
  library(shiny)
  library(DBI)
  library(duckdb)
  library(DT)
})

# -------------------------
# Helpers
# -------------------------
norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

# Parse hive-style partitions from folder paths
list_partitions_from_folders <- function(root) {
  root <- norm(root)
  if (!dir.exists(root)) return(list())
  
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  part_dirs <- dirs[
    grepl("Disease=", dirs) &
      grepl("BSource=", dirs) &
      grepl("BType=", dirs) &
      grepl("Isotype=", dirs)
  ]
  if (length(part_dirs) == 0) return(list())
  
  get_val <- function(x, key) sub(paste0(".*", key, "=?([^/]+).*"), "\\1", x)
  
  Disease <- unique(get_val(part_dirs, "Disease"))
  BSource <- unique(get_val(part_dirs, "BSource"))
  BType   <- unique(get_val(part_dirs, "BType"))
  Isotype <- unique(get_val(part_dirs, "Isotype"))
  
  list(
    Disease = sort(Disease),
    BSource = sort(BSource),
    BType   = sort(BType),
    Isotype = sort(Isotype),
    n_part_dirs = length(part_dirs)
  )
}

# Quote string for SQL safely (basic)
sql_quote <- function(x) {
  x <- gsub("'", "''", x, fixed = TRUE)
  paste0("'", x, "'")
}

# Build WHERE clause from inputs
build_where <- function(filters, col_uniqueness = NULL) {
  w <- character()
  
  # partition filters
  if (!identical(filters$Disease, "All")) w <- c(w, paste0("Disease = ", sql_quote(filters$Disease)))
  if (!identical(filters$BSource, "All")) w <- c(w, paste0("BSource = ", sql_quote(filters$BSource)))
  if (!identical(filters$BType, "All"))   w <- c(w, paste0("BType = ", sql_quote(filters$BType)))
  if (!identical(filters$Isotype, "All")) w <- c(w, paste0("Isotype = ", sql_quote(filters$Isotype)))
  
  # disease_specific
  if (filters$Disease_specific == "Yes") {
    w <- c(w, "COALESCE(Disease_presence_peptide, 0) = 1")
  } else if (filters$Disease_specific == "No") {
    w <- c(w, "COALESCE(Disease_presence_peptide, 0) <> 1")
  }
  
  # peptide uniqueness (n_distinct_peptide preferred; fallback to n_distinct_antibodies)
  if (!is.null(col_uniqueness)) {
    if (filters$Peptide_uniqueness == "Yes") {
      w <- c(w, sprintf("COALESCE(%s, 0) = 1", col_uniqueness))
    } else if (filters$Peptide_uniqueness == "No") {
      w <- c(w, sprintf("COALESCE(%s, 0) <> 1", col_uniqueness))
    }
  }
  
  # CDR3
  if (filters$CDR3 == "Yes") {
    w <- c(w, "COALESCE(peptide_in_cdr3, 0) = 1")
  } else if (filters$CDR3 == "No") {
    w <- c(w, "COALESCE(peptide_in_cdr3, 0) <> 1")
  }
  
  if (length(w) == 0) return("")
  paste("WHERE", paste(w, collapse = " AND "))
}

# Discover which uniqueness column exists in the dataset
discover_uniqueness_col <- function(con) {
  cols <- dbGetQuery(con, "DESCRIBE v;")$column_name
  if ("n_distinct_peptide" %in% cols) return("n_distinct_peptide")
  if ("n_distinct_antibodies" %in% cols) return("n_distinct_antibodies")
  return(NULL)
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  titlePanel("PeptideDB final_enriched Query (DuckDB + Parquet)"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("db_folder", "Parquet folder (final_enriched):", value = "D:/OAS/PeptideDB/final_enriched"),
      actionButton("load_db", "Load / Refresh partitions"),
      
      tags$hr(),
      
      uiOutput("partition_info"),
      
      selectInput("Disease",  "Disease",  choices = c("All")),
      selectInput("BSource",  "BSource",  choices = c("All")),
      selectInput("BType",    "BType",    choices = c("All")),
      selectInput("Isotype",  "Isotype",  choices = c("All")),
      
      tags$hr(),
      
      selectInput("Disease_specific", "Disease_specific (Disease_presence_peptide = 1?)",
                  choices = c("All", "Yes", "No"), selected = "All"),
      
      selectInput("Peptide_uniqueness", "Peptide_uniqueness (=1?)",
                  choices = c("All", "Yes", "No"), selected = "All"),
      
      selectInput("CDR3", "CDR3 (peptide_in_cdr3 = 1?)",
                  choices = c("All", "Yes", "No"), selected = "All"),
      
      tags$hr(),
      
      numericInput("max_preview", "Preview rows (table):", value = 200, min = 10, max = 5000, step = 10),
      numericInput("max_export",  "Max rows to export (0 = no limit):", value = 0, min = 0, max = 500000000, step = 10000),
      
      actionButton("run_preview", "Run Preview"),
      downloadButton("download_csv", "Download CSV")
    ),
    
    mainPanel(
      verbatimTextOutput("status"),
      DTOutput("preview")
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  # Keep a DuckDB connection alive for the app session
  con_rv <- reactiveVal(NULL)
  
  observeEvent(TRUE, {
    # create on startup
    con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
    dbExecute(con, "PRAGMA threads=6;")
    dbExecute(con, "PRAGMA enable_progress_bar=false;")
    con_rv(con)
    
    session$onSessionEnded(function() {
      try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    })
  }, once = TRUE)
  
  # Partition choices state
  parts_rv <- reactiveVal(NULL)
  
  output$partition_info <- renderUI({
    p <- parts_rv()
    if (is.null(p)) return(tags$small("Partitions not loaded yet. Click 'Load / Refresh partitions'."))
    tags$small(sprintf("Detected %d partition folders.", p$n_part_dirs))
  })
  
  # Load partitions from folder structure
  observeEvent(input$load_db, {
    root <- norm(input$db_folder)
    if (!dir.exists(root)) {
      parts_rv(NULL)
      showNotification("Folder does not exist. Check path.", type = "error")
      return()
    }
    
    p <- list_partitions_from_folders(root)
    if (length(p) == 0) {
      parts_rv(NULL)
      showNotification("No partition folders found (Disease=/BSource=/BType=/Isotype=).", type = "error")
      return()
    }
    parts_rv(p)
    
    updateSelectInput(session, "Disease", choices = c("All", p$Disease), selected = "All")
    updateSelectInput(session, "BSource", choices = c("All", p$BSource), selected = "All")
    updateSelectInput(session, "BType",   choices = c("All", p$BType),   selected = "All")
    updateSelectInput(session, "Isotype", choices = c("All", p$Isotype), selected = "All")
    
    showNotification("Partitions loaded.", type = "message")
  }, ignoreInit = TRUE)
  
  # Ensure we have a DuckDB view pointing to the selected folder
  ensure_view <- reactive({
    con <- con_rv()
    root <- norm(input$db_folder)
    if (is.null(con)) stop("DuckDB connection not ready.")
    if (!dir.exists(root)) stop("Folder does not exist.")
    
    glob <- paste0(gsub("\\\\", "/", root), "/**/*.parquet")
    
    # Recreate view each time folder changes (cheap)
    dbExecute(con, "DROP VIEW IF EXISTS v;")
    dbExecute(con, sprintf("CREATE VIEW v AS SELECT * FROM read_parquet('%s');", glob))
    
    # Determine uniqueness col available
    uniq_col <- discover_uniqueness_col(con)
    
    list(con = con, uniq_col = uniq_col)
  })
  
  # Build current filter list
  current_filters <- reactive({
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
  
  # Status text
  output$status <- renderText({
    root <- norm(input$db_folder)
    if (!dir.exists(root)) return("Set a valid folder path and click 'Load / Refresh partitions'.")
    p <- parts_rv()
    if (is.null(p)) return("Folder found. Click 'Load / Refresh partitions' to detect partition values.")
    paste0(
      "DB folder: ", root, "\n",
      "Partition folders detected: ", p$n_part_dirs, "\n",
      "Filters: ",
      "Disease=", input$Disease, "; ",
      "BSource=", input$BSource, "; ",
      "BType=", input$BType, "; ",
      "Isotype=", input$Isotype, "; ",
      "Disease_specific=", input$Disease_specific, "; ",
      "Peptide_uniqueness=", input$Peptide_uniqueness, "; ",
      "CDR3=", input$CDR3
    )
  })
  
  # Preview query
  observeEvent(input$run_preview, {
    v <- ensure_view()
    con <- v$con
    uniq_col <- v$uniq_col
    
    if (is.null(uniq_col) && input$Peptide_uniqueness != "All") {
      showNotification("Uniqueness column not found (expected n_distinct_peptide or n_distinct_antibodies). Uniqueness filter ignored.",
                       type = "warning")
    }
    
    where_sql <- build_where(current_filters(), col_uniqueness = uniq_col)
    
    limit_n <- as.integer(input$max_preview)
    sql <- sprintf("SELECT * FROM v %s LIMIT %d;", where_sql, limit_n)
    
    # For preview, fetch into R
    dat <- tryCatch(
      dbGetQuery(con, sql),
      error = function(e) {
        showNotification(paste("Query error:", e$message), type = "error")
        NULL
      }
    )
    
    output$preview <- renderDT({
      if (is.null(dat)) return(datatable(data.frame()))
      datatable(dat, options = list(scrollX = TRUE, pageLength = min(50, nrow(dat))))
    })
  }, ignoreInit = TRUE)
  
  # Download handler (fast export via DuckDB COPY)
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("final_enriched_query_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      v <- ensure_view()
      con <- v$con
      uniq_col <- v$uniq_col
      
      where_sql <- build_where(current_filters(), col_uniqueness = uniq_col)
      
      max_export <- as.numeric(input$max_export)
      limit_sql <- ""
      if (!is.na(max_export) && max_export > 0) {
        limit_sql <- paste0("LIMIT ", as.integer(max_export))
      }
      
      # Use COPY TO for fast CSV writing
      # Note: DuckDB wants forward slashes
      out_file <- gsub("\\\\", "/", file)
      
      sql_copy <- sprintf("
        COPY (
          SELECT * FROM v
          %s
          %s
        )
        TO '%s'
        (FORMAT CSV, HEADER, DELIMITER ',');
      ", where_sql, limit_sql, out_file)
      
      dbExecute(con, sql_copy)
    }
  )
}

shinyApp(ui, server)
