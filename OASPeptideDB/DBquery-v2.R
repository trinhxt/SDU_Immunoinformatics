# ==============================================================================
# DAT-DB Shiny Viewer (v2.1 - Performance Optimized)
# ------------------------------------------------------------------------------
# Optimizations:
# 1. CTE Query Structure: Forces filtering by Disease BEFORE joining.
# 2. Increased Memory Limit: Prevents disk spilling during Joins.
# 3. Asynchronous Filters: Metrics filters only apply after the join reduces data.
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
    .summary-stat { font-size: 14px; margin-bottom: 4px; }
  "))),
  
  titlePanel(HTML("<strong>DAT-DB</strong>: Disease-Specific Antibody Peptide Database")),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # --- 1. Data Selection ---
      div(class = "panel-container",
          div(class = "custom-header", "1. Database Source"),
          shinyDirButton("db_folder_btn", "Browse Folder", 
                         "Select the folder containing 'observations' and 'peptide_metrics'"),
          br(), br(),
          uiOutput("db_status_ui")
      ),
      
      # --- 2. Clinical Filters (Observations) ---
      div(class = "panel-container",
          div(class = "custom-header", "2. Clinical Filters"),
          # Note: Disease is the Primary Partition Key. Filtering this is instant.
          selectInput("filter_disease", "Disease (Required)", choices = c("None"), selected = "None"),
          selectInput("filter_isotype", "Isotype", choices = "All", selected = "All"),
          selectInput("filter_bsource", "B-Cell Source", choices = "All", selected = "All"),
          selectInput("filter_btype", "B-Cell Type", choices = "All", selected = "All")
      ),
      
      # --- 3. Metric Filters (Dimension) ---
      div(class = "panel-container",
          div(class = "custom-header", "3. Quality Metrics"),
          p(style="font-size:11px; color:#555;", "Filters applied AFTER Disease selection."),
          
          sliderInput("min_patients", "Min Patients", min = 1, max = 50, value = 1, step = 1),
          sliderInput("min_diseases", "Min Diseases", min = 1, max = 20, value = 1, step = 1),
          sliderInput("min_antibodies", "Min Antibodies", min = 1, max = 50, value = 1, step = 1)
      ),
      
      # --- 4. Downloads ---
      div(class = "panel-container",
          div(class = "custom-header", "4. Download Results"),
          uiOutput("download_ui")
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Data Preview", 
                 br(),
                 div(class = "panel-container",
                     div(class = "custom-header", "Live Summary"),
                     uiOutput("summary_stats")
                 ),
                 div(class = "panel-container",
                     DTOutput("preview_table")
                 )
        ),
        tabPanel("About",
                 br(),
                 div(class = "panel-container",
                     h4("Performance Note"),
                     p("This viewer uses a Star Schema (Fact + Dimension). To ensure fast queries:"),
                     tags$ul(
                       tags$li("Always select a ", strong("Disease"), " first."),
                       tags$li("DuckDB performs a 'Partition Pruning' scan on the Observations table first, reducing 2B rows to a manageable subset, before joining with Metrics.")
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
  
  # --- Setup DuckDB Connection ---
  con <- dbConnect(duckdb())
  
  # PERFORMANCE TUNING
  # 1. Increase Memory: 2GB is too small for Joins. 8GB+ is recommended if available.
  #    If user has 64GB RAM, we can safely use 16GB.
  dbExecute(con, "PRAGMA memory_limit='16GB'") 
  dbExecute(con, "PRAGMA threads=8")
  
  session$onSessionEnded(function() {
    dbDisconnect(con, shutdown = TRUE)
  })
  
  # --- Reactive State ---
  vals <- reactiveValues(
    db_loaded = FALSE,
    db_path = NULL
  )
  
  # --- Folder Selection Logic ---
  roots <- c(Computer = "/")
  if (.Platform$OS.type == "windows") roots <- shinyFiles::getVolumes()()
  shinyDirChoose(input, "db_folder_btn", roots = roots)
  
  observeEvent(input$db_folder_btn, {
    req(input$db_folder_btn)
    selected <- shinyFiles::parseDirPath(roots, input$db_folder_btn)
    if (length(selected) == 0) return()
    path <- norm_path(selected)
    
    # Check root or 'data' subfolder
    candidates <- c(path, file.path(path, "data"))
    valid_root <- NULL
    for (p in candidates) {
      if (dir.exists(file.path(p, "observations")) && dir.exists(file.path(p, "peptide_metrics"))) {
        valid_root <- p
        break
      }
    }
    
    if (!is.null(valid_root)) {
      tryCatch({
        dbExecute(con, "DROP VIEW IF EXISTS observations")
        dbExecute(con, "DROP VIEW IF EXISTS metrics")
        
        # Create Views
        obs_glob <- file.path(valid_root, "observations", "**", "*.parquet")
        met_glob <- file.path(valid_root, "peptide_metrics", "**", "*.parquet")
        
        dbExecute(con, sprintf("CREATE VIEW observations AS SELECT * FROM read_parquet('%s')", obs_glob))
        dbExecute(con, sprintf("CREATE VIEW metrics      AS SELECT * FROM read_parquet('%s')", met_glob))
        
        # Load Filters (Fast Distinct)
        # We assume Disease is the primary partition key
        dis <- dbGetQuery(con, "SELECT DISTINCT Disease FROM observations WHERE Disease IS NOT NULL ORDER BY Disease")
        updateSelectInput(session, "filter_disease", choices = c("None", "All", dis$Disease), selected = "None")
        
        iso <- dbGetQuery(con, "SELECT DISTINCT Isotype FROM observations WHERE Isotype IS NOT NULL ORDER BY Isotype")
        updateSelectInput(session, "filter_isotype", choices = c("All", iso$Isotype), selected = "All")
        
        bs <- dbGetQuery(con, "SELECT DISTINCT BSource FROM observations WHERE BSource IS NOT NULL ORDER BY BSource")
        updateSelectInput(session, "filter_bsource", choices = c("All", bs$BSource), selected = "All")
        
        bt <- dbGetQuery(con, "SELECT DISTINCT BType FROM observations WHERE BType IS NOT NULL ORDER BY BType")
        updateSelectInput(session, "filter_btype", choices = c("All", bt$BType), selected = "All")
        
        vals$db_loaded <- TRUE
        vals$db_path <- valid_root
        
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        vals$db_loaded <- FALSE
      })
    } else {
      showNotification("Invalid Folder! Missing 'observations' and 'peptide_metrics'.", type = "error")
    }
  })
  
  output$db_status_ui <- renderUI({
    if (vals$db_loaded) {
      tagList(div(class="status-ok", icon("check-circle"), " DB Connected"),
              div(style="font-size:11px; color:#666; word-break:break-all;", vals$db_path))
    } else {
      div(class="status-err", icon("exclamation-circle"), " No DB Loaded")
    }
  })
  
  # --- OPTIMIZED QUERY BUILDER (CTE Method) ---
  # This uses Common Table Expressions to force DuckDB to filter the Observations table
  # BEFORE performing the expensive join.
  get_sql_query <- function(count_only = FALSE, limit = NULL) {
    req(vals$db_loaded)
    req(input$filter_disease != "None")
    
    # 1. Base Filtered CTE for Observations (The Fact Table)
    #    Since this is partitioned by Disease, this is nearly instant.
    obs_where <- "1=1"
    if (input$filter_disease != "All") obs_where <- paste0(obs_where, sprintf(" AND Disease = '%s'", input$filter_disease))
    if (input$filter_isotype != "All") obs_where <- paste0(obs_where, sprintf(" AND Isotype = '%s'", input$filter_isotype))
    if (input$filter_bsource != "All") obs_where <- paste0(obs_where, sprintf(" AND BSource = '%s'", input$filter_bsource))
    if (input$filter_btype   != "All") obs_where <- paste0(obs_where, sprintf(" AND BType   = '%s'", input$filter_btype))
    
    # 2. Base Filtered CTE for Metrics (The Dimension Table)
    #    Filtering this table first reduces the size of the Join hash table.
    met_where <- "1=1"
    if (input$min_patients > 1)   met_where <- paste0(met_where, sprintf(" AND N_Patients >= %d", input$min_patients))
    if (input$min_diseases > 1)   met_where <- paste0(met_where, sprintf(" AND N_Diseases >= %d", input$min_diseases))
    if (input$min_antibodies > 1) met_where <- paste0(met_where, sprintf(" AND N_Antibodies >= %d", input$min_antibodies))
    
    # 3. Construct the CTE Query
    #    We select only necessary columns to reduce memory overhead.
    columns <- if (count_only) "COUNT(*) as n" else "
      o.Peptide,
      o.Antibody,
      o.Disease,
      o.Patient,
      o.Isotype,
      o.CDR3_spanning_count,
      m.N_Patients,
      m.N_Diseases,
      m.N_Antibodies
    "
    
    sql <- sprintf("
      WITH filtered_obs AS (
        SELECT * FROM observations WHERE %s
      ),
      filtered_met AS (
        SELECT * FROM metrics WHERE %s
      )
      SELECT %s
      FROM filtered_obs o
      JOIN filtered_met m ON o.Peptide = m.Peptide
    ", obs_where, met_where, columns)
    
    if (!is.null(limit) && !count_only) sql <- paste0(sql, " LIMIT ", limit)
    
    return(sql)
  }
  
  # --- Output: Summary Stats ---
  output$summary_stats <- renderUI({
    if (!vals$db_loaded) return(div("Please load a database."))
    if (input$filter_disease == "None") return(div("Select a Disease to begin."))
    
    withProgress(message = 'Calculating stats...', value = 0.5, {
      sql <- get_sql_query(count_only = TRUE)
      n_rows <- tryCatch(dbGetQuery(con, sql)$n, error = function(e) 0)
    })
    
    tagList(
      div(class="summary-stat", strong("Rows Matched:"), fmt_int(n_rows)),
      div(class="summary-stat", strong("Disease Filter:"), input$filter_disease)
    )
  })
  
  # --- Output: Preview Table ---
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
      downloadButton("dl_csv", "Download CSV"),
      downloadButton("dl_fasta", "Download FASTA")
    )
  })
  
  # --- Downloads ---
  output$dl_csv <- downloadHandler(
    filename = function() { paste0("OAS_Export_", Sys.Date(), ".csv") },
    content = function(file) {
      showNotification("Preparing CSV... (This may take a while)", duration = 5)
      sql <- get_sql_query()
      tmp_sql <- sprintf("COPY (%s) TO '%s' (FORMAT CSV, HEADER)", sql, norm_path(file))
      dbExecute(con, tmp_sql)
    }
  )
  
  output$dl_fasta <- downloadHandler(
    filename = function() { paste0("OAS_Export_", Sys.Date(), ".fasta") },
    content = function(file) {
      showNotification("Preparing FASTA...", duration = 5)
      sql <- get_sql_query()
      df <- dbGetQuery(con, sql) # Must load to R for Biostrings
      if (nrow(df) > 0) {
        headers <- paste(df$Antibody, df$Disease, df$Patient, sep="|")
        seqs <- AAStringSet(df$Peptide)
        names(seqs) <- headers
        writeXStringSet(seqs, file)
      } else {
        writeLines("No data found.", file)
      }
    }
  )
}

shinyApp(ui, server)