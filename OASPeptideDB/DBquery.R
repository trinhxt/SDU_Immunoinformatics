# ==============================================================================
# DAT-DB: Disease-exclusive Antibody Tryptic Peptide Database Viewer
# ==============================================================================
# This Shiny application acts as a GUI for exploring antibody peptide data.
# Key Features:
# 1. Connects to a local DuckDB instance to query Parquet files efficiently.
# 2. Provides dynamic filtering for Diseases, B-Cell types, and Specificity logic.
# 3. Generates live summary statistics and data previews.
# 4. Exports filtered data to Parquet, FASTA, and Summary Text reports.
# ==============================================================================

suppressPackageStartupMessages({
  # Core Shiny libraries for UI and Server logic
  library(shiny)
  library(shinythemes)   # For the "cerulean" visual theme
  library(shinyjs)       # For JavaScript operations (hide/show elements)
  library(shinyWidgets)  # For alerts and enhanced UI inputs
  library(shinyBS)       # For Bootstrap Tooltips/Popovers (Info Icons)
  
  # Data Handling & Database libraries
  library(DBI)           # Standard Database Interface
  library(duckdb)        # High-performance in-process SQL database (OLAP)
  library(DT)            # For interactive DataTables
  library(data.table)    # For fast data manipulation
  
  # System & Biology libraries
  library(shinyFiles)    # For OS-native file/folder selection dialogs
  library(Biostrings)    # For handling biological sequences (writing FASTA)
  library(svDialogs)     # For system save dialogs
  library(tools)         # For file path manipulation (file_path_sans_ext)
})

# --- Helper Functions ---

# Normalizes file paths to use forward slashes (Windows compatibility)
norm_path <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

# Formats integers with commas (e.g., 1,000,000)
fmt_int <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE)

# Custom NULL coalescing operator: if 'a' is not NULL return 'a', else return 'b'
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# 1. DATA DICTIONARY
# ==============================================================================
# Maps short internal disease codes (keys) to readable full names (values).
disease_code_to_full <- c(
  "AChR-MG"                         = "Acetylcholine receptor-positive myasthenia gravis",
  "Allergic-Rhinitis-In-Season"       = "Allergic rhinitis in season",
  "Allergic-Rhinitis-Out-Of-Season" = "Allergic rhinitis out of season",
  "Allergy%2FNoSIT"                    = "Allergy NoSIT",
  "Allergy/NoSIT"                      = "Allergy NoSIT",
  "Allergy%2FSIT"                      = "Allergy SIT",
  "Allergy/SIT"                        = "Allergy SIT",
  "Asthma"                             = "Asthma",
  "CLL"                                = "Chronic lymphocytic leukemia",
  "CMV"                                = "Cytomegalovirus",
  "CMV%2FEBV"                         = "Cytomegalovirus/Epstein–Barr virus",
  "CMV/EBV"                            = "Cytomegalovirus/Epstein–Barr virus",
  "Dengue"                             = "Dengue",
  "Ebola"                              = "Ebola",
  "EBV"                                = "Epstein–Barr virus",
  "HCV"                                = "Hepatitis C virus",
  "Healthy%2Fceliac-disease"           = "Celiac",
  "Healthy/celiac-disease"            = "Celiac",
  "HIV"                                = "HIV",
  "Light-Chain-Amyloidosis"            = "Light chain amyloidosis",
  "MS"                                 = "Multiple sclerosis",
  "MuSK-MG"                            = "Muscle-specific kinase myasthenia gravis",
  "Non-Dengue-Febrile-Illness"        = "Non-dengue febrile illness",
  "Obstructive-Sleep-Apnea"            = "Obstructive sleep apnea",
  "SARS-COV-2"                         = "SARS-COV-2",
  "SLE"                                = "Systemic lupus erythematosus",
  "Tonsillitis"                       = "Tonsillitis",
  "Tonsillitis%2FObstructive-Sleep-Apnea" = "Tonsillitis/Obstructive sleep apnea",
  "Tonsillitis/Obstructive-Sleep-Apnea"   = "Tonsillitis/Obstructive sleep apnea"
)

# Helper to look up disease names safely, handling URL encoding
get_full_disease_name <- function(code) {
  val <- if (length(code) > 0) code[1] else "None"
  if (is.null(val) || val %in% c("None", "All")) return(val)
  
  # Direct match
  if (val %in% names(disease_code_to_full)) return(disease_code_to_full[[val]])
  
  # Try decoding URL characters (e.g., %2F -> /)
  decoded <- tryCatch(utils::URLdecode(val), error = function(e) val)
  if (decoded %in% names(disease_code_to_full)) return(disease_code_to_full[[decoded]])
  
  return(decoded) # Return original if no match found
}

# ==============================================================================
# 2. USER INTERFACE (UI)
# ==============================================================================
ui <- fluidPage(
  useShinyjs(), # Initialize JavaScript capabilities
  theme = shinytheme("cerulean"), # Set global CSS theme
  
  # --- Custom CSS Styles ---
  tags$head(tags$style(HTML("
    /* General Panel Styling */
    .panel-container { 
      border: 1px solid #ddd; 
      padding: 15px; 
      border-radius: 5px; 
      background: #fcfcfc; 
      margin-bottom: 15px; 
      box-shadow: 0 1px 3px rgba(0,0,0,0.1); 
    }
    
    /* Section Headers */
    .custom-header { 
      font-size: 16px; 
      font-weight: bold; 
      border-left: 4px solid #007bc2; 
      padding-left: 10px; 
      margin-bottom: 15px; 
      color: #333; 
    }
    
    /* Summary Stats Formatting */
    .summary-stat { 
      font-size: 13px; 
      margin-bottom: 4px; 
      color: #444; 
      display: flex; 
      justify-content: space-between; 
      border-bottom: 1px dotted #eee; 
    }
    .summary-stat strong { 
      color: #222; 
      margin-right: 10px; 
      text-align: right; 
      max-width: 60%; 
    }
    
    /* Buttons */
    .btn-download { width: auto; margin-right: 10px; margin-bottom: 5px; }
    
    /* Status Messages */
    .status-ok { color: green; font-weight: bold; }
    .status-err { color: red; font-weight: bold; }
    
    /* Tooltip Icon Style */
    .info-icon { 
      color: #007bc2; 
      font-size: 0.9em; 
      margin-left: 5px; 
      cursor: help; 
    }
    
    /* [NEW] FIXED TOOLTIP WIDTH (Aggressive Override) */
    body .tooltip-inner {
      background-color: #ffffff !important;   /* White background */
      color: #333333 !important;              /* Dark text */
      border: 1px solid #cccccc;              /* Border */
      box-shadow: 0 4px 15px rgba(0,0,0,0.15); /* Shadow */
      
      /* WIDTH SETTINGS */
      min-width: 250px !important;            
      max-width: 600px !important;            
      
      font-size: 13px;
      padding: 12px;
      text-align: left;
    }
    
    body .tooltip.right .tooltip-arrow {
      border-right-color: #cccccc !important; 
    }
    body .tooltip {
      opacity: 1 !important; 
    }
  "))),
  
  # --- App Title ---
  titlePanel(HTML("<strong>DAT-DB</strong>: <strong>D</strong>isease-exclusive <strong>A</strong>ntibody <strong>T</strong>rypic peptide <strong>D</strong>atabase for <strong>B</strong>ottom-up proteomics")),
  
  sidebarLayout(
    # --- Sidebar Panel (Inputs) ---
    sidebarPanel(width = 3,
                 # Section 1: Database Connection
                 div(class = "panel-container",
                     div(class = "custom-header", "1. Database Source"),
                     uiOutput("db_status_ui"), # Shows green check or red error
                     hidden(actionButton("db_folder_btn", "Browse Folder", icon = icon("folder-open"))) # Hidden if DB loads automatically
                 ),
                 
                 # Section 2: Primary Filters (Metadata)
                 div(class = "panel-container",
                     div(class = "custom-header", "2. Clinical Filters"),
                     selectInput("filter_disease", "Disease", choices = c("None"), selected = "None"),
                     selectInput("filter_bsource", "B-Cell Source", choices = "All", selected = "All"),
                     selectInput("filter_btype", "B-Cell Type", choices = "All", selected = "All"),
                     selectInput("filter_isotype", "Isotype", choices = "All", selected = "All")
                 ),
                 
                 # Section 3: Logic Filters (Calculated fields) with Info Tooltips
                 div(class = "panel-container",
                     div(class = "custom-header", "3. Specificity & Quality"),
                     
                     # --- CDR3 Spanning with Tooltip [UPDATED] ---
                     radioButtons("logic_cdr3", 
                                  label = tagList(
                                    "CDR3 Spanning (Coverage %)",
                                    tags$span(id = "tip_cdr3", icon("info-circle"), class = "info-icon"),
                                    bsTooltip("tip_cdr3", "Percentage of the CDR3 region covered by the peptide overlap. Low (<10%), Medium (10-30%), High (>30%).", placement = "right", trigger = "hover")
                                  ),
                                  choices = c("All"="all", "Low (<10%)"="low", "Medium (10-30%)"="medium", "High (>30%)"="high"), selected = "high"),
                     
                     # --- Disease Specificity with Tooltip ---
                     radioButtons("logic_specificity", 
                                  label = tagList(
                                    "Disease Specificity",
                                    tags$span(id = "tip_spec", icon("info-circle"), class = "info-icon"),
                                    bsTooltip("tip_spec", "Disease-Exclusive: Peptide found ONLY in this disease. Shared: Peptide found in this disease AND others.", placement = "right", trigger = "hover")
                                  ), 
                                  choices = c("All"="all","Disease-Exclusive"="specific","Shared"="shared"), selected = "specific"),
                     
                     # --- Peptide Uniqueness with Tooltip ---
                     radioButtons("logic_uniqueness", 
                                  label = tagList(
                                    "Peptide Uniqueness",
                                    tags$span(id = "tip_uniq", icon("info-circle"), class = "info-icon"),
                                    bsTooltip("tip_uniq", "Unique: Peptide maps to exactly 1 antibody sequence. Common: Peptide maps to >1 antibody (conserved region).", placement = "right", trigger = "hover")
                                  ), 
                                  choices = c("All"="all","Unique"="unique","Common"="common"), selected = "all"),
                     
                     # --- Patient Frequency with Tooltip ---
                     radioButtons("logic_patients", 
                                  label = tagList(
                                    "Patient Frequency",
                                    tags$span(id = "tip_pat", icon("info-circle"), class = "info-icon"),
                                    bsTooltip("tip_pat", "Singleton: Found in only 1 patient. Recurring: Found in 2+ patients (potential public clonotype).", placement = "right", trigger = "hover")
                                  ), 
                                  choices = c("All"="all","Singleton"="single","Recurring"="multi"), selected = "multi")
                 )
    ),
    
    # --- Main Panel (Outputs) ---
    mainPanel(width = 9,
              tabsetPanel(
                # Tab 1: Live Data Preview
                tabPanel("Data Preview", br(),
                         div(class = "panel-container",
                             div(class = "custom-header", "4. Live Summary"),
                             uiOutput("summary_stats"), # Statistics block
                             
                             # Download buttons (Hidden until a disease is selected)
                             hidden(div(id = "download_section", hr(),
                                        div(class = "custom-header", "5. Download Results"),
                                        actionButton("btn_save_parquet", "Save Parquet", icon = icon("file-archive"), class="btn-download"),
                                        actionButton("btn_save_fasta", "Save FASTA", icon = icon("file-code"), class="btn-download"),
                                        actionButton("btn_save_summary", "Save Summary", icon = icon("file-text"), class="btn-download")
                             ))
                         ),
                         div(class = "panel-container", DTOutput("preview_table"))
                ),
                
                # Tab 2: About/Documentation
                tabPanel("About", br(),
                         div(class = "panel-container", h4("Welcome to DAT-DB app!"), 
                             p("The identification of human antibodies using bottom-up proteomics relies on database searches that match experimental peptide fragments to theoretical values derived from protein sequences in databases. Standard protein databases, such as UniProt and NCBI-RefSeq, contain a limited number of antibody sequences compared to the human body’s capacity to generate billions and currently lack disease-exclusive antibody sequences. This limitation could lead to the misidentification of antibodies in samples. To overcome this challenge, DAT-DB, a database of disease-exclusive antibodies, provides researchers with antibody tryptic peptides derived from next-generation sequencing of antibody repertoires.")),
                         div(class = "panel-container", h4("Database Status"), textOutput("debug_path"))
                )
              )
    )
  )
)

# ==============================================================================
# 3. SERVER LOGIC
# ==============================================================================
server <- function(input, output, session) {
  
  # Display current working directory for debugging path issues
  output$debug_path <- renderText({ paste("Current Working Directory:", getwd()) })
  
  # --- Database Initialization ---
  # Create an in-memory DuckDB connection
  con <- dbConnect(duckdb())
  
  # Configure DuckDB for performance
  dbExecute(con, "PRAGMA memory_limit='4GB'; PRAGMA threads=4")
  
  # Ensure cleanup when the user closes the app
  session$onSessionEnded(function() { dbDisconnect(con, shutdown = TRUE) })
  
  # Reactive values to track database state
  vals <- reactiveValues(db_loaded = FALSE, db_path = NULL)
  
  # Show/Hide download section based on if a disease is selected
  observe({
    if (vals$db_loaded && length(input$filter_disease) > 0 && input$filter_disease[1] != "None") {
      shinyjs::show("download_section")
    } else {
      shinyjs::hide("download_section")
    }
  })
  
  # --- Function: Load Database ---
  load_database <- function(path) {
    withProgress(message = 'Initializing Database', value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Resetting views...")
        dbExecute(con, "DROP VIEW IF EXISTS data")
        
        incProgress(0.2, detail = "Scanning file system...")
        clean_path <- gsub("\\\\", "/", path) # Fix Windows slashes
        glob <- file.path(clean_path, "**", "*.parquet") # Recursive glob pattern
        
        # Create a View directly on the files (Zero-copy load)
        incProgress(0.2, detail = "Indexing Parquet files...")
        dbExecute(con, sprintf("CREATE VIEW data AS SELECT * FROM read_parquet('%s')", glob))
        
        # --- Populate Metadata Dropdowns ---
        incProgress(0.2, detail = "Extracting Disease metadata...")
        dis <- dbGetQuery(con, "SELECT DISTINCT Disease FROM data WHERE Disease IS NOT NULL ORDER BY Disease")
        if (nrow(dis) == 0) stop("No data found.")
        updateSelectInput(session, "filter_disease", choices = c("None", dis$Disease), selected = "None")
        
        incProgress(0.2, detail = "Extracting Isotype/Source metadata...")
        iso <- dbGetQuery(con, "SELECT DISTINCT Isotype FROM data WHERE Isotype IS NOT NULL ORDER BY Isotype")
        updateSelectInput(session, "filter_isotype", choices = c("All", iso$Isotype), selected = "All")
        
        bs <- dbGetQuery(con, "SELECT DISTINCT BSource FROM data WHERE BSource IS NOT NULL ORDER BY BSource")
        updateSelectInput(session, "filter_bsource", choices = c("All", bs$BSource), selected = "All")
        
        bt <- dbGetQuery(con, "SELECT DISTINCT BType FROM data WHERE BType IS NOT NULL ORDER BY BType")
        updateSelectInput(session, "filter_btype", choices = c("All", bt$BType), selected = "All")
        
        incProgress(0.1, detail = "Finalizing...")
        vals$db_loaded <- TRUE
        vals$db_path <- norm_path(path)
        shinyjs::hide("db_folder_btn") # Hide manual browse button on success
        sendSweetAlert(session, "Success", paste0("Loaded:\n", vals$db_path), "success")
        
      }, error = function(e) {
        sendSweetAlert(session, "Error ", "No parquet files found in this folder.\n Please select a proper database folder.", "error")
        vals$db_loaded <- FALSE; shinyjs::show("db_folder_btn")
      })
    })
  }
  
  # --- Auto-Detection Logic ---
  observe({
    candidates <- c("db_stage3", "OASpepDB", "OASpeptideDB")
    found <- NULL
    for (p in candidates) if (dir.exists(p) && length(list.files(p, "\\.parquet$", recursive=T)) > 0) { found <- p; break }
    if (!is.null(found)) load_database(found) else shinyjs::show("db_folder_btn")
  })
  
  # Manual Browse Button Handler
  observeEvent(input$db_folder_btn, {
    res <- dlg_dir(default = getwd(), title = "Select DB Folder")$res
    if (length(res) == 1 && nzchar(res)) load_database(norm_path(res))
  })
  
  # Render Database Status Indicator
  output$db_status_ui <- renderUI({
    if (vals$db_loaded) tagList(div(class="status-ok", icon("check"), " Database Loaded"), div(style="font-size:11px; color:#666;", vals$db_path))
    else div(class="status-err", icon("exclamation"), " No Database Found")
  })
  
  # --- Query Builder Function ---
  get_sql_query <- function(count_only = FALSE, limit = NULL, custom_select = NULL) {
    req(vals$db_loaded, length(input$filter_disease) > 0, input$filter_disease[1] != "None")
    
    val <- function(x) if (length(x) > 0) x[1] else "All"
    
    # Determine columns to select
    sel <- if (!is.null(custom_select)) custom_select 
    else if (count_only) "COUNT(*) as n, COUNT(DISTINCT Peptide) as n_pep, COUNT(DISTINCT Patient) as n_pat, COUNT(DISTINCT Antibody) as n_ab, MIN(N_Patients) as min_p, MAX(N_Patients) as max_p" 
    else "Peptide, Disease, Patient, Isotype, CDR3_spanning_pct, CDR3_spanning_count, N_Patients, N_Diseases, N_Antibodies, cdr3_aa, Antibody, filename"
    
    q <- sprintf("SELECT %s FROM data WHERE 1=1", sel)
    
    # --- [NEW] UNDERGROUND FILTER ---
    # Strictly exclude any peptides found in the healthy background.
    q <- paste0(q, " AND is_healthy_background = FALSE")
    
    # Apply Clinical Filters
    f_dis <- val(input$filter_disease)
    if (f_dis != "All") q <- paste0(q, sprintf(" AND Disease = '%s'", f_dis))
    f_src <- val(input$filter_bsource)
    if (f_src != "All") q <- paste0(q, sprintf(" AND BSource = '%s'", f_src))
    f_typ <- val(input$filter_btype)
    if (f_typ != "All") q <- paste0(q, sprintf(" AND BType = '%s'", f_typ))
    f_iso <- val(input$filter_isotype)
    if (f_iso != "All") q <- paste0(q, sprintf(" AND Isotype = '%s'", f_iso))
    
    # Apply Logic Filters (Updated for Percentages)
    l_cdr3 <- val(input$logic_cdr3)
    if (l_cdr3 == "low")         q <- paste0(q, " AND CDR3_spanning_pct <= 0.1")
    else if (l_cdr3 == "medium") q <- paste0(q, " AND CDR3_spanning_pct > 0.1 AND CDR3_spanning_pct <= 0.3")
    else if (l_cdr3 == "high")   q <- paste0(q, " AND CDR3_spanning_pct > 0.3")
    # 'all' implies no extra filter on CDR3
    
    if (val(input$logic_specificity) == "specific") q <- paste0(q, " AND N_Diseases = 1")
    if (val(input$logic_specificity) == "shared")   q <- paste0(q, " AND N_Diseases > 1")
    if (val(input$logic_uniqueness) == "unique")    q <- paste0(q, " AND N_Antibodies = 1")
    if (val(input$logic_uniqueness) == "common")    q <- paste0(q, " AND N_Antibodies > 1")
    if (val(input$logic_patients) == "single")      q <- paste0(q, " AND N_Patients = 1")
    if (val(input$logic_patients) == "multi")       q <- paste0(q, " AND N_Patients > 1")
    
    # Apply Limit
    if (!is.null(limit) && !count_only && is.null(custom_select)) q <- paste0(q, " LIMIT ", limit)
    return(q)
  }
  
  # --- Summary Stats Output ---
  output$summary_stats <- renderUI({
    if (!vals$db_loaded) return(div("Please load a database."))
    if (length(input$filter_disease) == 0 || input$filter_disease[1] == "None") return(div("Select a Disease."))
    
    withProgress(message = 'Calculating Summary', value = 0, {
      # Get Basic Counts
      incProgress(0.2, detail = "Counting records...")
      res <- tryCatch(dbGetQuery(con, get_sql_query(count_only = TRUE)), error = function(e) data.frame(n=0, n_pep=0))
      
      # Get Metadata Summaries
      incProgress(0.5, detail = "Analyzing metadata categories...")
      res_meta <- tryCatch(dbGetQuery(con, get_sql_query(custom_select = "DISTINCT BSource, BType, Isotype")), error = function(e) data.frame())
      
      uniq_s <- if(nrow(res_meta)>0) paste(unique(res_meta$BSource), collapse=", ") else "-"
      uniq_t <- if(nrow(res_meta)>0) paste(unique(res_meta$BType), collapse=", ") else "-"
      uniq_i <- if(nrow(res_meta)>0) paste(unique(res_meta$Isotype), collapse=", ") else "-"
      incProgress(0.3, detail = "Rendering...")
    })
    
    # UI Label Map (Updated)
    map_dict <- list("all"="All", 
                     "low"="Low (<10%)", "medium"="Medium (10-30%)", "high"="High (>30%)", 
                     "specific"="Disease-Exclusive", "shared"="Shared", 
                     "unique"="Unique", "common"="Common", 
                     "single"="Singleton", "multi"="Recurring")
    lab <- function(x) map_dict[[x[1]]] %||% x[1]
    
    # Render HTML Block
    tagList(
      h4(strong("Distinct Peptides: "), fmt_int(res$n_pep), style="color: #007bc2;"),
      div(style="font-size:11px; margin-bottom:8px;", paste("Total Rows (Incidences):", fmt_int(res$n))),
      hr(style="margin:5px 0;"),
      div(class="summary-stat", span("Disease:"), strong(get_full_disease_name(input$filter_disease))),
      div(class="summary-stat", span("Sources:"), strong(uniq_s)),
      div(class="summary-stat", span("Types:"), strong(uniq_t)),
      div(class="summary-stat", span("Isotypes:"), strong(uniq_i)),
      div(class="summary-stat", span("CDR3 Coverage:"), strong(lab(input$logic_cdr3))),
      div(class="summary-stat", span("Specificity:"), strong(lab(input$logic_specificity))),
      div(class="summary-stat", span("Uniqueness:"), strong(lab(input$logic_uniqueness))),
      div(class="summary-stat", span("Patient Freq:"), strong(lab(input$logic_patients))),
      div(class="summary-stat", span("Healthy Excl:"), strong("Yes (Active)"))
    )
  })
  
  # --- Data Table Preview ---
  output$preview_table <- renderDT({
    if (!vals$db_loaded || length(input$filter_disease)==0 || input$filter_disease[1] == "None") return(NULL)
    withProgress(message='Previewing Table', value=0, { 
      incProgress(0.5, detail = "Fetching top 100 rows...")
      dbGetQuery(con, get_sql_query(limit=100)) 
    }) %>% datatable(options=list(scrollX=T, pageLength=10))
  })
  
  # ============================================================================
  # 4. DOWNLOAD HANDLERS
  # ============================================================================
  
  # --- Parquet Export Handler ---
  observeEvent(input$btn_save_parquet, {
    f <- dlg_save(default = paste0("OAS_", Sys.Date(), ".parquet"), title = "Save Parquet")$res
    if (length(f) == 1 && nzchar(f)) {
      tryCatch({
        withProgress(message="Exporting Parquet", value=0, {
          incProgress(0.2, detail = "Generating SQL...")
          sql <- get_sql_query(); 
          incProgress(0.4, detail = "Writing to disk (DuckDB COPY)...")
          dbExecute(con, sprintf("COPY (%s) TO '%s' (FORMAT PARQUET)", sql, norm_path(f)))
          incProgress(0.4, detail = "Done.")
        })
        sendSweetAlert(session, "Export Complete!", paste("Saved:", f), "success")
      }, error = function(e) sendSweetAlert(session, "Failed", e$message, "error"))
    }
  })
  
  # --- FASTA Export Handler ---
  observeEvent(input$btn_save_fasta, {
    f <- dlg_save(default = paste0("OAS_", Sys.Date(), ".fasta"), title = "Save FASTA")$res
    if (length(f) == 1 && nzchar(f)) {
      tryCatch({
        withProgress(message="Exporting FASTA", value=0, {
          incProgress(0.2, detail = "Querying peptide data...")
          df <- dbGetQuery(con, get_sql_query(custom_select = "DISTINCT Peptide"))
          if(nrow(df) == 0) stop("No data matches current filters.")
          incProgress(0.4, detail = "Formatting sequences...")
          seqs <- AAStringSet(df$Peptide)
          names(seqs) <- paste0("OAS|pep_", seq_along(df$Peptide))
          incProgress(0.3, detail = "Writing file...")
          writeXStringSet(seqs, f)
          incProgress(0.1, detail = "Done.")
        })
        sendSweetAlert(session, "Export Complete!", paste("FASTA saved to:", f), "success")
      }, error = function(e) sendSweetAlert(session, "Failed", e$message, "error"))
    }
  })
  
  # --- Summary Report Export Handler ---
  observeEvent(input$btn_save_summary, {
    f <- dlg_save(default = paste0("OAS_Summary_", Sys.Date(), ".txt"), title = "Save Summary")$res
    if (length(f) == 1 && nzchar(f)) {
      tryCatch({
        withProgress(message="Generating Summary Report", value=0, {
          incProgress(0.1, detail = "Counting totals...")
          cnt <- dbGetQuery(con, get_sql_query(count_only = T))
          tot <- cnt$n_pep
          
          stat <- function(col) {
            incProgress(0.1, detail = paste("Aggregating", col, "..."))
            df <- dbGetQuery(con, paste0(get_sql_query(custom_select=sprintf("%s, COUNT(DISTINCT Peptide) c", col)), sprintf(" GROUP BY %s ORDER BY c DESC", col)))
            if(nrow(df)==0) return("-")
            paste(sprintf("%s: %s (%.1f%%)", df[[col]], fmt_int(df$c), (df$c/tot)*100), collapse=", ")
          }
          
          map_dict <- list("all"="All", "low"="Low (<10%)", "medium"="Medium (10-30%)", "high"="High (>30%)", "specific"="Disease-Exclusive", "shared"="Shared", "unique"="Unique", "common"="Common", "single"="Singleton", "multi"="Recurring")
          
          incProgress(0.2, detail = "Formatting report...")
          lines <- c(
            "==================================================",
            "DAT-DB EXPORT SUMMARY",
            "==================================================",
            paste("Export Date:   ", Sys.time()),
            paste("Database Path: ", vals$db_path),
            "",
            "CORE METRICS:",
            "--------------------------------------------------",
            paste("Distinct Peptides: ", fmt_int(cnt$n_pep)),
            paste("Unique Patients:   ", fmt_int(cnt$n_pat)),
            paste("Unique Antibodies: ", fmt_int(cnt$n_ab)),
            paste("Total Rows/Hits:   ", fmt_int(cnt$n)),
            "",
            "ACTIVE FILTERS:",
            "--------------------------------------------------",
            paste("Disease:             ", get_full_disease_name(input$filter_disease)),
            paste("B-Cell Source:       ", input$filter_bsource[1]),
            paste("B-Cell Type:         ", input$filter_btype[1]),
            paste("Isotype:             ", input$filter_isotype[1]),
            paste("CDR3 Coverage:       ", map_dict[[input$logic_cdr3[1]]] %||% input$logic_cdr3[1]),
            paste("Specificity:         ", map_dict[[input$logic_specificity[1]]] %||% input$logic_specificity[1]),
            paste("Peptide Uniqueness:  ", map_dict[[input$logic_uniqueness[1]]] %||% input$logic_uniqueness[1]),
            paste("Patient Frequency:   ", map_dict[[input$logic_patients[1]]] %||% input$logic_patients[1]),
            "Healthy Exclusion:    YES (Enforced)",
            "",
            "DATA BREAKDOWN (by Peptide count):",
            "--------------------------------------------------",
            paste("Sources:             ", stat("BSource")),
            paste("Types:               ", stat("BType")),
            paste("Isotypes:            ", stat("Isotype")),
            paste("Patient Freq Range:  ", sprintf("Min: %s, Max: %s patients per peptide", cnt$min_p, cnt$max_p)),
            "=================================================="
          )
          
          incProgress(0.1, detail = "Writing to disk...")
          writeLines(lines, f)
        })
        sendSweetAlert(session, "Export Complete!", paste("Saved:", f), "success")
      }, error = function(e) sendSweetAlert(session, "Failed", e$message, "error"))
    }
  })
}

# Run the Application
shinyApp(ui, server)