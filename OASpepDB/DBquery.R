# ==============================================================================
# OASpepDB-v3: Interactive Disease-Exclusive Antibody Database Explorer
# Framework: R Shiny + DuckDB + Bslib + DT + Arrow
# Purpose: Sub-second Hive partition queries, Group-specific FDR FASTA export,
#          Parquet reverse-lookup, and multi-omics provenance tracking.
# ==============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(duckdb)
  library(DBI)
  library(arrow)
  library(ggplot2)
  library(plotly)
  library(dplyr)
  library(htmlwidgets)
  library(zip)
})

# ------------------------------------------------------------------------------
# Configuration & Constants
# ------------------------------------------------------------------------------
# Resolve application directory dynamically for cross-platform portability
APP_DIR <- tryCatch({
  normalizePath(dirname(sys.frame(1)$ofile))
}, error = function(e) {
  if (dir.exists("OASpepDB")) file.path(getwd(), "OASpepDB") else getwd()
})

resolve_data_path <- function(filename, pattern = NULL) {
  search_dirs <- c(
    file.path(APP_DIR, "Data"),
    file.path(getwd(), "OASpepDB", "Data"),
    file.path(getwd(), "Data"),
    file.path(dirname(APP_DIR), "Data")
  )
  for (d in search_dirs) {
    if (dir.exists(d)) {
      if (!is.null(pattern)) {
        matches <- list.files(d, pattern = pattern, full.names = TRUE)
        if (length(matches) > 0) return(matches[1])
      }
      p <- file.path(d, filename)
      if (file.exists(p)) return(p)
    }
  }
  return(file.path(search_dirs[1], filename))
}

CRAP_REF_PATH       <- resolve_data_path("cRAP.fasta")
UNIPROT_REF_PATH    <- resolve_data_path("Uniprot_SP_canonical_2026_09_14.fasta", pattern = "^Uni[pP]rot_SP_canonical.*\\.fasta$")
ENTRAPMENT_REF_PATH <- resolve_data_path("entrapment_cassettes.fasta")

detect_default_db_root <- function() {
  candidates <- c(
    "D:/OAS/unpaired/CDR3_db",
    file.path(getwd(), "CDR3_db"),
    file.path(APP_DIR, "CDR3_db"),
    file.path(dirname(APP_DIR), "CDR3_db"),
    "C:/OAS/unpaired/CDR3_db"
  )
  for (cand in candidates) {
    if (dir.exists(cand)) {
      dis_dirs <- list.dirs(cand, full.names = FALSE, recursive = FALSE)
      if (any(grepl("^Disease=", dis_dirs))) {
        return(normalizePath(cand, winslash = "/"))
      }
    }
  }
  return("")
}

get_db_cohort_count <- function(path) {
  if (!nzchar(path) || !dir.exists(path)) return(0L)
  dis_dirs <- list.dirs(path, full.names = FALSE, recursive = FALSE)
  cohorts <- grep("^Disease=", dis_dirs, value = TRUE)
  length(cohorts)
}

DEFAULT_DB_ROOT <- detect_default_db_root()
DB_ROOT <- DEFAULT_DB_ROOT

# Standardized Disease Dictionary (25 Cohorts)
DISEASE_MAP <- c(
  "COVID-19"              = "Severe Acute Respiratory Syndrome Coronavirus 2 (COVID-19)",
  "HIV"                   = "Human Immunodeficiency Virus Infection",
  "CMV-EBV"               = "Cytomegalovirus & Epstein-Barr Virus Co-infection",
  "Allergy-NoSIT"         = "Allergy without Specific Immunotherapy",
  "SLE"                   = "Systemic Lupus Erythematosus",
  "Ebola"                 = "Ebola Virus Disease",
  "Tonsillitis"           = "Chronic / Recurrent Tonsillitis",
  "Celiac-Disease"        = "Celiac Disease (Gluten-Sensitive Enteropathy)",
  "MS"                    = "Multiple Sclerosis",
  "OSA"                   = "Obstructive Sleep Apnea",
  "AChR-MG"               = "Acetylcholine Receptor Myasthenia Gravis",
  "AL-Amyloidosis"        = "Immunoglobulin Light Chain Amyloidosis",
  "Allergic-Rhinitis-In"  = "Allergic Rhinitis (Pollen Season In)",
  "Allergic-Rhinitis-Out" = "Allergic Rhinitis (Pollen Season Out)",
  "Allergy-SIT"           = "Allergy with Specific Immunotherapy",
  "Asthma"                = "Bronchial Asthma",
  "CLL"                   = "Chronic Lymphocytic Leukemia",
  "CMV"                   = "Cytomegalovirus Infection",
  "Dengue"                = "Dengue Virus Infection",
  "EBV"                   = "Epstein-Barr Virus Infection",
  "HCV"                   = "Hepatitis C Virus Infection",
  "MuSK-MG"               = "Muscle-Specific Kinase Myasthenia Gravis",
  "NDFI"                  = "Non-Dengue Acute Febrile Illness",
  "POEMS"                 = "POEMS Syndrome",
  "Tonsillitis-OSA"       = "Tonsillitis with Obstructive Sleep Apnea"
)

# Build UI selector labels (Short names as requested)
DISEASE_CHOICES <- sort(names(DISEASE_MAP))

# Categorized Clinical Cohorts for Landing Exploration
DISEASE_CATEGORIES <- list(
  "Infectious Diseases" = list(
    icon = "virus",
    color = "danger",
    diseases = c("COVID-19", "HIV", "HCV", "Dengue", "EBV", "CMV", "CMV-EBV", "Ebola", "NDFI")
  ),
  "Autoimmune & Neuro" = list(
    icon = "shield-halved",
    color = "primary",
    diseases = c("SLE", "Celiac-Disease", "MS", "AChR-MG", "MuSK-MG")
  ),
  "Hematology & Tumors" = list(
    icon = "droplet",
    color = "warning",
    diseases = c("CLL", "AL-Amyloidosis", "POEMS")
  ),
  "Allergy & Airways" = list(
    icon = "lungs",
    color = "success",
    diseases = c("Asthma", "Allergic-Rhinitis-In", "Allergic-Rhinitis-Out", "Allergy-SIT", "Allergy-NoSIT", "Tonsillitis", "Tonsillitis-OSA", "OSA")
  )
)

# Pre-aggregated 25 Cohort Scale Summary for Instant Bar Plot Rendering (0.00s latency)
COHORT_SUMMARY_DATA <- data.frame(
  Disease = c("COVID-19", "HIV", "CMV-EBV", "Allergy-NoSIT", "SLE", "Ebola", 
              "Tonsillitis", "Celiac-Disease", "MS", "OSA", "Tonsillitis-OSA", 
              "Asthma", "AL-Amyloidosis", "MuSK-MG", "EBV", "HCV", "AChR-MG", 
              "POEMS", "CMV", "Allergy-SIT", "Dengue", "Allergic-Rhinitis-In", 
              "CLL", "Allergic-Rhinitis-Out", "NDFI"),
  Category = c("Infectious Diseases", "Infectious Diseases", "Infectious Diseases", 
               "Allergy & Airways", "Autoimmune & Neuro", "Infectious Diseases", 
               "Allergy & Airways", "Autoimmune & Neuro", "Autoimmune & Neuro", 
               "Allergy & Airways", "Allergy & Airways", "Allergy & Airways", 
               "Hematology & Tumors", "Autoimmune & Neuro", "Infectious Diseases", 
               "Infectious Diseases", "Autoimmune & Neuro", "Hematology & Tumors", 
               "Infectious Diseases", "Allergy & Airways", "Infectious Diseases", 
               "Allergy & Airways", "Hematology & Tumors", "Allergy & Airways", 
               "Infectious Diseases"),
  Cassettes = c(21660594, 17600440, 4902218, 4029662, 2542357, 2487848, 
                2395769, 2320911, 2109543, 1347586, 1277457, 521877, 
                455763, 365541, 340276, 279566, 263982, 228038, 174440, 
                142322, 74481, 6611, 4490, 3903, 1877),
  Reads = c(94723857, 136352117, 20842914, 22466536, 7822296, 9144698, 
            13741829, 5213223, 17735211, 5649056, 7799262, 6872312, 
            2846474, 973822, 1199475, 654791, 805834, 865093, 186094, 
            277580, 113679, 10351, 11458, 5436, 2958),
  stringsAsFactors = FALSE
) %>% arrange(desc(Cassettes))

COHORT_SUMMARY_DATA$Code <- paste0("D", 1:nrow(COHORT_SUMMARY_DATA))

# Bi-directional mapping: Code (D1..D25) <-> Disease Name
CODE_TO_DISEASE <- setNames(COHORT_SUMMARY_DATA$Disease, COHORT_SUMMARY_DATA$Code)
DISEASE_TO_CODE <- setNames(COHORT_SUMMARY_DATA$Code, COHORT_SUMMARY_DATA$Disease)

# Helper to slugify names for consistent CSS selector mapping
slugify <- function(x) {
  s <- gsub("[^A-Za-z0-9_-]", "-", x)
  gsub("-+", "-", s)
}

# Helper to construct Liquid Glass responsive cohort capsule pills
build_glass_cohort_pills <- function(category_name, cat_slug = "default") {
  sub_df <- COHORT_SUMMARY_DATA %>% 
    filter(Category == category_name) %>% 
    arrange(desc(Cassettes))
  
  if (nrow(sub_df) == 0) return(NULL)
  
  tagList(
    lapply(1:nrow(sub_df), function(i) {
      row <- sub_df[i, ]
      d <- row$Disease
      cas <- row$Cassettes
      cas_str <- if (cas >= 1e6) sprintf("%.1fM", cas / 1e6) else if (cas >= 1e3) sprintf("%.0fK", cas / 1e3) else as.character(cas)
      
      tags$button(
        type = "button",
        class = paste0("cohort-pill-glass pill-cat-", cat_slug),
        `data-disease` = d,
        title = sprintf("%s (%s CDR3 peptides)", DISEASE_MAP[d], format(cas, big.mark = ",")),
        onclick = "Shiny.setInputValue('btn_cohort_click', this.getAttribute('data-disease'), {priority: 'event'});",
        div(class = "d-flex align-items-center text-truncate",
          tags$span(class = "cohort-color-dot"),
          tags$span(class = "pill-label text-truncate", d)
        ),
        tags$span(class = "pill-count ms-2 flex-shrink-0", cas_str)
      )
    })
  )
}

# Backward compatibility alias
build_cohort_buttons <- function(disease_vec, color_class) {
  sorted_diseases <- intersect(COHORT_SUMMARY_DATA$Disease, disease_vec)
  tagList(
    lapply(sorted_diseases, function(d) {
      cas <- COHORT_SUMMARY_DATA$Cassettes[COHORT_SUMMARY_DATA$Disease == d]
      cas_str <- if (cas >= 1e6) sprintf("%.1fM", cas / 1e6) else if (cas >= 1e3) sprintf("%.0fK", cas / 1e3) else as.character(cas)
      tags$button(
        type = "button",
        class = paste0("cohort-pill-glass pill-cat-", color_class),
        `data-disease` = d,
        title = sprintf("%s (%s CDR3 peptides)", DISEASE_MAP[d], format(cas, big.mark = ",")),
        onclick = "Shiny.setInputValue('btn_cohort_click', this.getAttribute('data-disease'), {priority: 'event'});",
        div(class = "d-flex align-items-center text-truncate",
          tags$span(class = "cohort-color-dot"),
          tags$span(class = "pill-label text-truncate", d)
        ),
        tags$span(class = "pill-count ms-2 flex-shrink-0", cas_str)
      )
    })
  )
}


# ------------------------------------------------------------------------------
# DuckDB Connection Helper
# ------------------------------------------------------------------------------
get_db_con <- function() {
  con <- dbConnect(duckdb::duckdb(dbdir = ":memory:"))
  n_cores <- max(1L, parallel::detectCores() - 1L)
  dbExecute(con, sprintf("PRAGMA threads=%d;", min(n_cores, 16L)))
  dbExecute(con, "PRAGMA memory_limit='16GB';")
  return(con)
}

# ------------------------------------------------------------------------------
# UI Definition
# ------------------------------------------------------------------------------
ui <- page_navbar(
  title = "Human Disease CDR3 Antibody Peptides Database",
  id = "main_nav",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#1b4965",
    secondary = "#62b6cb",
    success = "#2a9d8f",
    info = "#457b9d",
    base_font = font_google("Inter")
  ),
  header = tags$head(
    tags$style(HTML("
      /* OASpepDB Liquid Glass Design System */
      body {
        background-color: #f7f9fc !important;
        background-image: 
          radial-gradient(circle at 12% 18%, rgba(224, 238, 255, 0.5) 0%, transparent 45%),
          radial-gradient(circle at 88% 82%, rgba(220, 248, 240, 0.45) 0%, transparent 50%) !important;
        background-attachment: fixed !important;
        font-family: 'Inter', system-ui, -system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
        color: #1d1d1f;
        -webkit-font-smoothing: antialiased;
        -moz-osx-font-smoothing: grayscale;
      }
      
      /* Center Navbar Tabs: Query & Export, Reverse Lookup */
      .navbar {
        position: relative !important;
      }
      .navbar > .container-fluid {
        position: relative !important;
      }
      .navbar-brand {
        font-size: 1.05rem !important;
        font-weight: 700 !important;
        letter-spacing: -0.01em !important;
      }
      @media (min-width: 992px) {
        .navbar-collapse {
          display: flex !important;
          justify-content: center !important;
        }
        .navbar-nav#main_nav,
        .navbar .navbar-nav {
          position: absolute !important;
          left: 50% !important;
          top: 50% !important;
          transform: translate(-50%, -50%) !important;
          display: flex !important;
          flex-direction: row !important;
          align-items: center !important;
          justify-content: center !important;
          margin: 0 !important;
          padding: 0 !important;
          gap: 8px !important;
        }
      }
      @media (max-width: 991.98px) {
        .navbar-collapse {
          display: flex !important;
          justify-content: center !important;
          text-align: center !important;
        }
        .navbar-nav#main_nav,
        .navbar .navbar-nav {
          margin: 0 auto !important;
          display: flex !important;
          flex-direction: row !important;
          justify-content: center !important;
          align-items: center !important;
          gap: 8px !important;
        }
      }
      .navbar-nav .nav-link,
      .navbar-nav > li > a {
        font-size: 0.88rem !important;
        font-weight: 600 !important;
        letter-spacing: -0.01em !important;
        padding: 6px 14px !important;
        border-radius: 8px !important;
        transition: all 0.18s ease !important;
      }
      
      /* Liquid Glass Cards */
      .glass-card {
        background: rgba(255, 255, 255, 0.82) !important;
        backdrop-filter: blur(28px) saturate(200%) !important;
        -webkit-backdrop-filter: blur(28px) saturate(200%) !important;
        border: 1px solid rgba(255, 255, 255, 0.85) !important;
        border-radius: 16px !important;
        box-shadow: 0 8px 32px 0 rgba(31, 38, 135, 0.07), inset 0 1px 1.5px rgba(255, 255, 255, 0.9) !important;
        transition: box-shadow 0.25s ease, border-color 0.25s ease, transform 0.25s ease;
        position: relative;
        overflow: visible !important;
      }
      .glass-card > .card-header {
        position: relative !important;
        z-index: 1050 !important;
        overflow: visible !important;
      }
      .glass-card > .card-body {
        position: relative !important;
        z-index: 1 !important;
      }
      .glass-card:hover {
        box-shadow: 0 12px 36px 0 rgba(31, 38, 135, 0.10), inset 0 1px 1.5px rgba(255, 255, 255, 0.95) !important;
      }
      
      /* Liquid Glass Hero Canvas */
      .glass-hero-panel {
        background: rgba(255, 255, 255, 0.85);
        backdrop-filter: blur(28px) saturate(200%);
        -webkit-backdrop-filter: blur(28px) saturate(200%);
        border: 1px solid rgba(255, 255, 255, 0.9);
        border-radius: 16px;
        box-shadow: 0 8px 32px 0 rgba(31, 38, 135, 0.08), inset 0 1px 1.5px rgba(255, 255, 255, 0.95);
        padding: 16px 24px;
        position: relative;
        overflow: hidden;
      }
      .glass-hero-panel::before {
        content: '';
        position: absolute;
        top: 0; left: 0; right: 0;
        height: 3px;
        background: linear-gradient(90deg, #1b4965 0%, #0071e3 50%, #2a9d8f 100%);
      }
      .glass-hero-title {
        font-size: 1.25rem;
        font-weight: 700;
        letter-spacing: -0.025em;
        color: #1d1d1f;
      }
      .glass-hero-badge {
        font-size: 0.70rem;
        font-weight: 600;
        padding: 2px 8px;
        border-radius: 9999px;
        background: rgba(0, 113, 227, 0.08);
        color: #0071e3;
        border: 1px solid rgba(0, 113, 227, 0.18);
        font-family: ui-monospace, SFMono-Regular, monospace;
      }
      .glass-hero-sub {
        font-size: 0.83rem;
        color: #6e6e73;
        letter-spacing: -0.01em;
        margin-top: 3px;
      }
      
      /* Precision Liquid Glass KPI Chips */
      .glass-kpi-chip {
        display: inline-flex;
        align-items: center;
        gap: 5px;
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 0.76rem;
        font-weight: 600;
        letter-spacing: -0.01em;
        border: 1px solid rgba(255, 255, 255, 0.75);
        background: rgba(255, 255, 255, 0.65);
        backdrop-filter: blur(14px);
        -webkit-backdrop-filter: blur(14px);
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.03), inset 0 1px 1px rgba(255, 255, 255, 0.9);
        color: #1d1d1f;
        transition: transform 0.15s ease, background 0.15s ease;
        white-space: nowrap;
      }
      .glass-kpi-chip:hover {
        transform: translateY(-1px);
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.06), inset 0 1px 1px rgba(255, 255, 255, 0.95);
      }
      .glass-kpi-chip .kpi-num {
        font-family: ui-monospace, SFMono-Regular, monospace;
        font-weight: 700;
      }
      .glass-kpi-blue {
        background: rgba(0, 113, 227, 0.08);
        color: #0071e3;
        border-color: rgba(0, 113, 227, 0.20);
      }
      .glass-kpi-green {
        background: rgba(42, 157, 143, 0.10);
        color: #1b7a6e;
        border-color: rgba(42, 157, 143, 0.22);
      }
      .glass-kpi-amber {
        background: rgba(245, 158, 11, 0.10);
        color: #b45309;
        border-color: rgba(245, 158, 11, 0.24);
      }
      .glass-kpi-neutral {
        background: rgba(255, 255, 255, 0.75);
        color: #333336;
        border-color: rgba(255, 255, 255, 0.85);
      }
      
      /* Interactive Cohort Capsule Pills (Full-Width 1-Column Row) */
      .cohort-pill-glass {
        display: flex;
        justify-content: space-between;
        align-items: center;
        width: 100%;
        background: rgba(255, 255, 255, 0.75);
        backdrop-filter: blur(16px) saturate(180%);
        -webkit-backdrop-filter: blur(16px) saturate(180%);
        border: 1px solid rgba(255, 255, 255, 0.82);
        border-radius: 9px;
        padding: 6px 10px;
        margin: 0 0 3px 0;
        font-size: 0.76rem;
        font-weight: 500;
        color: #1d1d1f;
        cursor: pointer;
        user-select: none;
        transition: transform 0.14s cubic-bezier(0.16, 1, 0.3, 1), box-shadow 0.14s ease, border-color 0.14s ease, background-color 0.14s ease;
        box-shadow: 0 2px 6px rgba(0, 0, 0, 0.03), inset 0 1px 1px rgba(255, 255, 255, 0.9);
        text-decoration: none;
      }
      .cohort-pill-glass:hover {
        transform: translateX(2.5px);
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.07), inset 0 1px 1px rgba(255, 255, 255, 0.95);
        background-color: rgba(255, 255, 255, 0.92);
      }
      .cohort-pill-glass:active {
        transform: scale(0.99);
        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.04);
      }
      .cohort-pill-glass .cohort-color-dot {
        width: 7px;
        height: 7px;
        border-radius: 50%;
        display: inline-block;
        margin-right: 8px;
        flex-shrink: 0;
      }
      .cohort-pill-glass .pill-label {
        font-weight: 600;
        color: #1d1d1f;
        letter-spacing: -0.01em;
        font-size: 0.77rem;
      }
      .cohort-pill-glass .pill-count {
        font-size: 0.70rem;
        color: #6e6e73;
        font-weight: 600;
        font-family: ui-monospace, SFMono-Regular, monospace;
        background: rgba(245, 245, 247, 0.7);
        border: 1px solid rgba(0, 0, 0, 0.06);
        padding: 2px 7px;
        border-radius: 6px;
      }
      .pill-cat-infectious .cohort-color-dot { background-color: #D55E00; }
      .pill-cat-infectious:hover { border-color: rgba(213, 94, 0, 0.4); background-color: rgba(255, 255, 255, 0.95); }
      .pill-cat-allergy .cohort-color-dot { background-color: #009E73; }
      .pill-cat-allergy:hover { border-color: rgba(0, 158, 115, 0.4); background-color: rgba(255, 255, 255, 0.95); }
      .pill-cat-autoimmune .cohort-color-dot { background-color: #0072B2; }
      .pill-cat-autoimmune:hover { border-color: rgba(0, 114, 178, 0.4); background-color: rgba(255, 255, 255, 0.95); }
      .pill-cat-hematology .cohort-color-dot { background-color: #CC79A7; }
      .pill-cat-hematology:hover { border-color: rgba(204, 121, 167, 0.4); background-color: rgba(255, 255, 255, 0.95); }
      
      /* Liquid Glass Domain Header Badges */
      .glass-domain-badge {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        padding: 3px 9px;
        border-radius: 12px;
        font-size: 0.72rem;
        font-weight: 600;
        background: rgba(255, 255, 255, 0.7);
        backdrop-filter: blur(10px);
        -webkit-backdrop-filter: blur(10px);
        color: #333336;
        border: 1px solid rgba(255, 255, 255, 0.85);
        box-shadow: 0 1px 4px rgba(0, 0, 0, 0.03);
      }
      .glass-domain-badge .domain-dot {
        width: 7px;
        height: 7px;
        border-radius: 50%;
        background-color: var(--domain-color);
      }
      .glass-domain-group {
        border-bottom: 1px solid rgba(0, 0, 0, 0.05);
        padding-bottom: 8px;
      }
      .glass-domain-group:last-child {
        border-bottom: none;
      }
      .glass-count-meta {
        font-size: 0.72rem;
        color: #86868b;
        font-family: ui-monospace, SFMono-Regular, monospace;
      }
      .glass-guide-callout {
        background: rgba(255, 255, 255, 0.7);
        backdrop-filter: blur(14px);
        -webkit-backdrop-filter: blur(14px);
        border: 1px solid rgba(255, 255, 255, 0.85);
        border-radius: 12px;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.02);
      }
      .glass-interactive-hint {
        color: #0071e3;
        font-weight: 500;
        display: inline-flex;
        align-items: center;
      }
      
      /* Subtle Liquid Glass Scrollbar */
      .glass-scrollbar::-webkit-scrollbar {
        width: 5px;
      }
      .glass-scrollbar::-webkit-scrollbar-track {
        background: transparent;
      }
      .glass-scrollbar::-webkit-scrollbar-thumb {
        background: rgba(0, 0, 0, 0.12);
        border-radius: 10px;
      }
      .glass-scrollbar::-webkit-scrollbar-thumb:hover {
        background: rgba(0, 0, 0, 0.22);
      }
      
      /* Liquid Glass Sleek Filter Toolbar (Gray Capsule) */
      .glass-toolbar {
        display: inline-flex;
        align-items: center;
        background: #eaecf0 !important;
        backdrop-filter: blur(28px) saturate(200%);
        -webkit-backdrop-filter: blur(28px) saturate(200%);
        border: 1px solid #d5d8de !important;
        border-radius: 10px;
        padding: 3px 6px;
        gap: 6px;
        box-shadow: inset 0 1px 1.5px rgba(0, 0, 0, 0.04), 0 1px 3px rgba(0, 0, 0, 0.04) !important;
      }
      .glass-segment {
        display: inline-flex;
        align-items: center;
        gap: 4px;
      }
      .glass-label {
        font-size: 0.70rem;
        font-weight: 600;
        color: #6e6e73;
        letter-spacing: -0.01em;
        white-space: nowrap;
      }
      .glass-toolbar .shiny-input-container {
        margin-bottom: 0 !important;
        display: inline-block !important;
      }
      .glass-toolbar .shiny-input-container > label {
        display: none !important;
      }
      .glass-toolbar select.form-control,
      .glass-toolbar select.form-select,
      .glass-toolbar select.shiny-input-select {
        height: 25px !important;
        min-height: 25px !important;
        line-height: 23px !important;
        padding: 0 18px 0 7px !important;
        font-size: 0.72rem !important;
        font-weight: 500 !important;
        color: #1d1d1f !important;
        background-color: rgba(255, 255, 255, 0.85) !important;
        backdrop-filter: blur(12px) !important;
        -webkit-backdrop-filter: blur(12px) !important;
        border: 1px solid rgba(0, 0, 0, 0.12) !important;
        border-radius: 6px !important;
        cursor: pointer !important;
        background-image: url('data:image/svg+xml;utf8,<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"8\" height=\"5\" viewBox=\"0 0 8 5\"><path fill=\"none\" stroke=\"%2386868b\" stroke-width=\"1.5\" stroke-linecap=\"round\" stroke-linejoin=\"round\" d=\"M1 1l3 3 3-3\"/></svg>') !important;
        background-repeat: no-repeat !important;
        background-position: right 5px center !important;
        box-shadow: 0 1px 2px rgba(0, 0, 0, 0.04), inset 0 1px 1px rgba(255, 255, 255, 0.9) !important;
        transition: all 0.15s ease !important;
      }
      .glass-toolbar select.form-control:hover,
      .glass-toolbar select.shiny-input-select:hover {
        border-color: rgba(0, 113, 227, 0.4) !important;
        box-shadow: 0 1px 4px rgba(0, 0, 0, 0.08) !important;
      }
      .glass-toolbar select.form-control:focus,
      .glass-toolbar select.shiny-input-select:focus {
        outline: none !important;
        border-color: #0071e3 !important;
        box-shadow: 0 0 0 2px rgba(0, 113, 227, 0.20) !important;
      }
      .glass-divider {
        width: 1px;
        height: 16px;
        background-color: rgba(0, 0, 0, 0.12);
        margin: 0 2px;
      }
      .glass-action-btn {
        height: 25px !important;
        line-height: 23px !important;
        font-size: 0.72rem !important;
        font-weight: 600 !important;
        color: #ffffff !important;
        background: linear-gradient(180deg, #1b4965 0%, #153a51 100%) !important;
        border: 1px solid rgba(255, 255, 255, 0.35) !important;
        border-radius: 6px 0 0 6px !important;
        padding: 0 9px !important;
        display: inline-flex !important;
        align-items: center !important;
        gap: 5px !important;
        box-shadow: 0 2px 6px rgba(27, 73, 101, 0.25), inset 0 1px 1px rgba(255, 255, 255, 0.3) !important;
        text-decoration: none !important;
        white-space: nowrap !important;
      }
      .glass-action-btn:hover {
        background: linear-gradient(180deg, #225c80 0%, #1b4965 100%) !important;
        box-shadow: 0 3px 8px rgba(27, 73, 101, 0.35), inset 0 1px 1px rgba(255, 255, 255, 0.4) !important;
        color: #ffffff !important;
      }
      .glass-action-btn:active {
        background: #112e40 !important;
        box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.2) !important;
        color: #ffffff !important;
      }
      .glass-action-toggle {
        border-radius: 0 6px 6px 0 !important;
        padding: 0 5px !important;
        border-left: 1px solid rgba(255, 255, 255, 0.25) !important;
      }
      .glass-reset-btn {
        width: 26px;
        height: 26px;
        padding: 0;
        border-radius: 6px;
        border: 1px solid rgba(255, 255, 255, 0.85) !important;
        background: rgba(255, 255, 255, 0.85) !important;
        backdrop-filter: blur(14px) !important;
        -webkit-backdrop-filter: blur(14px) !important;
        color: #555559 !important;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        transition: all 0.2s cubic-bezier(0.16, 1, 0.3, 1);
        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05), inset 0 1px 1px rgba(255, 255, 255, 0.95);
      }
      .glass-reset-btn:hover {
        background: rgba(255, 255, 255, 0.95) !important;
        color: #1d1d1f !important;
        border-color: rgba(27, 73, 101, 0.3) !important;
        box-shadow: 0 2px 6px rgba(0, 0, 0, 0.08), inset 0 1px 1px rgba(255, 255, 255, 1);
        transform: rotate(-35deg);
      }
      .glass-reset-btn:active {
        background: #e5e5ea !important;
        transform: rotate(-180deg);
      }
      
      /* Liquid Glass Static Tooltips */
      .glass-tooltip-wrap {
        position: relative;
        display: inline-flex;
        align-items: center;
        gap: 3px;
      }
      .glass-tooltip-wrap:hover {
        z-index: 1060;
      }
      .glass-tooltip-wrap:has(.dropdown-menu.show) .glass-tooltip-box,
      .glass-tooltip-wrap:has(.dropdown-toggle.show) .glass-tooltip-box {
        display: none !important;
        visibility: hidden !important;
        opacity: 0 !important;
      }
      .glass-tooltip-wrap .glass-tooltip-box {
        visibility: hidden;
        opacity: 0;
        position: absolute;
        top: calc(100% + 8px);
        left: 50%;
        transform: translateX(-50%) translateY(-4px);
        background: rgba(255, 255, 255, 0.94) !important;
        backdrop-filter: blur(28px) saturate(200%) !important;
        -webkit-backdrop-filter: blur(28px) saturate(200%) !important;
        color: #1d1d1f !important;
        border: 1px solid rgba(255, 255, 255, 0.92) !important;
        box-shadow: 0 16px 40px -4px rgba(31, 38, 135, 0.18), 0 4px 12px rgba(0, 0, 0, 0.06), inset 0 1px 1.5px rgba(255, 255, 255, 0.95) !important;
        border-radius: 11px !important;
        padding: 9px 13px !important;
        font-size: 0.72rem !important;
        font-weight: 500 !important;
        line-height: 1.4 !important;
        white-space: normal !important;
        width: 250px;
        z-index: 999999 !important;
        pointer-events: none;
        transition: opacity 0.2s cubic-bezier(0.16, 1, 0.3, 1), transform 0.2s cubic-bezier(0.16, 1, 0.3, 1);
      }
      .glass-tooltip-wrap .glass-tooltip-box::before {
        content: '';
        position: absolute;
        bottom: 100%;
        left: 50%;
        transform: translateX(-50%);
        border-width: 6px;
        border-style: solid;
        border-color: transparent transparent rgba(255, 255, 255, 0.94) transparent;
      }
      .glass-tooltip-wrap:hover .glass-tooltip-box {
        visibility: visible;
        opacity: 1;
        transform: translateX(-50%) translateY(0);
      }
      .glass-tooltip-box.align-start {
        left: 0 !important;
        transform: translateY(-4px) !important;
      }
      .glass-tooltip-wrap:hover .glass-tooltip-box.align-start {
        transform: translateY(0) !important;
      }
      .glass-tooltip-box.align-start::before {
        left: 14px !important;
        transform: none !important;
      }
      .glass-tooltip-box.align-end {
        left: auto !important;
        right: 0 !important;
        transform: translateY(-4px) !important;
      }
      .glass-tooltip-wrap:hover .glass-tooltip-box.align-end {
        transform: translateY(0) !important;
      }
      .glass-tooltip-box.align-end::before {
        left: auto !important;
        right: 14px !important;
        transform: none !important;
      }
      .glass-tooltip-title {
        font-weight: 700;
        color: #1b4965;
        margin-bottom: 3px;
        display: flex;
        align-items: center;
        gap: 5px;
      }
      .glass-info-icon {
        font-size: 0.72rem;
        color: #86868b;
        transition: color 0.15s ease, transform 0.15s ease;
        cursor: help;
        margin-left: 2px;
      }
      .glass-tooltip-wrap:hover .glass-info-icon {
        color: #1b4965;
        transform: scale(1.15);
      }
      
      /* Global Bootstrap Tooltip Override to Liquid Glass Style */
      .tooltip {
        --bs-tooltip-bg: rgba(255, 255, 255, 0.88) !important;
        --bs-tooltip-color: #1d1d1f !important;
        --bs-tooltip-padding-x: 10px !important;
        --bs-tooltip-padding-y: 7px !important;
        --bs-tooltip-border-radius: 9px !important;
        --bs-tooltip-font-size: 0.72rem !important;
        z-index: 99999 !important;
      }
      .tooltip .tooltip-inner {
        background-color: rgba(255, 255, 255, 0.88) !important;
        backdrop-filter: blur(28px) saturate(200%) !important;
        -webkit-backdrop-filter: blur(28px) saturate(200%) !important;
        color: #1d1d1f !important;
        border: 1px solid rgba(255, 255, 255, 0.9) !important;
        box-shadow: 0 12px 36px -4px rgba(31, 38, 135, 0.15), inset 0 1px 1.5px rgba(255, 255, 255, 0.95) !important;
        font-weight: 500 !important;
        line-height: 1.4 !important;
        text-align: left !important;
        max-width: 250px !important;
      }
      .tooltip .tooltip-arrow::before {
        border-top-color: rgba(255, 255, 255, 0.88) !important;
        border-bottom-color: rgba(255, 255, 255, 0.88) !important;
      }

      /* Floating Interactive Liquid Glass Tooltip for Alluvial Diagram */
      #alluvial-glass-tooltip {
        position: fixed;
        z-index: 999999;
        pointer-events: none;
        opacity: 0;
        visibility: hidden;
        transition: opacity 0.12s cubic-bezier(0.16, 1, 0.3, 1), transform 0.08s ease-out;
        transform: translateY(4px);
        background: rgba(255, 255, 255, 0.88) !important;
        backdrop-filter: blur(28px) saturate(200%) !important;
        -webkit-backdrop-filter: blur(28px) saturate(200%) !important;
        border: 1px solid rgba(255, 255, 255, 0.92) !important;
        box-shadow: 0 16px 40px -4px rgba(31, 38, 135, 0.16), 0 4px 12px rgba(0, 0, 0, 0.06), inset 0 1px 1.5px rgba(255, 255, 255, 0.95) !important;
        border-radius: 12px;
        padding: 9px 13px;
        max-width: 320px;
        min-width: 170px;
        color: #1d1d1f;
        font-family: 'Inter', system-ui, -system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
      }
      #alluvial-glass-tooltip.active {
        opacity: 1;
        visibility: visible;
        transform: translateY(0);
      }
      #alluvial-glass-tooltip .tooltip-header {
        font-size: 0.78rem;
        font-weight: 700;
        color: #1b4965;
        margin-bottom: 3px;
        display: flex;
        align-items: center;
        gap: 6px;
      }
      #alluvial-glass-tooltip .tooltip-meta {
        font-size: 0.72rem;
        font-weight: 600;
        color: #2b2b2f;
        margin-bottom: 2px;
        font-family: ui-monospace, SFMono-Regular, monospace;
      }
      #alluvial-glass-tooltip .tooltip-desc {
        font-size: 0.69rem;
        color: #6e6e73;
        line-height: 1.35;
        margin-top: 3px;
        padding-top: 3px;
        border-top: 1px solid rgba(0, 0, 0, 0.06);
      }
    ")),
    tags$script(HTML("
      $(document).on('shiny:connected', function() {
        Shiny.addCustomMessageHandler('switch_tab', function(tabName) {
          const el = document.querySelector('a[data-value=\"' + tabName + '\"]');
          if (el) {
            if (window.bootstrap && bootstrap.Tab) {
              bootstrap.Tab.getOrCreateInstance(el).show();
            } else {
              el.click();
            }
          }
        });
      });

      // Liquid Glass Interactive Cursor-Tracking Tooltip for Alluvial SVG
      (function() {
        let tt = null;
        function getTooltip() {
          if (!tt) {
            tt = document.getElementById('alluvial-glass-tooltip');
            if (!tt) {
              tt = document.createElement('div');
              tt.id = 'alluvial-glass-tooltip';
              document.body.appendChild(tt);
            }
          }
          return tt;
        }

        function findGlassTarget(el) {
          while (el && el !== document) {
            if (el.getAttribute && el.getAttribute('data-glass-title')) {
              return el;
            }
            el = el.parentNode;
          }
          return null;
        }

        document.addEventListener('mouseover', function(e) {
          const target = findGlassTarget(e.target);
          if (target) {
            const el = getTooltip();
            const title = target.getAttribute('data-glass-title') || '';
            const meta = target.getAttribute('data-glass-meta') || '';
            const desc = target.getAttribute('data-glass-desc') || '';

            let html = '<div class=\"tooltip-header\">' + title + '</div>';
            if (meta) {
              html += '<div class=\"tooltip-meta\">' + meta + '</div>';
            }
            if (desc) {
              html += '<div class=\"tooltip-desc\">' + desc + '</div>';
            }
            el.innerHTML = html;
            el.classList.add('active');
          }
        });

        document.addEventListener('mousemove', function(e) {
          const el = getTooltip();
          if (el.classList.contains('active')) {
            const offset = 14;
            let left = e.clientX + offset;
            let top = e.clientY + offset;

            const rect = el.getBoundingClientRect();
            if (left + rect.width > window.innerWidth - 12) {
              left = e.clientX - rect.width - offset;
            }
            if (top + rect.height > window.innerHeight - 12) {
              top = e.clientY - rect.height - offset;
            }
            if (left < 10) left = 10;
            if (top < 10) top = 10;

            el.style.left = left + 'px';
            el.style.top = top + 'px';
          }
        });

        document.addEventListener('mouseout', function(e) {
          const target = findGlassTarget(e.target);
          if (target) {
            const related = findGlassTarget(e.relatedTarget);
            if (related !== target) {
              const el = getTooltip();
              el.classList.remove('active');
            }
          }
        });
      })();
    "))
  ),
  
  # Tab 1: Query & Export (Standalone Landing Page with Database Summary Flow)
  nav_panel(
    title = "Query & Export",
    value = "Query & Export",
    icon = icon("magnifying-glass-chart"),
    
    # Unified Master-Detail Alluvial Flow Canvas
    div(
      class = "row g-3 my-2 px-2",
      div(
        class = "col-12",
        card(
          class = "glass-card shadow-sm border",
          card_header(
            class = "bg-transparent border-bottom py-2 px-3 d-flex align-items-center justify-content-between flex-wrap gap-2",
            style = "position: relative; z-index: 1050; overflow: visible !important;",
            div(
              class = "d-flex align-items-center flex-nowrap",
              div(
                class = "glass-tooltip-wrap me-2",
                actionButton(
                  "btn_reset_landing",
                  label = NULL,
                  icon = icon("rotate-left"),
                  class = "btn glass-reset-btn"
                ),
                div(
                  class = "glass-tooltip-box align-start",
                  style = "width: 230px;",
                  div(class = "glass-tooltip-title", icon("rotate-left", class = "me-1"), "Reset Overview"),
                  "Unselect disease cohort, restore default filters, and free system memory."
                )
              ),
              tags$strong("Database overview", class = "text-dark small fw-semibold text-nowrap"),
              div(
                class = "glass-tooltip-wrap ms-2",
                actionButton(
                  "btn_open_db_modal_header",
                  label = tagList(icon("database", class = "me-1"), uiOutput("txt_db_status_badge", inline = TRUE)),
                  class = "btn glass-action-btn btn-sm",
                  style = "height: 28px !important; line-height: 26px !important; font-size: 0.74rem !important; border-radius: 8px !important; padding: 0 10px !important; display: inline-flex; align-items: center;"
                ),
                div(
                  class = "glass-tooltip-box align-start",
                  style = "width: 260px;",
                  div(class = "glass-tooltip-title", icon("database", class = "me-1"), "Connect / Switch Database"),
                  "Select or update the root folder of the Hive CDR3 Parquet database."
                )
              )
            ),
            div(
              class = "d-flex align-items-center justify-content-center mx-auto",
              style = "position: relative; z-index: 1050;",
              uiOutput("intro_center_toolbar")
            ),
            div(
              class = "d-flex align-items-center justify-content-end",
              style = "position: relative; z-index: 1050;",
              uiOutput("intro_export_controls")
            )
          ),
          card_body(
            class = "p-2 d-flex flex-column align-items-center justify-content-center position-relative",
            style = "min-height: 590px; overflow: hidden; position: relative; z-index: 1;",
            uiOutput("landing_unified_alluvial", style = "width: 100%; height: 100%;")
          )
        )
      )
    )
  ),
  
  # Tab 2: Reverse Lookup (Liquid Glass Design)
  nav_panel(
    title = "Reverse Lookup",
    value = "Reverse Lookup",
    icon = icon("magnifying-glass-location"),
    
    div(
      class = "row g-3 my-2 px-2",
      div(
        class = "col-12",
        card(
          class = "glass-card shadow-sm border",
          card_header(
            class = "bg-transparent border-bottom py-2 px-3 d-flex align-items-center justify-content-between flex-wrap gap-2",
            style = "position: relative; z-index: 1050;",
            div(
              class = "d-flex align-items-center gap-2",
              div(
                class = "glass-reset-btn me-1",
                style = "width: 28px; height: 28px; pointer-events: none;",
                icon("magnifying-glass-location", class = "text-primary")
              ),
              div(
                tags$strong("Tryptic Peptide Reverse Lookup", class = "text-dark small fw-bold d-block"),
                tags$span("Trace mass spectrometry peptides back to donor patients and full-length antibody clonotypes", class = "text-muted", style = "font-size: 0.72rem;")
              )
            ),
            div(
              class = "d-flex align-items-center gap-2",
              div(
                class = "glass-tooltip-wrap",
                actionButton(
                  "btn_open_db_modal_reverse",
                  label = tagList(icon("database", class = "me-1"), uiOutput("txt_db_status_badge_rev", inline = TRUE)),
                  class = "btn glass-action-btn btn-sm",
                  style = "height: 28px !important; line-height: 26px !important; font-size: 0.74rem !important; border-radius: 8px !important; padding: 0 10px !important; display: inline-flex; align-items: center;"
                ),
                div(
                  class = "glass-tooltip-box align-end",
                  style = "width: 260px;",
                  div(class = "glass-tooltip-title", icon("database", class = "me-1"), "Connect / Switch Database"),
                  "Select or update the root folder of the Hive CDR3 Parquet database."
                )
              ),
              div(
                class = "glass-tooltip-wrap",
                span(
                  class = "badge",
                  style = "background: rgba(0, 113, 227, 0.08); color: #0071e3; border: 1px solid rgba(0, 113, 227, 0.20); font-size: 0.72rem; font-weight: 600; padding: 4px 10px; border-radius: 20px; font-family: ui-monospace, SFMono-Regular, monospace; cursor: help;",
                  icon("dna", class = "me-1"), "Exact Fragment Match"
                ),
                div(
                  class = "glass-tooltip-box align-end",
                  style = "width: 260px;",
                  div(class = "glass-tooltip-title", icon("dna", class = "me-1"), "Tryptic Digestion Logic"),
                  "Performs high-speed exact matching against pre-computed in silico tryptic peptides (cleavage C-terminal to Lys/Arg, excluding Pro). Powered by DuckDB list_contains()."
                )
              )
            )
          ),
          card_body(
            class = "p-3",
            
            # Query Input Panel in Liquid Glass
            div(
              class = "glass-hero-panel p-3 mb-3",
              style = "border-radius: 12px; border: 1px solid rgba(255, 255, 255, 0.9);",
              div(
                class = "row g-3 align-items-center",
                div(
                  class = "col-lg-8 col-md-7 col-12",
                  tags$label(
                    class = "form-label fw-semibold text-dark d-flex align-items-center gap-1 mb-1",
                    style = "font-size: 0.78rem;",
                    icon("keyboard", class = "text-primary"),
                    "Tryptic Peptide Fragment(s)",
                    tags$span(class = "text-muted fw-normal ms-1", "(One per line, comma, semicolon, or space-separated):")
                  ),
                  tags$textarea(
                    id = "txt_reverse_pep",
                    class = "form-control",
                    rows = 3,
                    style = "background: rgba(255, 255, 255, 0.85); backdrop-filter: blur(12px); -webkit-backdrop-filter: blur(12px); border: 1px solid rgba(0, 0, 0, 0.12); border-radius: 8px; font-family: ui-monospace, SFMono-Regular, monospace; font-size: 0.82rem; font-weight: 600; letter-spacing: 0.5px; box-shadow: inset 0 1px 2px rgba(0, 0, 0, 0.04); resize: vertical;",
                    placeholder = "e.g. ARVRHLGRGMDV or GPSVFPLAPSSK",
                    "ARVRHLGRGMDV"
                  )
                ),
                div(
                  class = "col-lg-4 col-md-5 col-12 d-flex flex-column justify-content-between h-100",
                  div(
                    class = "glass-guide-callout p-2 mb-2",
                    style = "border-radius: 8px; font-size: 0.72rem; line-height: 1.35; color: #555559;",
                    div(class = "fw-semibold text-dark mb-1 d-flex align-items-center gap-1", icon("circle-info", class = "text-primary"), "MS/MS Identification"),
                    "Matches identified peptides from FragPipe / MaxQuant directly back to circulating OAS clonotypes, patient counts, and NGS files."
                  ),
                  div(
                    class = "d-flex gap-2",
                    actionButton(
                      "btn_run_lookup",
                      label = tagList(icon("magnifying-glass", class = "me-1"), "Search Provenance"),
                      class = "btn glass-action-btn flex-grow-1",
                      style = "height: 32px !important; line-height: 30px !important; font-size: 0.78rem !important; border-radius: 8px !important; justify-content: center !important;"
                    ),
                    actionButton(
                      "btn_clear_lookup",
                      label = tagList(icon("trash-can")),
                      class = "btn glass-reset-btn",
                      style = "width: 32px; height: 32px; border-radius: 8px; flex-shrink: 0;",
                      title = "Clear query input"
                    )
                  )
                )
              )
            ),
            
            # Results Table Container in Liquid Glass
            div(
              class = "glass-card p-3",
              style = "border-radius: 12px; background: rgba(255, 255, 255, 0.90) !important;",
              div(
                class = "d-flex align-items-center justify-content-between mb-2 pb-2 border-bottom",
                div(
                  class = "d-flex align-items-center gap-2",
                  tags$strong("Provenance Matches", class = "text-dark small fw-bold"),
                  uiOutput("lookup_results_count_badge", inline = TRUE)
                ),
                div(
                  class = "text-muted",
                  style = "font-size: 0.72rem;",
                  icon("table-cells", class = "me-1"),
                  "Sorted by N_Patients DESC, Redundancy DESC"
                )
              ),
              DTOutput("table_reverse_results")
            )
          )
        )
      )
    )
  )
)

# ------------------------------------------------------------------------------
# Server Logic
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Active in-memory DuckDB connection
  con <- get_db_con()
  onStop(function() {
    dbDisconnect(con, shutdown = TRUE)
  })
  
  # Automatically terminate the app and free resources when the browser window/tab is closed
  session$onSessionEnded(function() {
    message("Browser session closed. Server shutting down.")
    stopApp()
  })
  
  # Reactive state for query metrics (default empty until query runs)
  current_metrics <- reactiveValues(
    cassette_count = "—",
    read_count = "—",
    max_sharing = "—",
    query_time = "—",
    avg_cdr3_len = "—",
    avg_v_identity = "—"
  )

  # ----------------------------------------------------------------------------
  # Source Unified Master-Detail Alluvial Flow Canvas Generator
  # ----------------------------------------------------------------------------
  unified_candidates <- c(
    file.path(APP_DIR, "Scripts", "generate_unified_alluvial_svg.R"),
    file.path(getwd(), "OASpepDB", "Scripts", "generate_unified_alluvial_svg.R"),
    file.path(getwd(), "Scripts", "generate_unified_alluvial_svg.R"),
    file.path(dirname(APP_DIR), "Scripts", "generate_unified_alluvial_svg.R")
  )
  unified_gen_script <- unified_candidates[file.exists(unified_candidates)][1]
  if (!is.na(unified_gen_script) && file.exists(unified_gen_script)) {
    source(unified_gen_script, local = TRUE)
  }

  # ----------------------------------------------------------------------------
  # Interactive State & Observers for Unified Alluvial Canvas
  # ----------------------------------------------------------------------------
  selected_intro_disease <- reactiveVal("")
  default_export_dir <- if (dir.exists("D:/OAS/Fasta-export")) "D:/OAS/Fasta-export" else normalizePath(file.path(getwd(), "Exports"), winslash = "/", mustWork = FALSE)
  export_target_folder   <- reactiveVal(default_export_dir)
  db_root                <- reactiveVal(DEFAULT_DB_ROOT)

  # Render Database Status Badges
  render_db_status_content <- function() {
    path <- db_root()
    is_valid <- nzchar(path) && dir.exists(path)
    n_cohorts <- if (is_valid) get_db_cohort_count(path) else 0L
    if (is_valid && n_cohorts > 0) {
      tagList(
        tags$span(class = "badge bg-success-subtle text-success border border-success-subtle me-1", style = "font-size: 0.65rem; padding: 2px 5px;", sprintf("%d Cohorts", n_cohorts)),
        tags$span(class = "text-dark", "DB Loaded")
      )
    } else {
      tagList(
        tags$span(class = "badge bg-warning text-dark me-1", style = "font-size: 0.65rem; padding: 2px 5px;", "Not Loaded"),
        tags$span(class = "text-warning fw-semibold", "Load DB")
      )
    }
  }

  output$txt_db_status_badge     <- renderUI({ render_db_status_content() })
  output$txt_db_status_badge_rev <- renderUI({ render_db_status_content() })

  # Database Connection Modal Dialog
  open_db_modal <- function(err_msg = NULL) {
    cur_path <- isolate(db_root())
    if (!nzchar(cur_path)) {
      cur_path <- detect_default_db_root()
    }
    if (!nzchar(cur_path)) cur_path <- "D:/OAS/unpaired/CDR3_db"

    showModal(modalDialog(
      title = div(
        class = "d-flex align-items-center gap-2",
        icon("database", class = "text-primary"),
        tags$span("Load CDR3 Database (Hive Parquet)", class = "fw-bold", style = "font-size: 1.05rem;")
      ),
      size = "m",
      easyClose = TRUE,
      fade = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "btn_confirm_load_db",
          "Connect Database",
          icon = icon("plug"),
          class = "btn btn-primary fw-semibold",
          style = "background: #0071e3; border: none; border-radius: 8px; padding: 6px 18px;"
        )
      ),
      div(
        class = "p-1",
        tags$p(
          class = "text-muted small mb-3",
          "Select or enter the root directory containing the Hive-partitioned Parquet database (e.g. subfolders named ",
          tags$code("Disease=COVID-19"), ", ", tags$code("Disease=HIV"), ")."
        ),
        if (!is.null(err_msg) && nzchar(err_msg)) {
          div(
            class = "alert alert-danger py-2 px-3 small d-flex align-items-center gap-2 mb-3",
            icon("triangle-exclamation", class = "flex-shrink-0 text-danger"),
            tags$span(err_msg)
          )
        },
        tags$label(class = "form-label fw-bold small mb-1", "Database Directory:"),
        div(
          class = "d-flex gap-2 align-items-center mb-2",
          div(
            style = "flex-grow: 1;",
            textInput("txt_db_path_input", label = NULL, value = cur_path, width = "100%", placeholder = "e.g. D:/OAS/unpaired/CDR3_db")
          ),
          actionButton(
            "btn_browse_db_folder",
            "Browse...",
            icon = icon("folder-open"),
            class = "btn btn-outline-secondary",
            style = "height: 38px; padding: 0 14px; font-size: 0.85rem; display: inline-flex; align-items: center; gap: 6px; white-space: nowrap;"
          )
        ),
        div(
          class = "text-muted",
          style = "font-size: 0.76rem;",
          icon("circle-info", class = "me-1"),
          "Click 'Browse...' to select a folder on Windows, or paste/type any local or network directory path."
        )
      )
    ))
  }

  observeEvent(input$btn_open_db_modal_header, {
    open_db_modal()
  })

  observeEvent(input$btn_open_db_modal_reverse, {
    open_db_modal()
  })

  observeEvent(input$btn_browse_db_folder, {
    init_dir <- isolate(input$txt_db_path_input)
    if (is.null(init_dir) || !dir.exists(init_dir)) init_dir <- isolate(db_root())
    if (is.null(init_dir) || !dir.exists(init_dir)) init_dir <- getwd()
    
    new_dir <- tryCatch({
      if (exists("choose.dir", where = asNamespace("utils"))) {
        utils::choose.dir(default = init_dir, caption = "Select CDR3_db Hive Database Directory")
      } else {
        NA
      }
    }, error = function(e) NA)
    
    if (!is.na(new_dir) && nzchar(new_dir)) {
      new_dir <- gsub("\\\\", "/", new_dir)
      updateTextInput(session, "txt_db_path_input", value = new_dir)
    }
  })

  observeEvent(input$btn_confirm_load_db, {
    target_dir <- trimws(input$txt_db_path_input)
    if (is.null(target_dir) || !nzchar(target_dir)) {
      open_db_modal("Please enter or select a directory path.")
      return()
    }
    target_dir <- gsub("\\\\", "/", target_dir)
    if (!dir.exists(target_dir)) {
      open_db_modal(sprintf("Directory does not exist: %s", target_dir))
      return()
    }
    
    dis_dirs <- list.dirs(target_dir, full.names = FALSE, recursive = FALSE)
    cohorts <- grep("^Disease=", dis_dirs, value = TRUE)
    if (length(cohorts) == 0) {
      open_db_modal(sprintf("Directory exists but does not contain Hive partitions (e.g. 'Disease=COVID-19'): %s", target_dir))
      return()
    }
    
    db_root(target_dir)
    removeModal()
    showNotification(
      sprintf("Database loaded successfully! Connected to %d disease cohorts in %s.", length(cohorts), basename(target_dir)),
      type = "message",
      duration = 4
    )
  })
  
  observeEvent(input$btn_cohort_click, {
    req(input$btn_cohort_click)
    cur <- selected_intro_disease()
    if (identical(cur, input$btn_cohort_click)) {
      selected_intro_disease("")
    } else {
      selected_intro_disease(input$btn_cohort_click)
    }
  })
  
  observeEvent(input$sel_intro_quick, {
    if (nzchar(input$sel_intro_quick)) {
      selected_intro_disease(input$sel_intro_quick)
    }
  })
  
  observeEvent(input$btn_reset_landing, {
    selected_intro_disease("")
    updateSelectInput(session, "sel_intro_quick", selected = "")
    updateSelectInput(session, "intro_min_patients", selected = 1)
    updateSelectInput(session, "intro_min_redundancy", selected = 1)
    gc(verbose = FALSE, full = TRUE)
  })
  
  # ----------------------------------------------------------------------------
  # Reactive Micro Data for Landing Page (Shared between toolbar, alluvial, & export)
  # ----------------------------------------------------------------------------
  intro_micro_data <- reactive({
    dis <- selected_intro_disease()
    if (!nzchar(dis)) return(NULL)
    
    cur_db <- db_root()
    if (!nzchar(cur_db) || !dir.exists(cur_db)) {
      showNotification("Please load the CDR3 database directory first.", type = "warning", duration = 4)
      open_db_modal()
      return(NULL)
    }
    
    min_p <- as.numeric(input$intro_min_patients)
    if (is.null(min_p) || is.na(min_p)) min_p <- 1
    
    min_r <- as.numeric(input$intro_min_redundancy)
    if (is.null(min_r) || is.na(min_r)) min_r <- 1
    
    sql <- sprintf("
      SELECT BSource, BType, Isotype, count(*) as n_cassettes
      FROM read_parquet('%s/Disease=%s/**/*.parquet', hive_partitioning=true)
      WHERE N_Patients >= %d AND Redundancy >= %d
      GROUP BY BSource, BType, Isotype
    ", cur_db, dis, min_p, min_r)
    
    df <- tryCatch({
      dbGetQuery(con, sql)
    }, error = function(e) {
      data.frame(BSource = character(), BType = character(), Isotype = character(), n_cassettes = numeric())
    })
    
    df$BSource[is.na(df$BSource) | !nzchar(df$BSource)] <- "Unknown Source"
    df$BType[is.na(df$BType) | !nzchar(df$BType)] <- "Unsorted"
    df$Isotype[is.na(df$Isotype) | !nzchar(df$Isotype)] <- "Unknown Isotype"
    df
  })
  
  # ----------------------------------------------------------------------------
  # Center Unified Toolbar (Disease, Cassettes, Min Patients, Min Reads in 1 Frame)
  # ----------------------------------------------------------------------------
  output$intro_center_toolbar <- renderUI({
    dis <- selected_intro_disease()
    cur_p <- isolate(input$intro_min_patients)
    if (is.null(cur_p) || is.na(cur_p)) cur_p <- 1
    cur_r <- isolate(input$intro_min_redundancy)
    if (is.null(cur_r) || is.na(cur_r)) cur_r <- 1
    
    if (!nzchar(dis)) {
      # Baseline overview: Quick-select dropdown + 65.5M total count + Filters
      div(
        class = "glass-toolbar",
        # 1. Disease Quick Select
        div(
          class = "glass-segment",
          selectInput(
            "sel_intro_quick",
            label = NULL,
            choices = c("— Select Cohort —" = "", DISEASE_CHOICES),
            selected = "",
            selectize = FALSE,
            width = "165px"
          )
        ),
        # 2. Total Database CDR3 Peptide Count with Liquid Glass Tooltip
        div(
          class = "glass-segment",
          div(
            class = "glass-tooltip-wrap",
            span(
              class = "badge bg-white text-secondary border",
              style = "font-size: 0.72rem; font-weight: 600; padding: 4px 8px; border-radius: 6px; border-color: rgba(0,0,0,0.12) !important; font-family: ui-monospace, SFMono-Regular, monospace; cursor: help;",
              "65.5M CDR3 peptides"
            ),
            div(
              class = "glass-tooltip-box",
              style = "width: 240px;",
              div(class = "glass-tooltip-title", icon("database", class = "me-1"), "Database Summary"),
              "65,510,795 total non-redundant CDR3 peptides across all 25 disease cohorts."
            )
          )
        ),
        # Divider
        div(class = "glass-divider"),
        # 3. Min Patients Filter with Info Icon & Liquid Glass Tooltip
        div(
          class = "glass-segment",
          div(
            class = "glass-tooltip-wrap",
            tags$span("Min Patients:", class = "glass-label"),
            icon("circle-info", class = "glass-info-icon"),
            div(
              class = "glass-tooltip-box",
              div(class = "glass-tooltip-title", icon("users", class = "me-1"), "Min Patients (N)"),
              "Minimum number of individual patients sharing this CDR3 peptide (N ≥ 1, 2, 3, 5). Filters out private somatic mutations and singletons to reveal public, convergent antibodies."
            )
          ),
          selectInput(
            "intro_min_patients",
            label = NULL,
            choices = c("N ≥ 1" = 1, "N ≥ 2" = 2, "N ≥ 3" = 3, "N ≥ 5" = 5),
            selected = cur_p,
            selectize = FALSE,
            width = "78px"
          )
        ),
        # Divider
        div(class = "glass-divider"),
        # 4. Min Reads Filter with Info Icon & Liquid Glass Tooltip
        div(
          class = "glass-segment",
          div(
            class = "glass-tooltip-wrap",
            tags$span("Min Reads:", class = "glass-label"),
            icon("circle-info", class = "glass-info-icon"),
            div(
              class = "glass-tooltip-box",
              div(class = "glass-tooltip-title", icon("barcode", class = "me-1"), "Min Reads (Depth)"),
              "Minimum sequencing read depth / redundancy in OAS (≥ 1, 5, 10, 50 reads). Filters out sequencing noise to ensure high-confidence clonotype expression."
            )
          ),
          selectInput(
            "intro_min_redundancy",
            label = NULL,
            choices = c("≥ 1 read" = 1, "≥ 5 reads" = 5, "≥ 10 reads" = 10, "≥ 50 reads" = 50),
            selected = cur_r,
            selectize = FALSE,
            width = "92px"
          )
        )
      )
    } else {
      # Specific disease selected: Cohort Badge + Dynamic Filtered Count + Filters
      df <- intro_micro_data()
      row <- COHORT_SUMMARY_DATA %>% filter(Disease == dis)
      cat_name <- if (nrow(row) > 0) row$Category[1] else "Cohort"
      
      min_p_num <- as.numeric(cur_p)
      min_r_num <- as.numeric(cur_r)
      cas_num <- if (!is.null(df) && nrow(df) > 0) sum(df$n_cassettes) else (if (nrow(row) > 0) row$Cassettes[1] else 0)
      fmt_cas <- ifelse(cas_num >= 1e6, sprintf("%.1fM", cas_num / 1e6),
                 ifelse(cas_num >= 1e3, sprintf("%.0fK", cas_num / 1e3), as.character(cas_num)))
      
      badge_style <- switch(cat_name,
        "Infectious Diseases" = "background-color: rgba(213, 94, 0, 0.12); color: #D55E00; border: 1px solid rgba(213, 94, 0, 0.25);",
        "Allergy & Airways"   = "background-color: rgba(0, 158, 115, 0.12); color: #009E73; border: 1px solid rgba(0, 158, 115, 0.25);",
        "Autoimmune & Neuro"  = "background-color: rgba(0, 114, 178, 0.12); color: #0072B2; border: 1px solid rgba(0, 114, 178, 0.25);",
        "Hematology & Tumors" = "background-color: rgba(204, 121, 167, 0.12); color: #CC79A7; border: 1px solid rgba(204, 121, 167, 0.25);",
        "background-color: rgba(0, 113, 227, 0.10); color: #0071e3; border: 1px solid rgba(0, 113, 227, 0.20);"
      )
      
      div(
        class = "d-inline-flex align-items-center gap-2",
        tags$span("Selected data:", class = "text-dark small fw-semibold text-nowrap", style = "font-size: 0.78rem; letter-spacing: -0.01em;"),
        div(
          class = "glass-toolbar",
          # 1. Disease Name Badge with Liquid Glass Tooltip
          div(
            class = "glass-segment",
            div(
              class = "glass-tooltip-wrap",
              span(
                class = "badge",
                style = paste0(badge_style, " font-size: 0.74rem; font-weight: 700; padding: 4px 9px; border-radius: 6px; letter-spacing: -0.01em; cursor: help;"),
                dis
              ),
              div(
                class = "glass-tooltip-box",
                style = "width: 220px;",
                div(class = "glass-tooltip-title", icon("virus", class = "me-1"), "Selected Cohort"),
                sprintf("Cohort: %s (%s). Click Reset to restore all 25 cohorts.", dis, cat_name)
              )
            )
          ),
          # 2. Filtered CDR3 Peptides Count Badge with Liquid Glass Tooltip
          div(
            class = "glass-segment",
            div(
              class = "glass-tooltip-wrap",
              span(
                class = "badge bg-white text-dark border", 
                style = "font-size: 0.72rem; font-weight: 600; padding: 4px 8px; border-radius: 6px; border-color: rgba(0,0,0,0.12) !important; font-family: ui-monospace, SFMono-Regular, monospace; cursor: help;",
                sprintf("%s CDR3 peptides", fmt_cas)
              ),
              div(
                class = "glass-tooltip-box",
                style = "width: 240px;",
                div(class = "glass-tooltip-title", icon("layer-group", class = "me-1"), "Filtered CDR3 Peptides"),
                sprintf("Filtered non-redundant CDR3 peptides for %s (N ≥ %d, Reads ≥ %d): %s CDR3 peptides.", dis, min_p_num, min_r_num, format(cas_num, big.mark = ","))
              )
            )
          ),
          # Divider
          div(class = "glass-divider"),
          # 3. Min Patients Filter with Info Icon & Liquid Glass Tooltip
          div(
            class = "glass-segment",
            div(
              class = "glass-tooltip-wrap",
              tags$span("Min Patients:", class = "glass-label"),
              icon("circle-info", class = "glass-info-icon"),
              div(
                class = "glass-tooltip-box",
                div(class = "glass-tooltip-title", icon("users", class = "me-1"), "Min Patients (N)"),
                "Minimum number of individual patients sharing this CDR3 peptide (N ≥ 1, 2, 3, 5). Filters out private somatic mutations and singletons to reveal public, convergent antibodies."
              )
            ),
            selectInput(
              "intro_min_patients",
              label = NULL,
              choices = c("N ≥ 1" = 1, "N ≥ 2" = 2, "N ≥ 3" = 3, "N ≥ 5" = 5),
              selected = cur_p,
              selectize = FALSE,
              width = "78px"
            )
          ),
          # Divider
          div(class = "glass-divider"),
          # 4. Min Reads Filter with Info Icon & Liquid Glass Tooltip
          div(
            class = "glass-segment",
            div(
              class = "glass-tooltip-wrap",
              tags$span("Min Reads:", class = "glass-label"),
              icon("circle-info", class = "glass-info-icon"),
              div(
                class = "glass-tooltip-box",
                div(class = "glass-tooltip-title", icon("barcode", class = "me-1"), "Min Reads (Depth)"),
                "Minimum sequencing read depth / redundancy in OAS (≥ 1, 5, 10, 50 reads). Filters out sequencing noise to ensure high-confidence clonotype expression."
              )
            ),
            selectInput(
              "intro_min_redundancy",
              label = NULL,
              choices = c("≥ 1 read" = 1, "≥ 5 reads" = 5, "≥ 10 reads" = 10, "≥ 50 reads" = 50),
              selected = cur_r,
              selectize = FALSE,
              width = "92px"
            )
          )
        )
      )
    }
  })
  
  # ----------------------------------------------------------------------------
  # Top-Right Export Controls (Liquid Glass Split Button)
  # ----------------------------------------------------------------------------
  output$intro_export_controls <- renderUI({
    dis <- selected_intro_disease()
    if (!nzchar(dis)) {
      div(
        class = "glass-tooltip-wrap",
        tags$button(
          type = "button",
          class = "btn glass-action-btn btn-sm disabled",
          style = "opacity: 0.45; cursor: not-allowed;",
          icon("folder-open"),
          tags$span("Export", class = "ms-1")
        ),
        div(
          class = "glass-tooltip-box align-end",
          style = "width: 220px;",
          div(class = "glass-tooltip-title", icon("folder-open", class = "me-1"), "Export Dataset"),
          "Please select a disease cohort to enable dataset export."
        )
      )
    } else {
      div(
        class = "glass-tooltip-wrap",
        div(
          class = "btn-group",
          actionButton(
            "btn_intro_export_modal",
            label = "Export",
            icon = icon("folder-open"),
            class = "btn glass-action-btn btn-sm"
          ),
          tags$button(
            type = "button",
            class = "btn glass-action-btn glass-action-toggle btn-sm dropdown-toggle dropdown-toggle-split",
            `data-bs-toggle` = "dropdown",
            `aria-expanded` = "false",
            tags$span(class = "visually-hidden", "Toggle Dropdown")
          ),
          tags$ul(
            class = "dropdown-menu dropdown-menu-end shadow-sm",
            style = "font-size: 0.78rem; border-radius: 10px; border: 1px solid rgba(0,0,0,0.08); padding: 4px; backdrop-filter: blur(20px); background: rgba(255, 255, 255, 0.96);",
            tags$li(
              tags$h6(class = "dropdown-header px-2 py-1 text-uppercase text-muted", style = "font-size: 0.65rem; letter-spacing: 0.5px;", "Direct Folder Export")
            ),
            tags$li(
              actionLink(
                "btn_intro_export_modal_link",
                class = "dropdown-item py-1 px-2 rounded d-flex align-items-center fw-semibold text-primary",
                tagList(icon("folder-open", class = "me-2"), "Save 3 Files to Folder...")
              )
            ),
            tags$li(tags$hr(class = "dropdown-divider my-1")),
            tags$li(
              tags$h6(class = "dropdown-header px-2 py-1 text-uppercase text-muted", style = "font-size: 0.65rem; letter-spacing: 0.5px;", "Individual Formats")
            ),
            tags$li(
              downloadLink(
                "btn_intro_export_fasta",
                class = "dropdown-item py-1 px-2 rounded d-flex align-items-center",
                tagList(icon("file-lines", class = "text-primary me-2"), "FASTA Library (.fasta)")
              )
            ),
            tags$li(
              downloadLink(
                "btn_intro_export_parquet",
                class = "dropdown-item py-1 px-2 rounded d-flex align-items-center",
                tagList(icon("table", class = "text-success me-2"), "Parquet Table (.parquet)")
              )
            ),
            tags$li(
              downloadLink(
                "btn_intro_export_report",
                class = "dropdown-item py-1 px-2 rounded d-flex align-items-center",
                tagList(icon("file-alt", class = "text-secondary me-2"), "Audit Report (.txt)")
              )
            )
          )
        ),
        div(
          class = "glass-tooltip-box align-end",
          style = "width: 250px;",
          div(class = "glass-tooltip-title", icon("folder-open", class = "me-1"), "Export Dataset"),
          sprintf("Save 3 uncompressed files (FASTA, Parquet table, audit report) for %s to your chosen folder, or select individual formats.", dis)
        )
      )
    }
  })

  # ----------------------------------------------------------------------------
  # Render Unified Alluvial Flow SVG Output
  # ----------------------------------------------------------------------------
  output$landing_unified_alluvial <- renderUI({
    dis <- selected_intro_disease()
    if (!nzchar(dis)) {
      svg_code <- generate_unified_alluvial_svg(selected_dis = "", df_micro = NULL)
      return(HTML(svg_code))
    }
    
    df <- intro_micro_data()
    svg_code <- generate_unified_alluvial_svg(selected_dis = dis, df_micro = df)
    HTML(svg_code)
  })

  # Handle "Back to All 25 Cohorts" button
  observeEvent(input$btn_back_directory, {
    updateSelectInput(session, "sel_disease", selected = "")
    selected_intro_disease("")
    updateNavbarPage(session, "main_nav", selected = "Query & Export")
    session$sendCustomMessage("switch_tab", "Query & Export")
  })

  # Handle "Go to Query & Export" button from empty Query page
  observeEvent(input$btn_go_intro, {
    selected_intro_disease("")
    updateNavbarPage(session, "main_nav", selected = "Query & Export")
    session$sendCustomMessage("switch_tab", "Query & Export")
  })

  # ----------------------------------------------------------------------------
  # Dynamic Partition Filter Population (BSource, BType, Isotype)
  # ----------------------------------------------------------------------------
  observeEvent(input$sel_disease, {
    disease <- input$sel_disease
    if (!nzchar(disease)) {
      updateSelectInput(session, "sel_bsource", choices = c("All" = "ALL"), selected = "ALL")
      updateSelectInput(session, "sel_btype",   choices = c("All" = "ALL"), selected = "ALL")
      updateSelectInput(session, "sel_isotype", choices = c("All" = "ALL"), selected = "ALL")
      current_metrics$cassette_count <- "—"
      current_metrics$read_count     <- "—"
      current_metrics$max_sharing    <- "—"
      current_metrics$query_time     <- "—"
      current_metrics$avg_cdr3_len   <- "—"
      current_metrics$avg_v_identity <- "—"
      filtered_data(data.frame())
      return()
    }
    
    cur_db <- db_root()
    if (!nzchar(cur_db) || !dir.exists(cur_db)) {
      showNotification("Please load the CDR3 database directory first.", type = "warning", duration = 4)
      open_db_modal()
      return()
    }
    
    # Query partition directory names via DuckDB Hive metadata (takes ~15ms)
    meta_sql <- sprintf("
      SELECT DISTINCT BSource, BType, Isotype 
      FROM read_parquet('%s/Disease=%s/**/*.parquet', hive_partitioning=true);
    ", cur_db, disease)
    
    part_meta <- tryCatch({
      dbGetQuery(con, meta_sql)
    }, error = function(e) {
      data.frame(BSource = character(), BType = character(), Isotype = character())
    })
    
    if (nrow(part_meta) > 0) {
      bsources <- sort(unique(part_meta$BSource))
      btypes   <- sort(unique(part_meta$BType))
      isotypes <- sort(unique(part_meta$Isotype))
      
      updateSelectInput(session, "sel_bsource", choices = c("All" = "ALL", bsources), selected = "ALL")
      updateSelectInput(session, "sel_btype",   choices = c("All" = "ALL", btypes),   selected = "ALL")
      updateSelectInput(session, "sel_isotype", choices = c("All" = "ALL", isotypes), selected = "ALL")
    }
    
    # Auto-execute query to immediately load cohort profile and table
    execute_main_query()
  }, ignoreInit = TRUE)

  # ----------------------------------------------------------------------------
  # SQL Query Builder
  # ----------------------------------------------------------------------------
  get_query_components <- function() {
    req(input$sel_disease)
    disease <- input$sel_disease
    
    # Build partition path with pruning
    b_src <- if (input$sel_bsource != "ALL") paste0("BSource=", input$sel_bsource, "/") else "**/"
    b_typ <- if (input$sel_btype != "ALL")   paste0("BType=", input$sel_btype, "/")     else ""
    iso   <- if (input$sel_isotype != "ALL") paste0("Isotype=", input$sel_isotype, "/") else ""
    
    # Path globbing
    path_pattern <- sprintf("%s/Disease=%s/%s%s%s*.parquet", db_root(), disease, b_src, b_typ, iso)
    
    # Where clauses for non-partitioned metrics
    clauses <- c(
      sprintf("Disease = '%s'", disease),
      sprintf("N_Patients >= %d", input$num_min_patients),
      sprintf("Redundancy >= %d", input$num_min_redundancy)
    )
    
    if (input$sel_bsource != "ALL") clauses <- c(clauses, sprintf("BSource = '%s'", input$sel_bsource))
    if (input$sel_btype != "ALL")   clauses <- c(clauses, sprintf("BType = '%s'", input$sel_btype))
    if (input$sel_isotype != "ALL") clauses <- c(clauses, sprintf("Isotype = '%s'", input$sel_isotype))
    
    where_str <- paste(clauses, collapse = " AND ")
    list(path = path_pattern, where = where_str)
  }

  build_query_sql <- function(limit = NULL) {
    comp <- get_query_components()
    limit_clause <- if (!is.null(limit)) sprintf("LIMIT %d", limit) else ""
    
    sprintf("
      SELECT 
        clonotype_id,
        cdr3_cassette,
        cdr3_aa,
        cdr3_len,
        fr3_tail,
        fwr4_aa,
        isotype_tail,
        tryptic_peptides,
        v_call,
        d_call,
        j_call,
        v_identity,
        Redundancy,
        N_Patients,
        Data_file,
        Disease,
        BSource,
        BType,
        Isotype
      FROM read_parquet('%s', hive_partitioning=true)
      WHERE %s
      ORDER BY N_Patients DESC, Redundancy DESC
      %s
    ", comp$path, comp$where, limit_clause)
  }

  # ----------------------------------------------------------------------------
  # Reactive Data Fetching
  # ----------------------------------------------------------------------------
  filtered_data <- reactiveVal(data.frame())
  
  execute_main_query <- function() {
    if (!nzchar(input$sel_disease)) {
      showNotification("Please select a Disease Cohort from the sidebar first.", type = "warning", duration = 4)
      return()
    }
    withProgress(message = "Querying disease CDR3 peptides...", value = 0.5, {
      t0 <- Sys.time()
      comp <- get_query_components()
      
      # 1. Fetch Aggregated Metrics directly
      metric_sql <- sprintf("
        SELECT 
          COUNT(*) as total_cassettes,
          COALESCE(SUM(Redundancy), 0) as total_reads,
          COALESCE(MAX(N_Patients), 0) as max_sharing,
          COALESCE(ROUND(AVG(cdr3_len), 1), 0) as avg_cdr3_len,
          COALESCE(MIN(cdr3_len), 0) as min_cdr3_len,
          COALESCE(MAX(cdr3_len), 0) as max_cdr3_len,
          COALESCE(ROUND(AVG(v_identity), 1), 0) as avg_v_identity
        FROM read_parquet('%s', hive_partitioning=true)
        WHERE %s;
      ", comp$path, comp$where)
      
      m_df <- tryCatch({
        dbGetQuery(con, metric_sql)
      }, error = function(e) {
        data.frame(total_cassettes = 0, total_reads = 0, max_sharing = 0, 
                   avg_cdr3_len = 0, min_cdr3_len = 0, max_cdr3_len = 0, avg_v_identity = 0)
      })
      
      # 2. Fetch Top 10,000 for Interactive Datatable
      data_sql <- build_query_sql(limit = 10000)
      df <- tryCatch({
        dbGetQuery(con, data_sql)
      }, error = function(e) {
        data.frame()
      })
      
      t1 <- Sys.time()
      elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
      
      current_metrics$cassette_count <- format(m_df$total_cassettes, big.mark = ",")
      current_metrics$read_count     <- format(m_df$total_reads, big.mark = ",")
      current_metrics$max_sharing    <- paste0(m_df$max_sharing, " Patients")
      current_metrics$query_time     <- sprintf("%.3f sec", elapsed)
      current_metrics$avg_cdr3_len   <- sprintf("%.1f AA (Range: %d - %d AA)", m_df$avg_cdr3_len, m_df$min_cdr3_len, m_df$max_cdr3_len)
      current_metrics$avg_v_identity <- sprintf("%.1f%% V-Identity", m_df$avg_v_identity)
      
      filtered_data(df)
    })
  }
  
  # Trigger query on button press
  observeEvent(input$btn_run_query, {
    execute_main_query()
  })

  # ----------------------------------------------------------------------------
  # Value Box & Property Outputs
  # ----------------------------------------------------------------------------
  output$vb_cassette_count <- renderText({ current_metrics$cassette_count })
  output$vb_read_count     <- renderText({ current_metrics$read_count })
  output$vb_max_sharing    <- renderText({ current_metrics$max_sharing })
  output$vb_query_time     <- renderText({ current_metrics$query_time })

  # Cohort Header & Property Text Outputs
  output$txt_selected_disease_title <- renderText({
    req(input$sel_disease)
    input$sel_disease
  })
  
  output$txt_selected_disease_desc <- renderText({
    req(input$sel_disease)
    DISEASE_MAP[input$sel_disease]
  })

  output$txt_prop_cdr3_len <- renderText({
    current_metrics$avg_cdr3_len
  })
  
  output$txt_prop_shm <- renderText({
    current_metrics$avg_v_identity
  })
  
  output$txt_prop_partitions <- renderText({
    req(input$sel_disease)
    sprintf("Source: %s | Subtype: %s | Isotype: %s", input$sel_bsource, input$sel_btype, input$sel_isotype)
  })

  # ----------------------------------------------------------------------------
  # Main Interactive DataTable
  # ----------------------------------------------------------------------------
  output$table_clonotypes <- renderDT({
    df <- filtered_data()
    if (nrow(df) == 0) {
      msg <- if (!nzchar(input$sel_disease)) {
        "Please select a Disease Cohort from the sidebar and click 'Execute Query' to begin exploring."
      } else {
        "No matching disease CDR3 peptides found with the specified filter criteria."
      }
      return(datatable(
        data.frame(Status = msg), 
        options = list(dom = "t", ordering = FALSE),
        rownames = FALSE
      ))
    }
    
    display_df <- df %>%
      mutate(v_identity = round(as.numeric(v_identity), 2)) %>%
      select(
        `Clonotype ID` = clonotype_id,
        `CDR3 (AA)` = cdr3_aa,
        `Length` = cdr3_len,
        `Isotype` = Isotype,
        `Patients` = N_Patients,
        `Redundancy` = Redundancy,
        `V-Call` = v_call,
        `J-Call` = j_call,
        `V-Identity (%)` = v_identity,
        `Source Files` = Data_file
      )
    
    datatable(
      display_df,
      filter = "top",
      selection = "single",
      rownames = FALSE,
      extensions = c("Buttons"),
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 25, 50, 100),
        autoWidth = FALSE,
        scrollX = TRUE,
        order = list(list(4, "desc"), list(5, "desc")),
        columnDefs = list(
          list(
            targets = 9,
            render = JS(
              "function(data, type, row) {",
              "  if (typeof data !== 'string' || !data) return data || '';",
              "  var files = data.split(';');",
              "  if (files.length <= 1) return data;",
              "  return files[0] + ' ... (+' + (files.length - 1) + ' files)';",
              "}"
            )
          )
        )
      ),
      class = "display nowrap compact stripe hover"
    ) %>%
      formatRound("V-Identity (%)", 2)
  })

  # ----------------------------------------------------------------------------
  # Shared Export Logic (FASTA, Parquet, Text Audit)
  # Direct C++ DuckDB COPY for zero-memory sub-second streaming
  # ----------------------------------------------------------------------------
  export_fasta_worker <- function(file, disease = NULL, min_p = NULL, min_r = NULL) {
    if (is.null(disease) || !nzchar(disease)) disease <- selected_intro_disease()
    if (!nzchar(disease) && !is.null(input$sel_disease)) disease <- input$sel_disease
    req(nzchar(disease))
    
    if (is.null(min_p)) min_p <- as.numeric(input$intro_min_patients)
    if (is.null(min_p) || is.na(min_p)) min_p <- 1
    if (is.null(min_r)) min_r <- as.numeric(input$intro_min_redundancy)
    if (is.null(min_r) || is.na(min_r)) min_r <- 1
    
    out_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
    path_pattern <- sprintf("%s/Disease=%s/**/*.parquet", db_root(), disease)
    where_clause <- sprintf("Disease = '%s' AND N_Patients >= %d AND Redundancy >= %d", disease, min_p, min_r)
    
    sql <- sprintf("
      COPY (
        SELECT 
          '>sp|' || clonotype_id || '|' || clonotype_id || 
          ' OAS CDR3 Micro-cassette Disease=' || FIRST(Disease) || 
          ' N_Patients=' || CAST(MAX(N_Patients) AS VARCHAR) || 
          ' Redundancy=' || CAST(SUM(Redundancy) AS VARCHAR) || 
          ' OS=Homo sapiens OX=9606 GN=' || 
          CASE 
            WHEN FIRST(v_call) LIKE 'IGL%%' OR FIRST(j_call) LIKE 'IGL%%' THEN 'IGL'
            WHEN FIRST(v_call) LIKE 'IGK%%' OR FIRST(j_call) LIKE 'IGK%%' THEN 'IGK'
            ELSE 'IGH'
          END || 
          ' PE=4 SV=1\n' || replace(cdr3_cassette, 'I', 'L') AS fasta_record
        FROM read_parquet('%s', hive_partitioning=true)
        WHERE %s
        GROUP BY clonotype_id, cdr3_cassette
        ORDER BY MAX(N_Patients) DESC, SUM(Redundancy) DESC
      ) TO '%s' (FORMAT CSV, HEADER FALSE, QUOTE '', ESCAPE '');
    ", path_pattern, where_clause, out_file)
    
    tryCatch({
      dbExecute(con, sql)
      
      # 1. Append Common Contaminants (cRAP, PE=1)
      if (file.exists(CRAP_REF_PATH)) {
        file.append(out_file, CRAP_REF_PATH)
      }
      # 2. Append Canonical Human Swiss-Prot Reference (PE=1)
      if (file.exists(UNIPROT_REF_PATH)) {
        file.append(out_file, UNIPROT_REF_PATH)
      }
      # 3. Append Camelid Entrapment Controls (PE=4, Auto-scaled to ~1% in round hundreds)
      if (file.exists(ENTRAPMENT_REF_PATH)) {
        count_sql <- sprintf("
          SELECT COUNT(*) AS n FROM (
            SELECT clonotype_id 
            FROM read_parquet('%s', hive_partitioning=true)
            WHERE %s
            GROUP BY clonotype_id, cdr3_cassette
          );
        ", path_pattern, where_clause)
        n_cassettes <- tryCatch({
          as.integer(dbGetQuery(con, count_sql)$n[1])
        }, error = function(e) 0L)
        
        target_entrapment <- max(100L, as.integer(floor((n_cassettes * 0.01) / 100 + 0.5) * 100))
        
        entrap_lines <- readLines(ENTRAPMENT_REF_PATH, encoding = "UTF-8", warn = FALSE)
        header_indices <- grep("^>", entrap_lines)
        total_available <- length(header_indices)
        actual_entrapment <- min(target_entrapment, total_available)
        
        if (actual_entrapment < total_available) {
          end_line <- header_indices[actual_entrapment + 1] - 1
          entrap_subset <- entrap_lines[1:end_line]
        } else {
          entrap_subset <- entrap_lines
        }
        
        tmp_entrap <- tempfile(fileext = ".fasta")
        writeLines(entrap_subset, con = tmp_entrap, useBytes = TRUE)
        file.append(out_file, tmp_entrap)
        unlink(tmp_entrap)
      }
    }, error = function(e) {
      warning(paste("FASTA Export error:", e$message))
      writeLines(paste("Error generating FASTA:", e$message), con = file)
    })
  }

  export_parquet_worker <- function(file, disease = NULL, min_p = NULL, min_r = NULL) {
    if (is.null(disease) || !nzchar(disease)) disease <- selected_intro_disease()
    if (!nzchar(disease) && !is.null(input$sel_disease)) disease <- input$sel_disease
    req(nzchar(disease))
    
    if (is.null(min_p)) min_p <- as.numeric(input$intro_min_patients)
    if (is.null(min_p) || is.na(min_p)) min_p <- 1
    if (is.null(min_r)) min_r <- as.numeric(input$intro_min_redundancy)
    if (is.null(min_r) || is.na(min_r)) min_r <- 1
    
    out_file <- normalizePath(file, winslash = "/", mustWork = FALSE)
    path_pattern <- sprintf("%s/Disease=%s/**/*.parquet", db_root(), disease)
    where_clause <- sprintf("Disease = '%s' AND N_Patients >= %d AND Redundancy >= %d", disease, min_p, min_r)
    
    sql <- sprintf("
      COPY (
        SELECT 
          clonotype_id, cdr3_cassette, cdr3_aa, cdr3_len, fr3_tail, fwr4_aa, isotype_tail,
          tryptic_peptides, v_call, d_call, j_call, v_identity, Redundancy, N_Patients,
          Data_file, Disease, BSource, BType, Isotype
        FROM read_parquet('%s', hive_partitioning=true)
        WHERE %s
        ORDER BY N_Patients DESC, Redundancy DESC
      ) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD);
    ", path_pattern, where_clause, out_file)
    
    tryCatch({
      dbExecute(con, sql)
    }, error = function(e) {
      warning(paste("Parquet Export error:", e$message))
    })
  }

  export_report_worker <- function(file, disease = NULL, min_p = NULL, min_r = NULL) {
    if (is.null(disease) || !nzchar(disease)) disease <- selected_intro_disease()
    if (!nzchar(disease) && !is.null(input$sel_disease)) disease <- input$sel_disease
    req(nzchar(disease))
    
    if (is.null(min_p)) min_p <- as.numeric(input$intro_min_patients)
    if (is.null(min_p) || is.na(min_p)) min_p <- 1
    if (is.null(min_r)) min_r <- as.numeric(input$intro_min_redundancy)
    if (is.null(min_r) || is.na(min_r)) min_r <- 1
    
    path_pattern <- sprintf("%s/Disease=%s/**/*.parquet", db_root(), disease)
    where_clause <- sprintf("Disease = '%s' AND N_Patients >= %d AND Redundancy >= %d", disease, min_p, min_r)
    
    metric_sql <- sprintf("
      SELECT 
        COUNT(*) as total_cassettes,
        COALESCE(SUM(Redundancy), 0) as total_reads,
        COALESCE(MAX(N_Patients), 0) as max_sharing,
        COALESCE(ROUND(AVG(cdr3_len), 1), 0) as avg_cdr3_len,
        COALESCE(MIN(cdr3_len), 0) as min_cdr3_len,
        COALESCE(MAX(cdr3_len), 0) as max_cdr3_len,
        COALESCE(ROUND(AVG(v_identity), 1), 0) as avg_v_identity
      FROM read_parquet('%s', hive_partitioning=true)
      WHERE %s;
    ", path_pattern, where_clause)
    
    m_df <- tryCatch({
      dbGetQuery(con, metric_sql)
    }, error = function(e) {
      data.frame(total_cassettes = 0, total_reads = 0, max_sharing = 0, 
                 avg_cdr3_len = 0, min_cdr3_len = 0, max_cdr3_len = 0, avg_v_identity = 0)
    })
    
    top_sql <- sprintf("
      SELECT clonotype_id, cdr3_aa, cdr3_len, Isotype, N_Patients, Redundancy, Data_file
      FROM read_parquet('%s', hive_partitioning=true)
      WHERE %s
      ORDER BY N_Patients DESC, Redundancy DESC
      LIMIT 20;
    ", path_pattern, where_clause)
    top20 <- tryCatch({
      dbGetQuery(con, top_sql)
    }, error = function(e) {
      data.frame()
    })
    
    dis_desc <- if (!is.na(DISEASE_MAP[disease])) DISEASE_MAP[disease] else disease
    
    report_lines <- c(
      "================================================================================",
      "  OASpepDB: CLINICAL PROVENANCE & QUERY AUDIT REPORT",
      "================================================================================",
      sprintf("Generated On:       %s", Sys.time()),
      sprintf("Disease Target:     %s (%s)", disease, dis_desc),
      sprintf("Min Patients (N):   %d", min_p),
      sprintf("Min Redundancy:     %d", min_r),
      "--------------------------------------------------------------------------------",
      sprintf("Total CDR3 Peptides: %s", format(m_df$total_cassettes, big.mark = ",")),
      sprintf("Total Read Count:   %s", format(m_df$total_reads, big.mark = ",")),
      sprintf("Max Cohort Sharing: %d Patients", m_df$max_sharing),
      sprintf("Mean CDR3 Length:   %.1f AA (Range: %d - %d AA)", m_df$avg_cdr3_len, m_df$min_cdr3_len, m_df$max_cdr3_len),
      sprintf("Mean SHM Identity:  %.1f%% V-Identity", m_df$avg_v_identity),
      "================================================================================",
      "TOP CONVERGENT CLONOTYPES (First 20):",
      ""
    )
    
    if (nrow(top20) > 0) {
      for (i in seq_len(nrow(top20))) {
        row <- top20[i, ]
        report_lines <- c(
          report_lines,
          sprintf("[%02d] %s | CDR3: %s (%d aa) | Isotype: %s | Patients: %d | Redundancy: %d",
                  i, row$clonotype_id, row$cdr3_aa, row$cdr3_len, row$Isotype, row$N_Patients, row$Redundancy),
          sprintf("     Files: %s", row$Data_file)
        )
      }
    }
    
    writeLines(report_lines, con = file)
  }

  # ----------------------------------------------------------------------------
  # Direct Uncompressed Folder Export (FASTA, Parquet, Audit Report)
  # Saves 3 separate files directly to user-selected folder without zipping
  # ----------------------------------------------------------------------------
  open_export_modal <- function() {
    dis <- selected_intro_disease()
    if (!nzchar(dis) && !is.null(input$sel_disease)) dis <- input$sel_disease
    req(nzchar(dis))
    
    min_p <- as.numeric(input$intro_min_patients)
    if (is.null(min_p) || is.na(min_p)) min_p <- 1
    min_r <- as.numeric(input$intro_min_redundancy)
    if (is.null(min_r) || is.na(min_r)) min_r <- 1
    
    date_str <- format(Sys.Date(), "%Y%m%d")
    fasta_name   <- sprintf("OAS_%s_%s.fasta", dis, date_str)
    parquet_name <- sprintf("OAS_%s_%s.parquet", dis, date_str)
    report_name  <- sprintf("OAS_%s_%s.txt", dis, date_str)
    
    cur_folder <- export_target_folder()
    
    showModal(modalDialog(
      title = div(
        class = "d-flex align-items-center gap-2",
        icon("folder-open", class = "text-primary"),
        tags$span("Export Dataset (3 Uncompressed Files)", class = "fw-bold", style = "font-size: 1.05rem;")
      ),
      size = "m",
      easyClose = TRUE,
      fade = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton(
          "btn_confirm_export",
          "Save 3 Files",
          icon = icon("check"),
          class = "btn btn-primary fw-semibold",
          style = "background: #0071e3; border: none; border-radius: 8px; padding: 6px 18px;"
        )
      ),
      div(
        class = "p-1",
        div(
          class = "d-flex flex-wrap gap-2 mb-3 pb-2 border-bottom",
          span(class = "badge rounded-pill bg-light text-dark border px-3 py-2",
               icon("virus", class = "me-1 text-primary"), sprintf("Cohort: %s", dis)),
          span(class = "badge rounded-pill bg-light text-dark border px-3 py-2",
               icon("users", class = "me-1 text-info"), sprintf("Min Patients: N ≥ %d", min_p)),
          span(class = "badge rounded-pill bg-light text-dark border px-3 py-2",
               icon("barcode", class = "me-1 text-success"), sprintf("Min Reads: ≥ %d", min_r))
        ),
        tags$p(class = "text-muted small mb-2 fw-semibold", "The following 3 files will be generated directly into your selected folder without zip compression:"),
        tags$ul(
          class = "list-group list-group-flush mb-3 small",
          style = "border-radius: 8px; border: 1px solid #e5e5ea; background: #fafafa;",
          tags$li(class = "list-group-item bg-transparent d-flex align-items-center py-2",
                  icon("file-lines", class = "text-primary me-2"),
                  tags$code(fasta_name, style = "color: #1d1d1f; font-size: 0.82rem;"),
                  span(class = "ms-auto text-muted small", "FASTA Library")
          ),
          tags$li(class = "list-group-item bg-transparent d-flex align-items-center py-2",
                  icon("table", class = "text-success me-2"),
                  tags$code(parquet_name, style = "color: #1d1d1f; font-size: 0.82rem;"),
                  span(class = "ms-auto text-muted small", "Parquet Table")
          ),
          tags$li(class = "list-group-item bg-transparent d-flex align-items-center py-2",
                  icon("file-alt", class = "text-secondary me-2"),
                  tags$code(report_name, style = "color: #1d1d1f; font-size: 0.82rem;"),
                  span(class = "ms-auto text-muted small", "Audit Report")
          )
        ),
        tags$label(class = "form-label fw-bold small mb-1", "Destination Folder:"),
        div(
          class = "d-flex gap-2 align-items-center mb-2",
          div(
            style = "flex-grow: 1;",
            textInput("txt_export_folder", label = NULL, value = cur_folder, width = "100%")
          ),
          actionButton(
            "btn_browse_export_folder",
            "Browse...",
            icon = icon("folder-open"),
            class = "btn btn-outline-secondary",
            style = "height: 38px; padding: 0 14px; font-size: 0.85rem; display: inline-flex; align-items: center; gap: 6px; white-space: nowrap;"
          )
        ),
        div(class = "text-muted", style = "font-size: 0.76rem;",
            icon("circle-info", class = "me-1"),
            "Click 'Browse...' to select a folder via Windows Explorer, or edit the path directly.")
      )
    ))
  }

  observeEvent(input$btn_intro_export_modal, {
    open_export_modal()
  })

  observeEvent(input$btn_intro_export_modal_link, {
    open_export_modal()
  })

  observeEvent(input$btn_browse_export_folder, {
    init_dir <- isolate(input$txt_export_folder)
    if (is.null(init_dir) || !dir.exists(init_dir)) init_dir <- export_target_folder()
    
    new_dir <- tryCatch({
      if (exists("choose.dir", where = asNamespace("utils"))) {
        utils::choose.dir(default = init_dir, caption = "Select Destination Folder for Export")
      } else {
        NA
      }
    }, error = function(e) NA)
    
    if (!is.na(new_dir) && nzchar(new_dir)) {
      new_dir <- gsub("\\\\", "/", new_dir)
      export_target_folder(new_dir)
      updateTextInput(session, "txt_export_folder", value = new_dir)
    }
  })

  observeEvent(input$btn_confirm_export, {
    target_dir <- trimws(input$txt_export_folder)
    if (is.null(target_dir) || !nzchar(target_dir)) {
      showNotification("Please select or enter a destination folder.", type = "warning")
      return()
    }
    
    target_dir <- gsub("\\\\", "/", target_dir)
    if (!dir.exists(target_dir)) {
      ok <- tryCatch({
        dir.create(target_dir, recursive = TRUE)
        TRUE
      }, error = function(e) FALSE)
      if (!ok || !dir.exists(target_dir)) {
        showNotification(sprintf("Cannot create or write to folder: %s", target_dir), type = "error")
        return()
      }
    }
    
    export_target_folder(target_dir)
    removeModal()
    
    dis <- selected_intro_disease()
    if (!nzchar(dis) && !is.null(input$sel_disease)) dis <- input$sel_disease
    req(nzchar(dis))
    
    min_p <- as.numeric(input$intro_min_patients)
    if (is.null(min_p) || is.na(min_p)) min_p <- 1
    min_r <- as.numeric(input$intro_min_redundancy)
    if (is.null(min_r) || is.na(min_r)) min_r <- 1
    
    date_str <- format(Sys.Date(), "%Y%m%d")
    fasta_name   <- sprintf("OAS_%s_%s.fasta", dis, date_str)
    parquet_name <- sprintf("OAS_%s_%s.parquet", dis, date_str)
    report_name  <- sprintf("OAS_%s_%s.txt", dis, date_str)
    
    fasta_path   <- file.path(target_dir, fasta_name)
    parquet_path <- file.path(target_dir, parquet_name)
    report_path  <- file.path(target_dir, report_name)
    
    tryCatch({
      withProgress(message = sprintf("Exporting %s dataset...", dis), value = 0.1, {
        incProgress(0.2, detail = "Writing FASTA library...")
        export_fasta_worker(fasta_path, disease = dis, min_p = min_p, min_r = min_r)
        
        incProgress(0.4, detail = "Writing Parquet table...")
        export_parquet_worker(parquet_path, disease = dis, min_p = min_p, min_r = min_r)
        
        incProgress(0.3, detail = "Writing audit summary report...")
        export_report_worker(report_path, disease = dis, min_p = min_p, min_r = min_r)
      })
      
      showModal(modalDialog(
        title = div(class = "d-flex align-items-center gap-2 text-success",
                    icon("circle-check"), tags$span("Export Complete", class = "fw-bold")),
        size = "m",
        easyClose = TRUE,
        footer = modalButton("Close"),
        div(
          class = "p-2",
          tags$p("Successfully saved 3 uncompressed files directly to:"),
          tags$pre(class = "p-2 bg-light rounded text-break", style = "font-size: 0.82rem;", target_dir),
          tags$ul(
            class = "list-unstyled small mb-3",
            tags$li(icon("check", class = "text-success me-2"), fasta_name),
            tags$li(icon("check", class = "text-success me-2"), parquet_name),
            tags$li(icon("check", class = "text-success me-2"), report_name)
          ),
          actionButton("btn_open_export_dir", "Open Folder in Explorer",
                       icon = icon("folder-open"), class = "btn btn-outline-primary btn-sm")
        )
      ))
    }, error = function(e) {
      showNotification(sprintf("Export failed: %s", e$message), type = "error", duration = 10)
    })
  })

  open_folder_in_os <- function(path) {
    if (.Platform$OS.type == "windows") {
      shell.exec(path)
    } else if (Sys.info()["sysname"] == "Darwin") {
      system2("open", shQuote(path))
    } else {
      system2("xdg-open", shQuote(path))
    }
  }

  observeEvent(input$btn_open_export_dir, {
    d <- export_target_folder()
    if (dir.exists(d)) {
      open_folder_in_os(d)
    }
  })

  output$btn_intro_export_fasta <- downloadHandler(
    filename = function() {
      dis <- selected_intro_disease()
      sprintf("OAS_%s_%s.fasta", dis, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      export_fasta_worker(file)
    },
    contentType = "text/plain"
  )

  output$btn_intro_export_parquet <- downloadHandler(
    filename = function() {
      dis <- selected_intro_disease()
      sprintf("OAS_%s_%s.parquet", dis, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      export_parquet_worker(file)
    },
    contentType = "application/octet-stream"
  )

  output$btn_intro_export_report <- downloadHandler(
    filename = function() {
      dis <- selected_intro_disease()
      sprintf("OAS_%s_%s.txt", dis, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      export_report_worker(file)
    },
    contentType = "text/plain"
  )

  # Backward compatible handlers for standalone query
  output$btn_download_fasta <- downloadHandler(
    filename = function() {
      dis <- if (!is.null(input$sel_disease) && nzchar(input$sel_disease)) input$sel_disease else "query"
      sprintf("OAS_%s_%s.fasta", dis, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) export_fasta_worker(file),
    contentType = "text/plain"
  )
  output$btn_download_parquet <- downloadHandler(
    filename = function() {
      dis <- if (!is.null(input$sel_disease) && nzchar(input$sel_disease)) input$sel_disease else "query"
      sprintf("OAS_%s_%s.parquet", dis, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) export_parquet_worker(file),
    contentType = "application/octet-stream"
  )
  output$btn_download_report <- downloadHandler(
    filename = function() {
      dis <- if (!is.null(input$sel_disease) && nzchar(input$sel_disease)) input$sel_disease else "query"
      sprintf("OAS_%s_%s.txt", dis, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) export_report_worker(file),
    contentType = "text/plain"
  )

  # ----------------------------------------------------------------------------
  # Reverse Lookup Logic (Exact Tryptic Peptide Fragment Matching)
  # ----------------------------------------------------------------------------
  reverse_results <- reactiveVal(data.frame())
  
  observeEvent(input$btn_run_lookup, {
    req(input$txt_reverse_pep)
    raw_input <- input$txt_reverse_pep
    
    # Split by newline, comma, semicolon, space
    peptides <- unlist(strsplit(raw_input, "[\r\n,; ]+"))
    peptides <- toupper(trimws(peptides))
    peptides <- peptides[nzchar(peptides)]
    
    if (length(peptides) == 0) return()
    
    cur_db <- db_root()
    if (!nzchar(cur_db) || !dir.exists(cur_db)) {
      showNotification("Please load the CDR3 database directory first.", type = "warning", duration = 4)
      open_db_modal()
      return()
    }
    
    withProgress(message = "Searching database for identified peptide fragment(s)...", value = 0.5, {
      # Exact match with tryptic peptide fragment in semicolon-separated list
      # Escape single quotes to prevent SQL syntax errors / malformed queries
      peptides_clean <- gsub("'", "''", peptides)
      conds <- sprintf("list_contains(string_split(tryptic_peptides, ';'), '%s')", peptides_clean)
      where_clause <- paste(conds, collapse = " OR ")
      
      lookup_sql <- sprintf("
        SELECT 
          clonotype_id,
          Disease,
          BSource,
          BType,
          Isotype,
          cdr3_aa,
          cdr3_len,
          tryptic_peptides,
          N_Patients,
          Redundancy,
          sequence_alignment_aa AS Antibody_sequence,
          Data_file
        FROM read_parquet('%s/**/*.parquet', hive_partitioning=true)
        WHERE %s
        ORDER BY N_Patients DESC, Redundancy DESC
        LIMIT 500;
      ", cur_db, where_clause)
      
      res <- tryCatch({
        dbGetQuery(con, lookup_sql)
      }, error = function(e) {
        data.frame()
      })
      
      reverse_results(res)
    })
  })

  observeEvent(input$btn_clear_lookup, {
    updateTextAreaInput(session, "txt_reverse_pep", value = "")
    reverse_results(data.frame())
  })

  output$lookup_results_count_badge <- renderUI({
    df <- reverse_results()
    if (nrow(df) == 0) return(NULL)
    span(
      class = "badge",
      style = "background: rgba(0, 113, 227, 0.10); color: #0071e3; border: 1px solid rgba(0, 113, 227, 0.22); font-size: 0.70rem; font-weight: 600; padding: 2px 8px; border-radius: 12px; font-family: ui-monospace, SFMono-Regular, monospace;",
      sprintf("%d clonotypes found", nrow(df))
    )
  })
  
  output$table_reverse_results <- renderDT({
    df <- reverse_results()
    if (nrow(df) == 0) {
      return(datatable(data.frame(Message = "No matching disease CDR3 peptides found for specified peptide(s).")))
    }
    
    datatable(
      df,
      filter = "top",
      selection = "single",
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        columnDefs = list(
          list(
            targets = which(names(df) == "Antibody_sequence") - 1,
            render = JS(
              "function(data, type, row) {",
              "  if (type === 'display' && data && data.length > 35) {",
              "    return '<span title=\"' + data + '\" style=\"font-family: monospace; font-size: 0.78rem;\">' + data.substr(0, 32) + '...</span>';",
              "  }",
              "  return '<span style=\"font-family: monospace; font-size: 0.78rem;\">' + data + '</span>';",
              "}"
            )
          )
        )
      ),
      class = "display compact stripe"
    )
  })
}

# ------------------------------------------------------------------------------
# Launch Application
# ------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
