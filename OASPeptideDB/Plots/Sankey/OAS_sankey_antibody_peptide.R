#===============================================================================
# Draw 2 Sankey plots (Antibody + Peptide) with (possibly) different stage order
#
# Antibody input columns:
#   Species, Disease_presence_antibody, Disease, Isotype, BSource, BType, n_unique_antibodies
#
# Peptide input columns (your file):
#   Species, Disease_presence_peptide, Disease, CDR3, Uniqueness, n_unique_peptides
#
# Custom coloring:
#   - Only Disease nodes are colored (shared palette across both plots)
#   - All other nodes are auto-colored
#===============================================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsankey)
library(scales)

#------------------------------------------------------------------------------
# 1) Disease color palette (shared across antibody + peptide)
#------------------------------------------------------------------------------
disease_colors <- c(
  # Low-data diseases (soft / relaxed)
  "Allergic-Rhinitis-Out-Of-Season" = "#C6DBEF",
  "Allergic-Rhinitis-In-Season"     = "#F2C9C5",
  "Non-Dengue-Febrile-Illness"      = "#CDEAE3",
  "CLL"                             = "#D4C6E7",
  "Dengue"                          = "#F4D6A8",
  "Allergy/SIT"                     = "#C7E9C0",
  "CMV"                             = "#BFD7EA",
  "MuSK-MG"                         = "#E8B7B2",
  "Light-Chain-Amyloidosis"         = "#D7CBE8",
  "AChR-MG"                         = "#E4C1F9",
  "HCV"                             = "#F6C1CC",
  "EBV"                             = "#C3E6CB",
  "Asthma"                          = "#B3CDE3",
  "Tonsillitis/Obstructive-Sleep-Apnea" = "#E6C9A8",
  "Obstructive-Sleep-Apnea"         = "#C9E4CA",
  "Tonsillitis"                     = "#FFE0B2",
  
  # High-data diseases (Nature high-contrast)
  "SLE"                             = "#F781BF",
  "MS"                              = "#A65628",
  "Ebola"                           = "#666666",
  "Allergy/NoSIT"                   = "#FF7F00",
  "Healthy/celiac-disease"          = "#984EA3",
  "CMV/EBV"                         = "#4DAF4A",
  "HIV"                             = "#377EB8",
  "SARS-COV-2"                      = "#f8766d"
)

#------------------------------------------------------------------------------
# 2) Helper: order factor levels by decreasing total N_adj (big -> small)
#------------------------------------------------------------------------------
order_levels <- function(dt, col) {
  dt %>%
    group_by(.data[[col]]) %>%
    summarise(N_adj_sum = sum(N_adj), .groups = "drop") %>%
    arrange(desc(N_adj_sum)) %>%
    pull(.data[[col]]) %>%
    as.character()
}

#------------------------------------------------------------------------------
# 3) Flexible helper: build + plot + save sankey with disease-only fixed coloring
#------------------------------------------------------------------------------
plot_sankey <- function(dt,
                        stage_cols,     # vector of columns in the exact order to plot
                        presence_col,   # e.g. "Disease_presence_antibody" or "Disease_presence_peptide"
                        count_col,      # e.g. "n_unique_antibodies" or "n_unique_peptides"
                        title,
                        out_pdf,
                        disease_colors,
                        scale_div = 1000) {
  
  # ---- Validate columns ----
  needed <- unique(c(stage_cols, presence_col, count_col))
  miss <- setdiff(needed, colnames(dt))
  if (length(miss) > 0) stop("Missing columns in input: ", paste(miss, collapse = ", "))
  
  # ---- Add N_adj (replication weight) ----
  dt2 <- dt %>%
    mutate(N_adj = round(.data[[count_col]] / scale_div, 0)) %>%
    filter(N_adj > 0)
  
  if (nrow(dt2) == 0) stop("No rows left after scaling/filtering (all N_adj <= 0). Try smaller scale_div.")
  
  # ---- Ensure all diseases exist in the color map (auto-add fallback if not) ----
  diseases_present <- sort(unique(dt2$Disease))
  missing_colors <- setdiff(diseases_present, names(disease_colors))
  if (length(missing_colors) > 0) {
    fallback <- rep("#999999", length(missing_colors))
    names(fallback) <- missing_colors
    disease_colors <- c(disease_colors, fallback)
    warning("Added fallback grey colors for diseases not in disease_colors: ",
            paste(missing_colors, collapse = ", "))
  }
  
  # ---- Order levels for each stage (Species kept as-is; others ordered by size) ----
  stage_levels_map <- list()
  for (col in stage_cols) {
    if (col == "Species") {
      stage_levels_map[[col]] <- unique(dt2$Species)
    } else {
      stage_levels_map[[col]] <- order_levels(dt2, col)
    }
  }
  
  # ---- Apply factor levels for each stage col ----
  for (col in stage_cols) {
    dt2[[col]] <- factor(dt2[[col]], levels = stage_levels_map[[col]])
  }
  
  # ---- Replicate rows based on N_adj (controls flow thickness) ----
  dt2_rep <- dt2[rep(seq_len(nrow(dt2)), dt2$N_adj), ]
  dt2_rep <- as.data.frame(dt2_rep)
  
  # ---- Convert stage columns to character BEFORE make_long() ----
  for (col in stage_cols) dt2_rep[[col]] <- as.character(dt2_rep[[col]])
  
  # ---- Long format for ggsankey (uses stage_cols in the given order) ----
  long <- dt2_rep %>% make_long(!!!rlang::syms(stage_cols))
  setDT(long)
  
  # ---- Consistent node ordering across all stages ----
  stage_levels <- unlist(stage_levels_map, use.names = FALSE)
  stage_levels <- unique(as.character(stage_levels))
  
  long[, node := factor(node, levels = stage_levels)]
  long[, next_node := factor(next_node, levels = stage_levels)]
  
  # ---- Build palette: Disease fixed, everything else auto ----
  disease_nodes <- as.character(stage_levels_map[["Disease"]])
  all_nodes <- stage_levels
  other_nodes <- setdiff(all_nodes, disease_nodes)
  
  other_palette <- scales::hue_pal(l = 65, c = 100)(length(other_nodes))
  names(other_palette) <- other_nodes
  
  fill_values <- c(other_palette, disease_colors)
  
  # ---- Plot ----
  p <- ggplot(
    long,
    aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = node
    )
  ) +
    geom_sankey(flow.alpha = 0.55) +
    geom_sankey_label(aes(label = node), size = 3, color = "black") +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    theme_sankey() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = title, x = "", y = "")
  
  print(p)
  
  ggsave(
    filename = out_pdf,
    plot     = p,
    device   = cairo_pdf,
    width    = 12,
    height   = 8,
    units    = "in"
  )
  
  invisible(p)
}

#===============================================================================
# 4) Plot 1: Antibody (unchanged stage order)
#===============================================================================
sankey_data_ab <- fread("sankey_data_antibody.csv")

plot_sankey(
  dt         = sankey_data_ab,
  stage_cols = c("Species", "Disease_presence_antibody", "Disease", "Isotype", "BSource", "BType"),
  presence_col = "Disease_presence_antibody",
  count_col    = "n_unique_antibodies",
  title      = "Antibody Sankey Diagram",
  out_pdf    = "sankey_antibody.pdf",
  disease_colors = disease_colors,
  scale_div  = 1000
)

#===============================================================================
# 5) Plot 2: Peptide (your different columns)
#    Choose the order you want:
#      Option A: Species -> presence -> Disease -> Uniqueness -> CDR3
#===============================================================================
sankey_data_pep <- fread("sankey_data_peptide.csv")

plot_sankey(
  dt         = sankey_data_pep,
  stage_cols = c("Species", "Disease_presence_peptide", "Disease", "Uniqueness", "CDR3"),
  presence_col = "Disease_presence_peptide",
  count_col    = "n_unique_peptides",
  title      = "Peptide Sankey Diagram",
  out_pdf    = "sankey_peptide.pdf",
  disease_colors = disease_colors,
  scale_div  = 1000
)
