# ==============================================================================
# Unified Master-Detail Alluvial Flow Canvas SVG Generator for OASpepDB
# Stages:
# - Left Stage (Macro): OASpepDB Root -> 4 Clinical Domains -> 25 Disease Cohorts
# - Dynamic Visual Bridge: Connects selected disease row smoothly to Micro Stage
# - Right Stage (Micro): Selected Disease -> Isotype (Method 1) -> BSource -> BType
# ==============================================================================

suppressPackageStartupMessages(library(dplyr))

generate_unified_alluvial_svg <- function(selected_dis = "", df_micro = NULL, W = 1600, H = 570) {
  xml_esc <- function(s) {
    s <- gsub("&", "&amp;", s, fixed = TRUE)
    s <- gsub("<", "&lt;", s, fixed = TRUE)
    s <- gsub(">", "&gt;", s, fixed = TRUE)
    s <- gsub("\"", "&quot;", s, fixed = TRUE)
    s
  }
  
  fmt_n <- function(n) {
    ifelse(n >= 1e6, sprintf("%.1fM", n / 1e6),
    ifelse(n >= 1e3, sprintf("%.0fK", n / 1e3), as.character(n)))
  }
  
  slugify <- function(x) {
    s <- gsub("[^A-Za-z0-9_-]", "-", x)
    gsub("-+", "-", s)
  }
  
  tint_hex <- function(hex, factor = 0.5) {
    rgb_val <- col2rgb(hex) / 255
    new_rgb <- rgb_val + (1 - rgb_val) * (1 - factor)
    rgb(new_rgb[1], new_rgb[2], new_rgb[3])
  }
  
  calc_linear_heights <- function(counts, target_total = 380, h_small = 10, min_h_large = 16, max_col_h = 480) {
    k <- length(counts)
    if (k == 0) return(numeric(0))
    if (k == 1) return(min(target_total * 0.55, 180))
    
    is_small <- counts < 1000
    k_small  <- sum(is_small)
    k_large  <- k - k_small
    
    h <- numeric(k)
    
    if (k_large == 0) {
      h[] <- h_small
    } else {
      h[is_small] <- h_small
      total_small_h <- k_small * h_small
      rem_budget <- max(target_total - total_small_h, k_large * min_h_large)
      
      n_large <- counts[!is_small]
      prop <- n_large / sum(n_large)
      raw_h <- prop * rem_budget
      raw_h <- pmax(raw_h, min_h_large)
      h[!is_small] <- raw_h
    }
    
    # Ensure total height with minimum gaps (8px) does not overflow max_col_h
    gap_budget <- (k - 1) * 8
    if (sum(h) + gap_budget > max_col_h) {
      scale_fac <- (max_col_h - gap_budget) / sum(h)
      h <- pmax(h * scale_fac, if (k_large > 0) 6 else h_small)
    }
    
    round(h, 1)
  }
  
  layout_column_bars <- function(df_col, y_start = 32, y_end = 544, min_gap = 6, max_gap = 22) {
    k <- nrow(df_col)
    if (k == 0) return(df_col)
    total_span <- y_end - y_start
    sum_h <- sum(df_col$h)
    
    if (k == 1) {
      df_col$y_top <- y_start + (total_span - sum_h) / 2
      df_col$y_bot <- df_col$y_top + sum_h
    } else {
      raw_gap <- (total_span - sum_h) / (k - 1)
      gap <- max(min(raw_gap, max_gap), min_gap)
      eff_total <- sum_h + (k - 1) * gap
      curr_y <- y_start + max((total_span - eff_total) / 2, 0)
      
      df_col$y_top <- 0.0
      df_col$y_bot <- 0.0
      for (i in seq_len(k)) {
        df_col$y_top[i] <- curr_y
        df_col$y_bot[i] <- curr_y + df_col$h[i]
        curr_y <- curr_y + df_col$h[i] + gap
      }
    }
    df_col
  }
  
  make_ribbon <- function(xa, ya_top, ya_bot, xb, yb_top, yb_bot, col, alpha = 0.45, title = "", extra_cls = "") {
    dx <- (xb - xa) * 0.50
    path_d <- sprintf(
      "M %.1f,%.1f C %.1f,%.1f %.1f,%.1f %.1f,%.1f L %.1f,%.1f C %.1f,%.1f %.1f,%.1f %.1f,%.1f Z",
      xa, ya_top,
      xa + dx, ya_top,
      xb - dx, yb_top,
      xb, yb_top,
      xb, yb_bot,
      xb - dx, yb_bot,
      xa + dx, ya_bot,
      xa, ya_bot
    )
    tip_attr <- if (nzchar(title)) sprintf('data-glass-title="%s" data-glass-meta="Repertoire Flow"', xml_esc(title)) else ""
    sprintf('<path d="%s" fill="%s" fill-opacity="%.2f" stroke="none" class="flow-ribbon %s" %s></path>',
            path_d, col, alpha, extra_cls, tip_attr)
  }
  
  # Cohort Data Baseline (Order by Category, then by Cassettes)
  cat_order <- c("Infectious Diseases", "Allergy & Airways", "Autoimmune & Neuro", "Hematology & Tumors")
  cohort_data <- data.frame(
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
  )
  cohort_data$CatFactor <- factor(cohort_data$Category, levels = cat_order)
  cohort_data <- cohort_data[order(cohort_data$CatFactor, -cohort_data$Cassettes), ]
  rownames(cohort_data) <- NULL
  
  cat_colors <- c(
    "Infectious Diseases" = "#D55E00",
    "Allergy & Airways"   = "#009E73",
    "Autoimmune & Neuro"  = "#0072B2",
    "Hematology & Tumors" = "#CC79A7"
  )
  
  # Macro Coordinates
  x_root <- 55;   w_root <- 14
  x_dom  <- 275;  w_dom  <- 14
  x_dis  <- 525;  w_dis  <- 12
  
  N_dis <- nrow(cohort_data)
  y_dis_start <- 32
  y_dis_end   <- 544
  
  log_cas    <- log10(cohort_data$Cassettes)
  log_thresh <- 5.0
  log_max    <- max(log_cas)
  h_slender  <- 3.0
  h_max      <- 24.0
  
  cohort_data$h_bar <- ifelse(
    log_cas <= log_thresh,
    h_slender,
    round(h_slender + ((log_cas - log_thresh) / (log_max - log_thresh)) * (h_max - h_slender), 1)
  )
  
  total_bar_h <- sum(cohort_data$h_bar)
  total_span  <- y_dis_end - y_dis_start
  const_gap   <- (total_span - total_bar_h) / (N_dis - 1)
  
  cohort_data$y_top    <- 0.0
  cohort_data$y_bot    <- 0.0
  cohort_data$y_center <- 0.0
  
  curr_y <- y_dis_start
  for (i in 1:N_dis) {
    hb <- cohort_data$h_bar[i]
    cohort_data$y_top[i]    <- curr_y
    cohort_data$y_bot[i]    <- curr_y + hb
    cohort_data$y_center[i] <- curr_y + hb / 2
    curr_y <- curr_y + hb + const_gap
  }
  cohort_data$color <- unname(cat_colors[cohort_data$Category])
  
  domain_summary <- aggregate(Cassettes ~ Category, data = cohort_data, sum)
  domain_summary$CatFactor <- factor(domain_summary$Category, levels = cat_order)
  domain_summary <- domain_summary[order(domain_summary$CatFactor), ]
  rownames(domain_summary) <- NULL
  
  domain_info <- list()
  for (i in 1:nrow(domain_summary)) {
    cat_name <- domain_summary$Category[i]
    sub_df   <- cohort_data[cohort_data$Category == cat_name, ]
    y_mid    <- (min(sub_df$y_top) + max(sub_df$y_bot)) / 2
    dom_h    <- sum(sub_df$h_bar)
    
    domain_info[[cat_name]] <- list(
      Category = cat_name,
      Total    = domain_summary$Cassettes[i],
      Color    = unname(cat_colors[cat_name]),
      y_mid    = y_mid,
      Height   = dom_h,
      y_top    = y_mid - dom_h / 2,
      y_bot    = y_mid + dom_h / 2
    )
  }
  
  root_h    <- sum(cohort_data$h_bar)
  root_ymid <- (y_dis_start + y_dis_end) / 2
  root_ytop <- root_ymid - root_h / 2
  root_ybot <- root_ymid + root_h / 2
  
  curr_y <- root_ytop
  for (cat in cat_order) {
    dh <- domain_info[[cat]]$Height
    domain_info[[cat]]$root_slice_top <- curr_y
    domain_info[[cat]]$root_slice_bot <- curr_y + dh
    curr_y <- curr_y + dh
  }
  
  cohort_data$dom_slice_top <- 0.0
  cohort_data$dom_slice_bot <- 0.0
  for (cat in cat_order) {
    sub_df <- cohort_data[cohort_data$Category == cat, ]
    c_y <- domain_info[[cat]]$y_top
    for (k in 1:nrow(sub_df)) {
      dis_name <- sub_df$Disease[k]
      hb <- sub_df$h_bar[k]
      cohort_data$dom_slice_top[cohort_data$Disease == dis_name] <- c_y
      cohort_data$dom_slice_bot[cohort_data$Disease == dis_name] <- c_y + hb
      c_y <- c_y + hb
    }
  }
  
  is_selected <- nzchar(selected_dis) && (selected_dis %in% cohort_data$Disease)
  sel_row <- if (is_selected) cohort_data[cohort_data$Disease == selected_dis, ][1, ] else NULL
  sel_cat <- if (is_selected) sel_row$Category else ""
  sel_col <- if (is_selected) sel_row$color else "#0071e3"
  
  svg <- c(
    sprintf('<svg id="unified-alluvial-svg" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%" height="100%%" style="font-family: \'Inter\', system-ui, \'Segoe UI\', Roboto, sans-serif;">', W, H),
    '  <defs>',
    '    <style>',
    '      text { font-family: "Inter", system-ui, "Segoe UI", Roboto, sans-serif; -webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; }',
    '      .guide-header { font-size: 13px; font-weight: 700; fill: #6e6e73; text-transform: uppercase; letter-spacing: 0.6px; }',
    '      .alluvial-node { cursor: pointer; transition: opacity 0.25s ease, filter 0.2s ease; }',
    '      .alluvial-node rect { rx: 3px; }',
    '      .alluvial-node:hover rect { filter: drop-shadow(0 2px 8px rgba(0,0,0,0.25)); }',
    '      .node-text { font-size: 12px; font-weight: 600; fill: #1d1d1f; transition: fill 0.25s ease, opacity 0.25s ease; }',
    '      .meta-text { font-size: 11px; font-weight: 500; fill: #86868b; }',
    '      .flow-ribbon { transition: fill 0.25s ease, fill-opacity 0.25s ease, opacity 0.25s ease; pointer-events: auto; cursor: pointer; }',
    '      .macro-ribbon { pointer-events: auto; cursor: pointer; }',
    '      .bar-node { rx: 3px; cursor: pointer; transition: filter 0.15s ease; }',
    '      .bar-node:hover { filter: drop-shadow(0 2px 8px rgba(0,0,0,0.25)); }',
    '    </style>',
    '  </defs>',
    sprintf('  <rect width="%d" height="%d" fill="#ffffff" rx="14" />', W, H),
    
    # Column Headers
    sprintf('  <text x="%d" y="20" text-anchor="middle" class="guide-header">Root Database</text>', x_root + w_root/2),
    sprintf('  <text x="%d" y="20" text-anchor="middle" class="guide-header">Clinical Domains</text>', x_dom + w_dom/2),
    sprintf('  <text x="%d" y="20" text-anchor="start"  class="guide-header">25 Disease Cohorts</text>', x_dis)
  )
  
  if (is_selected) {
    svg <- c(svg,
      '  <text x="867" y="20" text-anchor="middle" class="guide-header">Isotype</text>',
      '  <text x="1137" y="20" text-anchor="middle" class="guide-header">B-Cell Source</text>',
      '  <text x="1407" y="20" text-anchor="middle" class="guide-header">B-Cell Subtype</text>'
    )
  }
  
  flow_active_alpha <- 0.45
  
  # Stage 1 Ribbons: Root -> Domains
  svg <- c(svg, '  <!-- Stage 1 Flows: Root -> 4 Domains -->', '  <g id="stage1-ribbons">')
  for (cat in cat_order) {
    info <- domain_info[[cat]]
    cat_slug <- slugify(cat)
    is_active_cat <- is_selected && (cat == sel_cat)
    fill_col <- if (is_selected) (if (is_active_cat) sel_col else "#d1d1d6") else "#d1d1d6"
    alpha <- if (is_selected) (if (is_active_cat) flow_active_alpha else 0.08) else 0.50
    
    ribbon <- make_ribbon(x_root + w_root, info$root_slice_top, info$root_slice_bot,
                          x_dom, info$y_top, info$y_bot,
                          fill_col, alpha = alpha, title = sprintf("Root → %s", cat),
                          extra_cls = sprintf("macro-ribbon ribbon-dom-%s", cat_slug))
    svg <- c(svg, paste0("    ", ribbon))
  }
  svg <- c(svg, '  </g>')
  
  # Stage 2 Ribbons: Domains -> 25 Diseases
  svg <- c(svg, '  <!-- Stage 2 Flows: Domains -> 25 Disease Cohorts -->', '  <g id="stage2-ribbons">')
  for (i in 1:nrow(cohort_data)) {
    row <- cohort_data[i, ]
    d_name <- row$Disease
    d_slug <- slugify(d_name)
    col    <- row$color
    is_cur_sel <- is_selected && (d_name == selected_dis)
    fill_col <- if (is_selected) (if (is_cur_sel) sel_col else "#d1d1d6") else "#d1d1d6"
    alpha <- if (is_selected) (if (is_cur_sel) flow_active_alpha else 0.08) else 0.50
    
    cas_ribbon <- if (is_cur_sel && !is.null(df_micro)) sum(df_micro$n_cassettes) else row$Cassettes
    ribbon <- make_ribbon(x_dom + w_dom, row$dom_slice_top, row$dom_slice_bot,
                          x_dis, row$y_top, row$y_bot,
                          fill_col, alpha = alpha, title = sprintf("%s → %s: %s CDR3 peptides", row$Category, d_name, fmt_n(cas_ribbon)),
                          extra_cls = sprintf("macro-ribbon ribbon-dis-%s", d_slug))
    svg <- c(svg, paste0("    ", ribbon))
  }
  svg <- c(svg, '  </g>')
  
  # Root Node
  root_alpha <- if (is_selected) 0.65 else 1.0
  svg <- c(svg,
    '  <!-- Column 1: Root Node -->',
    sprintf('  <g id="node-root" class="alluvial-node" style="opacity: %.2f;">', root_alpha),
    sprintf('    <rect x="%d" y="%.1f" width="%d" height="%.1f" rx="3" fill="#1d1d1f" />',
            x_root, root_ytop, w_root, root_h),
    sprintf('    <text x="%d" y="%.1f" text-anchor="middle" font-size="12.5px" font-weight="700" fill="#1d1d1f" transform="rotate(-90, %d, %.1f)">OASpepDB (65.5M)</text>',
            x_root - 16, root_ymid, x_root - 16, root_ymid),
    '  </g>'
  )
  
  # Domain Nodes
  svg <- c(svg, '  <!-- Column 2: 4 Clinical Domain Bars -->', '  <g id="nodes-domains">')
  for (cat in cat_order) {
    info <- domain_info[[cat]]
    cat_slug <- slugify(cat)
    is_active_cat <- is_selected && (cat == sel_cat)
    dom_alpha <- if (is_selected) (if (is_active_cat) 1.0 else 0.22) else 1.0
    stroke_style <- 'stroke="none"'
    
    svg <- c(svg,
      sprintf('    <g class="alluvial-node node-dom-%s" style="opacity: %.2f;" data-glass-title="%s" data-glass-meta="%s CDR3 peptides across %d cohorts" data-glass-desc="Clinical Disease Domain">',
              cat_slug, dom_alpha, xml_esc(cat), fmt_n(info$Total), sum(cohort_data$Category == cat)),
      sprintf('      <rect x="%d" y="%.1f" width="%d" height="%.1f" rx="3" fill="%s" %s />',
              x_dom, info$y_top, w_dom, info$Height, info$Color, stroke_style),
      sprintf('      <text x="%d" y="%.1f" text-anchor="end" dominant-baseline="central" font-size="12px" font-weight="700" fill="%s" style="paint-order: stroke fill; stroke: #ffffff; stroke-width: 3.5px; stroke-linejoin: round;">%s <tspan font-size="10.5px" font-weight="500" fill="#6e6e73" stroke="none">(%s)</tspan></text>',
              x_dom - 10, info$y_mid, info$Color, xml_esc(cat), fmt_n(info$Total)),
      '    </g>'
    )
  }
  svg <- c(svg, '  </g>')
  
  # Disease Nodes
  svg <- c(svg, '  <!-- Column 3: 25 Disease Bars -->', '  <g id="nodes-diseases">')
  for (i in 1:nrow(cohort_data)) {
    row <- cohort_data[i, ]
    d_name   <- row$Disease
    d_slug   <- slugify(d_name)
    col      <- row$color
    is_cur_sel <- is_selected && (d_name == selected_dis)
    
    node_alpha  <- if (is_selected) (if (is_cur_sel) 1.0 else 0.20) else 1.0
    text_color  <- if (is_selected) (if (is_cur_sel) "#000000" else "#8e8e93") else "#1d1d1f"
    text_weight <- if (is_cur_sel) "700" else "600"
    rect_stroke <- if (is_cur_sel) sprintf('stroke="%s" stroke-width="1.8" filter="drop-shadow(0 2px 8px rgba(0,0,0,0.25))"', col) else 'stroke="none"'
    
    # CDR3 peptide count: dynamic update for selected disease when filtered micro data is present
    cas_display <- if (is_cur_sel && !is.null(df_micro)) sum(df_micro$n_cassettes) else row$Cassettes
    
    if (is_cur_sel) {
      svg <- c(svg,
        sprintf('    <!-- Highlight Capsule with Bold Border for Selected Disease -->'),
        sprintf('    <rect x="520" y="%.1f" width="185" height="%.1f" rx="6" fill="%s" fill-opacity="0.12" stroke="%s" stroke-width="2.2" />',
                row$y_top - 3, row$h_bar + 6, col, col)
      )
    }
    
    svg <- c(svg,
      sprintf('    <g class="alluvial-node node-dis-%s" style="opacity: %.2f;" onclick="Shiny.setInputValue(\'btn_cohort_click\', \'%s\', {priority: \'event\'});" data-glass-title="%s" data-glass-meta="%s CDR3 peptides%s (%s)" data-glass-desc="Click cohort to expand single-cell resolution">',
              d_slug, node_alpha, xml_esc(d_name), xml_esc(d_name), fmt_n(cas_display), if (is_cur_sel) " [filtered]" else "", xml_esc(row$Category)),
      sprintf('      <rect x="%d" y="%.1f" width="%d" height="%.1f" rx="3" fill="%s" %s />',
              x_dis, row$y_top, w_dis, row$h_bar, col, rect_stroke),
      sprintf('      <text class="node-text" x="%d" y="%.1f" dominant-baseline="central" style="fill: %s; font-weight: %s;">%s <tspan class="meta-text">(%s)</tspan></text>',
              x_dis + w_dis + 6, row$y_center, text_color, text_weight, xml_esc(d_name), fmt_n(cas_display)),
      '    </g>'
    )
  }
  svg <- c(svg, '  </g>')
  
  # Transition Bridge & Right Stage
  if (!is_selected) {
    svg <- c(svg,
      '  <!-- Idle State: Liquid Glass Placeholder -->',
      '  <line x1="750" y1="30" x2="750" y2="545" stroke="#e5e5ea" stroke-width="1" stroke-dasharray="4,4" />',
      '  <rect x="780" y="32" width="780" height="512" rx="16" fill="#fcfcfd" stroke="#e5e5ea" stroke-dasharray="4,4" />',
      '  <g transform="translate(1170, 240)">',
      '    <circle cx="0" cy="0" r="38" fill="rgba(0,113,227,0.08)" />',
      '    <path d="M-12,-10 L12,-10 L6,10 L-6,10 Z" fill="none" stroke="#0071e3" stroke-width="2.5" stroke-linejoin="round" />',
      '    <circle cx="-12" cy="-10" r="4" fill="#0071e3" />',
      '    <circle cx="12" cy="-10" r="4" fill="#0071e3" />',
      '    <circle cx="0" cy="12" r="4" fill="#0071e3" />',
      '    <text x="0" y="65" text-anchor="middle" font-size="16px" font-weight="700" fill="#1d1d1f">Select a Disease Cohort to Deep-Dive</text>',
      '    <text x="0" y="90" text-anchor="middle" font-size="12px" fill="#86868b">Click any of the 25 disease bars on the left to expand single-cell resolution</text>',
      '    <text x="0" y="110" text-anchor="middle" font-size="12px" fill="#86868b">across Isotype distribution, B-cell tissue source, and cell subtypes.</text>',
      '  </g>'
    )
  } else {
    x_iso <- 860;  w_iso <- 14
    x_src <- 1130; w_src <- 14
    x_typ <- 1400; w_typ <- 14
    
    usable_H <- y_dis_end - y_dis_start
    dis_h <- usable_H * 0.50
    
    c_dis <- sel_col
    c_iso <- tint_hex(sel_col, 0.82)
    c_src <- tint_hex(sel_col, 0.58)
    c_typ <- tint_hex(sel_col, 0.38)
    
    if (is.null(df_micro) || nrow(df_micro) == 0 || sum(df_micro$n_cassettes) == 0) {
      svg <- c(svg, sprintf('  <text x="1170" y="280" text-anchor="middle" font-size="13px" font-weight="500" fill="#86868b">No CDR3 peptides match current filter criteria for %s.</text>', xml_esc(selected_dis)))
    } else {
      total_cassettes <- sum(df_micro$n_cassettes)
      iso_summary <- df_micro %>% group_by(Isotype) %>% summarize(n = sum(n_cassettes), .groups = "drop") %>% arrange(desc(n))
      src_summary <- df_micro %>% group_by(BSource) %>% summarize(n = sum(n_cassettes), .groups = "drop") %>% arrange(desc(n))
      typ_summary <- df_micro %>% group_by(BType) %>% summarize(n = sum(n_cassettes), .groups = "drop") %>% arrange(desc(n))
      
      # Calculate linear heights proportional to N (with shared small height for N < 1000)
      iso_df <- iso_summary
      iso_df$h <- calc_linear_heights(iso_df$n)
      iso_df <- layout_column_bars(iso_df, y_dis_start, y_dis_end)
      
      src_df <- src_summary
      src_df$h <- calc_linear_heights(src_df$n)
      src_df <- layout_column_bars(src_df, y_dis_start, y_dis_end)
      
      typ_df <- typ_summary
      typ_df$h <- calc_linear_heights(typ_df$n)
      typ_df <- layout_column_bars(typ_df, y_dis_start, y_dis_end)
      
      n_iso <- nrow(iso_df)
      n_src <- nrow(src_df)
      n_typ <- nrow(typ_df)
      
      # Links: Isotype -> BSource
      src_order_map <- setNames(seq_len(n_src), src_df$BSource)
      iso_order_map <- setNames(seq_len(n_iso), iso_df$Isotype)
      
      link_iso_src <- df_micro %>%
        group_by(Isotype, BSource) %>%
        summarize(val = sum(n_cassettes), .groups = "drop") %>%
        mutate(src_rank = src_order_map[BSource], iso_rank = iso_order_map[Isotype]) %>%
        arrange(iso_rank, src_rank)
      
      # Links: BSource -> BType
      typ_order_map <- setNames(seq_len(n_typ), typ_df$BType)
      link_src_typ <- df_micro %>%
        group_by(BSource, BType) %>%
        summarize(val = sum(n_cassettes), .groups = "drop") %>%
        mutate(src_rank = src_order_map[BSource], typ_rank = typ_order_map[BType]) %>%
        arrange(src_rank, typ_rank)
      
      # Disease slices at x=705 fanning directly into Isotype
      H_dis <- (sel_row$y_bot + 1) - (sel_row$y_top - 1)
      iso_df$dis_slice_top <- 0.0
      iso_df$dis_slice_bot <- 0.0
      c_y <- sel_row$y_top - 1
      for (i in seq_len(nrow(iso_df))) {
        slice_h <- (iso_df$n[i] / total_cassettes) * H_dis
        iso_df$dis_slice_top[i] <- c_y
        iso_df$dis_slice_bot[i] <- c_y + slice_h
        c_y <- c_y + slice_h
      }
      
      # Proportional slice heights for Isotype -> BSource ribbons
      link_iso_src$iso_slice_top <- 0.0
      link_iso_src$iso_slice_bot <- 0.0
      link_iso_src$src_slice_top <- 0.0
      link_iso_src$src_slice_bot <- 0.0
      
      for (iso_name in iso_df$Isotype) {
        idx <- which(link_iso_src$Isotype == iso_name)
        if (length(idx) > 0) {
          bar_y <- iso_df$y_top[iso_df$Isotype == iso_name]
          bar_h <- iso_df$h[iso_df$Isotype == iso_name]
          vals <- link_iso_src$val[idx]
          slice_heights <- (vals / sum(vals)) * bar_h
          c_y <- bar_y
          for (j in seq_along(idx)) {
            link_iso_src$iso_slice_top[idx[j]] <- c_y
            link_iso_src$iso_slice_bot[idx[j]] <- c_y + slice_heights[j]
            c_y <- c_y + slice_heights[j]
          }
        }
      }
      
      for (src_name in src_df$BSource) {
        idx <- which(link_iso_src$BSource == src_name)
        if (length(idx) > 0) {
          idx_sorted <- idx[order(link_iso_src$iso_rank[idx])]
          bar_y <- src_df$y_top[src_df$BSource == src_name]
          bar_h <- src_df$h[src_df$BSource == src_name]
          vals <- link_iso_src$val[idx_sorted]
          slice_heights <- (vals / sum(vals)) * bar_h
          c_y <- bar_y
          for (j in seq_along(idx_sorted)) {
            link_iso_src$src_slice_top[idx_sorted[j]] <- c_y
            link_iso_src$src_slice_bot[idx_sorted[j]] <- c_y + slice_heights[j]
            c_y <- c_y + slice_heights[j]
          }
        }
      }
      
      # Proportional slice heights for BSource -> BType ribbons
      link_src_typ$src_slice_top <- 0.0
      link_src_typ$src_slice_bot <- 0.0
      link_src_typ$typ_slice_top <- 0.0
      link_src_typ$typ_slice_bot <- 0.0
      
      for (src_name in src_df$BSource) {
        idx <- which(link_src_typ$BSource == src_name)
        if (length(idx) > 0) {
          idx_sorted <- idx[order(link_src_typ$typ_rank[idx])]
          bar_y <- src_df$y_top[src_df$BSource == src_name]
          bar_h <- src_df$h[src_df$BSource == src_name]
          vals <- link_src_typ$val[idx_sorted]
          slice_heights <- (vals / sum(vals)) * bar_h
          c_y <- bar_y
          for (j in seq_along(idx_sorted)) {
            link_src_typ$src_slice_top[idx_sorted[j]] <- c_y
            link_src_typ$src_slice_bot[idx_sorted[j]] <- c_y + slice_heights[j]
            c_y <- c_y + slice_heights[j]
          }
        }
      }
      
      for (typ_name in typ_df$BType) {
        idx <- which(link_src_typ$BType == typ_name)
        if (length(idx) > 0) {
          idx_sorted <- idx[order(link_src_typ$src_rank[idx])]
          bar_y <- typ_df$y_top[typ_df$BType == typ_name]
          bar_h <- typ_df$h[typ_df$BType == typ_name]
          vals <- link_src_typ$val[idx_sorted]
          slice_heights <- (vals / sum(vals)) * bar_h
          c_y <- bar_y
          for (j in seq_along(idx_sorted)) {
            link_src_typ$typ_slice_top[idx_sorted[j]] <- c_y
            link_src_typ$typ_slice_bot[idx_sorted[j]] <- c_y + slice_heights[j]
            c_y <- c_y + slice_heights[j]
          }
        }
      }
      
      # Ribbons: Selected Disease (x=705) -> Isotype (x_iso)
      svg <- c(svg, '  <!-- Micro Stage 1 Ribbons: Selected Disease -> Isotype -->')
      for (i in seq_len(nrow(iso_df))) {
        title_str <- sprintf("%s → %s: %s CDR3 peptides", selected_dis, iso_df$Isotype[i], fmt_n(iso_df$n[i]))
        ribbon <- make_ribbon(705, iso_df$dis_slice_top[i], iso_df$dis_slice_bot[i],
                              x_iso, iso_df$y_top[i], iso_df$y_bot[i],
                              c_iso, alpha = flow_active_alpha, title = title_str)
        svg <- c(svg, paste0("  ", ribbon))
      }
      
      # Ribbons: Isotype -> BSource
      svg <- c(svg, '  <!-- Micro Stage 2 Ribbons: Isotype -> BSource -->')
      for (i in seq_len(nrow(link_iso_src))) {
        r <- link_iso_src[i, ]
        title_str <- sprintf("%s → %s: %s CDR3 peptides", r$Isotype, r$BSource, fmt_n(r$val))
        ribbon <- make_ribbon(x_iso + w_iso, r$iso_slice_top, r$iso_slice_bot,
                              x_src, r$src_slice_top, r$src_slice_bot,
                              c_src, alpha = 0.45, title = title_str)
        svg <- c(svg, paste0("  ", ribbon))
      }
      
      # Ribbons: BSource -> BType
      svg <- c(svg, '  <!-- Micro Stage 3 Ribbons: BSource -> BType -->')
      for (i in seq_len(nrow(link_src_typ))) {
        r <- link_src_typ[i, ]
        title_str <- sprintf("%s → %s: %s CDR3 peptides", r$BSource, r$BType, fmt_n(r$val))
        ribbon <- make_ribbon(x_src + w_src, r$src_slice_top, r$src_slice_bot,
                              x_typ, r$typ_slice_top, r$typ_slice_bot,
                              c_typ, alpha = 0.40, title = title_str)
        svg <- c(svg, paste0("  ", ribbon))
      }
      
      # Col 4: Isotypes (WITH rx="3")
      for (i in seq_len(nrow(iso_df))) {
        r <- iso_df[i, ]
        svg <- c(svg,
          sprintf('  <rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="3" fill="%s" class="bar-node" data-glass-title="%s" data-glass-meta="%s CDR3 peptides (%.1f%% of cohort)" data-glass-desc="Antibody Isotype"></rect>',
                  x_iso, r$y_top, w_iso, r$h, c_iso, xml_esc(r$Isotype), fmt_n(r$n), (r$n/total_cassettes)*100)
        )
        if (r$h >= 22) {
          svg <- c(svg,
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="node-text" dominant-baseline="middle">%s</text>',
                    x_iso + w_iso + 6, r$y_top + r$h/2 - 6, xml_esc(r$Isotype)),
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="meta-text" dominant-baseline="middle">%s</text>',
                    x_iso + w_iso + 6, r$y_top + r$h/2 + 6, fmt_n(r$n))
          )
        } else {
          svg <- c(svg,
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="node-text" dominant-baseline="middle">%s <tspan class="meta-text">(%s)</tspan></text>',
                    x_iso + w_iso + 6, r$y_top + r$h/2, xml_esc(r$Isotype), fmt_n(r$n))
          )
        }
      }
      
      # Col 5: B-Cell Source (WITH rx="3")
      for (i in seq_len(nrow(src_df))) {
        r <- src_df[i, ]
        svg <- c(svg,
          sprintf('  <rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="3" fill="%s" class="bar-node" data-glass-title="%s" data-glass-meta="%s CDR3 peptides (%.1f%% of cohort)" data-glass-desc="B-Cell Tissue Source"></rect>',
                  x_src, r$y_top, w_src, r$h, c_src, xml_esc(r$BSource), fmt_n(r$n), (r$n/total_cassettes)*100)
        )
        if (r$h >= 22) {
          svg <- c(svg,
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="node-text" dominant-baseline="middle">%s</text>',
                    x_src + w_src + 6, r$y_top + r$h/2 - 6, xml_esc(r$BSource)),
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="meta-text" dominant-baseline="middle">%s</text>',
                    x_src + w_src + 6, r$y_top + r$h/2 + 6, fmt_n(r$n))
          )
        } else {
          svg <- c(svg,
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="node-text" dominant-baseline="middle">%s <tspan class="meta-text">(%s)</tspan></text>',
                    x_src + w_src + 6, r$y_top + r$h/2, xml_esc(r$BSource), fmt_n(r$n))
          )
        }
      }
      
      # Col 6: B-Cell Subtype (WITH rx="3")
      for (i in seq_len(nrow(typ_df))) {
        r <- typ_df[i, ]
        svg <- c(svg,
          sprintf('  <rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="3" fill="%s" class="bar-node" data-glass-title="%s" data-glass-meta="%s CDR3 peptides (%.1f%% of cohort)" data-glass-desc="B-Cell Subtype"></rect>',
                  x_typ, r$y_top, w_typ, r$h, c_typ, xml_esc(r$BType), fmt_n(r$n), (r$n/total_cassettes)*100)
        )
        if (r$h >= 22) {
          svg <- c(svg,
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="node-text" dominant-baseline="middle">%s</text>',
                    x_typ + w_typ + 6, r$y_top + r$h/2 - 6, xml_esc(r$BType)),
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="meta-text" dominant-baseline="middle">%s</text>',
                    x_typ + w_typ + 6, r$y_top + r$h/2 + 6, fmt_n(r$n))
          )
        } else {
          svg <- c(svg,
            sprintf('  <text x="%.1f" y="%.1f" text-anchor="start" class="node-text" dominant-baseline="middle">%s <tspan class="meta-text">(%s)</tspan></text>',
                    x_typ + w_typ + 6, r$y_top + r$h/2, xml_esc(r$BType), fmt_n(r$n))
          )
        }
      }
    }
  }
  
  svg <- c(svg, '</svg>')
  paste(svg, collapse = "\n")
}
