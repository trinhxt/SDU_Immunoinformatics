# R/helpers_paths.R
suppressPackageStartupMessages({
  library(yaml)
})

# Load config from YAML
load_config <- function(path = "config/paths.yml") {
  if (!file.exists(path)) {
    stop(
      "Config file not found: ", path, "\n",
      "Copy config/paths_example.yml to config/paths.yml and edit it."
    )
  }
  yaml::read_yaml(path)
}

# Normalize paths to absolute forward-slash (DuckDB safe on Windows)
norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

# Join paths using forward slashes (avoids backslashes from file.path on Windows)
pjoin <- function(...) norm(paste(..., sep = "/"))

# Optional: assert that a string contains no backslashes (useful for DuckDB globs)
assert_no_backslash <- function(x, name = deparse(substitute(x))) {
  if (grepl("\\\\", x)) stop("Backslash found in ", name, ":\n", x)
  invisible(TRUE)
}

# Convenience: normalize key paths from cfg (keeps scripts tidy)
get_paths <- function(cfg) {
  # allow missing optional keys
  get <- function(k, default = NULL) if (!is.null(cfg[[k]]) && nzchar(cfg[[k]])) cfg[[k]] else default
  
  list(
    oas_download_script = norm(get("oas_download_script")),
    oas_raw_dir         = {
      v <- get("oas_raw_dir", "")
      if (nzchar(v)) norm(v) else ""
    },
    oas_processed_dir   = norm(get("oas_processed_dir")),
    metadata_csv        = norm(get("metadata_csv")),
    
    work_root           = norm(get("work_root")),
    final_part2_dir      = norm(get("final_part2_dir")),
    
    part3_root          = norm(get("part3_root")),
    part4_root          = norm(get("part4_root")),
    part5_root          = norm(get("part5_root")),
    
    reference_tryptic_rdata = norm(get("reference_tryptic_rdata", "data/reference/UniProtNCBI_Tryptic.RData")),
    
    zenodo_out_root     = norm(get("zenodo_out_root"))
  )
}
