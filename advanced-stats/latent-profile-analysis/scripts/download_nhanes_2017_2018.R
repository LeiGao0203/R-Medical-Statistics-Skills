#!/usr/bin/env Rscript

# Download only public NHANES 2017-2018 XPT files used by the LPA example.
# Usage from repository root:
#   Rscript advanced-stats/latent-profile-analysis/scripts/download_nhanes_2017_2018.R

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else file.path("example", "nhanes-lpa", "data", "raw")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles"
files <- c("DEMO_J", "BMX_J", "BPX_J", "GHB_J", "TCHOL_J", "HDL_J")
for (file_stem in files) {
  destination <- file.path(output_dir, paste0(file_stem, ".XPT"))
  url <- paste0(base_url, "/", file_stem, ".XPT")
  message("Downloading ", url)
  utils::download.file(url, destination, mode = "wb", quiet = FALSE)
}
message("Downloaded ", length(files), " public NHANES files to ", normalizePath(output_dir))
