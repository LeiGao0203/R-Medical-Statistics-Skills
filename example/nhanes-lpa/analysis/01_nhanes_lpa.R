#!/usr/bin/env Rscript

# Reproducible NHANES 2017-2018 LPA example.
# Run from the repository root:
#   Rscript example/nhanes-lpa/analysis/01_nhanes_lpa.R

required_packages <- c("dplyr", "ggplot2", "haven", "mclust", "readr", "tibble", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(dplyr)
  library(mclust)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else "example/nhanes-lpa/analysis/01_nhanes_lpa.R"
project_root <- normalizePath(file.path(dirname(script_path), "../../.."), mustWork = TRUE)

raw_dir <- file.path(project_root, "example", "nhanes-lpa", "data", "raw")
derived_dir <- file.path(project_root, "example", "nhanes-lpa", "data", "derived")
results_dir <- file.path(project_root, "example", "nhanes-lpa", "results")
figures_dir <- file.path(project_root, "example", "nhanes-lpa", "figures")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

raw_files <- c("DEMO_J", "BMX_J", "BPX_J", "GHB_J", "TCHOL_J", "HDL_J")
raw_paths <- file.path(raw_dir, paste0(raw_files, ".XPT"))
if (!all(file.exists(raw_paths))) {
  stop("Missing raw files. Run advanced-stats/latent-profile-analysis/scripts/download_nhanes_2017_2018.R first.")
}

read_nhanes <- function(stem) {
  haven::read_xpt(file.path(raw_dir, paste0(stem, ".XPT")), .name_repair = "minimal")
}

demo <- read_nhanes("DEMO_J")
bmx <- read_nhanes("BMX_J")
bpx <- read_nhanes("BPX_J")
ghb <- read_nhanes("GHB_J")
tchol <- read_nhanes("TCHOL_J")
hdl <- read_nhanes("HDL_J")

bp_columns <- paste0("BPXSY", 1:4)
bpx_clean <- bpx %>%
  mutate(
    sbp_n = rowSums(!is.na(across(all_of(bp_columns)))),
    sbp_mean = if_else(
      sbp_n > 0,
      rowMeans(across(all_of(bp_columns)), na.rm = TRUE),
      NA_real_
    )
  ) %>%
  select(SEQN, sbp_mean)

nhanes_merged <- demo %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3, WTMEC2YR, SDMVPSU, SDMVSTRA) %>%
  left_join(bmx %>% select(SEQN, BMXBMI, BMXWAIST), by = "SEQN") %>%
  left_join(bpx_clean, by = "SEQN") %>%
  left_join(ghb %>% select(SEQN, LBXGH), by = "SEQN") %>%
  left_join(tchol %>% select(SEQN, LBXTC), by = "SEQN") %>%
  left_join(hdl %>% select(SEQN, LBDHDD), by = "SEQN") %>%
  mutate(
    age_years = as.numeric(RIDAGEYR),
    bmi = as.numeric(BMXBMI),
    waist_cm = as.numeric(BMXWAIST),
    systolic_bp = as.numeric(sbp_mean),
    hba1c = as.numeric(LBXGH),
    total_cholesterol = as.numeric(LBXTC),
    hdl_c = as.numeric(LBDHDD)
  )

indicators <- c("bmi", "waist_cm", "systolic_bp", "hba1c", "total_cholesterol", "hdl_c")

# Higher values indicate greater cardiometabolic burden after reversing HDL-C.
analysis_data <- nhanes_merged %>%
  filter(age_years >= 20) %>%
  mutate(hdl_c_risk = -hdl_c) %>%
  select(SEQN, age_years, RIAGENDR, RIDRETH3, WTMEC2YR, SDMVPSU, SDMVSTRA,
         all_of(indicators), hdl_c_risk) %>%
  filter(if_all(c(setdiff(indicators, "hdl_c"), "hdl_c_risk"), ~ is.finite(.x))) %>%
  mutate(
    complete_indicator_case = if_all(all_of(c(setdiff(indicators, "hdl_c"), "hdl_c_risk")), ~ !is.na(.x))
  ) %>%
  filter(complete_indicator_case) %>%
  select(-complete_indicator_case)

analysis_indicators <- c("bmi", "waist_cm", "systolic_bp", "hba1c", "total_cholesterol", "hdl_c_risk")
if (nrow(analysis_data) < 200) stop("Too few complete adult cases for this example: ", nrow(analysis_data))

readr::write_csv(analysis_data, file.path(derived_dir, "nhanes_lpa_analysis_data.csv"))

centers <- vapply(analysis_data[analysis_indicators], mean, numeric(1), na.rm = TRUE)
scales <- vapply(analysis_data[analysis_indicators], sd, numeric(1), na.rm = TRUE)
z_data <- as.data.frame(scale(analysis_data[analysis_indicators], center = centers, scale = scales))
readr::write_csv(tibble::tibble(indicator = names(centers), center = centers, scale = scales),
                 file.path(results_dir, "standardization_parameters.csv"))

candidate_models <- c("EEI", "VVI", "EEE", "VVV")
candidate_classes <- 1:5

fit_one <- function(model_name, classes) {
  error_message <- NA_character_
  fit <- tryCatch(
    if (classes == 1) {
      # mclust 6.1.1 has a model-specific one-class edge case for some
      # covariance names; one class is a common covariance baseline.
      mclust::Mclust(z_data, G = 1, verbose = FALSE)
    } else {
      mclust::Mclust(z_data, G = classes, modelNames = model_name, verbose = FALSE)
    },
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(fit)) {
    return(tibble::tibble(model = model_name, classes = classes, bic_mclust = NA_real_,
                          log_likelihood = NA_real_, n_parameters = NA_real_, error = TRUE,
                          error_message = error_message))
  }
  tibble::tibble(model = model_name, classes = classes, bic_mclust = fit$bic,
                 log_likelihood = fit$loglik, n_parameters = fit$df, error = FALSE,
                 error_message = error_message)
}

grid <- tidyr::expand_grid(model = candidate_models, classes = candidate_classes)
selection <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  fit_one(grid$model[[i]], grid$classes[[i]])
})) %>%
  mutate(model_description = dplyr::recode(model,
    EEI = "equal diagonal variance",
    VVI = "variable diagonal variance",
    EEE = "equal covariance matrix",
    VVV = "variable covariance matrix")) %>%
  arrange(desc(bic_mclust))
readr::write_csv(selection, file.path(results_dir, "model_selection.csv"))

# Choose the best BIC candidate, then retain the explicit choice in the output.
best_candidate <- selection %>% filter(!error, is.finite(bic_mclust)) %>% slice_max(bic_mclust, n = 1, with_ties = FALSE)
if (nrow(best_candidate) != 1) stop("No candidate model converged.")
chosen_model <- best_candidate$model[[1]]
chosen_classes <- best_candidate$classes[[1]]
set.seed(20260716)
final_fit <- if (chosen_classes == 1) {
  mclust::Mclust(z_data, G = 1, verbose = FALSE)
} else {
  mclust::Mclust(z_data, G = chosen_classes, modelNames = chosen_model, verbose = FALSE)
}

posterior <- as.data.frame(final_fit$z)
names(posterior) <- paste0("prob_class_", seq_len(ncol(posterior)))
classification <- tibble::tibble(
  SEQN = analysis_data$SEQN,
  assigned_class = final_fit$classification,
  maximum_posterior = apply(posterior, 1, max)
) %>%
  bind_cols(posterior) %>%
  mutate(uncertain_assignment = maximum_posterior < 0.70)
readr::write_csv(classification, file.path(results_dir, "posterior_classification.csv"))

class_sizes <- classification %>%
  count(assigned_class, name = "n") %>%
  mutate(proportion = n / sum(n))
mean_posteriors <- classification %>%
  group_by(assigned_class) %>%
  summarise(mean_assigned_posterior = mean(maximum_posterior), .groups = "drop")
class_sizes <- left_join(class_sizes, mean_posteriors, by = "assigned_class")
readr::write_csv(class_sizes, file.path(results_dir, "class_sizes.csv"))

mean_matrix <- t(as.matrix(final_fit$parameters$mean))
if (ncol(mean_matrix) != length(analysis_indicators)) {
  mean_matrix <- matrix(as.numeric(final_fit$parameters$mean), nrow = chosen_classes, byrow = TRUE)
}
colnames(mean_matrix) <- analysis_indicators
profile_z <- as.data.frame(mean_matrix) %>% mutate(class_id = seq_len(nrow(mean_matrix)))
profile_order <- profile_z %>% mutate(severity = rowMeans(across(all_of(analysis_indicators)))) %>% arrange(severity) %>% pull(class_id)
profile_labels <- setNames(paste0("Profile_", seq_along(profile_order)), profile_order)

profile_means <- profile_z %>%
  tidyr::pivot_longer(all_of(analysis_indicators), names_to = "indicator", values_to = "z_mean") %>%
  mutate(
    original_mean = z_mean * unname(scales[indicator]) + unname(centers[indicator]),
    profile = unname(profile_labels[as.character(class_id)])
  ) %>%
  arrange(match(profile, paste0("Profile_", seq_len(chosen_classes))), indicator)
readr::write_csv(profile_means, file.path(results_dir, "profile_means_original_scale.csv"))

plot_data <- profile_means %>%
  ggplot2::ggplot(ggplot2::aes(x = indicator, y = z_mean, group = profile, color = profile)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey75") +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(title = "NHANES 2017-2018 cardiometabolic profiles", x = NULL,
                y = "Standardized model-estimated mean", color = "Profile") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
ggplot2::ggsave(file.path(figures_dir, "nhanes_lpa_profile_plot.png"), plot_data, width = 9, height = 5.5, dpi = 160)

sink(file.path(results_dir, "session_info.txt"))
cat("Chosen model:", chosen_model, "\nChosen classes:", chosen_classes, "\n")
cat("Complete adult cases:", nrow(analysis_data), "\n")
print(best_candidate)
cat("\nSession information:\n")
print(sessionInfo())
sink()

message("LPA complete: ", nrow(analysis_data), " complete adult cases; selected ", chosen_model,
        " with ", chosen_classes, " profiles. Results written to ", results_dir)
