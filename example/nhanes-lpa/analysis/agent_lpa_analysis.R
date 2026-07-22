#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, warn = 1)

required <- c("haven", "mclust", "ggplot2")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Required packages are missing: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(haven)
  library(mclust)
  library(ggplot2)
})

# Resolve paths from the script location so that the command works from the
# repository root or from any working directory.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "analysis/agent_lpa_analysis.R"
target_dir <- normalizePath(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), ".."), mustWork = TRUE)
raw_dir <- file.path(target_dir, "data", "raw")
derived_dir <- file.path(target_dir, "data", "derived")
results_dir <- file.path(target_dir, "results")
figures_dir <- file.path(target_dir, "figures")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

read_selected_xpt <- function(file, columns) {
  path <- file.path(raw_dir, file)
  if (!file.exists(path)) stop("Missing NHANES file: ", path)
  x <- haven::read_xpt(path)
  if (!all(columns %in% names(x))) {
    stop("Missing required variables in ", file, ": ",
         paste(setdiff(columns, names(x)), collapse = ", "))
  }
  x[, columns, drop = FALSE]
}

format_md <- function(x, digits = 3) {
  y <- x
  y[] <- lapply(y, function(v) {
    if (is.numeric(v)) {
      format(round(v, digits), nsmall = digits, trim = TRUE, scientific = FALSE)
    } else {
      v <- as.character(v)
      v[is.na(v)] <- ""
      v
    }
  })
  lines <- c(
    paste0("| ", paste(names(y), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(y)), collapse = " | "), " |"),
    apply(y, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  )
  paste(lines, collapse = "\n")
}

seed <- 20260716L
set.seed(seed)
indicators <- c("BMXBMI", "BMXWAIST", "SBP_MEAN")
indicator_labels <- c(
  BMXBMI = "Body mass index",
  BMXWAIST = "Waist circumference",
  SBP_MEAN = "Mean systolic blood pressure"
)
indicator_units <- c(BMXBMI = "kg/m^2", BMXWAIST = "cm", SBP_MEAN = "mmHg")

# The six local XPT files include additional laboratory files. This minimal
# task only needs demographics, examination body measures, and blood pressure.
raw_files <- sort(list.files(raw_dir, pattern = "[.]XPT$", full.names = FALSE))
needed_files <- c("DEMO_J.XPT", "BMX_J.XPT", "BPX_J.XPT")
if (!all(needed_files %in% raw_files)) {
  stop("Expected local files are missing: ", paste(setdiff(needed_files, raw_files), collapse = ", "))
}

demo <- read_selected_xpt(
  "DEMO_J.XPT",
  c("SEQN", "RIDSTATR", "RIDAGEYR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
)
bmx <- read_selected_xpt("BMX_J.XPT", c("SEQN", "BMXBMI", "BMXWAIST"))
bpx <- read_selected_xpt("BPX_J.XPT", c("SEQN", paste0("BPXSY", 1:4)))

if (anyDuplicated(demo$SEQN) || anyDuplicated(bmx$SEQN) || anyDuplicated(bpx$SEQN)) {
  stop("SEQN is not unique in at least one input file")
}

sbp_columns <- paste0("BPXSY", 1:4)
bpx$SBP_MEAN <- rowMeans(bpx[, sbp_columns, drop = FALSE], na.rm = TRUE)
bpx$SBP_MEAN[is.nan(bpx$SBP_MEAN)] <- NA_real_
bpx$SBP_N_READINGS <- rowSums(!is.na(bpx[, sbp_columns, drop = FALSE]))

dat <- Reduce(
  function(x, y) merge(x, y, by = "SEQN", all = FALSE, sort = FALSE),
  list(demo, bmx, bpx[, c("SEQN", "SBP_MEAN", "SBP_N_READINGS")])
)

# Adult MEC-examined participants are used because BMI, waist, and BP are
# examination measurements. Complete cases are required by ordinary Mclust.
adult <- dat[!is.na(dat$RIDSTATR) & dat$RIDSTATR == 2 &
             !is.na(dat$RIDAGEYR) & dat$RIDAGEYR >= 20, , drop = FALSE]
complete <- adult[complete.cases(adult[, indicators, drop = FALSE]), , drop = FALSE]

flow <- data.frame(
  stage = c("DEMO rows", "Inner merge of DEMO/BMX/BPX", "MEC examined", "Adults age >=20", "Complete cases"),
  n = c(nrow(demo), nrow(dat), sum(dat$RIDSTATR == 2, na.rm = TRUE), nrow(adult), nrow(complete)),
  rule = c("DEMO_J.XPT", "Merge by SEQN", "RIDSTATR == 2", "RIDAGEYR >= 20", "BMI, waist, and mean SBP all observed")
)

missingness <- data.frame(
  indicator = indicators,
  label = unname(indicator_labels[indicators]),
  unit = unname(indicator_units[indicators]),
  adult_n = nrow(adult),
  missing_n = vapply(adult[indicators], function(x) sum(is.na(x)), integer(1))
)
missingness$missing_pct <- 100 * missingness$missing_n / missingness$adult_n

raw_values <- as.matrix(complete[, indicators, drop = FALSE])
centers <- colMeans(raw_values)
scales <- apply(raw_values, 2, sd)
if (any(!is.finite(scales) | scales <= 0)) stop("An indicator has no positive complete-case SD")
z_values <- scale(raw_values, center = centers, scale = scales)
colnames(z_values) <- paste0("z_", indicators)

analysis_data <- cbind(
  complete[, c("SEQN", "RIDAGEYR", "RIDSTATR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
               "SBP_N_READINGS", indicators)],
  as.data.frame(z_values)
)
write_csv(analysis_data, file.path(derived_dir, "analysis_data.csv"))

# mclust is fit to standardized continuous indicators. WTMEC2YR, SDMVPSU,
# and SDMVSTRA are retained for future survey-aware sensitivity analyses, but
# ordinary Mclust() does not use the NHANES complex sampling design.
model_names <- c("EEI", "VVI", "EEE")
candidate_fits <- list()
candidate_rows <- list()

for (model_name in model_names) {
  for (g in 1:4) {
    key <- paste(g, model_name, sep = "_")
    set.seed(seed + g + 100L * match(model_name, model_names))
    fit <- tryCatch(
      mclust::Mclust(z_values, G = g, modelNames = model_name, verbose = FALSE),
      error = function(e) e
    )
    if (inherits(fit, "error") || is.null(fit)) {
      candidate_rows[[key]] <- data.frame(
        G = g, model = model_name, status = "error", logLik = NA_real_, df = NA_real_,
        AIC = NA_real_, BIC = NA_real_, SABIC = NA_real_, ICL = NA_real_,
        min_class_n = NA_integer_, min_class_prop = NA_real_,
        avg_max_posterior = NA_real_, uncertain_lt_0_70 = NA_real_,
        relative_entropy = NA_real_, quality_ok = FALSE, selected = FALSE,
        error_message = conditionMessage(fit)
      )
      next
    }

    candidate_fits[[key]] <- fit
    posterior <- as.matrix(fit$z)
    class_n <- tabulate(fit$classification, nbins = g)
    max_posterior <- apply(posterior, 1, max)
    entropy <- -sum(posterior * log(pmax(posterior, .Machine$double.eps)))
    relative_entropy <- if (g == 1) 1 else 1 - entropy / (nrow(z_values) * log(g))
    bic <- 2 * fit$loglik - fit$df * log(nrow(z_values))
    sabic <- 2 * fit$loglik - fit$df * log((nrow(z_values) + 2) / 24)

    candidate_rows[[key]] <- data.frame(
      G = g, model = model_name, status = "ok", logLik = fit$loglik, df = fit$df,
      AIC = 2 * fit$loglik - 2 * fit$df, BIC = bic, SABIC = sabic,
      ICL = bic - 2 * entropy, min_class_n = min(class_n),
      min_class_prop = min(class_n) / nrow(z_values),
      avg_max_posterior = mean(max_posterior),
      uncertain_lt_0_70 = mean(max_posterior < 0.70),
      relative_entropy = relative_entropy, quality_ok = FALSE, selected = FALSE,
      error_message = ""
    )
  }
}

model_selection <- do.call(rbind, candidate_rows)
model_selection$quality_ok <- with(
  model_selection,
  status == "ok" & G >= 2 & min_class_prop >= 0.05 &
    avg_max_posterior >= 0.80 & uncertain_lt_0_70 <= 0.20
)

eligible <- which(model_selection$quality_ok & is.finite(model_selection$BIC))
fallback <- FALSE
if (!length(eligible)) {
  eligible <- which(model_selection$status == "ok" & model_selection$G >= 2)
  fallback <- TRUE
}
if (!length(eligible)) stop("No usable G >= 2 candidate model converged")

# mclust BIC is 2 logLik - p log(n), so larger BIC is better. ICL and
# classification diagnostics remain visible in model_selection.csv/report.
selected_index <- eligible[which.max(model_selection$BIC[eligible])]
model_selection$selected[selected_index] <- TRUE
selected_row <- model_selection[selected_index, , drop = FALSE]
selected_key <- paste(selected_row$G, selected_row$model, sep = "_")
final_fit <- candidate_fits[[selected_key]]
if (is.null(final_fit)) stop("Selected fit is unavailable")
model_selection <- model_selection[order(-model_selection$BIC, model_selection$G, model_selection$model), ]
write_csv(model_selection, file.path(results_dir, "model_selection.csv"))

G <- as.integer(selected_row$G)
classification <- as.integer(final_fit$classification)
posterior <- as.matrix(final_fit$z)
colnames(posterior) <- paste0("posterior_class_", seq_len(G))
max_posterior <- apply(posterior, 1, max)

posterior_classification <- data.frame(
  SEQN = complete$SEQN,
  assigned_class = classification,
  posterior,
  max_posterior = max_posterior,
  uncertain_lt_0_70 = max_posterior < 0.70
)
write_csv(posterior_classification, file.path(results_dir, "posterior_classification.csv"))

average_posterior <- do.call(rbind, lapply(seq_len(G), function(k) {
  data.frame(assigned_class = k, n_assigned = sum(classification == k),
             as.list(colMeans(posterior[classification == k, , drop = FALSE])))
}))

class_n <- tabulate(classification, nbins = G)
profile_z <- do.call(rbind, lapply(seq_len(G), function(k) {
  colMeans(z_values[classification == k, , drop = FALSE])
}))
profile_names <- paste0("Profile ", seq_len(G))
if (G >= 2) {
  burden <- rowMeans(profile_z)
  profile_names[which.min(burden)] <- "Lower cardiometabolic burden"
  profile_names[which.max(burden)] <- "Higher cardiometabolic burden"
}

class_sizes <- data.frame(
  profile = seq_len(G), label = profile_names, n = class_n,
  proportion = class_n / nrow(complete),
  average_assigned_posterior = vapply(seq_len(G), function(k) mean(posterior[classification == k, k]), numeric(1)),
  n_below_0_70 = vapply(seq_len(G), function(k) sum(max_posterior[classification == k] < 0.70), integer(1))
)
write_csv(class_sizes, file.path(results_dir, "class_sizes.csv"))

profile_means <- do.call(rbind, lapply(seq_len(G), function(k) {
  data.frame(
    profile = k, label = profile_names[k], indicator = indicators,
    indicator_label = unname(indicator_labels[indicators]),
    unit = unname(indicator_units[indicators]),
    mean = colMeans(raw_values[classification == k, , drop = FALSE]),
    n_hard_class = class_n[k]
  )
}))
write_csv(profile_means, file.path(results_dir, "profile_means_original_scale.csv"))

plot_data <- do.call(rbind, lapply(seq_len(G), function(k) {
  data.frame(profile = profile_names[k], indicator = indicators, z_mean = profile_z[k, ])
}))
plot_data$indicator <- factor(plot_data$indicator, levels = indicators,
                              labels = unname(indicator_labels[indicators]))
plot_data$profile <- factor(plot_data$profile, levels = profile_names)
profile_plot <- ggplot(plot_data, aes(indicator, z_mean, group = profile, color = profile)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  labs(
    title = "NHANES 2017–2018 LPA profiles",
    subtitle = paste0("G=", G, " ", selected_row$model,
                      "; complete-case, unweighted exploratory analysis"),
    x = NULL, y = "Mean standardized indicator", color = "Profile"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
ggplot2::ggsave(file.path(figures_dir, "nhanes_lpa_profile.png"), profile_plot,
                width = 8.2, height = 5.0, dpi = 300)

model_plot_data <- model_selection[model_selection$status == "ok", , drop = FALSE]
model_plot <- ggplot(model_plot_data, aes(x = G, y = BIC, color = model, group = model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_point(
    data = model_plot_data[model_plot_data$selected, , drop = FALSE],
    shape = 21, size = 4, fill = "#f6c85f", color = "#333333", stroke = 1.1
  ) +
  scale_x_continuous(breaks = sort(unique(model_plot_data$G))) +
  labs(
    title = "Candidate model comparison",
    subtitle = "Higher BIC is preferred for this mclust specification; highlighted point is selected",
    x = "Number of profiles (G)", y = "BIC", color = "Covariance model"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggplot2::ggsave(file.path(figures_dir, "nhanes_lpa_model_selection.png"), model_plot,
                width = 8.2, height = 5.0, dpi = 300)

class_plot_data <- class_sizes
class_plot_data$label <- factor(class_plot_data$label, levels = class_plot_data$label)
class_plot <- ggplot(class_plot_data, aes(x = label, y = proportion, fill = label)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(
    aes(label = paste0(round(100 * proportion, 1), "%")),
    vjust = -0.35, size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, max(class_plot_data$proportion) * 1.18),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    title = "Estimated profile sizes",
    subtitle = "Hard classification counts shown as a proportion of complete cases",
    x = NULL, y = "Complete-case proportion"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
ggplot2::ggsave(file.path(figures_dir, "nhanes_lpa_class_sizes.png"), class_plot,
                width = 8.2, height = 5.0, dpi = 300)

session_lines <- c(
  paste0("Seed: ", seed),
  paste0("Input XPT files found: ", paste(raw_files, collapse = ", ")),
  paste0("Selected model: G=", G, " ", selected_row$model),
  paste0("Analysis N: ", nrow(complete)),
  capture.output(sessionInfo())
)
session_lines <- sub("[[:space:]]+$", "", session_lines)
writeLines(session_lines, file.path(results_dir, "session_info.txt"), useBytes = TRUE)

selected_profile_text <- paste(
  apply(profile_means[profile_means$profile == selected_row$G[1] * 0 + 1, c("indicator_label", "mean")], 1,
        function(r) paste0(r[[1]], "=", format(round(as.numeric(r[[2]]), 2), nsmall = 2))),
  collapse = "; "
)
report <- c(
  "# NHANES 2017–2018 minimal exploratory LPA",
  "",
  paste0("Run time: ", format(Sys.time(), tz = "Asia/Shanghai")),
  paste0("Analysis N: ", nrow(complete), "; random seed: ", seed, "."),
  "",
  "## Skill and analysis specification",
  "",
  "This run followed `advanced-stats/latent-profile-analysis/SKILL.md`, `references/model-selection.md`, and `references/datasets.md`. The three pre-specified continuous indicators were BMI (BMXBMI), waist circumference (BMXWAIST), and the mean of available BPXSY1–BPXSY4 systolic readings (SBP_MEAN).",
  "",
  paste0("The local raw directory contained: ", paste(raw_files, collapse = ", "), ". Only DEMO_J.XPT, BMX_J.XPT, and BPX_J.XPT were needed for these three indicators; the other local files were not used."),
  "",
  "## Sample selection and missingness",
  "",
  "Files were merged by SEQN. Participants were restricted to MEC-examined adults (RIDSTATR == 2 and RIDAGEYR >= 20). Ordinary `mclust` requires observed indicators, so complete cases on all three indicators were analyzed; no imputation was performed.",
  "",
  format_md(flow, 0),
  "",
  format_md(missingness, 2),
  "",
  "## Standardization and model selection",
  "",
  "The three indicators were standardized using the complete-case mean and SD. The derived file keeps the original-scale variables, z-score variables, MEC weight (WTMEC2YR), PSU (SDMVPSU), and strata (SDMVSTRA). Candidate models were G=1–4 with `EEI`, `VVI`, and `EEE` covariance structures. `mclust` BIC/AIC/SABIC/ICL were calculated with the larger BIC/ICL values indicating better relative fit.",
  "",
  "For a non-singleton candidate to be considered eligible, the minimum class proportion had to be at least 5%, mean maximum posterior probability at least 0.80, and the proportion with maximum posterior below 0.70 no more than 20%. The selected model was the eligible candidate with the largest BIC; ICL and classification diagnostics are retained for joint interpretation.",
  "",
  paste0("Selected model: G=", G, " ", selected_row$model,
         if (fallback) ". No candidate met all guardrails, so the best available G>=2 model by BIC was used." else "."),
  "",
  format_md(model_selection[, c("G", "model", "status", "BIC", "SABIC", "ICL", "min_class_n",
                                "min_class_prop", "avg_max_posterior", "uncertain_lt_0_70",
                                "relative_entropy", "quality_ok", "selected")], 3),
  "",
  "## Classification quality and key results",
  "",
  format_md(class_sizes, 3),
  "",
  format_md(average_posterior, 3),
  "",
  paste0("Overall mean maximum posterior probability was ", round(mean(max_posterior), 3),
         "; ", sum(max_posterior < 0.70), " of ", nrow(complete), " participants (",
         round(100 * mean(max_posterior < 0.70), 1),
         "%) had maximum posterior probability below 0.70 and were retained."),
  "",
  "Original-scale hard-class profile means:",
  "",
  format_md(profile_means[, c("profile", "label", "indicator_label", "unit", "mean", "n_hard_class")], 3),
  "",
  paste0("The accompanying figure plots standardized profile means. For orientation, Profile 1 has: ", selected_profile_text, "."),
  "",
  "## Limitations",
  "",
  "- Complete-case analysis can change the target population and may be biased if missingness is informative.",
  "- Ordinary `mclust::Mclust()` does not use NHANES WTMEC2YR, SDMVPSU, or SDMVSTRA. This is an unweighted, sample-internal exploratory LPA, not a nationally representative estimate; a formal study needs a survey-aware or sensitivity-analysis strategy.",
  "- Gaussian mixtures are sensitive to skewness, outliers, covariance restrictions, initialization, and local solutions. The candidate comparison is not a substitute for replication or substantive validation.",
  "- The three indicators define the profiles and should not be reused as independent external validation variables. Downstream analyses should account for classification uncertainty.",
  "- The data are cross-sectional; profiles are descriptive and do not establish causality or clinical subtypes.",
  "",
  "## Output files",
  "",
  "The script is `analysis/agent_lpa_analysis.R`. Outputs are `data/derived/analysis_data.csv`, `results/model_selection.csv`, `results/class_sizes.csv`, `results/profile_means_original_scale.csv`, `results/posterior_classification.csv`, `results/analysis_report.md`, `results/session_info.txt`, and the three figures in `figures/`."
)
writeLines(report, file.path(results_dir, "analysis_report.md"), useBytes = TRUE)

cat("Completed: N=", nrow(complete), "; selected G=", G, " ", selected_row$model,
    "; mean max posterior=", round(mean(max_posterior), 3), "\n", sep = "")
