# PUBLIC SOURCE-CODE COPY — SCRIPT 16
# This file is safe to share as SOURCE CODE after project review. It does NOT
# include clinical records, patient IDs, methylation arrays, or fitted NMF data.
# The expected input filenames are neutral examples; rename local files or edit
# the CONFIG paths to match authorized inputs. Runtime output is patient-level.
# The NMF and UMAP calculations are unchanged from the original analysis script.
# The chosen clinical covariates reflect a configurable research example.

#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT 16 — NMF H-MATRIX UMAP AND INTERACTIVE EXPLORATION
# ==============================================================================
# PROJECT: Pediatric AML DNA methylation / combined-sex NMF (k = 23)
# GOAL:
#   Project the 23-component NMF H representation of 854 patients into two
#   dimensions using uwot, then explore sample structure against clinical and
#   molecular annotations in both static ggplot2 and interactive Plotly plots.
#
# KEY DESIGN PRINCIPLES:
#   * Never refit NMF and never alter the underlying fitted H matrix.
#   * Retain all 23 components, including components with zero winning samples.
#   * Generate separate RAW and PATIENT-SUM-NORMALIZED H UMAPs.
#   * Use identical sample coordinates for every covariate plot within an
#     embedding; annotations only change color, not the UMAP geometry.
#   * Never interpret missing molecular lesion calls as lesion-negative.
#   * Preserve original NMF Margin and additionally calculate a relative margin.
#   * Provide a user-editable configuration block and timestamped run log.
#   * Retain software versions, parameters, data, QC, figures, and HTML outputs.
#
# METHODS NOTE:
#   This analysis is exploratory. UMAP does not test group separation or
#   survival association. The 2-D coordinates have no mechanistic meaning, and
#   inter-island distances and apparent gaps can be misleading. Statistical
#   inference belongs in the project's existing clinical/survival scripts.
#
# EXECUTION:
#   Rscript 16_NMF_H_UMAP_Interactive_Exploration.R
#   Adjust SECTION 01 below to customize paths, covariates and plot aesthetics.
# ==============================================================================

# ==============================================================================
# SECTION 01 — USER CONFIGURATION: EDIT THIS SECTION
# ==============================================================================
# Paths may be absolute or relative to ROOT_DIR. Input filenames are generic and can be adapted by users. The fitted NMF
# model is used only if the standalone H matrix is unavailable.
CONFIG <- list(
  # GitHub-ready portability: configure NMF_PROJECT_DIR on your own computer.
  # For example, in ~/.Renviron (NOT in this repository):
  # NMF_PROJECT_DIR=/path/to/private_project
  # Keep all controlled-access input files and generated outputs under this
  # private project root, outside the public GitHub checkout.
  ROOT_DIR = Sys.getenv("NMF_PROJECT_DIR", unset = ""),
  OUTPUT_DIR = "16_NMF_H_UMAP_Interactive_Exploration",
  H_FILE = "inputs/h_matrix.rds",
  MODEL_FILE = "inputs/nmf_model.rds",
  ASSIGNMENTS_FILE = "inputs/nmf_assignments.csv",
  CLINICAL_FILE = "inputs/clinical_annotations.csv",
  FALLBACK_CLINICAL_FILE = "inputs/clinical_annotations_fallback.csv",
  PATIENT_ID_COL = "USI",
  EXPECTED_PATIENTS = 854L,
  EXPECTED_COMPONENTS = 23L,
  SEED = 190L,

  # Each embedding is saved separately. Patients whose total H is zero cannot
  # be sum-normalized; fail loudly rather than silently dropping such patients.
  REPRESENTATIONS = c("raw", "normalized"),
  # IMPORTANT: cosine compares profile directions, Euclidean uses magnitudes.
  # We use cosine for normalized profiles, Euclidean for raw H initially.
  METRIC_BY_REPRESENTATION = c(raw = "euclidean", normalized = "cosine"),
  PRIMARY_N_NEIGHBORS = 30L,
  PRIMARY_MIN_DIST = 0.30,
  UMAP_SPREAD = 1.0,
  UMAP_N_EPOCHS = 500L,
  UMAP_N_THREADS = 1L,  # one thread for more reproducible results
  # Optional parameter grid: compared by neighborhood preservation only.
  RUN_SENSITIVITY_GRID = TRUE,
  GRID_N_NEIGHBORS = c(15L, 30L, 50L),
  GRID_MIN_DIST = c(0.1, 0.3, 0.5),
  GRID_MAX_RUNS_PER_REPRESENTATION = 9L,
  QC_K = 15L,

  # SELECT THE COVARIATES YOU WANT PLOTTED HERE. Missing columns are reported
  # and skipped (not silently substituted). Any numeric or categorical column
  # from the merged annotation table is valid. Column names must match exactly.
  # Example selection from the analyzed cohort. Missing optional fields are
  # reported and skipped; adapt exact names to your source clinical table.
  PLOT_COVARIATES = c(
    "NMF_Group", "Sex", "Risk.group.classification..by.protocol.1031",
    "Margin", "Margin_Relative", "Age.Category", "AgeGroup", "M6_AML",
    "M7_AML", "FLT3.ITD.positive.", "WT1.mutation.", "NPM.mutation.",
    "CEBPA.mutation.", "RAS.Gene", "Mutations.Category", "Primary.CNV",
    "Primary.Fusion", "Cytogenetic.Category.1", "Cytogenetic.Category.2",
    "t.9.11..p22.q23.", "t.8.21", "del5q", "inv.16.", "trisomy.21",
    "MLL", "Minus.Y", "Minus.X", "Primary.Cytogenetic.Code"
  ),
  # Optional: additional clinical fields to display in hover tooltips.
  TOOLTIP_COVARIATES = c(
    "Sex", "Risk.group.classification..by.protocol.1031", "NMF_Group",
    "Margin", "Margin_Relative", "OS.event.ID", "AgeGroup", "Age.Category"
  ),
  # Optional variables to interpret as binary lesion fields. If listed here,
  # e.g. "NUP98.NSD1" or "FLT3.ITD", 0/1, TRUE/FALSE, yes/no, +/- and
  # positive/negative are mapped to Present/Absent; unknown stays Unknown.
  # Include only columns truly coded as presence/absence. Check values first.
  BINARY_LESION_COLS = c(
    "M6_AML", "M7_AML", "FLT3.ITD.positive.", "WT1.mutation.",
    "NPM.mutation.", "CEBPA.mutation.", "t.9.11..p22.q23.", "t.8.21",
    "del5q", "inv.16.", "trisomy.21", "MLL", "Minus.Y", "Minus.X"
  ),
  # Optional more patient-level annotation files. Each must contain the patient ID and
  # one row per patient; other columns will be joined by patient ID.
  EXTRA_ANNOTATION_FILES = character(0),

  # Named override colors. The key is the EXACT column name (not a regex).
  # Add additional named palettes below as needed. Any unspecified categories
  # receive a consistent discrete color automatically.
  CATEGORICAL_PALETTES = list(
    Sex = c(Male = "#2B6CB0", Female = "#DD6B98", Unknown = "#999999"),
    Risk.group.classification..by.protocol.1031 = c(Low = "#239A3B",
                         High = "#D33F49", Unknown = "#999999")
  ),
  LESION_COLORS = c(Present = "#D62728", Absent = "#8C8C8C",
                    Unknown = "#D8D8D8"),
  CONTINUOUS_LOW_COLOR = "#2B6CB0",
  CONTINUOUS_MID_COLOR = "#F7F7F7",
  CONTINUOUS_HIGH_COLOR = "#C53030",
  # Additional continuous fields here are never converted to categorical.
  FORCE_CONTINUOUS = c("Margin", "Margin_Relative"),
  FORCE_CATEGORICAL = c(
    "NMF_Group", "Sex", "Risk.group.classification..by.protocol.1031",
    "Age.Category", "AgeGroup", "RAS.Gene", "Mutations.Category",
    "Primary.CNV", "Primary.Fusion", "Cytogenetic.Category.1",
    "Cytogenetic.Category.2", "Primary.Cytogenetic.Code"
  ),

  # VISUAL OPTIONS — applied to all figures unless edited in a plot function.
  PLOT_TITLE_PREFIX = "Pediatric AML NMF H UMAP",
  POINT_SIZE = 1.8,
  POINT_ALPHA = 0.85,
  FIGURE_WIDTH_IN = 9,
  FIGURE_HEIGHT_IN = 7,
  PNG_DPI = 350L,
  FONT_SIZE = 12,
  LEGEND_POSITION = "right",
  HTML_SELF_CONTAINED = TRUE,
  SAVE_PDF = TRUE,
  SAVE_PNG = TRUE,
  SAVE_HTML = TRUE,
  # Optional visualization-only patient filter: NULL means ALL patients. The
  # function receives the complete annotated data frame and returns TRUE/FALSE
  # of identical length. UMAP STILL FITS ALL PATIENTS before this filtering.
  PLOT_FILTER_FUNCTION = NULL,
  # Optional highlight IDs: when nonempty, other patients remain visible but
  # are shown faintly; the highlighted cohort uses the selected covariate.
  HIGHLIGHT_PATIENT_IDS = character(0)
)

# ==============================================================================
# SECTION 02 — PATHS, PACKAGE VALIDATION AND LOGGING
# ==============================================================================
# Console and log are both written directly
# An error handler records the failure and session information.
resolve_path <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("Invalid path")
  if (grepl("^([A-Za-z]:[/\\\\]|/|~)", path)) return(path.expand(path))
  file.path(CONFIG$ROOT_DIR, path)
}
if (!nzchar(CONFIG$ROOT_DIR)) {
  stop("Set NMF_PROJECT_DIR to the absolute private project directory before running Script 16")
}
if (!dir.exists(CONFIG$ROOT_DIR)) {
  stop("CONFIG$ROOT_DIR does not exist: ", CONFIG$ROOT_DIR,
       "\nUpdate SECTION 01 to point to the analysis directory.")
}
OUTPUT_DIR <- resolve_path(CONFIG$OUTPUT_DIR)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_PATH <- file.path(OUTPUT_DIR, "16_NMF_H_UMAP_Interactive_Exploration_Log.txt")
LOG_CON <- file(LOG_PATH, open = "wt")
START_TIME <- Sys.time()
log_line <- function(..., level = "INFO") {
  value <- paste(..., collapse = "")
  line <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  level, value)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = LOG_CON)
  flush(LOG_CON)
  invisible(line)
}
log_section <- function(section, title) {
  log_line(paste(rep("=", 76), collapse = ""))
  log_line(sprintf("SECTION %02d — %s", section, title))
  log_line(paste(rep("=", 76), collapse = ""))
}
output_record <- character(0)
record_output <- function(path) {
  output_record <<- unique(c(output_record, normalizePath(path, mustWork = FALSE)))
  log_line("Saved: ", path)
}
add_warning <- function(w) {
  log_line(conditionMessage(w), level = "WARNING")
  invokeRestart("muffleWarning")
}
close_logging <- function() {
  if (isOpen(LOG_CON)) {
    log_line("Completed logging; total runtime (minutes): ",
             round(as.numeric(difftime(Sys.time(), START_TIME, units = "mins")), 2))
    close(LOG_CON)
  }
}

# ==============================================================================
# SECTION 03 — GENERIC INPUT AND VALIDATION HELPERS
# ==============================================================================
required_packages <- c("uwot", "ggplot2", "plotly", "htmlwidgets")
check_packages <- function() {
  missing <- required_packages[!vapply(required_packages, requireNamespace,
                                       quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) stop("Install the missing packages first: ",
                            paste(missing, collapse = ", "),
                            ". For CRAN packages: install.packages(c('",
                            paste(missing, collapse = "','"), "'))")
  invisible(TRUE)
}
read_csv_checked <- function(path, label) {
  if (!file.exists(path)) stop(label, " file not found: ", path)
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  log_line(label, ": ", nrow(df), " rows × ", ncol(df),
           " columns; file=", path)
  df
}
assert_unique_ids <- function(ids, what) {
  ids <- as.character(ids)
  if (anyNA(ids) || any(!nzchar(trimws(ids))))
    stop(what, " has missing or empty patient IDs")
  if (anyDuplicated(ids))
    stop(what, " has duplicate patient IDs; examples: ",
         paste(head(unique(ids[duplicated(ids)]), 5), collapse = ", "))
  invisible(TRUE)
}
extract_H <- function() {
  h_path <- resolve_path(CONFIG$H_FILE)
  model_path <- resolve_path(CONFIG$MODEL_FILE)
  if (file.exists(h_path)) {
    log_line("Loading H from Script 01: ", h_path)
    H <- readRDS(h_path)
  } else if (file.exists(model_path)) {
    log_line("Script 01 H unavailable; reading fitted NMF object: ", model_path)
    fit <- readRDS(model_path)
    H <- attr(fit, "h")
    if (is.null(H) && is.list(fit)) H <- fit$h
    if (is.null(H)) stop("Cannot extract H using attr(model, 'h') or model$h")
    rm(fit)
  } else stop("Neither H_FILE nor MODEL_FILE exists. Edit SECTION 01.")
  if (!(is.matrix(H) || inherits(H, "Matrix") || is.data.frame(H)))
    stop("H must be a numeric matrix or Matrix object; got: ", paste(class(H), collapse = ", "))
  H <- as.matrix(H)
  storage.mode(H) <- "double"
  if (nrow(H) != CONFIG$EXPECTED_COMPONENTS &&
      ncol(H) == CONFIG$EXPECTED_COMPONENTS) {
    log_line("H appears patient-by-component; transposing to component-by-patient")
    H <- t(H)
  }
  if (nrow(H) != CONFIG$EXPECTED_COMPONENTS)
    stop("H component count mismatch: expected ", CONFIG$EXPECTED_COMPONENTS,
         ", found ", nrow(H), ". Do NOT drop components with no winning patient.")
  if (ncol(H) != CONFIG$EXPECTED_PATIENTS)
    stop("H patient count mismatch: expected ", CONFIG$EXPECTED_PATIENTS,
         ", found ", ncol(H))
  if (anyNA(H) || any(!is.finite(H))) stop("H contains NA/NaN/Inf")
  if (any(H < -1e-10)) stop("H contains substantially negative values: invalid NMF H")
  H[H < 0] <- 0
  if (is.null(rownames(H))) rownames(H) <- paste0("nmf", seq_len(nrow(H)))
  assert_unique_ids(rownames(H), "H component names")
  # Without trustworthy patient identifiers, a by-position join is unsafe.
  if (is.null(colnames(H)) || any(!nzchar(trimws(colnames(H)))))
    stop("H lacks complete patient-ID column names. Please recover IDs from the fitted NMF object or Script 01 H; refusing positional matching to clinical rows.")
  assert_unique_ids(colnames(H), "H column names")
  log_line("Validated H: ", nrow(H), " components × ", ncol(H), " patients")
  log_line("H value range: ", paste(signif(range(H), 5), collapse = " to "))
  H
}

# ==============================================================================
# SECTION 04 — LOAD ANNOTATIONS AND ALIGN BY EXPLICIT PATIENT IDS
# ==============================================================================
# CSVs from other scripts may have a synthetic row-number column; these are
# removed. We use match(USI), NEVER assume column/row order represents patients.
remove_csv_index <- function(df) {
  idx <- intersect(names(df), c("X", "...1", "row.names"))
  for (nm in idx) {
    if (identical(as.character(df[[nm]]), as.character(seq_len(nrow(df)))))
      df[[nm]] <- NULL
  }
  df
}
read_annotation_inputs <- function() {
  ap <- resolve_path(CONFIG$ASSIGNMENTS_FILE)
  cp <- resolve_path(CONFIG$CLINICAL_FILE)
  if (!file.exists(cp)) {
    fallback <- resolve_path(CONFIG$FALLBACK_CLINICAL_FILE)
    log_line("Merged clinical CSV unavailable; trying fallback ", fallback, level = "WARNING")
    cp <- fallback
  }
  clin <- remove_csv_index(read_csv_checked(cp, "Clinical annotation"))
  key <- CONFIG$PATIENT_ID_COL
  if (!key %in% names(clin)) stop("Clinical file lacks patient ID column: ", key)
  clin[[key]] <- as.character(clin[[key]])
  assert_unique_ids(clin[[key]], "Clinical data")
  if (file.exists(ap)) {
    assign <- remove_csv_index(read_csv_checked(ap, "NMF assignments"))
    if (!key %in% names(assign)) stop("Assignment file lacks column: ", key)
    assign[[key]] <- as.character(assign[[key]])
    assert_unique_ids(assign[[key]], "Assignments")
    # Preserve clinical fields; overlay canonical NMF assignment attributes.
    for (nm in intersect(c("NMF_Group", "NMF_Score", "Second_Score", "Margin"),
                         names(assign))) {
      newvals <- assign[[nm]][match(clin[[key]], assign[[key]])]
      if (nm %in% names(clin)) {
        conflict <- !is.na(newvals) & !is.na(clin[[nm]]) &
          as.character(newvals) != as.character(clin[[nm]])
        if (any(conflict)) log_line("Replacing ", nm, " for ", sum(conflict),
                                   " discrepant rows using canonical assignment file", level = "WARNING")
      }
      clin[[nm]] <- newvals
    }
  } else log_line("Assignment CSV absent: Margin and NMF_Group will be derived from H when needed", level = "WARNING")
  if (length(CONFIG$EXTRA_ANNOTATION_FILES)) {
    for (f in CONFIG$EXTRA_ANNOTATION_FILES) {
      extra <- read_csv_checked(resolve_path(f), paste0("Additional annotations: ", f))
      if (!key %in% names(extra)) stop("Extra annotation missing ", key, ": ", f)
      extra[[key]] <- as.character(extra[[key]])
      assert_unique_ids(extra[[key]], paste0("Extra annotation: ", f))
      duplicates <- intersect(setdiff(names(extra), key), names(clin))
      if (length(duplicates)) stop("Extra file has duplicate clinical columns: ",
                                    paste(duplicates, collapse = ", "),
                                    ". Rename them in the input first.")
      clin <- merge(clin, extra, by = key, all.x = TRUE, sort = FALSE)
    }
  }
  clin
}

# ==============================================================================
# SECTION 05 — NMF PROFILE PREPARATION AND MARGIN CALCULATION
# ==============================================================================
# Raw H keeps the fitted numerical component weights. Sum-normalized H divides
# each patient's component vector by the sum of all 23 weights; rows therefore
# sum to one, and represent relative contributions (not posterior probabilities).
# For patient i with ordered weights h_(1)>=h_(2), original Margin is compared
# to calculated raw margin h_(1)-h_(2); relative margin is that difference
# divided by the patient's total H weight. Neither constitutes a calibrated
# confidence probability, and NMF magnitude is scale-dependent.
prepare_h <- function(H) {
  X_raw <- t(H)
  rownames(X_raw) <- colnames(H)
  colnames(X_raw) <- rownames(H)
  total <- rowSums(X_raw)
  if (any(total <= 0)) stop("Zero-sum H profiles detected; unable to normalize. IDs: ",
                            paste(head(rownames(X_raw)[total <= 0], 8), collapse = ", "))
  X_norm <- sweep(X_raw, 1, total, "/")
  if (max(abs(rowSums(X_norm) - 1)) > 1e-10)
    stop("Normalized H rows do not sum to one")
  ordered <- t(apply(X_raw, 1, function(z) sort(z, decreasing = TRUE)[1:2]))
  # NMF rank k >= 2 is required by the above top-two operation.
  top <- ordered[, 1]
  second <- ordered[, 2]
  best_index <- max.col(X_raw, ties.method = "first")
  derived <- data.frame(
    USI = rownames(X_raw),
    NMF_Group_Derived = colnames(X_raw)[best_index],
    NMF_Score_Derived = top,
    Second_Score_Derived = second,
    Margin_Raw_Derived = top - second,
    Margin_Relative = (top - second) / total,
    H_Total = total,
    stringsAsFactors = FALSE
  )
  list(raw = X_raw, normalized = X_norm, calculated = derived)
}
combine_annotations <- function(prepared, clinical) {
  key <- CONFIG$PATIENT_ID_COL
  derived <- prepared$calculated
  if (!identical(key, "USI")) names(derived)[names(derived) == "USI"] <- key
  i <- match(derived[[key]], clinical[[key]])
  out <- cbind(derived, clinical[i, setdiff(names(clinical), key), drop = FALSE])
  out[[key]] <- as.character(out[[key]])
  matched <- sum(!is.na(i))
  log_line("Annotation match: ", matched, "/", nrow(derived),
           " H patients matched to clinical data")
  if (matched != nrow(derived))
    log_line("Patients missing clinical rows: ", paste(head(out[[key]][is.na(i)], 12), collapse = ", "), level = "WARNING")
  if (!"NMF_Group" %in% names(out)) out$NMF_Group <- out$NMF_Group_Derived
  if (!"NMF_Score" %in% names(out)) out$NMF_Score <- out$NMF_Score_Derived
  if (!"Second_Score" %in% names(out)) out$Second_Score <- out$Second_Score_Derived
  if (!"Margin" %in% names(out)) out$Margin <- out$Margin_Raw_Derived
  # Fill assignment-derived values only when assignments are missing. All
  # deviations between provided and calculated values are reported for review.
  check <- !is.na(out$Margin) & is.finite(suppressWarnings(as.numeric(out$Margin)))
  difference <- abs(suppressWarnings(as.numeric(out$Margin[check])) - out$Margin_Raw_Derived[check])
  log_line("Provided vs derived raw Margin maximum absolute discrepancy: ",
           if (length(difference)) signif(max(difference), 6) else "not comparable")
  if (length(difference) && any(difference > 1e-5))
    log_line("Provided Margin differs from derived raw Margin for ", sum(difference > 1e-5),
             " samples. Original Margin retained. Check original definition/scaling.", level = "WARNING")
  for (nm in c("NMF_Group", "NMF_Score", "Second_Score", "Margin")) {
    replacement <- switch(nm, NMF_Group = "NMF_Group_Derived",
                          NMF_Score = "NMF_Score_Derived",
                          Second_Score = "Second_Score_Derived",
                          Margin = "Margin_Raw_Derived")
    missing <- is.na(out[[nm]]) | trimws(as.character(out[[nm]])) == ""
    out[[nm]][missing] <- out[[replacement]][missing]
  }
  mismatch <- as.character(out$NMF_Group) != as.character(out$NMF_Group_Derived)
  mismatch[is.na(mismatch)] <- FALSE
  if (any(mismatch)) log_line("NMF group assignment differs from max.col(H) for ",
                              sum(mismatch), " patients (possible tie/scaling issue).", level = "WARNING")
  log_line("Original NMF groups observed: ", length(unique(out$NMF_Group)),
           "; H components retained: ", ncol(prepared$raw))
  out
}

# ==============================================================================
# SECTION 06 — ANNOTATION CODING: CATEGORICAL, CONTINUOUS, BINARY LESIONS
# ==============================================================================
normalize_binary <- function(z) {
  y <- tolower(trimws(as.character(z)))
  out <- rep("Unknown", length(y))
  out[y %in% c("1", "true", "t", "yes", "y", "positive", "pos", "+", "present")] <- "Present"
  out[y %in% c("0", "false", "f", "no", "n", "negative", "neg", "-", "absent")] <- "Absent"
  unresolved <- !is.na(y) & nzchar(y) & !(y %in% c(
    "1", "true", "t", "yes", "y", "positive", "pos", "+", "present",
    "0", "false", "f", "no", "n", "negative", "neg", "-", "absent",
    "na", "n/a", "unknown", "not tested", "not available", "missing"))
  if (any(unresolved)) log_line("Unrecognized binary values retained as Unknown: ",
                                paste(unique(head(y[unresolved], 12)), collapse = ", "), level = "WARNING")
  factor(out, levels = c("Absent", "Present", "Unknown"))
}
make_plot_field <- function(data, variable) {
  x <- data[[variable]]
  if (variable %in% CONFIG$BINARY_LESION_COLS) return(normalize_binary(x))
  if (variable %in% CONFIG$FORCE_CATEGORICAL) {
    vals <- as.character(x)
    vals[is.na(vals) | !nzchar(trimws(vals))] <- "Unknown"
    return(factor(vals))
  }
  if (variable %in% CONFIG$FORCE_CONTINUOUS || is.numeric(x)) {
    nx <- suppressWarnings(as.numeric(as.character(x)))
    if (variable %in% CONFIG$FORCE_CONTINUOUS &&
        any(!is.na(x) & is.na(nx)))
      stop("FORCE_CONTINUOUS column contains nonnumeric nonmissing data: ", variable)
    return(nx)
  }
  vals <- as.character(x)
  vals[is.na(vals) | !nzchar(trimws(vals))] <- "Unknown"
  factor(vals)
}

# ==============================================================================
# SECTION 07 — REPRODUCIBLE UMAP FITTING
# ==============================================================================
# UMAP is run on patient × component matrices. With cosine metric, relative
# profiles are compared by direction; with Euclidean, distances depend on H
# magnitude. 
fit_umap <- function(X, metric, n_neighbors, min_dist, seed, label) {
  if (n_neighbors >= nrow(X)) stop("n_neighbors must be less than patient count")
  log_line("UMAP ", label, ": n_neighbors=", n_neighbors,
           ", min_dist=", min_dist, ", metric=", metric,
           ", seed=", seed, ", threads=", CONFIG$UMAP_N_THREADS)
  set.seed(seed)
  fit <- uwot::umap(
    X = X,
    n_neighbors = as.integer(n_neighbors),
    n_components = 2L,
    metric = metric,
    min_dist = min_dist,
    spread = CONFIG$UMAP_SPREAD,
    n_epochs = as.integer(CONFIG$UMAP_N_EPOCHS),
    n_threads = as.integer(CONFIG$UMAP_N_THREADS),
    init = "spectral",
    verbose = FALSE,
    ret_model = FALSE
  )
  fit <- as.matrix(fit)
  if (!identical(dim(fit), c(nrow(X), 2L))) stop("Unexpected UMAP dimensions")
  rownames(fit) <- rownames(X)
  colnames(fit) <- c("UMAP1", "UMAP2")
  fit
}

# ==============================================================================
# SECTION 08 — NEIGHBORHOOD-PRESERVATION QC (EXPLORATORY, NOT INFERENCE)
# ==============================================================================
# Compute the mean fraction of each patient's k nearest high-dimensional
# neighbors also found among its k nearest 2-D UMAP neighbors. Random baseline
# is approximately k/(n-1). This is a descriptive diagnostic, not a test of
# biological clustering or a substitute for trustworthiness/continuity metrics.
# For 854 patients, all-pairs distance calculations are tractable (~729k pairs).
cosine_dist_matrix <- function(X) {
  row_norm <- sqrt(rowSums(X * X))
  if (any(row_norm <= 0)) stop("Cosine distance undefined for all-zero row")
  Z <- sweep(X, 1, row_norm, "/")
  D <- 1 - tcrossprod(Z)
  D <- pmax(D, 0)
  diag(D) <- 0
  D
}
neighbor_overlap <- function(X, U, metric, k) {
  n <- nrow(X)
  k <- min(as.integer(k), n - 1L)
  D_high <- if (metric == "cosine") cosine_dist_matrix(X) else
    as.matrix(stats::dist(X, method = "euclidean"))
  D_low <- as.matrix(stats::dist(U, method = "euclidean"))
  diag(D_high) <- Inf
  diag(D_low) <- Inf
  a <- t(apply(D_high, 1, order))[, seq_len(k), drop = FALSE]
  b <- t(apply(D_low, 1, order))[, seq_len(k), drop = FALSE]
  mean(vapply(seq_len(n), function(i) sum(a[i, ] %in% b[i, ]) / k, numeric(1)))
}

# ==============================================================================
# SECTION 09 — PLOT CUSTOMIZATION AND CUSTOM HOVER INFORMATION
# ==============================================================================
# The user-editable PLOT_COVARIATES vector controls the main color variable.
# Hover fields are independent of it, and can contain multiple clinical values.
# The named palette list controls specific category colors. Unknown covariate
# columns are skipped with explicit log messages instead of silently changing
# a biological category definition.
slug <- function(x) {
  x <- gsub("[^[:alnum:]_-]+", "_", x)
  gsub("^_+|_+$", "", x)
}
make_hover <- function(df, variable) {
  keep <- unique(c(CONFIG$PATIENT_ID_COL, variable, CONFIG$TOOLTIP_COVARIATES))
  keep <- intersect(keep, names(df))
  vapply(seq_len(nrow(df)), function(i) {
    paste(vapply(keep, function(nm) {
      val <- as.character(df[[nm]][i])
      if (is.na(val) || !nzchar(trimws(val))) val <- "Unknown"
      paste0(nm, ": ", val)
    }, character(1)), collapse = "<br>")
  }, character(1))
}
get_palette <- function(variable, values) {
  if (variable %in% CONFIG$BINARY_LESION_COLS) return(CONFIG$LESION_COLORS)
  candidate <- CONFIG$CATEGORICAL_PALETTES[[variable]]
  values <- unique(as.character(values))
  if (is.null(candidate)) candidate <- character(0)
  missing <- setdiff(values, names(candidate))
  if (length(missing)) {
    # hues are deterministic and avoid requiring another plotting package.
    hues <- grDevices::hcl(h = seq(15, 375, length.out = length(missing) + 1L)[
      seq_len(length(missing))], c = 70, l = 55)
    names(hues) <- missing
    candidate <- c(candidate, hues)
  }
  candidate
}
apply_plot_filter <- function(df) {
  f <- CONFIG$PLOT_FILTER_FUNCTION
  if (is.null(f)) return(df)
  if (!is.function(f)) stop("PLOT_FILTER_FUNCTION must be NULL or a function")
  keep <- f(df)
  if (!is.logical(keep) || length(keep) != nrow(df) || anyNA(keep))
    stop("PLOT_FILTER_FUNCTION must return one non-missing TRUE/FALSE per patient")
  df[keep, , drop = FALSE]
}
make_ggplot <- function(df, variable, representation) {
  plot_data <- apply_plot_filter(df)
  if (!nrow(plot_data)) stop("PLOT_FILTER_FUNCTION leaves no patients")
  plot_data$PlotValue <- make_plot_field(plot_data, variable)
  plot_data$HoverText <- make_hover(plot_data, variable)
  is_numeric <- is.numeric(plot_data$PlotValue)
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = UMAP1, y = UMAP2)
  )
  highlighted <- CONFIG$HIGHLIGHT_PATIENT_IDS
  if (length(highlighted)) {
    other <- plot_data[!(plot_data[[CONFIG$PATIENT_ID_COL]] %in% highlighted), , drop = FALSE]
    selected <- plot_data[plot_data[[CONFIG$PATIENT_ID_COL]] %in% highlighted, , drop = FALSE]
    if (nrow(other)) p <- p + ggplot2::geom_point(data = other,
                                                  color = "#D9D9D9", alpha = 0.20,
                                                  size = CONFIG$POINT_SIZE)
    plot_data <- selected
    if (!nrow(plot_data)) stop("No highlighted patient IDs remain after filtering")
  }
  p <- p + ggplot2::geom_point(
    data = plot_data,
    ggplot2::aes(color = PlotValue, text = HoverText),
    size = CONFIG$POINT_SIZE, alpha = CONFIG$POINT_ALPHA,
    show.legend = TRUE
  )
  if (is_numeric) {
    p <- p + ggplot2::scale_color_gradient2(
      low = CONFIG$CONTINUOUS_LOW_COLOR,
      mid = CONFIG$CONTINUOUS_MID_COLOR,
      high = CONFIG$CONTINUOUS_HIGH_COLOR,
      midpoint = stats::median(plot_data$PlotValue, na.rm = TRUE),
      na.value = "#BDBDBD", name = variable
    )
  } else {
    p <- p + ggplot2::scale_color_manual(
      values = get_palette(variable, plot_data$PlotValue),
      drop = FALSE, name = variable
    )
  }
  p + ggplot2::labs(
    title = paste(CONFIG$PLOT_TITLE_PREFIX, "—", toupper(representation)),
    subtitle = paste0("Colored by ", variable, " | n = ", nrow(plot_data)),
    x = "UMAP 1", y = "UMAP 2",
    caption = "UMAP of NMF H profiles; visualization is exploratory, not a significance test"
  ) + ggplot2::theme_classic(base_size = CONFIG$FONT_SIZE) +
    ggplot2::theme(
      legend.position = CONFIG$LEGEND_POSITION,
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "#555555")
    ) + ggplot2::coord_equal()
}

# ==============================================================================
# SECTION 10 — FIGURE EXPORT: PUBLICATION AND INTERACTIVE OUTPUTS
# ==============================================================================
# Plotly::ggplotly() preserves the mapping of UMAP coordinates. Plotly hover
# shows patient-level fields in HoverText. HTML widgets may include dependencies;
# selfcontained = TRUE creates a portable HTML, if the environment has pandoc.
export_plots <- function(df, representation, target_dir) {
  plot_dir <- file.path(target_dir, "Plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  requested <- unique(CONFIG$PLOT_COVARIATES)
  available <- intersect(requested, names(df))
  unavailable <- setdiff(requested, names(df))
  if (length(unavailable)) log_line("Requested plot variables absent; skipped: ",
                                    paste(unavailable, collapse = ", "), level = "WARNING")
  if (!length(available)) stop("None of PLOT_COVARIATES is available")
  for (variable in available) {
    log_line("Plotting ", representation, " / ", variable)
    p <- make_ggplot(df, variable, representation)
    stem <- file.path(plot_dir, paste0("UMAP_", representation, "_", slug(variable)))
    if (isTRUE(CONFIG$SAVE_PDF)) {
      dest <- paste0(stem, ".pdf")
      ggplot2::ggsave(dest, p, width = CONFIG$FIGURE_WIDTH_IN,
                      height = CONFIG$FIGURE_HEIGHT_IN,
                      units = "in", device = grDevices::cairo_pdf)
      record_output(dest)
    }
    if (isTRUE(CONFIG$SAVE_PNG)) {
      dest <- paste0(stem, ".png")
      ggplot2::ggsave(dest, p, width = CONFIG$FIGURE_WIDTH_IN,
                      height = CONFIG$FIGURE_HEIGHT_IN, units = "in",
                      dpi = CONFIG$PNG_DPI, bg = "white")
      record_output(dest)
    }
    if (isTRUE(CONFIG$SAVE_HTML)) {
      dest <- paste0(stem, ".html")
      interactive_plot <- plotly::ggplotly(p, tooltip = "text")
      # Portable HTML requires pandoc; fallback to dependency-directory mode
      # rather than dropping the interactive output.
      tryCatch({
        htmlwidgets::saveWidget(interactive_plot, dest,
                                selfcontained = isTRUE(CONFIG$HTML_SELF_CONTAINED))
      }, error = function(e) {
        log_line("Self-contained HTML failed: ", conditionMessage(e),
                 "; retrying non-self-contained", level = "WARNING")
        htmlwidgets::saveWidget(interactive_plot, dest, selfcontained = FALSE)
      })
      record_output(dest)
    }
  }
  invisible(available)
}

# ==============================================================================
# SECTION 11 — MAIN ANALYSIS RUNNER
# ==============================================================================
run_script16 <- function() {
  log_section(1, "Analysis initiation and software environment")
  log_line("Started: ", format(START_TIME, "%Y-%m-%d %H:%M:%S %Z"))
  log_line("ROOT_DIR: ", CONFIG$ROOT_DIR)
  log_line("OUTPUT_DIR: ", OUTPUT_DIR)
  log_line("R version: ", R.version.string)
  check_packages()
  for (pkg in required_packages)
    log_line(pkg, " version: ", as.character(utils::packageVersion(pkg)))
  log_line("Selected plot covariates: ", paste(CONFIG$PLOT_COVARIATES, collapse = ", "))
  log_line("Random seed: ", CONFIG$SEED)

  log_section(2, "H matrix validation")
  H <- extract_H()
  prepared <- prepare_h(H)
  for (r in CONFIG$REPRESENTATIONS) {
    if (!r %in% c("raw", "normalized")) stop("Unknown representation: ", r)
    log_line("Representation ", r, ": ", nrow(prepared[[r]]),
             " patients × ", ncol(prepared[[r]]), " components")
  }

  log_section(3, "Clinical and molecular annotation merging")
  clinical <- read_annotation_inputs()
  annotations <- combine_annotations(prepared, clinical)
  # Report annotation missingness for each requested covariate.
  for (v in intersect(CONFIG$PLOT_COVARIATES, names(annotations)))
    log_line("Missing in ", v, ": ", sum(is.na(annotations[[v]])),
             "/", nrow(annotations))
  write.csv(annotations, file.path(OUTPUT_DIR, "Patient_Annotations_With_Margins.csv"),
            row.names = FALSE, na = "")
  record_output(file.path(OUTPUT_DIR, "Patient_Annotations_With_Margins.csv"))

  log_section(4, "Save exact H matrices used for embedding")
  saveRDS(list(raw = prepared$raw, normalized = prepared$normalized,
               component_names = colnames(prepared$raw)),
          file.path(OUTPUT_DIR, "H_Representations_Used.rds"))
  record_output(file.path(OUTPUT_DIR, "H_Representations_Used.rds"))
  all_qc <- list()
  all_embeddings <- list()

  for (rep_name in CONFIG$REPRESENTATIONS) {
    rep_dir <- file.path(OUTPUT_DIR, paste0("H_", rep_name))
    dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
    X <- prepared[[rep_name]]
    metric <- unname(CONFIG$METRIC_BY_REPRESENTATION[[rep_name]])
    if (is.null(metric) || !metric %in% c("euclidean", "cosine"))
      stop("Metric must be euclidean or cosine for ", rep_name)

    log_section(5, paste0("Fit primary ", rep_name, " H UMAP"))
    embedding <- fit_umap(X, metric, CONFIG$PRIMARY_N_NEIGHBORS,
                          CONFIG$PRIMARY_MIN_DIST, CONFIG$SEED, rep_name)
    coords <- data.frame(USI = rownames(embedding), UMAP1 = embedding[, 1],
                         UMAP2 = embedding[, 2], stringsAsFactors = FALSE)
    if (!identical(CONFIG$PATIENT_ID_COL, "USI")) names(coords)[1] <- CONFIG$PATIENT_ID_COL
    # Preserve the embedding's exact order and join clinical metadata by ID.
    idx <- match(coords[[CONFIG$PATIENT_ID_COL]], annotations[[CONFIG$PATIENT_ID_COL]])
    if (anyNA(idx)) stop("Missing patient annotation row for an embedded sample")
    coords <- cbind(coords, annotations[idx, setdiff(names(annotations), CONFIG$PATIENT_ID_COL),
                                        drop = FALSE])
    stopifnot(nrow(coords) == CONFIG$EXPECTED_PATIENTS)
    coords_path <- file.path(rep_dir, "UMAP_Coordinates_Annotated.csv")
    write.csv(coords, coords_path, row.names = FALSE, na = "")
    record_output(coords_path)
    saveRDS(list(umap = embedding, annotations = coords, representation = rep_name,
                 metric = metric, n_neighbors = CONFIG$PRIMARY_N_NEIGHBORS,
                 min_dist = CONFIG$PRIMARY_MIN_DIST, seed = CONFIG$SEED),
            file.path(rep_dir, "UMAP_Primary_Objects.rds"))
    record_output(file.path(rep_dir, "UMAP_Primary_Objects.rds"))
    all_embeddings[[rep_name]] <- coords

    log_section(6, paste0("Neighborhood preservation QC: ", rep_name))
    primary_score <- neighbor_overlap(X, embedding, metric, CONFIG$QC_K)
    log_line("Primary mean ", CONFIG$QC_K, "-NN overlap: ", signif(primary_score, 4),
             "; approximate random expectation: ",
             signif(CONFIG$QC_K / (nrow(X) - 1), 4))
    qc <- data.frame(representation = rep_name,
                     n_neighbors = CONFIG$PRIMARY_N_NEIGHBORS,
                     min_dist = CONFIG$PRIMARY_MIN_DIST,
                     metric = metric, seed = CONFIG$SEED,
                     qc_k = CONFIG$QC_K, mean_neighbor_overlap = primary_score,
                     run_type = "primary", stringsAsFactors = FALSE)

    # Sensitivity embeddings are saved for review and QC. They do NOT overwrite
    # the predefined primary coordinates or silently change final figures.
    if (isTRUE(CONFIG$RUN_SENSITIVITY_GRID)) {
      grid <- expand.grid(n_neighbors = CONFIG$GRID_N_NEIGHBORS,
                          min_dist = CONFIG$GRID_MIN_DIST,
                          KEEP.OUT.ATTRS = FALSE)
      grid <- grid[seq_len(min(nrow(grid),
                               CONFIG$GRID_MAX_RUNS_PER_REPRESENTATION)), , drop = FALSE]
      sensitivity_dir <- file.path(rep_dir, "Sensitivity_Embeddings")
      dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)
      for (j in seq_len(nrow(grid))) {
        nn <- as.integer(grid$n_neighbors[j]); md <- as.numeric(grid$min_dist[j])
        if (nn == CONFIG$PRIMARY_N_NEIGHBORS &&
            isTRUE(all.equal(md, CONFIG$PRIMARY_MIN_DIST))) next
        # Same seed across settings permits consistent reproducibility.
        U <- fit_umap(X, metric, nn, md, CONFIG$SEED,
                      paste(rep_name, "sensitivity", j))
        score <- neighbor_overlap(X, U, metric, CONFIG$QC_K)
        log_line("Sensitivity result: k=", nn, ", min_dist=", md,
                 ", overlap=", signif(score, 4))
        qc <- rbind(qc, data.frame(representation = rep_name,
                                   n_neighbors = nn, min_dist = md,
                                   metric = metric, seed = CONFIG$SEED,
                                   qc_k = CONFIG$QC_K,
                                   mean_neighbor_overlap = score,
                                   run_type = "sensitivity"))
        dest <- file.path(sensitivity_dir,
                          sprintf("UMAP_neighbors_%s_min_dist_%s.rds", nn,
                                  gsub("\\.", "p", as.character(md))))
        saveRDS(U, dest)
        record_output(dest)
      }
    }
    qc_path <- file.path(rep_dir, "UMAP_Neighborhood_QC.csv")
    write.csv(qc, qc_path, row.names = FALSE)
    record_output(qc_path)
    all_qc[[rep_name]] <- qc

    log_section(7, paste0("Save static and interactive ", rep_name, " H plots"))
    plotted <- export_plots(coords, rep_name, rep_dir)
    log_line("Successfully plotted variables: ", paste(plotted, collapse = ", "))
  }

  log_section(8, "Provenance, final QC and output manifest")
  saveRDS(list(config = CONFIG, annotations = annotations, embeddings = all_embeddings,
               qc = all_qc, start_time = START_TIME,
               session_info = utils::sessionInfo()),
          file.path(OUTPUT_DIR, "Script_16_Complete_Analysis_Objects.rds"))
  record_output(file.path(OUTPUT_DIR, "Script_16_Complete_Analysis_Objects.rds"))
  write.csv(do.call(rbind, all_qc),
            file.path(OUTPUT_DIR, "All_UMAP_Parameter_QC.csv"), row.names = FALSE)
  record_output(file.path(OUTPUT_DIR, "All_UMAP_Parameter_QC.csv"))
  writeLines(capture.output(utils::sessionInfo()),
             con = file.path(OUTPUT_DIR, "R_Session_Info.txt"))
  record_output(file.path(OUTPUT_DIR, "R_Session_Info.txt"))
  log_line("Output files created: ", length(output_record) + 1L,
           " (including manifest)")
  manifest <- file.path(OUTPUT_DIR, "Output_Manifest.txt")
  writeLines(c(output_record, normalizePath(manifest, mustWork = FALSE),
               normalizePath(LOG_PATH, mustWork = FALSE)), con = manifest)
  record_output(manifest)
  log_line("SUCCESS: Script 16 finished without fatal errors")
  invisible(list(annotations = annotations, embeddings = all_embeddings,
                 qc = all_qc, outputs = output_record))
}

# ==============================================================================
# SECTION 12 — TOP-LEVEL ERROR REPORTING AND LOG CLEANUP
# ==============================================================================
# withCallingHandlers captures nonfatal warnings. Fatal errors are written to
# the log and rethrown after cleanup, ensuring a nonzero Rscript exit status.
script16_result <- tryCatch(
  withCallingHandlers(run_script16(), warning = add_warning),
  error = function(e) {
    log_line("FATAL ERROR: ", conditionMessage(e), level = "ERROR")
    log_line("Traceback follows when available:", level = "ERROR")
    tb <- capture.output(traceback())
    for (line in tb) log_line(line, level = "ERROR")
    structure(list(error = e), class = "script16_failure")
  },
  finally = close_logging()
)
if (inherits(script16_result, "script16_failure")) stop(script16_result$error)
# ==============================================================================
# END OF SCRIPT 16
# ==============================================================================
