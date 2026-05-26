# R/00_utils.R
# ------------------------------------------------------------
# Project utilities
# - Robust root detection
# - Package checks
# - Logging and session info
# - Small helpers used by scripts/
# ------------------------------------------------------------

.detect_root <- function(max_up = 6L) {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  is_root <- function(p) {
    file.exists(file.path(p, "scripts")) && file.exists(file.path(p, "R"))
  }
  if (is_root(wd)) return(wd)
  for (k in seq_len(max_up)) {
    p <- normalizePath(file.path(wd, paste(rep("..", k), collapse = "/")), winslash = "/", mustWork = FALSE)
    if (is_root(p)) return(p)
  }
  wd
}

init_project <- function(n_cores = 1L, seed = 2026L, out_dir = "output") {
  root <- .detect_root()
  try(setwd(root), silent = TRUE)

  out_dir <- as.character(out_dir)
  if (!nzchar(out_dir)) out_dir <- "output"

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  log_dir <- file.path(out_dir, "logs")
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

  obj <- list(
    root = root,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    out_dir = out_dir,
    log_dir = log_dir
  )

  obj$require_pkgs <- function(pkgs) {
    pkgs <- unique(as.character(pkgs))
    missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing) > 0) {
      stop(
        "Missing required packages: ", paste(missing, collapse = ", "), "\n",
        "Install them (example):\n",
        "  install.packages(c(\"", paste(missing, collapse = "\", \""), "\"))\n",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  obj$log_line <- function(path, text) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    cat(text, file = path, append = TRUE)
    invisible(TRUE)
  }

  obj$write_session_info <- function(prefix = "run") {
    f <- file.path(obj$log_dir, paste0("sessionInfo_", prefix, ".txt"))
    txt <- capture.output(utils::sessionInfo())
    writeLines(txt, con = f)
    invisible(f)
  }

  obj
}

# Deterministic seed helper
seed_for_job <- function(base_seed, scenario_id, truth = "", rep_id = 1L) {
  s <- sum(utf8ToInt(as.character(scenario_id))) %% 100000L
  t <- sum(utf8ToInt(as.character(truth))) %% 1000L
  as.integer(base_seed + 1000L * s + 10L * t + as.integer(rep_id))
}

# Safe file writer helper
write_dt_csv <- function(dt, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(dt, path)
  invisible(path)
}
