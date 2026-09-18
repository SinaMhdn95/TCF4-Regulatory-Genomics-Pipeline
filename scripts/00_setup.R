# Shared project paths -----------------------------------------------------

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Project root not found. Run this code from inside the repository.")
    }
    current <- parent
  }
}

project_root <- find_project_root()
data_root <- Sys.getenv(
  "TCF4_DATA_DIR",
  unset = file.path(project_root, "data", "raw")
)
interim_dir <- file.path(project_root, "data", "interim")
processed_dir <- file.path(project_root, "data", "processed")
results_dir <- file.path(project_root, "results")

dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", project_root)
message("Data root: ", data_root)
