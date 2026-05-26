# install_packages.R
# Install CRAN packages required to reproduce the simulation and figures.

pkgs <- c("data.table", "ggplot2", "survival")
inst <- rownames(installed.packages())
need <- setdiff(pkgs, inst)
if (length(need) > 0) {
  install.packages(need, repos = "https://cloud.r-project.org")
}

invisible(lapply(pkgs, library, character.only = TRUE))
cat("Required packages are available.\n")
sessionInfo()
