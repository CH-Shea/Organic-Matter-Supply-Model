# Unique packages to install
pkgs <- c(
  "knitr", "readxl", "MASS", "FactoMineR", "factoextra", "runjags", "coda",
  "DT", "ggplot2", "ggpubr", "gridExtra", "reshape2", "magrittr",
  "compositions", "DirichletReg", "Ternary", "vegan", "kableExtra",
  "ggforce", "plotly", "graphics", "psych", "openxlsx"
)

# Check which packages are not installed
not_installed <- pkgs[!sapply(pkgs, require, character.only = TRUE, quietly = TRUE)]

# Install missing packages
if (length(not_installed) > 0) {
  install.packages(not_installed, repos = "https://cran.rstudio.com/", dependencies = TRUE)
}

# Load all packages quietly
invisible(sapply(pkgs, library, character.only = TRUE, quietly = TRUE))