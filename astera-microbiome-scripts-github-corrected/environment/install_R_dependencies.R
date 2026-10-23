# Dependency recipe, NOT a historical version lock. Run in a dedicated R environment.
cran <- c("vegan", "ggplot2", "dplyr", "tidyr", "tibble", "stringr", "readr", "purrr", "forcats", "igraph", "Hmisc", "RColorBrewer", "scales", "BiocManager")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) install.packages(missing, repos="https://cloud.r-project.org")
for(p in c("phyloseq", "Maaslin2")) if(!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,ask=FALSE,update=FALSE)
writeLines(capture.output(sessionInfo()), "R_sessionInfo_installed.txt")
