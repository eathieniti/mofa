install.packages(data.table)
install.packages(ggplot2)
install.packages(tidyverse)

library(MOFA2)
install.packages(randomForest)
install.packages(survival)
install.packages(survminer)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MOFA2")