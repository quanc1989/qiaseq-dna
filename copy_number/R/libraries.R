# no warnings here (production)
#options(warn=-1)

# load libraries 
library(grid)
library(ggplot2)
library(gridExtra)
library(extrafont)
library(MASS)
library(naturalsort)
library(scales)

# turn warnings back on
#options(warn=1)

# load own sources 
source("config.R")
source("multiplot.R")
source("weighted_ttest.R")
source("cnv_calc_functions.R")
source("import.R")
source("segment.R")
source("themes.R")
source("vcf_header.R")
