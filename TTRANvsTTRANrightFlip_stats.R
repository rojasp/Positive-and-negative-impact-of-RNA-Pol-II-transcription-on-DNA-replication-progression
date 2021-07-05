################################################################################
#Project: genetics breakdown
#File: TTRANvsTTRANrightFlip_stats.r
#Author: Giovanni Guglielmi
#Description: perform the analysis on the right reads
################################################################################

#set the environment
setwd("C:/Users/Gianni/OneDrive - The University of Melbourne/uniPhD/projects/computational_genomics")

#library
source("my_functionR/boot2samples.R")
source("my_functionR/perm2samples.R")
source("my_functionR/clean_data.R")
source("my_functionR/zeroToTrue.R")
source("my_functionR/histBootstrap.R")
source("my_functionR/histPermutation.R")

################################################################################
##########################READ DATA#############################################
################################################################################

TTRAN3 <- read.csv(
  "2.datasets_panel/ES3_flip/ES3_peaks_transcribed_genes_tran_distance_breakpoint_r_flip",
  header = F, sep = "\t")

TTRAN4 <- read.csv(
  "2.datasets_panel/ES4_flip/ES4_peaks_tran_genes_distance_breakpoint_r_flip",
  header = F, sep = "\t")

TTRAN5 <- read.csv(
  "2.datasets_panel/ES5_flip/ES5_peaks_tran_genes_distance_breakpoint_r_flip",
  header = F, sep = "\t")

TTRANLS <- read.csv(
  "2.datasets_panel/LS_flip/LS_peaks_tran_genes_distance_breakpoint_r_flip",
  header = F, sep = "\t")

TTRANMS <- read.csv(
  "2.datasets_panel/MS_flip/MS_peaks_tran_genes_distance_breakpoint_r_flip",
  header = F, sep = "\t")

###############################################################################
############################CLEAN DATA#########################################
###############################################################################
ttran3_clean <- clean_data(TTRAN3[, ncol(TTRAN3)])
ttran3 <- zeroToTrue(ttran3_clean)

ttran4_clean <- clean_data(TTRAN4[, ncol(TTRAN4)])
ttran4 <- zeroToTrue(ttran4_clean)

ttran5_clean <- clean_data(TTRAN5[, ncol(TTRAN5)])
ttran5 <- zeroToTrue(ttran5_clean)

ttranls_clean <- clean_data(TTRANLS[, ncol(TTRANLS)])
ttranls <- zeroToTrue(ttranls_clean)

ttranms_clean <- clean_data(TTRANMS[, ncol(TTRANMS)])
ttranms <- zeroToTrue(ttranms_clean)


#remove the litter from the workspace
rm(TTRAN3, TTRAN4, TTRAN5, TTRANLS, TTRANMS)
rm(ttran3_clean, ttran4_clean, ttran5_clean, ttranms_clean, ttranls_clean)

################################################################################
#######################BOOTSTRAP DISTRIBUTION###################################
################################################################################
nSim <- 10^4

#0. general scheme
#set.seed(123)
#b_name <- boot2samples(mean, sample1, sample2, nSim)

#1.  ttran3 vs ttran4
set.seed(123)
b_ttran3_ttran4 <- boot2samples(mean, ttran3, ttran4, nSim)

#2. ttran3 vs ttran5
set.seed(123)
b_ttran3_ttran5 <- boot2samples(mean, ttran3, ttran5, nSim)

#3. ttran3 vs ttranls
set.seed(123)
b_ttran3_ttranls <- boot2samples(mean, ttran3, ttranls, nSim)

#4. ttran3 vs ttranms
set.seed(123)
b_ttran3_ttranms <- boot2samples(mean, ttran3, ttranms, nSim)


################################################################################
#######################PERMUTATION TEST#########################################
################################################################################
nSim <- 10^4

#0. general scheme
#set.seed(123)
#p_name <- perm2samples(mean, sample1, sample2, nSim)

#1. ttran3 vs ttran4
set.seed(123)
p_ttran3_ttran4 <- perm2samples(mean, ttran3, ttran4, nSim)

#2. ttran3 vs ttran5
set.seed(123)
p_ttran3_ttran5 <- perm2samples(mean, ttran3, ttran5, nSim)

#3. ttran3 vs ttranls
set.seed(123)
p_ttran3_ttranls <- perm2samples(mean, ttran3, ttranls, nSim)

#4. ttran3 vs ttranms
set.seed(123)
p_ttran3_ttranms <- perm2samples(mean, ttran3, ttranms, nSim)

################################################################################
#######################GRAPHICAL ANALYSIS#######################################
################################################################################

my_side <- "RIGHT FLIP"

#########################BOOTSTRAP##############################################
histBootstrap(b_ttran3_ttran4, "TTRAN3", "TTRAN4", my_side)
histBootstrap(b_ttran3_ttran5, "TTRAN3", "TTRAN5", my_side)
histBootstrap(b_ttran3_ttranls, "TTRAN3", "TTRANLS", my_side)
histBootstrap(b_ttran3_ttranms, "TTRAN3", "TTRANMS", my_side)

#########################PERMUTATION TEST#######################################
histPermutation(p_ttran3_ttran4, "TTRAN3", "TTRAN4", my_side)
histPermutation(p_ttran3_ttran5, "TTRAN3", "TTRAN5", my_side)
histPermutation(p_ttran3_ttranls, "TTRAN3", "TTRANLS", my_side)
histPermutation(p_ttran3_ttranms, "TTRAN3", "TTRANMS", my_side)