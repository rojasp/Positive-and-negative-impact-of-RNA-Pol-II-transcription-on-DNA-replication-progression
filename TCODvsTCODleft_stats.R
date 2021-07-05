################################################################################
#Project: genetics breakdown
#File: TCODvsTCODleft_stats.r
#Author: Giovanni Guglielmi
#Description: perform the analysis on the left reads
################################################################################

#set the environment
setwd("C:/Users/gguglielmi/OneDrive - The University of Melbourne/uniPhD/projects/computational_genomics")

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

TCOD3 <- read.csv(
  "2.datasets_panel/ES3/ES3_peaks_transcribed_genes_cod_distance_breakpoint_l",
  header = F, sep = "\t")

TCOD4 <- read.csv(
  "2.datasets_panel/ES4/ES4_peaks_transcribed_genes_cod_distance_breakpoint_l",
  header = F, sep = "\t")

TCOD5 <- read.csv(
  "2.datasets_panel/ES5/ES5_peaks_transcribed_genes_cod_distance_breakpoint_l",
  header = F, sep = "\t")

TCODLS <- read.csv(
  "2.datasets_panel/LS/LS_peaks_transcribed_genes_cod_distance_breakpoint_l",
  header = F, sep = "\t")

TCODMS <- read.csv(
  "2.datasets_panel/MS/MS_peaks_transcribed_genes_cod_distance_breakpoint_l",
  header = F, sep = "\t")

###############################################################################
############################CLEAN DATA#########################################
###############################################################################
tcod3_clean <- clean_data(TCOD3[, ncol(TCOD3)])
tcod3 <- zeroToTrue(tcod3_clean)

tcod4_clean <- clean_data(TCOD4[, ncol(TCOD4)])
tcod4 <- zeroToTrue(tcod4_clean)

tcod5_clean <- clean_data(TCOD5[, ncol(TCOD5)])
tcod5 <- zeroToTrue(tcod5_clean)

tcodls_clean <- clean_data(TCODLS[, ncol(TCODLS)])
tcodls <- zeroToTrue(tcodls_clean)

tcodms_clean <- clean_data(TCODMS[, ncol(TCODMS)])
tcodms <- zeroToTrue(tcodms_clean)


#remove the litter from the workspace
rm(TCOD3, TCOD4, TCOD5, TCODLS, TCODMS)
rm(tcod3_clean, tcod4_clean, tcod5_clean, tcodms_clean, tcodls_clean)

################################################################################
#######################BOOTSTRAP DISTRIBUTION###################################
################################################################################
nSim <- 10^5

#0. general scheme
#set.seed(123)
#b_name <- boot2samples(mean, sample1, sample2, nSim)

#1.  tcod3 vs tcod4
set.seed(123)
b_tcod3_tcod4 <- boot2samples(mean, tcod3, tcod4, nSim)

#2. tcod3 vs tcod5
set.seed(123)
b_tcod3_tcod5 <- boot2samples(mean, tcod3, tcod5, nSim)

#3. tcod3 vs tcodls
set.seed(123)
b_tcod3_tcodls <- boot2samples(mean, tcod3, tcodls, nSim)

#4. tcod3 vs tcodms
set.seed(123)
b_tcod3_tcodms <- boot2samples(mean, tcod3, tcodms, nSim)


################################################################################
#######################PERMUTATION TEST#########################################
################################################################################
nSim <- 10^5

#0. general scheme
#set.seed(123)
#p_name <- perm2samples(mean, sample1, sample2, nSim)

#1. tcod3 vs tcod4
set.seed(123)
p_tcod3_tcod4 <- perm2samples(mean, tcod3, tcod4, nSim)

#2. tcod3 vs tcod5
set.seed(123)
p_tcod3_tcod5 <- perm2samples(mean, tcod3, tcod5, nSim)

#3. tcod3 vs tcodls
set.seed(123)
p_tcod3_tcodls <- perm2samples(mean, tcod3, tcodls, nSim)

#4. tcod3 vs tcodms
set.seed(123)
p_tcod3_tcodms <- perm2samples(mean, tcod3, tcodms, nSim)

################################################################################
#######################GRAPHICAL ANALYSIS#######################################
################################################################################

my_side <- "LEFT"

#########################BOOTSTRAP##############################################
histBootstrap(b_tcod3_tcod4, "TCOD3", "TCOD4", my_side)
histBootstrap(b_tcod3_tcod5, "TCOD3", "TCOD5", my_side)
histBootstrap(b_tcod3_tcodls, "TCOD3", "TCODLS", my_side)
histBootstrap(b_tcod3_tcodms, "TCOD3", "TCODMS", my_side)

#########################PERMUTATION TEST#######################################
histPermutation(p_tcod3_tcod4, "TCOD3", "TCOD4", my_side)
histPermutation(p_tcod3_tcod5, "TCOD3", "TCOD5", my_side)
histPermutation(p_tcod3_tcodls, "TCOD3", "TCODLS", my_side)
histPermutation(p_tcod3_tcodms, "TCOD3", "TCODMS", my_side)