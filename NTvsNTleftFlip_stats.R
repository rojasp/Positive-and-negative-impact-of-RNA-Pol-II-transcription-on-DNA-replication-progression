################################################################################
#Project: genetics breakdown
#File: NTvsNTleftFlip_stats.r
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

NT3 <- read.csv(
  "2.datasets_panel/ES3_flip/ES3_peaks_notranscribed_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

NT4 <- read.csv(
  "2.datasets_panel/ES4_flip/ES4_peaks_notranscribed_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

NT5 <- read.csv(
  "2.datasets_panel/ES5_flip/ES5_peaks_nontrancribed_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

NTLS <- read.csv(
  "2.datasets_panel/LS_flip/LS_peaks_nontrancribed_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

NTMS <- read.csv(
  "2.datasets_panel/MS_flip/MS_peaks_nontrancribed_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

###############################################################################
############################CLEAN DATA#########################################
###############################################################################
nt3_clean <- clean_data(NT3[, ncol(NT3)])
nt3 <- zeroToTrue(nt3_clean)

nt4_clean <- clean_data(NT4[, ncol(NT4)])
nt4 <- zeroToTrue(nt4_clean)

nt5_clean <- clean_data(NT5[, ncol(NT5)])
nt5 <- zeroToTrue(nt5_clean)

ntls_clean <- clean_data(NTLS[, ncol(NTLS)])
ntls <- zeroToTrue(ntls_clean)

ntms_clean <- clean_data(NTMS[, ncol(NTMS)])
ntms <- zeroToTrue(ntms_clean)


#remove the litter from the workspace
rm(NT3, NT4, NT5, NTLS, NTMS)
rm(nt3_clean, nt4_clean, nt5_clean, ntms_clean, ntls_clean)

################################################################################
#######################BOOTSTRAP DISTRIBUTION###################################
################################################################################
nSim <- 10^4

#0. general scheme
#set.seed(123)
#b_name <- boot2samples(mean, sample1, sample2, nSim)

#1.  nt3 vs nt4
set.seed(123)
b_nt3_nt4 <- boot2samples(mean, nt3, nt4, nSim)

#2. nt3 vs nt5
set.seed(123)
b_nt3_nt5 <- boot2samples(mean, nt3, nt5, nSim)

#3. nt3 vs ntls
set.seed(123)
b_nt3_ntls <- boot2samples(mean, nt3, ntls, nSim)

#4. nt3 vs ntms
set.seed(123)
b_nt3_ntms <- boot2samples(mean, nt3, ntms, nSim)


################################################################################
#######################PERMUTATION TEST#########################################
################################################################################
nSim <- 10^4

#0. general scheme
#set.seed(123)
#p_name <- perm2samples(mean, sample1, sample2, nSim)

#1. nt3 vs nt4
set.seed(123)
p_nt3_nt4 <- perm2samples(mean, nt3, nt4, nSim)

#2. nt3 vs nt5
set.seed(123)
p_nt3_nt5 <- perm2samples(mean, nt3, nt5, nSim)

#3. nt3 vs ntls
set.seed(123)
p_nt3_ntls <- perm2samples(mean, nt3, ntls, nSim)

#4. nt3 vs ntms
set.seed(123)
p_nt3_ntms <- perm2samples(mean, nt3, ntms, nSim)

################################################################################
#######################GRAPHICAL ANALYSIS#######################################
################################################################################

my_side <- "LEFT FLIP"

#########################BOOTSTRAP##############################################
histBootstrap(b_nt3_nt4, "NT3", "NT4", my_side)
histBootstrap(b_nt3_nt5, "NT3", "NT5", my_side)
histBootstrap(b_nt3_ntls, "NT3", "NTLS", my_side)
histBootstrap(b_nt3_ntms, "NT3", "NTMS", my_side)

#########################PERMUTATION TEST#######################################
histPermutation(p_nt3_nt4, "NT3", "NT4", my_side)
histPermutation(p_nt3_nt5, "NT3", "NT5", my_side)
histPermutation(p_nt3_ntls, "NT3", "NTLS", my_side)
histPermutation(p_nt3_ntms, "NT3", "NTMS", my_side)