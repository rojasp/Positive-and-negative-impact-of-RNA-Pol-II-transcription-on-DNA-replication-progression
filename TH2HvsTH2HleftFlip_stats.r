################################################################################
#Project: genetics breakdown
#File: TH2HvsTH2HleftFlip_stats.r
#Author: Giovanni Guglielmi
#Description: perform the analysis on the left reads
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

TH2H3 <- read.csv(
  "2.datasets_panel/ES3_flip/ES3_peaks_transcribed_genes_h2h_distance_breakpoint_l_flip",
  header = F, sep = "\t")

TH2H4 <- read.csv(
  "2.datasets_panel/ES4_flip/ES4_peaks_h2h_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

TH2H5 <- read.csv(
  "2.datasets_panel/ES5_flip/ES5_peaks_h2h_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

TH2HLS <- read.csv(
  "2.datasets_panel/LS_flip/LS_peaks_h2h_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

TH2HMS <- read.csv(
  "2.datasets_panel/MS_flip/MS_peaks_h2h_genes_distance_breakpoint_l_flip",
  header = F, sep = "\t")

###############################################################################
############################CLEAN DATA#########################################
###############################################################################
th2h3_clean <- clean_data(TH2H3[, ncol(TH2H3)])
th2h3 <- zeroToTrue(th2h3_clean)

th2h4_clean <- clean_data(TH2H4[, ncol(TH2H4)])
th2h4 <- zeroToTrue(th2h4_clean)

th2h5_clean <- clean_data(TH2H5[, ncol(TH2H5)])
th2h5 <- zeroToTrue(th2h5_clean)

th2hls_clean <- clean_data(TH2HLS[, ncol(TH2HLS)])
th2hls <- zeroToTrue(th2hls_clean)

th2hms_clean <- clean_data(TH2HMS[, ncol(TH2HMS)])
th2hms <- zeroToTrue(th2hms_clean)


#remove the litter from the workspace
rm(TH2H3, TH2H4, TH2H5, TH2HLS, TH2HMS)
rm(th2h3_clean, th2h4_clean, th2h5_clean, th2hms_clean, th2hls_clean)

################################################################################
#######################BOOTSTRAP DISTRIBUTION###################################
################################################################################
nSim <- 10^4

#0. general scheme
#set.seed(123)
#b_name <- boot2samples(mean, sample1, sample2, nSim)

#1.  th2h3 vs th2h4
set.seed(123)
b_th2h3_th2h4 <- boot2samples(mean, th2h3, th2h4, nSim)

#2. th2h3 vs th2h5
set.seed(123)
b_th2h3_th2h5 <- boot2samples(mean, th2h3, th2h5, nSim)

#3. th2h3 vs th2hls
set.seed(123)
b_th2h3_th2hls <- boot2samples(mean, th2h3, th2hls, nSim)

#4. th2h3 vs th2hms
set.seed(123)
b_th2h3_th2hms <- boot2samples(mean, th2h3, th2hms, nSim)


################################################################################
#######################PERMUTATION TEST#########################################
################################################################################
nSim <- 10^4

#0. general scheme
#set.seed(123)
#p_name <- perm2samples(mean, sample1, sample2, nSim)

#1. th2h3 vs th2h4
set.seed(123)
p_th2h3_th2h4 <- perm2samples(mean, th2h3, th2h4, nSim)

#2. th2h3 vs th2h5
set.seed(123)
p_th2h3_th2h5 <- perm2samples(mean, th2h3, th2h5, nSim)

#3. th2h3 vs th2hls
set.seed(123)
p_th2h3_th2hls <- perm2samples(mean, th2h3, th2hls, nSim)

#4. th2h3 vs th2hms
set.seed(123)
p_th2h3_th2hms <- perm2samples(mean, th2h3, th2hms, nSim)

################################################################################
#######################GRAPHICAL ANALYSIS#######################################
################################################################################

my_side <- "LEFT FLIP"

#########################BOOTSTRAP##############################################
histBootstrap(b_th2h3_th2h4, "TH2H3", "TH2H4", my_side)
histBootstrap(b_th2h3_th2h5, "TH2H3", "TH2H5", my_side)
histBootstrap(b_th2h3_th2hls, "TH2H3", "TH2HLS", my_side)
histBootstrap(b_th2h3_th2hms, "TH2H3", "TH2HMS", my_side)

#########################PERMUTATION TEST#######################################
histPermutation(p_th2h3_th2h4, "TH2H3", "TH2H4", my_side)
histPermutation(p_th2h3_th2h5, "TH2H3", "TH2H5", my_side)
histPermutation(p_th2h3_th2hls, "TH2H3", "TH2HLS", my_side)
histPermutation(p_th2h3_th2hms, "TH2H3", "TH2HMS", my_side)