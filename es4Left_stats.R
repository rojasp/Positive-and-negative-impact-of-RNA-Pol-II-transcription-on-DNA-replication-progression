################################################################################
#Project: genetics breakdown
#File: es4Left_stats.r
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

NT <- read.csv(
  "2.datasets_panel/ES4/ES4_peaks_nontranscribed_genes_distance_breakpoint_l",
  header = F, sep = "\t")

TCOD <- read.csv(
  "2.datasets_panel/ES4/ES4_peaks_transcribed_genes_cod_distance_breakpoint_l",
  header = F, sep = "\t")

TH2H <- read.csv(
  "2.datasets_panel/ES4/ES4_peaks_transcribed_genes_h2h_distance_breakpoint_l",
  header = F, sep = "\t")


TTRAN <- read.csv(
  "2.datasets_panel/ES4/ES4_peaks_transcribed_genes_tran_distance_breakpoint_l",
  header = F, sep = "\t")


###############################################################################
############################CLEAN DATA#########################################
###############################################################################
nt_clean <- clean_data(NT[, ncol(NT)])
nt <- zeroToTrue(nt_clean)

tcod_clean <- clean_data(TCOD[, ncol(TCOD)])
tcod <- zeroToTrue(tcod_clean)

th2h_clean <- clean_data(TH2H[, ncol(TH2H)])
th2h <- zeroToTrue(th2h_clean)

ttran_clean <- clean_data(TTRAN[, ncol(TTRAN)])
ttran <- zeroToTrue(ttran_clean)

#remove the litter from the workspace
rm(NT, TCOD, TH2H, TTRAN)
rm(nt_clean, tcod_clean, th2h_clean, ttran_clean)

################################################################################
#######################BOOTSTRAP DISTRIBUTION###################################
################################################################################

#0. general scheme
#set.seed(123)
#b_name <- boot2samples(mean, sample1, sample2)

#1.  nt vs tcod
set.seed(123)
b_nt_tcod <- boot2samples(mean, nt, tcod)

#2. nt vs th2h
set.seed(123)
b_nt_th2h <- boot2samples(mean, nt, th2h)

#3. nt vs ttran
set.seed(123)
b_nt_ttran <- boot2samples(mean, nt, ttran)

#4. tcod vs th2h
set.seed(123)
b_tcod_th2h <- boot2samples(mean, tcod, th2h)

#5. tcod vs ttran
set.seed(123)
b_tcod_ttran <- boot2samples(mean, tcod, ttran)

#6. th2h vs ttran
set.seed(123)
b_th2h_ttran <- boot2samples(mean, th2h, ttran)

################################################################################
#######################PERMUTATION TEST#########################################
################################################################################

#0. general scheme
#set.seed(123)
#p_name <- perm2samples(mean, sample1, sample2)


#1. nt vs tcod
set.seed(123)
p_nt_tcod <- perm2samples(mean, nt, tcod)


#2. nt vs th2h
set.seed(123)
p_nt_th2h <- perm2samples(mean, nt, th2h)

#3. nt vs ttran
set.seed(123)
p_nt_ttran <- perm2samples(mean, nt, ttran)

#4. tcod vs th2h
set.seed(123)
p_tcod_th2h <- perm2samples(mean, tcod, th2h)

#5. tcod vs ttran
set.seed(123)
p_tcod_ttran <- perm2samples(mean, tcod, ttran)

#6. th2h vs ttran
set.seed(123)
p_th2h_ttran <- perm2samples(mean, th2h, ttran)

################################################################################
#######################GRAPHICAL ANALYSIS#######################################
################################################################################

my_side <- "LEFT 4"

#########################BOOTSTRAP##############################################
histBootstrap(b_nt_tcod, "NT", "TCOD", my_side)
histBootstrap(b_nt_th2h, "NT", "TH2H", my_side)
histBootstrap(b_nt_ttran, "NT", "TTRAN", my_side)

histBootstrap(b_tcod_th2h, "TCOD", "TH2H", my_side)
histBootstrap(b_tcod_ttran, "TCOD", "TTRAN", my_side)

histBootstrap(b_th2h_ttran, "TH2H", "TTRAN", my_side)
#########################PERMUTATION TEST#######################################
histPermutation(p_nt_tcod, "NT", "TCOD", my_side)
histPermutation(p_nt_th2h, "NT", "TH2H", my_side)
histPermutation(p_nt_ttran, "NT", "TTRAN", my_side)

histPermutation(p_tcod_th2h, "TCOD", "TH2H", my_side)
histPermutation(p_tcod_ttran, "TCOD", "TTRAN", my_side)

histPermutation(p_th2h_ttran, "TH2H", "TTRAN", my_side)