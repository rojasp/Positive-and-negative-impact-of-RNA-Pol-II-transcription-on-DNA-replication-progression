##############################
# 0 - Load libraries
##############################
library(tidyverse)

############################## 
# 1 - Source file 
##############################
# load breakpoints file (Cosmic.csv)
cosmic = read.csv("Cosmic.csv",header = T, sep = "\t")

#####################################################
# 2- Create two windows (left and right)
#####################################################

break_point_l <- cosmic
break_point_l  <- break_point_l  %>% 
  select (Sample.name, Mutation.Type , Mutation.ID,  Chrom.From, Location.From.min, Location.From.max,length )   %>%  
  mutate(
    Location.From.min = Location.From.min -1, 
    Location.From.max= Location.From.max+1,  
    Chrom.From = paste("chr",  Chrom.From, sep=""))              
write.table(break_point_l , "break_point_left", sep = "\t")

break_point_r <- cosmic
break_point_r <- break_point_r  %>% 
  select (Sample.name,  Mutation.Type ,Mutation.ID, Chrom.To, Location.To.min, Location.To.max,length )   %>%  
  mutate(
    Location.To.min = Location.To.min -1, 
    Location.To.max= Location.To.max+1,  
    Chrom.To = paste("chr",  Chrom.To, sep=""))              
write.table(break_point_r , "break_point_rigth", sep = "\t")


