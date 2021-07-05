#!/bin/bash  

#run closest-features from bedops 
date
echo --------------------------START closest-features -----------------


bedops_tools/./closest-features --closest --dist --delim '\t' *.breakpoints
