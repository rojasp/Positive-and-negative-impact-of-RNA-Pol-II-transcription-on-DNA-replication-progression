********************************************************
********************************************************
Project:“Positive and negative impact of RNA Pol II transcription on DNA replication progression”

File: README.txt

Author: 
Patricia Rojas (Leader) - University of Birmingham

Giovanni Guglielmi      - University of Melbourne, University of Birmingham
		
Description: This documentation is principally written to support the M&M published in the paper 
             “Positive and negative impact of RNA Pol II transcription on DNA replication  progression”. 
	     This spcripts and pipeline users should be able to reporcue the results. 
	     
       

**************
*Requirements*
**************
    I.   Unix-like operating system (Linux, macOS, etc) 
    II.  The R Project for Statistical Computing
    III. Python 
    IV.  All the file should be contained in the same directory
    
    
**************
*Workflow*
**************
   Please, download the complete layout (workflow.pdf) to follow the whole process 

***********************	
*Workflow for mRNA seq*
***********************
    i.   FastQC 0.11.9
    ii.  BBTools (bbduk)
    iii. STAR  v.020201
    iv.  SAMtools v.1.4 
    v.   featureCounts
    vi.  DESeq2
**********************************************************  
*Workflow for peak calling and breakpoints identification*
**********************************************************
    a. deepTools 
    b. MACS2 v.2.1.0
    c. Bedtool 2.28.0
    d. BEDOPS v2.4.39

*************************************************
*Workflow for comparing the different conditions*
*************************************************
	1. In the directory, please create a new folder named "my_functionR", which will contain:
	   boot2sample.r, clean_data.r, histBootstrap.r, histPermutation.r, perm2Samples.r, resultsPlotMain.r, volcano.r, zeroToTrue.r
	
	2. Please, deposit the following file in the directory. The different comparisons can be run by using all the following scripts. The file should be run through R software.
       
	   ***ANALYSES OF DIFFERENT CONDITIONS***
	   -es3Left_stats.r
	   -es3Right_stats.r
	   -es3LeftFlip_stats.r
	   -es3RightFlip_stats.r
	   -es4Left_stats.r
	   -es4Right_stats.r
	   -es4LeftFlip_stats.r
	   -es4RightFlip_stats.r
	   -es5Left_stats.r
	   -es5Right_stats.r
	   -es5LeftFlip_stats.r
	   -es5RightFlip_stats.r
	   -LSLeft_stats.r
	   -LSRight_stats.r
	   -LSLeftFlip_stats.r
	   -lsRightFlip_stats.r
	   -MSLeft_stats.r
	   -MSRight_stats.r
	   -MSLeftFlip_stats.r
	   -MSRightFlip_stats.r
	   
	   ***ANALYSES OF SAME CONDITION AT DIFFERENT TIME***
	   -NTvsNTleft_stats.r
	   -NTvsNTright_stats.r
	   -NTvsNTleftFlip_stats.r
	   -NTvsNTrightFilp_stats.r
	   -TCODvsTCODleft_stats.r
	   -TCODvsTCODright_stats.r
	   -TCODvsTCODleftFlip_stats.r
	   -TCODvsTCODrightFlip_stats.r
	   -TH2HvsTH2Hleft_stats.r
	   -TH2HvsTH2Hright_stats.r
	   -TH2HvsTH2HleftFlip_stats.r
	   -TH2HvsTH2HrightFlip_stats.r
	   -TTRANvsTTRANleft_stats.r
	   -TTRANvsTTRANright_stats.r
	   -TTRANvsTTRANleftFlip_stats.r
	   -TTRANvsTTRANrightFlip_stats.r
