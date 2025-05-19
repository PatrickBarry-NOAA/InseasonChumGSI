#### SW Analysis Control File ####
library(stringr)
library(tidyverse)

#Set the Current Statistical Week
CrntSW <- 32

#Check to see if we are combining SW
if(nchar(CrntSW)>2){
       CmbWk<-T
       }else { CmbWk<-F
       }

if(CmbWk==T){SW2Cmb<-str_extract_all(CrntSW, "[0-9]{2}") %>% unlist()}


CrntYr <- format(Sys.Date(),"%Y") #Get the current year

#Run the GSI analysis on the stat week
source('SW_GSI.R')

#Render the report
rmarkdown::render('InseasonGSIAnalysis.Rmd',
                  output_file = paste('./Output/',CrntYr,'/',CrntSW,'/',"SW",CrntSW,"_", Sys.Date(), 
                                      '.pdf', sep=''))