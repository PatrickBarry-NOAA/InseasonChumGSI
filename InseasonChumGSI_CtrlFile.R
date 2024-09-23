#### SW Analysis Control File ####
library(stringr)
library(tidyverse)
#Set the Current Statistical Week
CrntSW <- 36 
CrntYr <- format(Sys.Date(),"%Y") #Get the current year

#Run the GSI analysis on the stat week
source('SW_GSI.R')

#Render the report
rmarkdown::render('InseasonGSIAnalysis.Rmd',
                  output_file = paste('./Output/',CrntYr,'/',CrntSW,'/',"SW",CrntSW,"_", Sys.Date(), 
                                      '.pdf', sep=''))