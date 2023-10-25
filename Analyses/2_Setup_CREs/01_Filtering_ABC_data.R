# ---
# title: "Filtering ABC data"
# author: "Martin Loza"
# data: "2022/11/11"
# ---

#On this notebook I aim to code filter the ABC data using the ABC.score. 
#I aim to reduce the memory needed when data from a large number of cell types are used in future codes

###############
## Init
###############
library(dplyr)
library(here)

# Global variables
date = "221111"

#Dowload the ABC dataset from the Engreitz lab resources website: https://www.engreitzlab.org/resources

#change/add your Data directory to the ABC datasets.
in_dir <- "Data/Nasser/ABC/All_predictions/Predictions/" 
out_dir <- "Data/raw_enhancers/"

###############
## Load data
###############

#Get the folders names
folders <- list.files(path = in_dir)
cat("Number of available cell types: \n", length(folders))

###############
## Loop
###############

#for each dataset
for(i in seq_along(folders)){
  
  ###############
  # Load dataset 
  ###############
  d <- paste0(in_dir,"/",folders[i],"/Predictions_AvgHiC/EnhancerPredictionsAllPutative.txt")
  data <- as.data.frame(read.table(file = d, header = TRUE, sep = '\t'))
  
  ###############
  # Filtering 
  ###############
  
  ### Filter enhancers which have an ABCscore lower than 0.015 (score used in Nasser)
  data <- data %>% filter(ABC.Score >= 0.015)
  
  ###Starting sites are always before the end sites (see Test 1 below.). Then, we can check the elements that overlaps the TSS using the starting and ending site . e.g. overlapping -> start_site < TSS_site < end_site
  data["hg19_overlap_TSS"] <- FALSE 
  idx <- which(data$start < data$TargetGeneTSS & data$end > data$TargetGeneTSS)
  data$hg19_overlap_TSS[idx] <- TRUE
  
  ###############
  # Save data
  ###############
  
  saveRDS(object = data, file = paste0(out_dir,"/",folders[i],"_filtered_",date,".rds" ) )
}
