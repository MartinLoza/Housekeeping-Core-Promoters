# ---
# title: "Get core enhancers and MD"
# author: "Martin loza"
# date: "2022/11/27"
# ---

print("This a code is to get the others (not CTSE nor HKE) using a cloud computing service.")
print("This is the second step in the process to merge the enhancers")

library(dplyr)
library(GenomicRanges)
library(bedr)

## Global variables
date = "221127"
seed = 777
in_dir = "/Data/merged"
out_dir = in_dir

## Functions

### Global functions
source("Analysis/Functions.R")

### Local functions
source("Analysis/Functions_merge.R")

print("Load raw_filtered from HKE")
raw_filtered <- readRDS(file = paste0(out_dir,"/raw_filtered_hke_221127.rds"))

print("Number of enhancers in raw_filtered: \n")
nrow(raw_filtered)

print("Step 3. Loop over the rest of enhancers") 

###TEST TEST
# print("Test. Sampling raw_filtered data")
# set.seed(seed)
# raw_filtered <- raw_filtered %>% sample_n(size = 30000)

print("LOOP")
# Now, we want to get just the "core" enhancer for each number of cell types. 
# For this, we would like to use the strict-overlap.
# The problem is that this overlap can create a large number of short regions that would be difficult to filter later.
# Then, I would like to loop over a decreasing number of celltypes,
#   get the strict overlap, select the core enhancers for the max number of cell types,
#   and filter overlapping regions from the "raw_filtered" dataset.

#loop for data of less than 90% of cell types
for(n in round(x = 50*0.9, digits = 0):1){
  
  print("Iteration: ")
  print(n)
  
  ########
  ### Strict overlap
  ########
  # Get strict overlap of filtered raw regions 
  me = MergeRegions_strict(regions_df = raw_filtered, start = "start", end = "end", verbose = FALSE)
  
  ########
  ### Get core enhancers of the maximun number of cell types
  ########
  max_cts <- max(me$n_celltypes)
  print(paste0("Maximum number of cell types: ",max_cts))
  core <- me %>% filter(n_celltypes == max_cts)
  
  ########
  ### Merge core enhancers
  ########
  # We would like core enhancer closer than a given threshold to be merged.
  # For now I will use 20 bps.
  core <- bedr.merge.region(x = core, distance = 20, verbose = FALSE)
  colnames(core) <- c("chr", "start", "end", "length")
  core <- core %>% mutate(n_celltypes = max_cts)
  
  ###TEST TEST
  # print("Test. Sampling the strict merge enhancers.")
  # set.seed(seed)
  # core <- core %>% sample_n(size = min(nrow(core), 100), replace = FALSE)
  
  #set up region
  core <- SetupRegions(regions = core)
  
  ########
  ## Add MD
  ########
  core_expanded <- TransferMD(core_enhancers = core, raw_enhancers = raw_filtered)
  
  ########
  ## Update merged_expanded data frame
  ########
  if(!exists("merged_expanded")){
    merged_expanded <- core_expanded
  }else{
    merged_expanded <- rbind(merged_expanded, core_expanded)
  }
  
  ########
  ### Find core and raw_filtered overlapping regions and filter them from the raw_filtered dataset
  ########
  rmv_idx <- WhichOverlapEnhancer(query_enhancer = core, subject_enhancer = raw_filtered, min_overlap = 0)
  raw_filtered <- raw_filtered[-rmv_idx,]
  
  #Save tmp data
  print("Saving temporal data.")
  saveRDS(object = raw_filtered, file = paste0(out_dir, "/tmp/raw_filtered_", n,"_",date,".rds"))
  saveRDS(object = merged_expanded, file = paste0(out_dir, "/tmp/me_expanded_", n,"_",date,".rds"))
}


