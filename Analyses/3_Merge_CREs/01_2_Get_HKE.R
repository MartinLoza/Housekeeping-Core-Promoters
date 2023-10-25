# ---
# title: "Get core enhancers and MD"
# author: "Martin loza"
# date: "2022/11/27"
# ---

print("This a code is to get the HKE using the Shirokane server.")
print("This is the second step in the process to merge the enhancers")

library(dplyr)
library(GenomicRanges)
library(here)
library(bedr)

## Global variables
date = "221127"
seed = 777
out_dir = "/Data/merged"
in_dir = "/Data/raw_enhancers/merged/"

## Functions

### Global functions
source("Analysis/Functions.R")

### Local functions
source("Analysis/Functions_merge.R")

print("Load raw_filtered from CTSE")
raw_filtered <- readRDS(file = paste0(out_dir,"/raw_filtered_ctse_221126.rds"))

cat("\nNumber of enhancers in raw_filtered: \n", nrow(raw_filtered))

print("Step 2. HKE") 
print("We would like to get HKE are defined using a double-merge")

print("First overlap")
# In the first overlap, we strictly define HKE as those regions that completely overlap across cell types.
# We do it in the "raw_filtered" object, e.g. CTSE related regions were already filtered.

#Get merged enhancers
me = MergeRegions_strict(regions_df = raw_filtered, start = "start", end = "end", verbose = FALSE)

print("Second overlap") 
print("Relax the definition of HKE with 90%. Merge overlapping regions at least in 90 % of the total cell types.") 
hke <- RelaxHKE(merged_enhancers = me, per = 0.9, distance = 20)
#add n_celltype information
hke <- hke %>% mutate(n_celltypes = length(unique(raw_filtered$CellType)))

### save hke 
saveRDS(object = hke, file = paste0(out_dir,"/backup/hke_50cts_before_transfer_md_",date,".rds"))

###TEST TEST TEST
# print("Sampling HKE for test")
# set.seed(seed)
# hke <- hke %>% sample_n(size = 50, replace = FALSE)

print("Setup regions")
hke <- SetupRegions(regions = hke)

print("Get the meta data of the HKE (using 10% of overlap).")
# This is an iterative process not paralelizable.
hke_expanded <- TransferMD(core_enhancers = hke, raw_enhancers = raw_filtered)

print("Add raw_CTSE information")
hke_expanded$raw_CTSE <- FALSE

print("Find HKE and raw overlapping regions and filter them from the raw dataset")
rmv_idx <- WhichOverlapEnhancer(query_enhancer = hke, subject_enhancer = raw_filtered, min_overlap = 0)
raw_filtered <- raw_filtered[-rmv_idx,]

print("Save hke and raw_filtered")
saveRDS(object = hke_expanded, file = paste0(out_dir,"/hke_expanded_",date,".rds"))
saveRDS(object = raw_filtered, file = paste0(out_dir,"/raw_filtered_hke_",date,".rds"))
