# ---
# title: "Get core enhancers and MD"
# author: "Martin loza"
# date: "2022/11/26"
# ---

# This a code to get the CTSE using in a clud server.
#This is the first step in the process to merge the enhancers

library(dplyr)
library(GenomicRanges)
library(here)
library(bedr)

## Global variables
date = "221126"
seed = 777
out_dir = "/Data/merged"
in_dir = "/Data/raw_enhancers/merged"

## Functions

### Global functions
source(here("Analysis/Functions.R"))

### Local functions
source("Analysis/Functions_merge.R")

## Load enhancers data 

#Load raw enhancers from 50 celltypes
raw <- readRDS(file = paste0(in_dir, "/hg38_annotated_simple_221126.rds"))


#Check raw data
CheckEnhancers(enhancers = raw)
length(unique(raw$CellType))
unique(raw$CellType)
colnames(raw)


#Everything looks good.

## Set up

#I want to have a region id in the raw data for downstream analysis and comparisons
raw <- raw %>% mutate(region_id = paste0(chr,":",start,"-",end))

## Get core enhancers 

## Init

#Init raw_filtered. I will filter raw enhancer overlapping with core enhancers. I could use the original raw df, but I prefer to keep it in case I want to verify something from the original raw enhancers.
raw_filtered = raw


## Step 1. CTSE
# We need to:
# - Get CTSE using the non-strict merge of enhancers.
# - Get the meta data of the CTSE (using 10% of overlap).
# - Then, filter the raw regions overlapping with the CTSE (using 10% of overlap).

#TEST TEST TEST
# raw <- raw[1:1000,]

### Get non-strict merged enhancers, select CTSE
merged_ns <- MergeRegions(regions_df = raw, start = "start", end = "end")
ctse <- merged_ns %>% filter(n_celltypes == 1)

#save ctse and merged data
saveRDS(object = ctse, file = paste0(out_dir,"/backup/ctse_before_transfer_MD_",date,".rds"))
saveRDS(object = merged_ns, file = paste0(out_dir,"/backup/merged_ns_",date,".rds"))

#load ctse and merged data
merged_ns <- readRDS(file = paste0(out_dir,"/backup/merged_ns_",date,".rds"))
ctse <- readRDS(file = paste0(out_dir,"/backup/ctse_before_transfer_MD_",date,".rds"))

table(merged_ns$n_celltypes)
### Setup regions
ctse <- SetupRegions(regions = ctse)

#TEST TEST TES
# ctse <- ctse %>% sample_n(size = 50)

### Get the meta data of the CTSE (using 10% of overlap).
#This is an iterative process not paralelizable.
ctse_expanded = TransferMD(core_enhancers = ctse, raw_enhancers = raw_filtered)

#Add raw_CTSE information
ctse_expanded$raw_CTSE <- TRUE

### Find CTSE and raw overlapping regions and filter them from the raw dataset
rmv_idx <- WhichOverlapEnhancer(query_enhancer = ctse, subject_enhancer = raw_filtered, min_overlap = 0)
raw_filtered <- raw_filtered[-rmv_idx,]

### Save ctse and raw_filtered
saveRDS(object = ctse_expanded, file = paste0(out_dir,"/ctse_expanded_",date,".rds"))
saveRDS(object = raw_filtered, file = paste0(out_dir,"/raw_filtered_ctse_",date,".rds"))
