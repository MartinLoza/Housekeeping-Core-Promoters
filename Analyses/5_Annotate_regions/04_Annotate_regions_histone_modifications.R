# ---
# title: "Overlap ENCODE histone with ENCODE+NS"
# author: "Martin loza"
# date: "2023/02/15"
# ---

# Following a similar normalization method as SCREEN I normalized the data as log1p+scaling without adding an offset
# In this workflow I would like to annotate the enhancer using this histone data. 

library(dplyr)
library(GenomicRanges)
library(stringr)

## Global variables
date = "230215"
seed = 777
in_dir = "Data/ENCODE_histones_CHIP/"
out_dir = "Data/annotated/"

## Global functions
source("Analysis/Functions.R")

#Local functions
#get overlaps over a list
WhichRegionsOverlapList <- function(subject_data_list = NULL, query_enhancer = NULL, min_overlap = 0){
  subject_data_list <- lapply(X = subject_data_list, FUN = function(data){
    ovl <- WhichRegionsOverlap(query_enhancer = query_enhancer, subject_enhancer = data, min_overlap = min_overlap)
    ovl_idx <- subjectHits(ovl)
    #we remove duplicated rows, this doesn't mean we are eliminating duplicated regions from different cell types
    ovl_idx <- unique(ovl_idx)
    data <- data[ovl_idx,]
    return(data)
  })
  return(subject_data_list)
}

## Load data 

# Load histones data 
histone_data <- readRDS(file = paste0(in_dir, "/histone_data_normalized_230215.rds"))

cat("Number of cell types:\n ")
sapply(X = names(histone_data), FUN = function(histone){
  return(length(unique(histone_data[[histone]]$celltype)))
})

#Load ABC+NS regions 
regions <- readRDS(file = paste0(out_dir,"/merged_u_regions_NS_230210.rds"))
table(regions$type)

#short number of regions  
test <- table(regions$type) == c(225140,10000)
if((!test[1]) | (!test[2])){
  stop(call. = TRUE, "Number of regions do not match!!!!")
}
rm(test)

cat("Expanding regions")
#Following SCREEN workflow we would like to average data from an expanded version of the regions (e.g. =-500 bps)
#to capture the histone marks flanking the cis-Regulatory Elements (CRE)
regions_expanded <- regions %>% mutate(start = start-500, end = end+500)
rm(regions)
#note that we didn't modify the regions_id, so we can always link this annotations with the original data

#short test.
#no overlapping regions after expansion
abc_expanded <- regions_expanded %>% filter(type == "ABC")
ns_expanded <- regions_expanded %>% filter(type == "NS")
ovl <- WhichRegionsOverlap(query_enhancer = abc_expanded,
                           subject_enhancer = ns_expanded,min_overlap = 0)
cat("Number of overlaps: \n", length(subjectHits(ovl)))
if(length(subjectHits(ovl)) !=0){
  stop(call. = TRUE, "NS and ABC regions overlapping. Check the definition of NS.")
}
rm(ovl, abc_expanded, ns_expanded)

#sort regions before overlapping. This step doesn't seem to affect, but just in case...
histone_data <- lapply(X = histone_data, FUN = function(histone){
  #arrange by chr and start
  histone <- histone %>% arrange(chr, start)
  return(histone)
})

#sort regions before overlapping
regions_expanded <- regions_expanded %>% arrange(chr, start)


cat("Subset overlapping histone's peaks")
#to speed up computation we will focus on the histone's peaks overlapping with at least 250bp (half of the expansion)
histone_data <- WhichRegionsOverlapList(subject_data_list = histone_data,
                                        query_enhancer = regions_expanded,
                                        min_overlap = 250) 

cat("Overlap of ABC enhancers with histones data")

#get the minimum on each histone
min_histone <- sapply(X = histone_data, function(histone){
  #get the minimum value of the normalized data
  min_value <- min(histone$normalized)
  return(min_value)
})

histone_names <- names(histone_data)
chromosomes <- levels(regions_expanded$chr)

### TEST TEST TEST
# cat("TEST TEST TEST sampling regions.")
# set.seed(seed)
# regions_expanded <- regions_expanded %>% sample_n(size = 150)
# chromosomes = "chr18"

#for each choromosome
regions_annotated <- lapply(X = chromosomes, FUN = function(ch){
  
  #current regions
  current_regions <- regions_expanded %>% filter(chr == ch) %>% arrange(chr, start)
  #set the rownames as the region_ids to later subsets
  rownames(current_regions) <- current_regions$region_id
  #transform to matrix to speed up access to the data
  #select only numeric numbers to speed up calculations
  current_regions <- as.matrix(current_regions[,c("start","end")])
  #create matrix to store the mean and number of celltypes
  new_data <- matrix(data = 0, nrow = nrow(current_regions), ncol = 2*length(histone_names))
  colnames(new_data) <- paste0(rep(x = histone_names, each = 2), c("", "_n_cts"))
  current_regions <- cbind(current_regions,new_data)
  rm(new_data)
  #create matrix to store the celltypes (characters)
  current_regions_cts <- matrix(data = NA, nrow = nrow(current_regions), ncol = length(histone_names))
  colnames(current_regions_cts) <- paste0(histone_names,"_cts")
  rownames(current_regions_cts) <- rownames(current_regions)
  
  #assing the minimum value for each histone
  for(h in histone_names){
    current_regions[,h] <- min_histone[h]
  }
  
  ## setup histone data
  #to speed up annotation we would like to separate the celltypes from the integer data
  #get the current chromosome data
  current_histone <- lapply(histone_data, FUN = function(h){
    h <- h %>% filter(chr == ch)
    return(h)
  })
  #get current celltypes
  current_histone_cts <- lapply(X = current_histone, FUN = function(h){
    h <- h %>% pull(celltype)
    return(h)
  })
  #get current integer values: start, end, offset_10
  current_histone <- lapply(current_histone, FUN = function(h){
    sel_cols <- c("start", "end", "normalized")
    h <- h[,sel_cols, drop = FALSE]
    h <- as.matrix(h)
    return(h)
  })
  
  ## We have separated integers from characters, this speeds up the computations !!!
  
  #for each histone
  for(histone in histone_names){
    
    cat(paste0("\nCurrent chromosome: ", ch, "\n"))
    print(paste0("Current histone: ", histone))
    
    #We can get the whole set of overlaps first and then just index them
    overlaps <- WhichRegionsOverlap(chr = ch,
                                    query_enhancer = current_regions,
                                    subject_enhancer = current_histone[[histone]],
                                    min_overlap = 250)
    #transform to matrix
    overlaps <- as.matrix(overlaps)
    #change colnames for easy their indexing
    #the original order is "queryHits", and "subjectHits"
    colnames(overlaps) <- c("regions", "histones")
    
    #init progress bar
    pb <- txtProgressBar(min = 0, max = nrow(current_regions), style = 3,  width = 50,  char = "=")
    for(r in 1:nrow(current_regions)){
      
      #subset the overlaps from current regions
      current_ovl <- FilterMatrix(m = overlaps, col = "regions", wt = r)
      
      #if there are overlaps..
      if(length(current_ovl) != 0){
        #get overlapping histone data
        ovl_histone <- current_histone[[histone]][current_ovl[,"histones"], , drop = FALSE]
        #get overlapping histone cts data
        ovl_histone_cts <- current_histone_cts[[histone]][current_ovl[,"histones"]]
        #average of overlapping histone data
        current_regions[r,histone] <- mean(ovl_histone[,"normalized"])
        #number of celltypes in overlapping histone data
        current_regions[r,paste0(histone,"_n_cts")] <- length(unique(ovl_histone_cts))
        #celltypes in overlapping histone data
        current_regions_cts[r,paste0(histone,"_cts")] <- paste(unique(ovl_histone_cts),collapse = "/")
      }
      #progress bar
      setTxtProgressBar(pb, r) 
    }
  }
  
  #short test
  t <- identical(rownames(current_regions), rownames(current_regions_cts))
  if(!t){
    stop(call. = TRUE, "Rowname are not the same... check for some errors.")
  }
  rm(t)
  
  current_regions_df <- cbind(as.data.frame(current_regions), as.data.frame(current_regions_cts))
  
  saveRDS(object = current_regions_df,
          file = paste0(out_dir, "/Histones_intermediate/histones_annotated_",ch,"_", date, ".rds"))
  
  return(current_regions_df)
})

regions_annotated <- Reduce(f = rbind, regions_annotated)
regions_annotated <- regions_annotated %>% mutate(region_id = rownames(regions_annotated))
 
#save final annotated regions
saveRDS(object = regions_annotated, file = paste0(out_dir, "/histones_annotated_", date, ".rds")) 


























