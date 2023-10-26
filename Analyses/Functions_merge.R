##Check enhancers
CheckEnhancers <- function(enhancers = NULL){
  #check length
  cat("\nLength OK: \n", 
      sum(enhancers$length == (enhancers %>% mutate(t = end - start) %>% pull(t))) == nrow(enhancers))
  #check center
  cat("\nCenter OK: \n", 
      sum(enhancers$center == (enhancers %>% mutate(t = (start + floor(length/2))) %>% pull(t))) == nrow(enhancers))
  #check chromosome
  cat("\nChromosome OK: \n", 
      sum(enhancers$chr == factor(x = enhancers$chr, levels = paste0("chr", c(1:22,"X")))) == nrow(enhancers))
  #check distance
  cat("\nDistance OK: \n", 
      sum(enhancers$distance == (enhancers %>% mutate(t = abs(tss - center)) %>% pull(t))) == nrow(enhancers))
  #check overlapping_TSS
  cat("\nEnhancer-TSS overlapping OK: \n", 
      sum(enhancers %>% 
            mutate(t = if_else(condition = (tss >= start & tss <= end), 
                               true = TRUE, false = FALSE)) %>% pull(t)) == 0)
}

#setup regions to transfer MD
SetupRegions <- function(regions = NULL){
  
  # Get the length (sometimes we miss one base)
  regions <- regions %>% mutate(length = end - start) 
  # Get center of region
  regions <- regions %>% mutate(center = start + floor(length/2))
  # Create region id for regions
  regions <- regions %>% mutate(region_id = paste0(chr,":",start,"-",end))
  # Set chromosome information as factor
  chr_levels <- paste0("chr", c(1:22,"X"))
  regions <- regions %>% mutate(chr = factor(x = chr, levels = chr_levels))
  
  return(regions)
}

## Transfer MD from raw enhancers to core enhancers
AddMD <- function(overlapping_ids = NULL, raw_data = NULL, current_region = NULL){
  
  # initialize meta data in current_merged_data
  new_colnames <- c("n_genes", "gene", "tss", "distance", "signed_distance", "genes_all", "celltypes",
                    "test_n_celltypes","n_regions", "mean_ABCScores",
                    "Median_ABCScores","Fst_quantile_ABCScores", "Trd_quantile_ABCScores",
                    "min_ABCScores", "max_ABCScores", "raw_CTSE",
                    "filtered_10prc_enhancers", "filtered_10prc_ids")
  new_columns <- data.frame(matrix(data = NA, ncol = length(new_colnames))) 
  colnames(new_columns) <- new_colnames
  
  #Get the meta data and fill out the data frame
  current_region <- cbind(current_region, new_columns)
  
  #Get the overlapping regions
  overlapping_raw <- raw_data[overlapping_ids,]
  #statistics for ABCscores
  stat_abc <- quantile(x = overlapping_raw$ABC.Score, probs = c(0.25,0.5,0.75), na.rm = TRUE)
  
  #get and set the missing values in current_merged_data
  ## genes info
  current_region["n_genes"] <- length(unique(overlapping_raw$TargetGene))
  current_region["genes_all"] <-  paste(unique(overlapping_raw$TargetGene), collapse = "/")
  ## celltypes info
  current_region["celltypes"] <- paste(unique(overlapping_raw$CellType), collapse = "/")
  current_region["test_n_celltypes"] <- length(unique(overlapping_raw$CellType))
  ## regions info
  current_region["n_regions"] <- length(overlapping_ids)
  ## ABCscores info
  current_region["mean_ABCScores"] <- mean(overlapping_raw$ABC.Score, na.rm = TRUE)
  current_region["Median_ABCScores"] <- stat_abc["50%"]
  current_region["Fst_quantile_ABCScores"] <- stat_abc["25%"]
  current_region["Trd_quantile_ABCScores"] <- stat_abc["75%"]
  current_region["min_ABCScores"] <- min(overlapping_raw$ABC.Score, na.rm = TRUE)
  current_region["max_ABCScores"] <- max(overlapping_raw$ABC.Score, na.rm = TRUE)
  
  ################
  #Gene loop
  ################
  
  # I want to get the information for each target gene
  genes = unique(overlapping_raw$TargetGene)
  # Prepare the current region data frame. Duplicate the info by the number of genes
  current_region <- current_region[rep(x = 1, times = length(genes)),]
  # For each target gene
  for(i in seq_len(length(genes))){
    #get current gene
    g = genes[i]
    #get transcription start site
    tss <- overlapping_raw %>% filter(TargetGene == g) %>% pull(tss) %>% unique()
    current_region[i,"gene"] <- g
    current_region[i,"tss"] <- tss
    current_region[i,"distance"] <- abs(tss - current_region[i,"center"])
    current_region[i,"signed_distance"] <- (tss - current_region[i,"center"])
  }
  
  return(current_region)
}

## Transfer MD from raw enhancers to core enhancers
TransferMD <- function(core_enhancers = NULL, raw_enhancers = NULL){
  #for each region in the core enhancers 
  for(r in seq_len(nrow(core_enhancers))){
    
    # print(r)
    #get current region
    current_region = core_enhancers[r,]
    #get raw enhancers overlapping with the current region
    ## overlapping with at least 10%
    overlapping_ids <- WhichOverlapEnhancer(query_enhancer = current_region,
                                            subject_enhancer = raw_enhancers,
                                            min_overlap = floor(current_region$length*0.1))
    
    #Transfer meta data
    current_region <- AddMD(overlapping_ids = overlapping_ids, 
                            raw_data = raw_enhancers, current_region = current_region )
    
    #get test idx with minimum of 0 bases overlapp (even 1 base of overlap is OK)
    t_idx <- WhichOverlapEnhancer(query_enhancer = current_region,
                                  subject_enhancer = raw_enhancers,
                                  min_overlap = 0)
    
    if(length(t_idx) != length(overlapping_ids)){
      print(r)
      #save info for debug
      current_region["filtered_10prc_enhancers"] = TRUE
      #get the idx which are different
      diff_idx <- t_idx[which(! overlapping_ids %in% t_idx)]
      current_region["filtered_10prc_ids"] = paste(raw_filtered[diff_idx,"region_id"], collapse = "/")
    }
    
    #merge expanded CTSE enhancers
    if(r == 1){
      region_expanded = current_region
    }else{
      region_expanded <- rbind(region_expanded,current_region)
    }
  }
  return(region_expanded)
}

# Merge regions, non strict
MergeRegions <- function(regions_df = NULL, start = "start", end = "end"){
  #Reduce regions using Granges
  gr <- ReduceRegions(regions_df = regions_df)
  #Setup reduced regions to add MD
  merged_data <- data.frame("chr" = factor(seqnames(gr), levels = paste0("chr",c(1:22,"X"))),
                            GenomicRanges::ranges(gr))
  colnames(merged_data) <- c("chr", "start", "end", "length")
  #Pass the MD from the original regions to the merged ones
  merged_data <- ReduceRegion_addMD(from_df = regions_df, to_df = merged_data)
  #Return the merged data
  return(merged_data)
}

#Merged regions using Granges. The output already contains the MD
MergeRegions_strict <- function(regions_df = NULL, start = "start", end = "end", verbose = FALSE, build = "hg38"){
  #Setup the regions to use bedr
  #Set the regions in the bedr format
  regions_df <- regions_df %>% mutate(CellType = factor(CellType))
  regions_ls <- regions_df %>% select(all_of(c("chr", start, end))) %>%  split(f =  regions_df$CellType)
  
  #In previous test I found that there are problems if we used the cell types names. Then, for now I will just create dummy names
  #Create dummy names
  dummy_names <- paste0("ct_", 1:length(regions_ls))
  names(dummy_names) <- names(regions_ls)
  #Assign dummy names
  names(regions_ls) <- dummy_names
  
  ##Make sure chromosome information is a character vector. It will give an error if the chromosome information is a factor.
  regions_ls <- lapply(X = regions_ls, FUN = function(r){
    r[["chr"]] <- as.character(r[["chr"]])
    return(r)
  })
  
  #Merge and sort region across cell types to avoid errors when joining with bedr
  regions_ls <- lapply(X = regions_ls, FUN = bedr.merge.region, verbose = verbose)
  regions_ls <- lapply(X = regions_ls, FUN = bedr.sort.region, verbose = verbose)
  #Join multile using bedr. This is equivalent to the strict merge we want
  merged_regions <- bedr.join.multiple.region(x = regions_ls,
                                              check.sort = FALSE,
                                              check.chr = FALSE, 
                                              check.valid = FALSE, 
                                              check.merge = FALSE,
                                              build = build,
                                              verbose = verbose)
  
  #Setup the output to match the standard output used until now
  merged_regions <- merged_regions %>% 
    mutate(chr = V1, start = V2, end = V3) %>% 
    mutate(length = end - start,
           n_celltypes = merged_regions$n.overlaps,
           dummy_celltypes = merged_regions$names) %>% 
    select(chr, start, end, length, n_celltypes, dummy_celltypes)
  
  #Get the original cell types
  merged_regions <- merged_regions %>% mutate(celltypes = dummy_celltypes)
  # we go inverse to avoid incorrect detection, e.g. detedt ct11 as ct1
  for (i in length(dummy_names):1){
    merged_regions$celltypes <- stringr::str_replace(string = merged_regions$celltypes,
                                                     pattern = dummy_names[i], 
                                                     replacement = names(dummy_names)[i])
  }
  
  #Contiguous regions sometimes the overlap in the start and end, then I will add one base
  #to the start of each merged regions
  merged_regions <- merged_regions %>%
    mutate(start = (merged_regions$start + 1)) %>% #add one to the start 
    mutate(length = end - start) # update length
  
  #Change  n_celltypes from character to integer
  merged_regions <- merged_regions %>% mutate(n_celltypes = as.numeric(n_celltypes))
  
  return(merged_regions)
}

#Merge HKE after strict overlap "MergeRegions_strict" function
RelaxHKE <- function(merged_enhancers = NULL, per = 0.8, distance = 20){
  #get the max number of cell types
  n_cts = max(merged_enhancers$n_celltypes)
  cat("Number of cell types used as HKE:\n", round(x = n_cts*per, digits = 0), "-", n_cts)
  #setup for merge using bedr
  hke <- merged_enhancers %>%
    filter(n_celltypes >= round(x = n_cts*per, digits = 0)) %>% 
    select(chr, start, end)
  #merge with bedr
  hke <- bedr.merge.region(x = hke, verbose = FALSE, distance = distance)
  hke <- hke %>% mutate(length = end - start)
  return(hke)
}

