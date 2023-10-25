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
                    "filtered_10prc_enhancers", "filtered_10prc_ids", 
                    "n_abc_classes", "abc_classes",
                    "n_classes_genic", "classes_genic",
                    "n_classes", "classes",
                    "n_classes_advanced", "classes_advanced")
  new_columns <- data.frame(matrix(data = NA, ncol = length(new_colnames))) 
  colnames(new_columns) <- new_colnames
  
  #Get the meta data and fill out the data frame
  current_region <- cbind(current_region, new_columns)
  
  #Get the overlapping regions
  overlapping_raw <- raw_data[overlapping_ids,]
  
  ################
  # I DONT NEED THIS FILTER ANYMORE
  ################
  # #Filter genes overlapping with the core enhancers
  # ovl_ids <- which(overlapping_raw$tss >= current_region$start &
  #                    overlapping_raw$tss <= current_region$end)
  # #if tss overlapping with current enhancer
  # if(length(ovl_ids) != 0){
  #   #save info for debug
  #   current_region["Enhancer_TSS_overlap"] = TRUE
  #   current_region["ovl_ids"] = paste(overlapping_raw[ovl_ids,"region_id"], collapse = "/")
  #   #filter enhancers overlapping with TSS
  #   overlapping_raw <- overlapping_raw[-ovl_ids,]
  # }
  
  #statistics for ABCscores
  stat_abc <- quantile(x = overlapping_raw$ABC.Score, probs = c(0.25,0.5,0.75), na.rm = TRUE)
  
  #get and set the missing values in current_merged_data
  ## genes info
  current_region["n_genes"] <- length(unique(overlapping_raw$TargetGene))
  current_region["genes_all"] <-  paste(unique(overlapping_raw$TargetGene), collapse = "/")
  ## celltypes info
  current_region["celltypes"] <- paste(unique(overlapping_raw$CellType), collapse = "/")
  current_region["test_n_celltypes"] <- length(unique(overlapping_raw$CellType))
  ## classes info
  current_region["n_abc_classes"] <- length(unique(overlapping_raw$abc_class))
  current_region["abc_classes"] <- paste(unique(overlapping_raw$abc_class), collapse = "/")
  
  current_region["n_classes_genic"] <- length(unique(overlapping_raw$class_genic))
  current_region["classes_genic"] <- paste(unique(overlapping_raw$class_genic), collapse = "/")
  
  current_region["n_classes"] <- length(unique(overlapping_raw$class))
  current_region["classes"] <- paste(unique(overlapping_raw$class), collapse = "/")
  
  current_region["n_classes_advanced"] <- length(unique(overlapping_raw$class_advanced))
  current_region["classes_advanced"] <- paste(unique(overlapping_raw$class_advanced), collapse = "/")
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
    #add missing meta data
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
    
    #Test test test 
    ## overlapping with 0 bases. I am wondering if the non-matching regions are due to the min_overlap condition.
    #Then, I will save this info for future debug.
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