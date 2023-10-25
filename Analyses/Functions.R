# Functions used on this project

## Function name: SummarizeData
## input:     df,   Data's data frame
## output:    Printed summary of different features of the input data
# SummarizeData <- function(df = NULL){
#   cat("Number of E-G regions:", nrow(df))
#   cat("\nNumber of cell types: ", length(table(df$celltype)))
#   cat("\nAnalysis by cell type: \n")
#   SummarizeDataByCelltype(df = df)
# }

## Function name: SummarizeDataByCelltype
## input:     df,   Data's data frame
## output:    Printed summary of different features of the input data by celltype
# SummarizeDataByCelltype <- function(df = NULL){
#   tmp_data <- sapply(X = unique(df$celltype), FUN = function(ct){
#     tmp_idx <- which(df$celltype == ct)
#     tmp_df <- df[tmp_idx,]
#     n_genes <- length(unique(tmp_df$TargetGene))
#     n_chr <- length(unique(tmp_df$chr))
#     n_regions <- nrow(tmp_df)
#     rv <- c(n_chr, n_genes, n_regions)
#     names(rv) <- c("chr", "genes", "regions")
#     return(rv)
#   })
#   
#   cat("\n\tE-G regions summary: \n")
#   print(summary(tmp_data["regions",]))
#   
#   cat("\n\tChromosomes summary: \n")
#   print(summary(tmp_data["chr",]))
#   
#   cat("\n\tGenes summary: \n")
#   print(summary(tmp_data["genes",]))
# }

## Function name: LoadFulcoABCData
## input:     dir,   Data's directory. The directory must point to the folder with the celltype data fron Fulco ABC predictions.
## output:    List with the ABC data tables for different celltypes
# LoadFulcoABCData <- function(dir = NULL){
#   #List to store the tables in the input directory
#   table_list <- list()
#   
#   #Folders in the directory
#   folder_names <- list.files(path = dir)
#   
#   #Remove README file, if exists
#   tmp_idx <- which(str_detect(folder_names, pattern = "README"))
#   if(!is_empty(tmp_idx)){
#     folder_names <- folder_names[-tmp_idx]  
#   }
#   
#   # For each folder, extract files
#   for(folder in folder_names){
#     
#     #Get file names inside the folder
#     tmp_file_names <- list.files(path = paste(dir, folder, sep = "/"))
#     #Remove .gz files from the list
#     tmp_idx <- which(str_detect(string = tmp_file_names, pattern = ".gz"))
#     if(!is_empty(tmp_idx)){
#       tmp_file_names <- tmp_file_names[-tmp_idx]  
#     }
#     #File list directories
#     tmp_files <- paste(dir, folder, tmp_file_names, sep = "/")
#     
#     #Load the tables in the directories
#     tmp_file_list <- lapply(X = tmp_files, FUN = read.table, header = TRUE)
#     names(tmp_file_list) <- tmp_file_names
#     
#     #Add the data to the main data list
#     table_list[[folder]] <- tmp_file_list
#   }
#   
#   #Return data list
#   return(table_list)
# }

## Function name: GetFulcoABCDataFrame
## input:     data_ls,   Data list containing the data fron Fulco ABC predictions.
##            query_table,   query data to get from Fulco ABC predictions. Options are: 
##                          - AllPredictions
##                          - CandidateRegions
##                          - GeneSummary
##                          - PositivePredictions
## output:    Data frame with the query tables merged.
# GetFulcoABCDataFrame <- function(data_ls = NULL, query_table = "AllPredictions"){
#   
#   #Get the query table from each element in the list
#   data_ls <- lapply(X = data_ls, FUN = function(d){
#     names <- names(d)
#     idx <- which(str_detect(names, query_table))
#     return(d[[idx]])
#   })
#   
#   #Set as data frames
#   data_ls <- lapply(X = data_ls, FUN = as.data.frame)
#   
#   #Add cell type labels to each element in the list
#   celltypes_names <- names(data_ls)
#   data_ls <- lapply(celltypes_names, FUN = function(cell){
#     cell_table <- data_ls[[cell]]
#     cell_table["celltype"] <- cell
#     return(cell_table)
#   })
#   
#   #Merge the elements in the list
#   data_ls <- Reduce(f = rbind, x = data_ls)
#   
#   return(data_ls)
# }

## Function name: MySummary
## output:    Main statistics across `summarize_by` column. Includes violin plots of data and statistics.
# MySummary <- function(.data = NULL, col2summarize = NULL, summarize_by = NULL, plots = TRUE, box_width = 0.2, ...){
#   
#   names <- c(col2summarize, summarize_by)
#   .data <- .data[,names]
#   colnames(.data) <- c("w1", "type")
#   
#   summary <- .data %>%
#     group_by(type) %>%    #Group by "summarize_by"
#     summarise(n_elements = length(w1), #elements
#               min = min(w1, na.rm = TRUE),    #Group w1 by celltype
#               FirstQtl = quantile(w1, 0.25, na.rm = TRUE),
#               Median = quantile(w1, 0.5, na.rm = TRUE),
#               Mean = mean(w1, na.rm = TRUE),
#               ThirdQtl = quantile(w1, 0.75, na.rm = TRUE),
#               max = max(w1, na.rm = TRUE))
#   
#   return_summary = summary
#   
#   if(plots){
#     
#     ### Violin plots by type
#     p1 <- .data %>%
#       group_by(type) %>%     #Group by "summarize_by"
#       MyViolinPlot(x = "type", y = "w1",
#                    color = "type",box_width = box_width, ...) +
#       xlab(names[2]) + ylab(names[1])
#     
#     ### Violing plots by statistic
#     #Remove max and min columns. For now they are quite balanced.
#     summary <- summary[,-which(colnames(summary) %in% c("n_elements","min", "max"))]
#     #We first setup the data frame. We need the statistics names, values and types.
#     summary_list <- list()    #Create empty list to store the data frames.
#     # For each statistics in summary
#     for(i in seq_len(ncol(summary)-1)){
#       tmp_name <- colnames(summary)[i+1]
#       #create df with the celltypes and the values 
#       summary_list[[tmp_name]] <- data.frame("type" = summary[["type"]], "val" = summary[[tmp_name]]) 
#       #create statistics label
#       summary_list[[tmp_name]][["statistics"]] = tmp_name
#     }
#     
#     #Merge the datasets in the list and create violin plots
#     p2 <- Reduce(f = rbind, x = summary_list) %>%
#       ggplot( aes(x = statistics, y = val)) +  
#       geom_violin(trim = FALSE, width = 1.0) + 
#       theme_bw() + 
#       geom_boxplot(width=box_width, outlier.alpha = 0) +
#       geom_jitter(shape=16, position=position_jitter(0.0), mapping = aes(color = type), size = 2, alpha = 0.8)
#   }
#   
#   return(list(summary = return_summary, type_vln = p1, stats_vln = p2))
# }

## Function name: NoLegend
## output:    Removes the legend box and info from a ggplot object.
NoLegend <- function(){
  return(theme(legend.position = "none"))
}

## Function name: ExpandValue
## output:    Expands merged labels
ExpandValue <- function(.data = NULL, column = NULL, row = NULL, split = NULL){
  expanded <-  .data[row,] %>% 
    select(all_of(column)) %>%   #select column
    unlist() %>% #prepare for strsplit function
    strsplit(split = split) %>% #split the row by split
    unlist() #unlist the results
  
  return(expanded)
}

## Function name: GetRegionsCenter
## output:    Get the center of a given region. It needs the start and the length of the regions
GetRegionsCenter <- function(.data){
  #if length feature doesn't exist
  if(!"length" %in% colnames(.data)){
    .data <- .data %>% mutate("length" = end - start)
  }
  #get the center
  center = .data$start + floor(.data$length/2)
  .data <- .data %>% mutate("center" = center)
  return(.data)
}

## Function name: ExpandGenes
## output:    Expands the gene's information. Each Region-Gene appears once. After expanding, each region relates with only one gene/
ExpandGenes <- function(.data = NULL){
  data_expanded <- lapply(X = seq_len(nrow(.data)), FUN = function(region){
    #Get current region info
    current_region = .data[region,]
    #Get unique genes
    genes <- ExpandValue(.data = current_region, column = "genes", row = 1, split = "/")
    #Get the unique genes TSS (this gives a data frame with genes and TSS)
    current_genes_df <- genes_data %>% filter(TargetGene %in% genes)
    #Rename colums
    colnames(current_genes_df) <- c("Target_gene", "TSS")
    # Get the distance from center of regions to TSS
    current_genes_df <- current_genes_df %>% mutate(distance = abs(current_genes_df$TSS - current_region$region_center))
    #Set region info
    out_df <- data.frame(region_id = paste0(current_region$chr,":",current_region$start,"-", current_region$end,"|", current_genes_df$Target_gene))
    out_df <- out_df %>% mutate(chr = current_region$chr,
                                    start = current_region$start, 
                                    end = current_region$end,
                                    length = current_region$length,
                                    target_gene = current_genes_df$Target_gene,
                                    TSS = current_genes_df$TSS,
                                    distance = current_genes_df$distance,
                                    genes = current_region$genes,
                                    n_genes = current_region$n_genes,
                                    celltypes = current_region$celltypes,
                                    n_celltypes = current_region$n_celltypes,
                                    classes = current_region$classes,
                                    n_classes = current_region$n_classes,
                                    regions = current_region$regions,
                                    n_regions = current_region$n_regions)
    
    # genes_df <- genes_df %>% mutate(Region_center = current_region$region_center,
    #                                 Region_start = current_region$start, 
    #                                 Region_end = current_region$end,
    #                                 Region_chr = current_region$chr,
    #                                 Ngenes_in_region = current_region$n_genes)
    
    return(out_df)
  })
  
  data_expanded <- Reduce(f = rbind, x = data_expanded)
  return(data_expanded)
}

## Function name: RemoveDuplicatedRows
## output:    Remove duplicated rows in a df
RemoveDuplicatedRows <- function(.data){
  is_duplicated <- duplicated(x = .data)
  return(.data[!is_duplicated,])
}


## Function name: GetUniqueRegions
## output:    Df containing only unique regions
GetUniqueRegions <- function(.data = NULL){
  .data <- .data %>% filter(!duplicated(region_id))
  return(.data)
}

## Reduce regions using Granges
ReduceRegions <- function(regions_df = NULL, seqnames = "chr", start = "start", end = "end"){
  gr <- GRanges(seqnames = regions_df[[seqnames]], 
                ranges = IRanges(start = regions_df[[start]], end = regions_df[[end]]))
  gr <- reduce(gr)
  return(gr)
}

## Add metadata to regions reduced with Granges
ReduceRegion_addMD <- function(from_df = NULL, to_df = NULL){
  #we need to do this for each chromosome
  names_chromosomes <- levels(to_df$chr)
  
  merged_data <- lapply(X = names_chromosomes, FUN = function(c){ 
    current_data = from_df %>% filter(chr == c)
    #transform current_data to matrix to speed up computation
    current_data <- as.matrix(current_data, drop = FALSE)
    current_merged_data <- to_df %>% filter(chr == c)
    # initialize missing values in current_merged_data
    current_merged_data <- current_merged_data %>% mutate("n_celltypes" = NA)
    # for each merged_region in current_merged_data
    for(region in seq_len(nrow(current_merged_data))){
      
      #get the current region in current_merged_data
      current_region <- as.matrix(current_merged_data[region,])
      #find the overlapping regions from the original data, e.g. current_data
      overlap_idx <- which(as.numeric(current_data[,"start"]) >= as.numeric(current_region[,"start"]) &
                             as.numeric(current_data[,"end"]) <= as.numeric(current_region[,"end"]))
      current_overlap_data <- current_data[overlap_idx,,drop = FALSE]
      
      ## Set celltypes info
      current_merged_data[region, "n_celltypes"] <- length(unique(current_overlap_data[,"CellType"]))
    }
    
    return(current_merged_data)
  })
  merged_data <- Reduce(f = rbind, x = merged_data)
  return(merged_data)
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
  
  #It's weird that contiguous regions sometimes the overlap in the start and end, then I will add one base
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

# Which raw and merged regions overlaps
WhichOverlapEnhancer <- function(query_enhancer = NULL, subject_enhancer = NULL, min_overlap = 0){
  #get the query and subject genomic ranges objects
  gr_query <- GRanges(seqnames = query_enhancer$chr,
                      ranges = IRanges(start = query_enhancer$start,
                                       end = query_enhancer$end))
  gr_subject <- GRanges(seqnames = subject_enhancer$chr,
                        ranges = IRanges(start = subject_enhancer$start,
                                         end = subject_enhancer$end))
  #Find the overlaps
  overlaps <- findOverlaps(query = gr_query, subject = gr_subject,
                           minoverlap = min_overlap)
  #Get the overlapping region's ids
  overlaps_ids <- unique(subjectHits(overlaps))
  #return the ids
  return(overlaps_ids)
}

#Check enhancers data frame
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
  # cat("\nEnhancer-TSS overlapping OK: \n", 
  #     sum(enhancers %>% 
  #           mutate(t = if_else(condition = (tss >= start & tss <= end), 
                               # true = TRUE, false = FALSE)) %>% pull(t)) == 0)
  #check number of cell mismatch
  # cat("\nNumber of cell types mismatch OK: \n", 
  #     sum(enhancers %>% 
  #           mutate(t = test_n_celltypes - n_celltypes) %>% filter(t >0) %>% pull(t)) == 0)
  #check region ids
  cat("\nRegion ids OK:\n",
      sum(enhancers %>% 
            mutate(t = paste0(chr,":",start,"-",end)) %>%
            mutate (tt = (t == region_id)) %>% 
            pull(tt)) == nrow(enhancers))
  #check unique regions ids
  cat("\nUnique region ids OK:\n",
      enhancers %>% 
        pull(region_id) %>% 
        unique() %>% length() == nrow(enhancers))
}

#function to merge columns preserving md
MergeColumns <- function(.data = NULL, columns = NULL){
  merged_columns <- lapply(X = columns, FUN = function(c){
    tmp_df <- data.frame(value = .data[[c]], type = c)
    return(tmp_df)
  })
  merged_columns <- Reduce(f = rbind, x = merged_columns)
  return(merged_columns)
}

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

# Get regions overlapping
WhichRegionsOverlap <- function(chr = NULL, query_enhancer = NULL, subject_enhancer = NULL, min_overlap = 0){
  
  #get the query and subject genomic ranges objects
  if(is.null(chr)){
    gr_query <- GRanges(seqnames = query_enhancer[,"chr"],
                        ranges = IRanges(start = query_enhancer[,"start"],
                                         end = query_enhancer[,"end"]))
    gr_subject <- GRanges(seqnames = subject_enhancer[,"chr"],
                          ranges = IRanges(start = subject_enhancer[,"start"],
                                           end = subject_enhancer[,"end"]))
  }else{
    gr_query <- GRanges(seqnames = chr,
                        ranges = IRanges(start = query_enhancer[,"start"],
                                         end = query_enhancer[,"end"]))
    gr_subject <- GRanges(seqnames = chr,
                          ranges = IRanges(start = subject_enhancer[,"start"],
                                           end = subject_enhancer[,"end"]))
  }
  
  #Find the overlaps
  overlaps <- findOverlaps(query = gr_query, subject = gr_subject,
                           minoverlap = min_overlap, )
  #return the ids
  return(overlaps)
}


#Create region ids using the chr, start and end
GetRegionIds <- function(.data){
  .data <- .data %>% mutate(region_id = paste0(chr,":",start,"-",end))
  return(.data)
}

#Filter matrix (similar to dplry one)
FilterMatrix <- function(m = NULL, col = NULL, wt = NULL){
  idx <- which(m[,col] == wt)
  return(m[idx, , drop = FALSE])
}

#Transfer columns between dfs
TransferColumns <- function(reference = NULL, query = NULL, anchors = NULL, cols_to_tranfer = NULL){
  
  #verify if the columns to transfer doesn't exist in the reference data. If exist give a warning and erase the column
  if(cols_to_tranfer %in% colnames(reference)){
    warning(call. = TRUE, "The column to transfer already exists in the reference. The old column will be replaced.")
    reference <- reference %>% RemoveColumn(column = cols_to_tranfer)
  }
  
  #set the anchors as rownames
  rownames(query) <- query[,anchors]
  #arrange the query as the reference
  query <- query[reference[,anchors],]
  #short test
  if(identical(query[,anchors], reference[,anchors])){
    #transfer the columns 
    reference <- cbind(reference, query[,cols_to_tranfer, drop = FALSE]) 
  }else{
    stop(call. = TRUE, "ERROR: the anchors in reference and query are different.")
  }
  return(reference)
}

#clustering with louvain algorithm
ClusterLouvain <- function(x, k = 10, resolution = 0.5) {
  
  g <- bluster::makeSNNGraph(x, k = k)
  res <- igraph::cluster_louvain(g, resolution = resolution)
  
  memberships <- igraph::membership(res)
  
  return(memberships)
}

#Rename columns
RenameColumn <- function(.data = NULL, old_column_name = NULL, new_column_name = NULL){
  #get column index
  col_idx <- which(colnames(.data) == old_column_name)
  #replace the name
  colnames(.data)[col_idx] <- new_column_name
  return(.data)
}

RelabelClusters <- function(.data = NULL, clusters = NULL, old_labels = NULL, new_labels = NULL, new_clusters = NULL){
  #check that clusters are in .data
  if(!clusters %in% colnames(.data)){
    stop(call. = TRUE, "ERROR: The clusters labels were not found.")
  }
  
  #check that the size of old and new labels match
  if(length(old_labels) != length(new_labels)){
    stop(call. = TRUE, "ERROR: The number of labels doesn't match.")
  }
  
  #check that the clusters are a factor
  if(!"factor" %in% is(.data[[clusters]])){
    stop(call. = TRUE, "ERROR: The clusters should be a factor.")
  }
  
  #factor's levels are stored as characters. For comparisons let's transform the old/new labels to character
  old_labels <- as.character(old_labels)
  new_labels <- as.character(new_labels)
  
  #if the new_clusters is NULL, we assign the new labels in the previous factor
  if(is.null(new_clusters)){
    new_clusters <- clusters
  }
  
  #get the old clusters for comparison
  old_clusters <- as.character(.data[[clusters]])
  #init the new clusters
  ncl <- as.character(.data[[clusters]])
  
  for(i in seq_along(old_labels)){
    
    #get the index of the current element
    idx <- which(old_clusters == old_labels[i])
    #assign the new clusters
    ncl[idx] <- new_labels[i]
  }
  
  #if the new and old cluster has the same name, we need to remove the previous clusters
  if(new_clusters == clusters){
    .data <- .data[,-which(colnames(.data) == clusters)]
  }
  
  #new colnames
  new_colnames <- c(colnames(.data), new_clusters)
  #assign the new clusters
  .data <- cbind(.data,factor(x = ncl))
  #assisng colnames
  colnames(.data) <- new_colnames
  
  return(.data)
}

RemoveColumn <- function(.data = NULL, column = NULL){
  #ut. Columns in colnames
  if(sum(!column %in% colnames(.data)) == length(column)){
    stop(call. = TRUE, "ERROR: Columns not found in the data.")
  }
  idx <- which(colnames(.data) %in% column)
  .data <- .data[,-idx]
  return(.data)
}

ScaleFeatures <- function(.data = NULL, features = NULL, center = TRUE, scale = TRUE){
  
  new_colnames <- paste0(features,"_scaled")
  
  if(new_colnames %in% colnames(.data)){
    warning(call. = TRUE, "The column to transfer already exists in the reference. The old column will be replaced.")
    .data <- .data %>% RemoveColumn(column = features)
  }
  
  fea_data <- .data %>% dplyr::select(all_of(features)) %>% as.matrix()
  fea_data <- scale(x = fea_data, center = center, scale = scale)
  colnames(fea_data) <- new_colnames
  .data <- cbind(.data, fea_data)
  return(.data)
}

MapGeneLabels <- function(genes_labels = NULL, from_type = NULL, to_type = NULL){
  mapping <- AnnotationDbi::select(x = org.Hs.eg.db, 
                                   keys = genes_labels,
                                   keytype=from_type,
                                   columns = to_type)
  colnames(mapping) = c("original_name", "mapped_name")
  #remove unmapped labels
  mapping <- mapping %>% filter(!is.na(mapped_name))
  return(mapping)
}
 

