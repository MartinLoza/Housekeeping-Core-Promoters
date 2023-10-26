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
                                    regions = current_region$regions,
                                    n_regions = current_region$n_regions)
    
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
 

