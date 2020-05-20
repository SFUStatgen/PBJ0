#
# get_unique_part_and_weight() function takes the partitions list and sample_data and returns the 
# unique reconstructed partitions with their weights.
#
get_unique_part_and_weight = function(partition_list, sample_data){
  
  pos           = sample_data$Posn$SNV_Position
  pos           = c(0,pos)
  
  dupIndx       = which( duplicated(partition_list) == FALSE )
  
  weight.vector = numeric( length = length(dupIndx) )
  for(i in 1:(length(dupIndx)-1)){
    weight.vector[i] =  (pos[dupIndx[i+1]] - pos[dupIndx[i]]) / pos[length(pos)]
  }
  
  
  weight.vector[length(dupIndx)] = 1 - sum(weight.vector)
  partition_list                 = partition_list[ dupIndx ]
  
  
  return( list(partition_list = partition_list, weights = weight.vector) )
}






# 
# dCorIG is the key function to calculate the GNN statistic and return the dCor profile along the genome.
# It takes the following inputs:
#
# 1. partition_list = A list of reconstructed partitions
#
# 2. w = A numeric vector recording the weight of each reconstructed partition.
#
# 3. distance_matrix_list = A list of distance matrices. 
#
# 4. CClabels: A dataframe recording the sequence labels and their case/control status.
#
dCorIG                     = function(partition_list, w, distance_matrix_list, CClabels, q_cut_off){
  
  # 1. Calculate the GNN statistic for all sequences in CClabels
  GNN_Statistic = get_GNN_statistic(part = partition_list, w, CClabels)
  
  
  # 2. Reclassify the case sequences according to the 80-th percentile of the 
  #   GNN statistic in the control sequences
  rCClabels = reclassify_cases_gnn(CClabels = CClabels, GNN_Statistic = GNN_Statistic, q_cut_off) 
  
  
  # 3. Calculating the dCor profile accross the genome
  # We pass the first and fourth column of rCClabels dataframe (shown below).
  # 
  #   SeqID   Status  GNN_Statistic       GNN_implied_Status
  #   5202       1      0.6809889                  1
  #   3945       0      0.4152137                  0
  #   2012       1      0.8785675                  1 
  #    ...      
  # 
  # The first column is the sequence labels. The second column is the GNN implied status.
  # 
  dCorIG = dCor_Profile(cc_sample = rCClabels[,c(1,4)], distance_matrix_list = distance_matrix_list)
  
  
  # return dCor profile along the genome.
  return(dCorIG)
}






# get_GNN_statistic() is an internal function which is called in dCorIG() function.
get_GNN_statistic          = function(part, w, CClabels){
  
  
  # 
  # Initialize the GNN as a matrix of NA values. 
  # This matrix has 200 rows, each of them representing one of the sequences in the sample.
  # The number of columns in this matrix equal to the number of unique reconstructed partitions.
  # 
  # We go over each reconstructed partition. We pass every reconstructed partition and the list of 
  # sequences to GNN_Calculation() function. This function returns a numeric vector with
  # the same size as the given list of the sequences. We multiply this vector by the corresponding 
  # weight of each reconstructed partition. We record this vector at one of the columns of GNN matrix. 
  # 
  # In summary, each column records 200 GNN values for the corresponding reconstructed partition.
  #
  GNN = matrix(NA, nrow = nrow(CClabels), ncol = length(part) )
  
  for(i in 1:length(part) ){
    value   = GNN_Calculation( reconstructed_partition = part[[i]], CClabels = CClabels )
    GNN[,i] = value * w[i]
  }
  
  
  
  #
  # Calculating GNN statistic for each sequence.
  # The GNN statistic for each haplotype is the weighted average across the region.
  # We have already multiplied the GNN proportion by the weight of each partition in the previous stage.
  # So, for each sequence, we only need to sum all values across the genomic region.
  # 
  GNN_Statistic = numeric( length = nrow(CClabels) )
  names(GNN_Statistic) = CClabels$SeqID
  for( i in 1:nrow(GNN) ){
    GNN_Statistic[i] = sum( GNN[i,] )
  }
  
  return(GNN_Statistic)
}






# GNN_Calculation() is an internal function which is called in get_GNN_statistic() function.
GNN_Calculation            = function(reconstructed_partition, CClabels){
  
  
  # Define the Prop_Vector to save the proportion of case haplotypes in the 
  # nearest genealogical of any sequences/haplotypes. 
  # A vector of size 200 to save the proportion for each sequence/haplotype
  # Also add the names attribute to make it easier to find the value for each haplotype
  Prop_Vector           = numeric(length = nrow(CClabels) )
  hap                   = CClabels$SeqID
  names(Prop_Vector)    = hap
  
  
  # Identify the labels of case sequences
  case_sequences = as.character( CClabels[ which(CClabels$Status == 1) , 1] )
  
  
  
  #  Extracting the tip labels of the i-th reconstructed partition
  tip.label.list = strsplit( reconstructed_partition$tip.label , '-')
  
  
  
  # Create a data frame to map the sequence/haplotype id to the tip label in the 
  # reconstructed partition.
  # 
  # First column of this dataframe is haplo_id_1 which saves the tip label on 
  # the reconstructed partition
  # Second column of this dataframe is haplo_id_2 which saves the sequence/haplotype IDs.
  #
  # We need to have this dataframe to map the sequence/haplotype label to the tip label on
  # the reconstructed partition.
  #
  haplo_id_1 = c()
  haplo_id_2 = c()
  for( k in 1:length(tip.label.list) ){
    haplo_id_1 = c( haplo_id_1 , rep(k, length(tip.label.list[[k]])) )
    haplo_id_2 = c( haplo_id_2 , tip.label.list[[k]] )
  }
  mapping.haplo = data.frame( haplo_id_1, haplo_id_2 )
  
  
  
  # Iterate over all 200 sequence/haplotype in the sample and find the proportion of case 
  # haplotypes that are in the genealogical nearest neighbor of them.
  for( hap_iteration in 1:length(hap) ){
    
    
    # Find the tip label of each sequence/haplotype on the reconstructed partition
    # Save this tip label as target_node
    indx = which( mapping.haplo$haplo_id_2 == hap[hap_iteration] ) 
    target_node = mapping.haplo[indx, "haplo_id_1" ]
    
    
    # Find the immediate node in the path to the root for this target node
    # and save it as parent_node.
    indx = which(reconstructed_partition$edge == target_node, TRUE)
    parent_node = reconstructed_partition$edge[indx[1,1]]
    
    
    # Extract the clade descending from this parent_node and calculating the 
    # proportion of case hapltoypes in this clade, excluding the target haplotype
    subtree = ape::extract.clade(reconstructed_partition, parent_node)
    nodes   = unlist( strsplit( subtree$tip.label , '-' ) )
    nodes   = nodes[ nodes != hap[hap_iteration] ]
    prop    = mean( nodes %in% case_sequences )
    
    Prop_Vector[hap_iteration] = prop
  }
  
  return(Prop_Vector)
}





#
# reclassify_cases_gnn() function takes the following inputs:
# 1. CClabels: A dataframe recording the sequence labels and their case/control status.
#
# 2. GNN_Statistic: A vector recording the GNN statistic for each sequence.
#
# 3. q_cut_off: The quantile value defined by user. This is the quantile of the GNN statistic 
#    in the control group that we would like to consider to reclassify case sequences.
#
# This function returns back CClabels dataframe with two additional columns. 
#  1)GNN_Statistic: recording the GNN statistic of each sequence, 
#  2) GNN_implied_Status: indicates the status of sequence after reclassification 0 = Control, 1 = Case
# 
reclassify_cases_gnn       = function(CClabels, GNN_Statistic, q_cut_off){
  
  CClabels$"GNN_Statistic" = GNN_Statistic 
  
  # Finding the user requested percentile (q_cut_off) in the control group.
  cut_off                = quantile( x = CClabels[CClabels$Status == 0, "GNN_Statistic"], probs = q_cut_off )
  
  # Find which case sequences should be reclassified as carrier.
  case_carrier_indx      = which( CClabels$Status == 1 & CClabels$GNN_Statistic >= cut_off )
  
  # Find which case sequences shuld be reclassified as non-carrier
  case_non_carrier_indx  = which( CClabels$Status == 1 & CClabels$GNN_Statistic  < cut_off )
  
  # Find where are the control sequences
  control_indx           = which( CClabels$Status == 0 )
  
  # Apply the necessary changes and reclassify case sequences.
  CClabels[case_carrier_indx,     "GNN_implied_Status"] = 1
  CClabels[case_non_carrier_indx, "GNN_implied_Status"] = 0
  CClabels[control_indx,          "GNN_implied_Status"] = 0

  return(CClabels)
  
}






#
# pheno_dist_matrix_calc() calculates the phenotypic distance matrix. It takes 2 inputs:
#
# 1. A numeric vector y of length n, representing the phenotype of each sequence.
#    In the case of dichomotous trait it is coded as 1 = Case and 0 = Control
# 
# 2. The probability of disease in the population as prob. This must be a number 
#    between 0 and 1.
# 
# It returns a nxn phenotypic distance matrix according to calculation explained in:
# Burkett et al. 2014.(Link: https://www.karger.com/Article/Abstract/363443 )
# 
pheno_dist_matrix_calc     = function(cc_sample, prob){
  
  seq_names = cc_sample[,1]
  y         = cc_sample[,2]
  
  d = matrix(NA, nrow = length(y), ncol = length(y) )
  for( i in 1:length(y) ){
    for( j in 1:length(y) ){
      s_ij = (y[i] - prob) * (y[j] - prob)
      d[i,j] = 1 - s_ij
    }
  }
  
  diag(d) = 0
  rownames(d) = seq_names
  colnames(d) = seq_names
  
  return(d)
}






#
# dCor_Profile() is the key function to calculate the dCor profile accross the region.
# It takes the following arguments as input:
# 
# 1. cc_sample: A dataframe with two columns: SeqID and Status
# The SeqID shows the label of the sequences.
# The Status represents the case/control status of each sequence. 0 = control, 1 = case
# 
dCor_Profile           = function(cc_sample, distance_matrix_list){
  
  phenotye_dist_matrix  = pheno_dist_matrix_calc(cc_sample, prob = 0.05)
  
  n         = length(distance_matrix_list)
  dCor_Stat = numeric( length = n )
  for(i in 1:n){
    partition_dist  = distance_matrix_list[[i]][rownames(phenotye_dist_matrix),colnames(phenotye_dist_matrix)]
    dCor_Stat[i]    = perfectphyloR::dCorTest(Dx = partition_dist, Dy = phenotye_dist_matrix, nperm = 0)$Stat
    #print( paste0( 'Position ',i , ' out of ', n ) )
  }
  
  return(dCor_Stat)
  
}








# dCorIG_permute() can be called to perform the permutation test for dCorIG.
# The arguments of this function are as follows:
# 1. nperm = The desired number of permutation
# 2. sample_data = the sample object generated by Sample_from_population() function.
# 3. distance_matrix_list = A list of distance matrices.
#
dCorIG_permute          = function(nperm, sample_data, partition_list, w, distance_matrix_list, q_cut_off){
  
  # A matrix for recording the results of permutation.
  # Each row represents one permutation. We add one extra row to save the results of 
  # our observed statistic later.
  #
  # Columns are the SNVs along the genomic region.
  dCorIG_permutation  = matrix(NA, nrow = nperm + 1, ncol = length(sample_data$Posn$SNV_Position) )
  
  
  # ccLabel is the original case/control labeling of the individuals.
  ccLabel = sample_data$Genos$ccStatus
  
  
  
  # Load the permutation indices
  load("permute_indx.RData")
  
  

  # Running the permutation
  
  # pCCLabel = c() 
  # for(i in 1:nperm){
  #   pCCLabel$SeqID         = c(sample_data$CaseHapID, sample_data$ControlHapID) 
  #   pCCLabel$Status        = ccLabel[ permute_indx[i,] ]
  #   dCorIG_permutation[i,] = dCorIG(partition_list, w, distance_matrix_list, CClabels = pCCLabel )
  # }

  
  # Running the permutation in parallel
  library(foreach)
  library(doParallel)
  
  # To run parallel on compute canada:
  # (Ref: https://docs.computecanada.ca/wiki/R)
  #
  # Create an array from the NODESLIST environnement variable
  nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
  
  
  # Create the cluster with the nodes name. One process per count of node name.
  # nodeslist = node1 node1 node2 node2, means we are starting 2 processes on node1, likewise on node2.
  cl = makeCluster(nodeslist, type = "PSOCK") 
  registerDoParallel(cl)
  
  
  
  # To run parallel on your own PC/Laptop:#
  # cl = makeCluster(detectCores() - 1)
  # registerDoParallel(cl)
  
  
  #
  # List of the functions that need to export on each worker/slave when we run the code in parallel.
  # dCorIG() is the key function to calculate the GNN statistic and obtain the dCor profile across 
  # region every time we permute the sample. We call various functions in side this function to 
  # run the calculations. Hence, we need to use .export argument in foreach to copy every function
  # that is needed to the workers/slaves that are running the job in parallel.
  # We create a character vector export.list to hold the name of these functions.
  # 
  export.list = c("get_GNN_statistic",
                  "reclassify_cases_gnn", 
                  "dCor_Profile", 
                  "GNN_Calculation", 
                  "pheno_dist_matrix_calc", 
                  "dCorIG")
  
  

  res = foreach(i = 1:nperm, .export = export.list) %dopar%{
    pCCLabel = data.frame(SeqID  = c(sample_data$CaseHapID, sample_data$ControlHapID),
                          Status = unlist( lapply( ccLabel[ permute_indx[i,] ], FUN = function(x){rep(x,2)} ) )
                          )
    dCorIG(partition_list, w, distance_matrix_list, CClabels = pCCLabel, q_cut_off = q_cut_off )
  }
  
  stopCluster(cl)
  
  for( i in 1:length(res) ){
    dCorIG_permutation[i,] = res[[i]]
  }
  
  
  
  colnames(dCorIG_permutation) = sample_data$Posn$SNV_Names
  return( dCorIG_permutation )
  
}
