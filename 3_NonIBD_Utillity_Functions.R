# Function to run Fishers exact test along the genome.
#
# The Fishers.Exact.Test() takes the following inputs and run the Fisher's exact test:
# 1. sample_data which is a list. We generated this in 1_SimulateData.R script.
# 2. A vector indicating the case/control status of individuals. 1 = Case , 0 = Control
#
Fishers.Exact.Test = function(sample_data, CCLabel ){
  
  sample_geno = sample_data$Genos$sample_geno
  
  pvalue      = numeric( length = nrow(sample_geno) )
  for( i in 1:nrow(sample_geno) ){
    pvalue[i] = fisher.test(x = sample_geno[i,], CCLabel)$p.value
  }
 
  return( list(pvalue = pvalue, pos = sample_data$Posn$SNV_Position) )
}



#
# The permute_FET function to run the permutataion.
# 
# Parameters of this function:
# 1. nperm is the desired number of permutations
# 2. sample_data is the sample object generated in 1_SimulateData.R script. 
permute_FET = function(nperm, sample_data){
  
  # A matrix to record the pvalue from each permutation. 
  # The rows of this matrix represents the permutation number.
  # The last row of this matrix, records the pvalue for the original case/control labeling.
  # The columns are the SNV positions.
  FET_permutations  = matrix(NA, nrow = nperm + 1, ncol = length(sample_data$Posn$SNV_Position) )
  
  
  # ccLabel is the original case/control labeling of the individuals.
  ccLabel = sample_data$Genos$ccStatus
  
  
  # permute_indx is a matrix of indices for each permutation. 
  permute_indx = matrix(NA, nrow = nperm, ncol = length(ccLabel) )
  for(i in 1:nperm){
    permute_indx[i,] = sample( x = 1:length(ccLabel), replace = FALSE)
  }
  
  

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
  
  
  # To run parallel on your own PC/Laptop:
  # cl = makeCluster(detectCores() - 1)
  # registerDoParallel(cl)
  
  
  #
  # We need to export a copy of functions to each node/worker that runs the job in 
  # parallel. We can do this by using .export argument in foreach.
  #
  
  res = foreach(i = 1:nperm, .export = "Fishers.Exact.Test" ) %dopar%{
    pCCLabel              = ccLabel[ permute_indx[i,] ]
    x                     = Fishers.Exact.Test(sample_data = sample_data, CCLabel = pCCLabel)
    x$pvalue 
  }
  
  stopCluster(cl)
  
  for(i in 1:length(res)){
     FET_permutations[i,]  = res[[i]]
  }


  colnames(FET_permutations) = sample_data$Posn$SNV_Names
  return(FET_permutations = FET_permutations)
  
}




#
# Function to run SKATO test along the genome.
#
# The SKATO_TEST() takes the following inputs and run the SKATO test along the genome:
# 1. sample_data which is a list. We generated this in 1_SimulateData.R script.
# 2. A vector indicating the case/control status of individuals. 1 = Case , 0 = Control
# 3. type: Indicates the type of the phenotype. D = Dichotomous, C = Continous
SKATO_TEST = function(sample_data, CCLabel, type = "D" ){
  
  library(SKAT)
  
  sample_geno = sample_data$Genos$sample_geno 

  if( type == "D"){
    obj = SKAT::SKAT_Null_Model( CCLabel ~ 1, out_type = "D")
    
    
    pvalue = numeric(length = nrow(sample_geno))
    for( i in 1:nrow(sample_geno) ){
      pvalue[i] = SKAT::SKAT( Z = as.matrix(sample_geno[i,]), obj, method = "SKATO")$p.value
    }
    
    return( list( pvalue = pvalue, pos = sample_data$Posn$SNV_Position ) )
  }
  
  
  
  # This section will be developed later when we start working on continous trait.
  # 
  # if( type == "C" ){
  #   obj = SKAT_Null_Model( phenos ~ 1, out_type = "C")
  #   
  #   
  #   skat.pvalue = c()
  #   for( i in 1:nrow(Genomat) ){
  #     skat.pvalue[i] = SKAT( Z = as.matrix(Genomat[i,]), obj, method = "SKATO"  )$p.value
  #   }
  #   
  #   return(skat.pvalue)
  # }
  
  
}



#
# The permute_SKATO function to run the permutataion.
# 
# Parameters of this function:
# 1. nperm is the desired number of permutations
# 2. sample_data is the sample object generated in 1_SimulateData.R script. 
permute_SKATO = function(nperm, sample_data){
  
  
  # A matrix to record the pvalue from each permutation. 
  # The rows of this matrix represents the permutation number.
  # The last row of this matrix, records the pvalue for the original case/control labeling.
  # The columns are the SNV positions.
  SKATO_permutations = matrix(NA, nrow = nperm + 1, ncol = length(sample_data$Posn$SNV_Position) )
  
  
  
  # ccLabel is the original case/control labeling of the individuals.
  ccLabel = sample_data$Genos$ccStatus
  
  
  # permute_indx is a matrix of indices for each permutation. 
  permute_indx = matrix(NA, nrow = nperm, ncol = length(ccLabel) )
  for(i in 1:nperm){
    permute_indx[i,] = sample( x = 1:length(ccLabel), replace = FALSE)
  }
  
  
  # Running the permutation
  # for(i in 1:nperm){
  #   
  #   pCCLabel               = ccLabel[ permute_indx[i,] ]
  #   x                      = SKATO_TEST(sample_data = sample_data, CCLabel = pCCLabel)
  #   SKATO_permutations[i,] = x$pvalue 
  #   print( paste0("Currently at permutation ",i , " out of ",nperm ) )
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
  # We need to export a copy of functions to each node/worker that runs the job in 
  # parallel. We can do this by using .export argument in foreach.
  #
  
  res = foreach(i = 1:nperm, .export = "SKATO_TEST") %dopar%{
     pCCLabel               = ccLabel[ permute_indx[i,] ]
     x                      = SKATO_TEST(sample_data = sample_data, CCLabel = pCCLabel)
     x$pvalue
  }
  
  
  stopCluster(cl)
  
  
  for(i in 1:length(res)){
    SKATO_permutations[i,] = res[[i]]
  }
  
  colnames(SKATO_permutations) = sample_data$Posn$SNV_Names
  return(SKATO_permutations = SKATO_permutations)

}

