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
  
  
  load('permute_indx.RData')
  

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
# 2. CCLabel is a vector indicating the case/control status of individuals. 1 = Case , 0 = Control
# 3. type: Indicates the type of the phenotype. D = Dichotomous, C = Continous
# 4. window.size: The size of the window (including target SNV) in base pair
SKATO_TEST         = function(sample_data, CCLabel, window.size, type = "D"){
  
  if( (window.size %% 2) == 0 ){
    stop('window.size must be an even number and greater than 1')
  }
  
  library(SKAT)
  
  sample_geno = t(sample_data$Genos$sample_geno) 
  
  if( type == "D"){
    
    obj = SKAT_Null_Model( CCLabel ~ 1, out_type = "D")
    
    iter = seq(1, dim(sample_geno)[2], by = 1)

    #
    # SKAT function uses random sampling methods while calculating p-value.
    # The values in each run may be slightly different. 
    # In order to be able to reproduce the same p-values (did not happen in our case but if needed in the future), 
    # we use set.seed() function as recommended in SKAT package document.
    # 
    # Direct from SKATO documentation in R (https://cran.r-project.org/web/packages/SKAT/SKAT.pdf) at page 44:
    # Since small sample adjustment uses random sampling to estimate the kurtosis of the test statistics,
    # SKAT with the (kurtosis-based) small sample adjustment can yield slightly different p-values for
    # each run. If you want to reproduce p-values, please set a seed number using set.seed function in R
    # 
    # 
    set.seed(123)
    
    UB <- (window.size - 1) / 2
    LB <- (window.size - 1) / 2
    
    pvalue  = c()
    posn    = c()
    for( k in 1:length(iter) ){
      
       snvnos = seq(iter[k] - LB, iter[k] + UB, by = 1)
       snvnos <- snvnos[snvnos>0]
       snvnos <- snvnos[snvnos<ncol(sample_geno) + 1]

       Z = sample_geno[,snvnos]

       pvalue[k] = SKAT(Z, obj, kernel = "linear.weighted", method = "optimal.adj")$p.value

       posn[k]   = sample_data$Posn$SNV_Position[k]

    }

    return(list( pvalue = pvalue, pos = posn ))
    
  }

}




#
# The permute_SKATO function to run the permutataion.
# 
# Parameters of this function:
# 1. nperm is the desired number of permutations
# 2. sample_data is the sample object generated in 1_SimulateData.R script. 
# 3. window.size: The size of the sliding window in base pair
permute_SKATO = function(nperm, sample_data, window.size){
  
  
  # Determine the number of columns for SKATO_permutations
  # sample_geno = t(sample_data$Genos$sample_geno) 
  # iter = seq(1, dim(sample_geno)[2], by = (window.size - overlap))
    
  
  # A matrix to record the pvalue from each permutation. 
  # The rows of this matrix represents the permutation number.
  # The last row of this matrix, records the pvalue for the original case/control labeling.
  # The columns are the SNV positions.
  SKATO_permutations = matrix(NA, nrow = nperm + 1, ncol = nrow(sample_data$Posn) )
  
  
  
  # ccLabel is the original case/control labeling of the individuals.
  ccLabel = sample_data$Genos$ccStatus
  
  
  load('permute_indx.RData')
  
  
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
  #cl = makeCluster(detectCores() - 1)
  #registerDoParallel(cl)
  
  
  
  #
  # We need to export a copy of functions to each node/worker that runs the job in 
  # parallel. We can do this by using .export argument in foreach.
  #
  
  res = foreach(i = 1:nperm, .export = "SKATO_TEST") %dopar%{
     pCCLabel               = ccLabel[ permute_indx[i,] ]
     x                      = SKATO_TEST(sample_data = sample_data, CCLabel = pCCLabel, window.size = window.size)
     x$pvalue
  }
  
  
  stopCluster(cl)
  
  
  for(i in 1:length(res)){
    SKATO_permutations[i,] = res[[i]]
  }
  
  # colnames(SKATO_permutations) = sample_data$Posn$SNV_Names
  return(SKATO_permutations = SKATO_permutations)

}

