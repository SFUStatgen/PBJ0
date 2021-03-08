# 
# This function takes the path to the folder of any simulated dataset 
# and calculates the results of detection for each method.
# The output is a numeric vector where ach element is the calculated pvalue
# based on the permutation of case/control status of individuals.
# 
# The output vector reports the pvalues in the following order:
# "FET", "SKATO", "dCor", "Mantel"
# 
get_detection_results    = function(path){
  

  # Create an empty numeric vector to record the pvalue obtained by permutation
  pvalue = numeric( length = 4)
  

  if( file.exists( paste0(path, "FET_perm.RData") ) ){
    load( paste0(path, "FET_perm.RData") )  
    n_perm = nrow(FET_permutations) - 1
    perms = apply(X = -log(FET_permutations,base=10), MARGIN = 1, FUN = max)
    pvalue[1] = mean( perms >= perms[ n_perm + 1 ] )
  }else{
    pvalue[1] = NA
  }

  
  if( file.exists(paste0(path, "SKATO_perm.RData")) ){
    load( paste0(path, "SKATO_perm.RData") )
    n_perm = nrow(SKATO_permutations) - 1
    perms = apply(X = -log(SKATO_permutations,base=10), MARGIN = 1, FUN = max)
    pvalue[2] = mean( perms >= perms[ n_perm + 1 ] )
  }else{
    pvalue[2] = NA
  }
  
  
  
  if( file.exists(paste0(path, "dCorN_perm.RData")) ){
    load( paste0(path, "dCorN_perm.RData") )
    n_perm = nrow(dCorN_permutation) - 1
    perms = apply(X = dCorN_permutation , MARGIN = 1, FUN = max)
    pvalue[3] = mean( perms >= perms[ n_perm + 1 ] )
  }else{
    pvalue[3] = NA
  }

  
  if( file.exists(paste0(path, "Mantel_perm.RData")) ){
    load( paste0(path, "Mantel_perm.RData") )
    n_perm = nrow(Mantel_permutation) - 1
    perms = apply(X = Mantel_permutation^2 , MARGIN = 1, FUN = max)
    pvalue[4] = mean( perms >= perms[ n_perm + 1 ] )
  }else{
    pvalue[4] = NA
  }


  names(pvalue) = c("FET", "SKATO", "dCor", "Mantel")
  return(pvalue)
}


