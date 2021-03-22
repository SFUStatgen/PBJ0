# This function takes 2 arguments:
# 1. peakposn: The position on the genome where the maximum value of any association profile has occured.
#
# 2. reegion: A numeric vector of size two indicating the lower and upper bound of the causal region.
#
# This function returns the distance from the risk region (passed as region to function).
# In the scenario where there are multiple peaks, it returns the average of the distances from the region.
# 
avgPeakDist = function(peakposn, region){
  
  LowerBound = region[1]
  UpperBound = region[2]
  
  
  if(any(peakposn < LowerBound)){
    
    dist1 = LowerBound - peakposn[peakposn < LowerBound]
    
  }else{dist1 = NA}
  
  if(any(peakposn > UpperBound)){
    
    dist2 = peakposn[peakposn > UpperBound] - UpperBound
    
  }else{dist2 = NA}
  
  if(any( peakposn >= LowerBound & peakposn <= UpperBound )){
    
    
    dist3 = rep(0,length(peakposn >= LowerBound & peakposn <= UpperBound))
    
  }else{dist3 = NA}
  
  dist = c(dist1, dist2, dist3)
  
  if(any(is.na(dist))){
    
    avgDist = mean(dist[!is.na(dist)])
    
  }else{
    
    avgDist = mean(dist)
  }
  
  return(peakAvgDist = avgDist)
}



# 
# This function takes the path to the folder of any simulated dataset 
# and calculates the distance from the risk causal region for every method.
# The output is a numeric vector. Each element of this vector, denotes the 
# distance from the risk causal region. 
# 
# The output vector reports the distances in the following order:
# "FET", "SKATO", "dCorN", "dCorIT", "dCorIG
# 
get_localization_results = function(path){
  

  # Create an empty numeric vector to record the average distance 
  # from the risk region for each method
  average_distance = numeric( length = 4)
  

  
  # Loading the localization results for non IBD methods and IBD based methods
  # A numeric vector of positions on the genome.
  # We can use this for all methods except SKATO.
  load( paste0(path, "FET.RData")   )
  pos = FET$pos
  
  
  if( file.exists(paste0(path, "FET.RData") ) ){
    load( paste0(path, "FET.RData")   )
    max.indx            = which.max( -log( FET$pvalue, base = 10 ) )
    average_distance[1] = avgPeakDist( peakposn = pos[max.indx], region = c(900000,1100000) )
  }else{
    average_distance[1] = NA
  }
  
  
  
  if( file.exists( paste0(path, "SKATO.RData" )) ){
    load( paste0(path, "SKATO.RData") )
    pos.skato           = SKATO$pos
    max.indx            = which.max( -log( SKATO$pvalue, base = 10 ) )
    average_distance[2] = avgPeakDist( peakposn = pos.skato[max.indx], region = c(900000,1100000) )
  }else{
    average_distance[2] = NA
  }
  
  
  
  if( file.exists(paste0(path, "dCorN.RData")) ){
    load( paste0(path, "dCorN.RData")  )
    max.indx            = which.max( dCorN )
    average_distance[3] = avgPeakDist( peakposn = pos[max.indx], region = c(900000,1100000) )
  }else{
    average_distance[3] = NA
  }


  if( file.exists(paste0(path, "Mantel.RData")) ){
    load( paste0(path, "Mantel.RData") )
    max.indx            = which.max( Mantel^2 )
    average_distance[4] = avgPeakDist( peakposn = pos[max.indx], region = c(900000,1100000) )
  }else{
    average_distance[4] = NA
  }
  

  names(average_distance) = c("FET", "SKATO", "dCorN", "Mantel")
  return(average_distance)
}





# 
# This function takes the path to the folder of any simulated dataset 
# and calculates the results of detection for each method.
#
# The output is a numeric vector where each element is the calculated pvalue
# based on the permutation of case/control status of individuals.
# 
# The output vector reports the pvalues in the following order:
# "FET", "SKATO", "dCorN", "Mantel"
# 
get_detection_results    = function(path){
  

  # Create an empty numeric vector to record the pvalue obtained by permutation.
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
  

  
  names(pvalue) = c("FET", "SKATO", "dCorN", "Mantel")
  return(pvalue)
}


