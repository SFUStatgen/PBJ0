#
# This function takes the partitions along the genome as a list and returns 
# the pairwise distances between the sequences at each partition. 
#
# The output of this function is a list with the same length as the partitions.
#

get_dist_matrices = function(partition_list){
  
  n = length(partition_list)
  
  dists = list()
  for( i in 1:n ){
    dists[[i]] = perfectphyloR::rdistMatrix( dend = partition_list[[i]] )
  }
  
  return(dists)
}