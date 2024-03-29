#
# This function takes the following arguments.
# 1. pop_data: A list of population elements, generated by Simulate_Population_Data()
# 2. n_pop_case_ind: The desired number of case individuals in the null population.
#
Get_the_null_population = function(pop_data, n_pop_case_ind){
  
  #
  # Get the null population
  #
  population_individuals               <- names( pop_data$DISCRETE$BinaryTrait )
  null_case_individuals                <- sample(  x = population_individuals, size = n_pop_case_ind, replace = FALSE)
  null_control_individuals             <- setdiff( x = population_individuals,  y = null_case_individuals ) 
  
  
  #
  # Update the pop_data object
  #
  pop_data$DISCRETE$CaseIndividuals    <- null_case_individuals
  pop_data$DISCRETE$ControlIndividuals <- null_control_individuals
  pop_data$DISCRETE$BinaryTrait        <- c( rep(1,n_pop_case_ind), rep(0,length(population_individuals)-n_pop_case_ind) )
  names(pop_data$DISCRETE$BinaryTrait) <- c(null_case_individuals, null_control_individuals)
  
  
  #
  # return the null population
  #
  return(pop_data)
}



# This function takes a random sample from the pop_data list generated by Simulate_Population_Data() function 
# The parameters of this function:
#
# 1. pop_data: A list of population elements, generated by Simulate_Population_Data()
# 2. nCase: The desired number of case individuals in the sample.
# 3. nCon: The desired number of control individuals in the sample.
Sample_from_population   = function(pop_data, nCase, nCon){
  
  
  Variants                  = pop_data$Variants
  Population.Mapping        = pop_data$Population.Mapping
  Genotype.Matrix           = pop_data$Genotype.Matrix
  Positions                 = pop_data$Positions
  cSNV                      = pop_data$cSNV
  

    
  # Randomly select case and control from the population, respectively.
  D_case.sample    = sample( x = pop_data$DISCRETE[[1]], size = nCase ) 
  D_control.sample = sample( x = pop_data$DISCRETE[[2]], size = nCon  ) 
    
    
  # We find the sequences/haplotypes of each individual in the case and control sample.
  # Find_Sequences_of_Individuals() function takes the list of individuals label and population mapping.
  # It returns the sequence/haplotype labels of the case/control sample.
  D_case_haplotypes    = Find_Sequences_of_Individuals(ind_list = D_case.sample ,    PM = Population.Mapping)
  D_control_haplotypes = Find_Sequences_of_Individuals(ind_list = D_control.sample , PM = Population.Mapping)
    
  

  
  # Remove any SNV with no variation in the sample
  sample_variants = Variants[ ,  c(D_case_haplotypes,D_control_haplotypes) ]
  remove_indx     = which( rowSums(sample_variants) == 0 | rowSums(sample_variants) == (2*nCase + 2*nCon) )
  sample_haps     = sample_variants[-remove_indx,]
  hap_ccStatus    = c( rep(1,2*nCase) , rep(0,2*nCon) )
  Haps            = list(sample_haps = sample_haps, ccStatus = hap_ccStatus )
  
  
  
  sample_genotype_matrix = Genotype.Matrix[ , c(D_case.sample,D_control.sample) ]
  remove_indx            = which( rowSums(sample_genotype_matrix) == 0 | rowSums(sample_genotype_matrix) == (2*nCase + 2*nCon) )
  sample_geno            = sample_genotype_matrix[-remove_indx,]
  geno_ccStatus          = c( rep(1,nCase), rep(0,nCon) )
  Genos                  = list(sample_geno = sample_geno, ccStatus = geno_ccStatus) 
  
  
  Posn                   = Positions[-remove_indx,]
  poly_cSNV              = cSNV[ cSNV %in% rownames(sample_haps) ]
  
  return( list(Haps         = Haps, 
               Genos        = Genos, 
               Posn         = Posn, 
               poly_cSNV    = poly_cSNV, 
               CaseIND      = D_case.sample, 
               ControlIND   = D_control.sample,
               CaseHapID    = D_case_haplotypes,
               ControlHapID = D_control_haplotypes) )
  
}




#
# This function takes the case/control individuals labeling and returns their sequences/haplotypes.
# Each individual is diploid and has 2 sequences.
#
Find_Sequences_of_Individuals = function(ind_list, PM){
  
  a = PM[ gsub(x = ind_list,  pattern = "IND" ,replacement =  "") , 2:3 ] 
  Sequences_List = c()
  for( i in 1:nrow(a) ){
    Sequences_List = c( Sequences_List, a[i,1] , a[i,2] )
  }
  
  return(Sequences_List)
}