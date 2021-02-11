#
# This function simulates the population. Parameters of this function:
#
# 1. prob_disease:  A numeric number between 0 and 1. 
#                   The probability of disease in the entire population of 3100 individuals.
#
# 2. causal_region: A numeric vector of size 2. Entered values have to be in base pairs.
#    First element of this vector is the start of the causal region.
#    Second element of this vector is the end of the causal region
#    
# 3. N1: A positive integer number indicating the desired number of individuals that 
#        are carrying only one copy of the cSNV  
#
# 4. equifrequent_cSNV: A positive integer number.
#    This number defines the number of equifrequent cSNVs that user specifies in the simulation.
#
# 5. Beta0 and Beta: Two real numbers. beta0 is the sporadic disease. beta is the peneterance parameter.
#    ( Refer to N0_N1_Calcl_Formula.pdf file for more details )
#
#
Simulate_Population_Data = function(prob_disease, causal_region, N1, equifrequent_cSNV, Beta0, Beta ){
  
  
  # Initialize the parameters
  start_causal_region  = causal_region[1]
  end_causal_region    = causal_region[2]
  clade_size = round(N1 / equifrequent_cSNV)
  number_individuals_in_population = 3100
  
  
  
  simulation.stop.cond = FALSE
  while(simulation.stop.cond == FALSE){
    
    # Call python to simulate tree seqeuence file according to the parameters defined in simulation_in_msprime.py
    system("python 1_simulation_in_msprime.py")
    
    # Load the msprime package functionality into R
    msprime       = reticulate::import("msprime")
    
    # Load the .trees file into tree_sequence variable
    tree_sequence = msprime$load('TreeSequence.trees')
    
    
    # Use ExtractVariants() and ExtractPositions() functions to extract variants and positions from 
    # TreeSequence.trees into R
    # Variants is a matrix of 0/1 with rows representing each SNV on the simulated genomic region and 
    # columns represents the sequences/haplotypes labels/names. 
    # 0 = ancestral allele
    # 1 = derived allele
    Variants  = ExtractVariants(  ts = tree_sequence )
    Positions = ExtractPositions( ts = tree_sequence )
    
    
    
    # Randomly pair 6200 sequences to build the population of 3100 individuals.
    # population_mapping() fucntion takes the Variants objects and return a list of 2:
    # 1. Population.Mapping is a dataframe holding the haplotype labels of each individual.
    # 2. Genotype.Matrix counts the number of copies of the minor allele for each individual in the population.
    # We temporarily save the returned list in tt. After saving the objects, we remove tt as we do not need it anymore.
    tt                 = population_mapping( var = Variants )
    Population.Mapping = tt[[1]]
    Genotype.Matrix    = tt[[2]]
    rm(tt)
    
    
    
    
    # Find the SNVs in the causal region and save their name in SNV_in_causal_region
    # Calculate their derived allele frequency (daf) and derived allele percentage.
    indx                      = which( Positions$SNV_Position >= start_causal_region & 
                                       Positions$SNV_Position <= end_causal_region)
    SNV_in_causal_region      = Positions[indx,'SNV_Names']
    SNV_in_causal_region      = as.character(SNV_in_causal_region)
    derived.allele.frequency  = rowSums(Variants[SNV_in_causal_region,]) / ncol(Variants) 
    derived.allele.percentage = derived.allele.frequency * 100 
    
    
    
    # Number of sequences that are carrying any copy of the SNVs in the causal region.
    number_of_carrier_seq = sort( rowSums( Variants[SNV_in_causal_region, ]) )
    
    
    # We call the causal_SNV_selection() function to choose equifrequent cSNVs with a defined
    # minimum clade size.  
    cSNV = causal_SNV_selection(x = number_of_carrier_seq, CladeSize = clade_size, n_cSNV = equifrequent_cSNV)

    
    # Updating the values in the logistic regression model
    S_j = colSums( Genotype.Matrix[cSNV, ] ) 
    l_j = Beta0 + (Beta * S_j)
    p_j = exp(l_j) / ( 1 + exp(l_j) )
    p   = sum(p_j) / number_individuals_in_population
    
    
    # Check if p is close enough to prob_disease (population disease prevelance)
    # If not, add more cSNVs by calling fill_gap function to fill the gap.
    check = abs(p - prob_disease) <= 0.003
    if( !check ){
      cand  = names(number_of_carrier_seq)[!(names(number_of_carrier_seq) %in% cSNV)]
      newcSNV = fill_gap(current = cSNV, candidate = cand, GM = Genotype.Matrix, b0 = Beta0, b = Beta, N = number_individuals_in_population, prob = prob_disease, ep = 0.003)
      cSNV    = c(cSNV, newcSNV)
      S_j = colSums( Genotype.Matrix[cSNV, ] ) 
      l_j = Beta0 + (Beta * S_j)
      p_j = exp(l_j) / ( 1 + exp(l_j) )
      p =  sum(p_j) / number_individuals_in_population
    }
    
    
    
    # Distribution of the number of cSNVs that individuals are carrying in the population
    cSNV_Carrier_Frequency = table(S_j)
    print(cSNV_Carrier_Frequency)
    
    # We stop the simulation as long as the number of individuals carrying one copy of cSNV are close enough to N1
    cSNV_Carrier_Frequency <- data.frame( cSNV_Carrier_Frequency, stringsAsFactors = FALSE )    
    cSNV_Carrier_Frequency$S_j <- as.numeric( levels(cSNV_Carrier_Frequency$S_j) )[cSNV_Carrier_Frequency$S_j]

        
    # Check if cSNV_Carrier_Frequency only contains {0, 1} or {0,1,and 2}
    # It is almost impossible to see only 0 and 2. So the possible outcomes are: {0,1} or {0,1,2} or {0,1,2,x}
    # where x is greater than 3. So checking if the sum is 1 or 3 gaurantees that we have {0,1} or {0,1,2}
    
    if( sum(cSNV_Carrier_Frequency[,1]) == 1 ){  
      
      Freq_1 <- cSNV_Carrier_Frequency[ which( cSNV_Carrier_Frequency[,1] == 1 ) , 2 ]
      
      if( Freq_1 %in% c(N1 - 3, N1 - 2, N1 - 1, N1, N1 + 1, N1 + 2, N1 + 3) ){ 
        simulation.stop.cond <- TRUE 
      }

    }

    if( sum(cSNV_Carrier_Frequency[,1]) == 3 ){  
    
      Freq_1 <- cSNV_Carrier_Frequency[ which( cSNV_Carrier_Frequency[,1] == 1 ) , 2 ]
      Freq_2 <- cSNV_Carrier_Frequency[ which( cSNV_Carrier_Frequency[,1] == 2 ) , 2 ]

      if( ( Freq_1 %in% c(N1 - 4, N1 - 3, N1 - 2, N1 - 1, N1, N1 + 1, N1 + 2, N1 + 3) ) & (Freq_2 == 1) ){ 
        simulation.stop.cond <- TRUE 
      }
      
      if( ( Freq_1 %in% c(N1 - 4, N1 - 3, N1 - 2, N1 - 1, N1, N1 + 1, N1 + 2, N1 + 3) ) & (Freq_2 == 2) ){ 
        simulation.stop.cond <- TRUE 
      }
      
      
    }

    
  }
  
  
  
  # Randomly assign the case/control status to the individuals in the population.
  # The vector p_j is a numeric vector of population size. Each element of this vector reoresents the 
  # probability of getting the disease for each individual in the population.
  BinaryTrait        = rbinom(n = ncol(Genotype.Matrix), size = 1, prob = p_j)
  names(BinaryTrait) = colnames(Genotype.Matrix)
  CaseIndividuals    = names( BinaryTrait[ BinaryTrait == 1 ] )
  ControlIndividuals = names( BinaryTrait[ BinaryTrait == 0 ] )
  DISCRETE           = list(CaseIndividuals = CaseIndividuals, ControlIndividuals = ControlIndividuals, BinaryTrait = BinaryTrait)
  
  
  
  # Return the simulated population as a list of 7.
  pop = list( Variants  = Variants,
              Positions = Positions,
              Population.Mapping = Population.Mapping,
              Genotype.Matrix    = Genotype.Matrix,
              causal_region      = causal_region,
              cSNV = cSNV,
              DISCRETE = DISCRETE)
  
  
  return( SimulatedPopulation = pop )
  
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






# causal_SNV_selection() chooses causal SNVs from the causal region. 
# It takes the followings as inputs:
# 
# 1. x = Must be a numeric vector with names attribute showing the number of sequences that are carrying any SNV 
#        in the causal region. Example:  
#
#  SNV_1090  SNV_1102  SNV_1113  SNV_1117  SNV_1124  SNV_1129  .... 
#      1         4         9         14         26      30     ....
# 
# 
# 2. CladeSize: An integer number defining the size of each clade.
# 
#
# 3. n_cSNV: The number of equifrequent cSNVs that should be selected.
#
causal_SNV_selection = function(x, CladeSize, n_cSNV){
  
  
  # Find the positions in x that are as close as possible to CladeSize. It could return more than one position.
  minimum_positions = which( abs(x - CladeSize) == min( abs(x - CladeSize) ) )
  
  # The position of x vector that is close to the CladeSize
  w = minimum_positions[ length(minimum_positions) ]
  
  
  # Select the initial position of the moving pointer
  check = average_out( x[w], x[w+1] )
  if( check >= CladeSize ){
    pointer = w
  }else{
    pointer = w + 1
  }

  
  # Defining the indices of left vector (LV) and right vector (RV) from the pointer.  
  LV                = 1:(pointer - 1)
  RV                = (pointer + 1):length(x)
  
  
  # Choosing the first cSNV
  cSNV              = names( x[pointer] )
  
  
  # Assigning the pointer to the current position and start selecting the cSNVs based on the 
  # idea of moving average. We keep selecting until we have n_cSNV in the cSNV vector.
  current.position  = pointer
  
  
  while( length(cSNV) < n_cSNV ){
    
    # Getting the average 
    check = average_out( x[ LV[length(LV)] ], x[ RV[1] ])
    
    # Check if the average is greater or less than the CladeSize and 
    # update LV, RV, and the cSNVs accordingly.
    if( check >= CladeSize ){

      cSNV_indx         = which( names(x[LV])  %in% cSNV)
      LV                = LV[ !(LV %in% cSNV_indx) ]
      indx              = LV[length(LV)]
      cSNV              = c(cSNV, names( x[indx] ) )
      current.position  = indx
      LV                = 1:(current.position - 1)
      RV                = (current.position + 1):length(x)
      
    }else{
      
      cSNV_indx         = which( names(x) %in% cSNV )
      RV                = RV[ !(RV %in% cSNV_indx) ]
      current.position  = RV[1]
      cSNV              = c(cSNV, names(x[current.position]) )
      LV                = 1:(current.position - 1)
      RV                = (current.position + 1):length(x)
      
    }
    

  }

  # Return the selected cSNV
  return(cSNV)
}



#
# This function takes two numbers and returns their average. 
# This function is an internal one and is called inside causal_SNV_selection function.
#
average_out = function(x,y){
  return( mean( c(x,y) ) )
}



#
# fill_gap function is an internal function in Simulate_Population_Data(). The purpose of this function is
# to randomly choose additioanl SNVs as causal SNVs to reach the desired prob_disease in simulated population, 
# if it has not already been reached by causal_SNV_selection(). The following are the inputs to this function:
# 
# 1. current = A character vector of cSNVs returned by causal_SNV_selection function
# 2. candidate = A character vector of possible SNVs to choose from randomly.
# 3. GM = Genotype matrix with rows as SNV and columns as individuals
# 4. b0 and b: intercept and slope in logistic regression model
# 5. N = Number of individuals in population
# 6. prob = The desired probability of disease to reach in population
# 7. ep = Tolerance to reach the probability of disease (default = 0.003)
#
fill_gap = function(current, candidate, GM, b0, b, N, prob, ep = 0.003){
  
  # Setup keep and trash cSNV object
  keep.cSNV  = c()
  trash.cSNV = c()
  backup = candidate
  
  cond = FALSE
  while( cond == FALSE){
    
    # Randomly choose one SNV from candidate causal SNVs
    temp = sample(x = candidate, size = 1, replace = TRUE)
    
    # Calculate p
    S_j  = colSums( GM[c(current,keep.cSNV,temp), ] ) 
    l_j  = b0 + (b * S_j)
    p_j  = exp(l_j) / ( 1 + exp(l_j) )
    p    = sum(p_j) / N
    
    if( (p - prob) < ep ){
      candidate = candidate[ candidate != temp ]
      keep.cSNV = c(keep.cSNV, temp)
    }    
    
    if( (p - prob) > ep ){
      candidate = candidate[ candidate != temp ]
      trash.cSNV = c(trash.cSNV, temp)
    }    
    
    if( abs(p-prob) <= ep ){
      keep.cSNV = c(keep.cSNV, temp)
      cond = TRUE
    }
    
    if( length(candidate) == 0 ){
      candidate = backup
    }
  }
  
  keep.cSNV <- unique(keep.cSNV)
  return(keep.cSNV)
  
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




# This function takes ts (tree sequence), extract the variation, and returns it as a matrix.
# Rows in this matrix represents the simulated SNVs along the genome.
# Columns in this matrix represents the simulated sequences.
# The elements of this matrix are 0 and 1 where:
# 0 = Ancestral allele
# 1 = Derived allele
ExtractVariants       = function(ts){
  
  Variants = ts$genotype_matrix()
  
  if( dim(Variants)[1] != 0 ){
    colnames(Variants) = 1:ncol(Variants)
    rownames(Variants) = paste0('SNV_', 1:nrow(Variants) )
  }
  
  return(Variants)
}






#
# This function takes ts (tree sequence), extracts the position of SNVs on the genome, 
# and returns it as a data frame.
# The first column of this dataframe is "SNV_Names" which is the assigned labeles to each SNV to distinguish
# them by their names. 
# The "SNV_Position" column represents the position of each SNV along the genome in base pairs.
#
ExtractPositions      = function(ts){
  
  tables              = ts$dump_tables()
  SNVPosition         = tables$sites$position 
  SNVNames            = paste0('SNV_' , 1:length(SNVPosition) )
  Positions           = data.frame(SNVNames, SNVPosition)
  colnames(Positions) = c('SNV_Names', 'SNV_Position')
  
  return(Positions = Positions)
  
}






# This function takes the variants matrix, randomly pair haplotypes/sequences to build the population.
# This function returns a list of two: 
# 
# 1) pm or population mapping: This is a dataframe representing the population mapping. 
#    Each row belongs to one unique individual. The two columns holds the label of sequences 
#    that belongs to that individual.
# 
# 2) gm or genotype matrix: Counts the number of copies of minor allele for each individual in the population.
#    This is a matrix where rows represents simulated SNVs on the genome and columns represents the individuals.
# 
population_mapping    = function(var){
  
  sample_size = ncol(var)
  temp = split( sample(colnames(var) , replace = FALSE), rep(1:(sample_size/2),each=2) )
  col1 = unlist( lapply(temp, '[[', 1) )
  col2 = unlist( lapply(temp, '[[', 2) )
  indx = 1:(sample_size/2)
  
  # pm is the population mapping
  pm = data.frame(indx, col1, col2, stringsAsFactors = FALSE)
  colnames(pm) = c("Index", "Haplotype1", "Haplotype2")
  
  # gm is the genotype matrix.
  # The genotype matrix is prepared by calling the Genotype() function.
  gm = Genotype(var, pm = pm)
  
  return( list(pm, gm) )
}





#
# This function takes var (variants) and pm (population mapping) and returns the genotype matrix. 
# This is an internal function called in population_mapping() function.
#
Genotype              = function(var, pm){
  
  new.column.order = apply( X = pm , MARGIN = 1, FUN = function(x){ c(x[2] , x[3]) } ) 
  new.column.order = as.character(new.column.order)
  Variants = var[,new.column.order]
  
  Genotype.Matrix = matrix(NA, nrow = nrow(Variants), ncol = (ncol(Variants)/2), TRUE)
  loop = c( seq(1,ncol(Variants)-1, 2) )
  j = 0
  for( i in loop ){
    j = j + 1
    Genotype.Matrix[,j] = as.numeric(Variants[,i]) + as.numeric(Variants[,(i+1)])
  }
  
  rownames(Genotype.Matrix) = rownames(Variants)
  colnames(Genotype.Matrix) = paste0('IND',1:(ncol(Genotype.Matrix)))
  
  return( Genotype.Matrix ) 
}





#
# A utility function to report the number of cSNV that each individual is carrying.
#
cSNV.count           = function(x, POP){
  
  hap = unlist( POP[ gsub(x = x, pattern = "IND" ,replacement =  "") , 2:3 ] )
  hap1 = hap[1]
  hap2 = hap[2]
  names(hap1) = names(hap2) = NULL
  a = t( Variants[cSNV, c(hap1,hap2)] )
  return( sum( rowSums(a) ) )
  
}



