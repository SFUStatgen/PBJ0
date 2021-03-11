source('GetNullSample_Utillity_Functions.R')


#
# Iteration through all datasets
#
for( dataset_index in 1:500 ){

    #
    # Load the original pop data (simulated under alternative)
    #
    load( paste0('Dataset', dataset_index, '/pop_data.RData') )
    
    
    
    #
    # Get the null population
    #
    pop_data <- Get_the_null_population(pop_data = pop_data, n_pop_case_ind = 155)
    
    
    
    
    #
    # Get the null sample
    #
    sample_data <- Sample_from_population(pop_data = pop_data, nCase = 50, nCon = 50)
    
    
    
    
    #
    # Generate the permutation index for null sample
    #
    nperm = 1000
    ccLabel = sample_data$Genos$ccStatus
    
    permute_indx = matrix(NA, nrow = nperm, ncol = length(ccLabel) )
    for(i in 1:nperm){
      permute_indx[i,] = sample( x = 1:length(ccLabel), replace = FALSE)
    }
    
    
    
    
    #
    # Make null directory and save the null population and sample
    #
    system( paste0('mkdir Dataset' , dataset_index, '/', 'null') )
    save( pop_data,     file = paste0('Dataset', dataset_index, '/null/pop_data.RData')     )
    save( sample_data,  file = paste0('Dataset', dataset_index, '/null/sample_data.RData')  )
    save( permute_indx, file = paste0('Dataset', dataset_index, '/null/permute_indx.RData') ) 
    
    
    
    
    # Optional step
    #
    # We only run this stage as it saves us time later when we pick a dataset as example dataset and we want 
    # to present it. So this stage is totally optional.
    # 
    # We can save two tables here:
    #
    # Table1. Case/Control and cSNV information in the sample
    # Table2. Case/Control and cSNV information in the population
    # 
    # 
    ###########################################################################################################################
    # Table1. Case/Control and cSNV information
    Variants                  = pop_data$Variants
    cSNV                      = pop_data$cSNV
    D_case_haplotypes         = Find_Sequences_of_Individuals(sample_data$CaseIND,    PM = pop_data$Population.Mapping)
    D_control_haplotypes      = Find_Sequences_of_Individuals(sample_data$ControlIND, PM = pop_data$Population.Mapping)
    derived.allele.frequency  = rowSums( Variants ) / ncol(Variants)
    derived.allele.percentage = derived.allele.frequency * 100
    cSNV_HapMat = Variants[cSNV, c(D_case_haplotypes , D_control_haplotypes) ]
    dim(cSNV_HapMat)
    
    t = data.frame(
      cSNV = cSNV,
      daf = derived.allele.percentage[cSNV],
      Carrier_Case_Haplo    = rowSums(cSNV_HapMat[,1:100]),  
      Carrier_Control_Haplo = rowSums(cSNV_HapMat[,101:200])  
    )
    rownames(t) = NULL
    t = t[order(t$daf , decreasing = TRUE),]
    write.table(t, file = paste0('Dataset', dataset_index, '/null/Table1_cSNV_null_sample_info.txt'), sep = ",", quote = FALSE, row.names = F)
    rm(t)
    ##############################################################################################################################
    
    
    
    
    
    
    ##############################################################################################################################
    # Table2. Case/Control and cSNV information in the population
    pop_D_case_haplotypes    = Find_Sequences_of_Individuals(ind_list = pop_data$DISCRETE[[1]], PM = pop_data$Population.Mapping)
    pop_D_control_haplotypes = Find_Sequences_of_Individuals(ind_list = pop_data$DISCRETE[[2]], PM = pop_data$Population.Mapping)
    pop_cSNV_HapMat          = Variants[cSNV, c(pop_D_case_haplotypes , pop_D_control_haplotypes) ]
    
    t = data.frame(
      cSNV = cSNV,
      daf = derived.allele.percentage[cSNV],
      Carrier_Case_Haplo    = rowSums(pop_cSNV_HapMat[,1:length(pop_D_case_haplotypes) ]),  
      Carrier_Control_Haplo = rowSums(pop_cSNV_HapMat[,(length(pop_D_case_haplotypes)+1):ncol(pop_cSNV_HapMat) ])  
    )
    rownames(t) = NULL
    t = t[order(t$daf , decreasing = TRUE),]
    write.table(t, file = paste0('Dataset', dataset_index, '/null/Table2_cSNV_null_population_info.txt'), sep = ",", quote = FALSE, row.names = F)
    rm(t)
    ##############################################################################################################################
    
    print( paste0("Currently at dataset ", dataset_index , " out of 500") )
}