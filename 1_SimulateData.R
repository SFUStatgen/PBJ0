# Parameters used in msprime simulation
# The following parameters have been set in simulation_in_msprime.py script.
# Ne (effective diploid population size) = 3100, European population, Tenesa et al paper
# mutation and recombination rate = 2e-8 and 1e-8  per base per generation, respectively.
# Sequence length = 2 Mbp  (2e6 bp)
# Number of sequences = 6200
# demographic_events=[msprime.SimulationModelChange(time=5000, model="hudson")]
# In 5000 generations back in time, we switch from Wright-Fisher model to Hudson CwR model.





# Step 1. Load all the necessary functions for this script
source('1_SimulateData_Utillity_Functions.R')





# Step 2. Simulate the population data
# 
# Criteria for the simulated dataset:
# 
# 1. The probability of disease in the population should be around prob_disease (as close as possible)
#    (In our current simulations it is set to 0.05)
#
# 2. Number of individuals in the population that carry one copy of cSNV should be around N1.
#    (In our current simulations this value is 155)
#    (Refer to N0_N1_Calcl_Formula.pdf file for more details)
#
# 3. Aim for 15 equifrequent causal SNVs. These cSNVs need to have clade sizes of around 10 (= 155 / 15) in the population.
#
# 4. If needed, randomly choose more causal SNVs to gaurantee that the disease probability in population is close enough to 5%.
# (0.047 - 0.053)
#
pop_data = Simulate_Population_Data(prob_disease = 0.05,
                                    causal_region = c(900000,1100000),
                                    N1 = 155, 
                                    equifrequent_cSNV = 15, 
                                    Beta0 = -10,
                                    Beta  = 16
                                    )





# Step 3. Randomly sample 50 case and 50 control individuals from the population.
#
# We use the Sample_from_population() function. Parameters of this function:
# 
# 1. pop_data is a list resulting from Simulate_Population_Data() function in the previous step.
# 
# 2. nCase is the desired number of case individuals we want to sample from this population
#
# 3. nCon is the desired number of control individuals we want to sample from this population.
# 
sample_data = Sample_from_population(pop_data, nCase = 50, nCon = 50)



# The sample_data is a list with the following objects 
# 
# Elements of this list:
# 1. Haps: This is a list of two recording all the necessary information for sequences/haplotypes:
#    1.1: sample_haps: The case/control haplotype matrix. A matrix of 0/1 where 0 = Ancestral allele and 1 = Derrived allele.
#    1.2: ccStatus: A numeric vector indicating the case/control status of the sequences. 1 = Case and 0 = Control
#
# 2. Genos: This is a list of two recording all the necessary information for individuals:
#    2.1. sample_geno: Sample genotype matrix. A matrix with values 0/1/2 which counts the number of copies of the minor allele
#                      in each individual.
#    2.2. ccStatus: A numeric vector indicating the case/control status of each individual. 1 = Case and 0 = Control
#
# 3. Posn: A data.frame with two columns recording the label of SNVs and also their position on the genome.
#
# 4. poly_cSNV: A character vector recording the name of the polymorphic causal SNV in the sample.
# 
# 5. CaseIND: A character vector recording the label of each case individual in the sample.
# 
# 6. ControlIND: A character vector recording the label of each control individual in the sample.
# 
# 7. CaseHapID: A character vector recording the label of case sequences in the sample.
# 
# 8. ControlHapID: A character vector recording the label of control sequences in the sample.
#
# We can extract each element of this list by $ sign. 
# For example, sample_data$genos extracts the sample genotype matrix.










# 
# Step4. Create the permutation indices for all methods:
# 
# permute_indx is a matrix of indices for each permutation. 
#
nperm = 1000
ccLabel = sample_data$Genos$ccStatus


permute_indx = matrix(NA, nrow = nperm, ncol = length(ccLabel) )
for(i in 1:nperm){
  permute_indx[i,] = sample( x = 1:length(ccLabel), replace = FALSE)
}






# Step 5. Save the population data as pop_data.RData.
#         Save the sample data as sample_data.RData.
#         Save the permutation indices.
#         We need this data in the next scripts where we do IBD and none-IBD based methods.
#
#
# Important note:
# Make sure that all the files and scripts for each simulated dataset are at the same directory. 
#

save(pop_data,         file = "pop_data.RData")
save(sample_data,      file = "sample_data.RData")
save(permute_indx,     file = 'permute_indx.RData')


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
write.table(t, file = "Table1_cSNV_sample_info.txt", sep = ",", quote = FALSE, row.names = F)
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
write.table(t, file = "Table2_cSNV_population_info.txt", sep = ",", quote = FALSE, row.names = F)
rm(t)
##############################################################################################################################



