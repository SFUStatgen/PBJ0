# Load perfectphyloR for reconstructing partitions and also finding the partition distance matrices.
library(perfectphyloR)


# Load all the required functions for this script.
source("2_part_dist_Utillity_Functions.R")


# Load the sample data
load("sample_data.RData")




#
# Step 1.
# Initialize the variables for reconstructing the partitions.
#
#
# D_sample_hap_data is the sample haplotype matrix. Rows in this matrix, represents the case/control sample.
# Columns represents the SNVs along the genome.
# Entries of this matrix are 0 = ancestral allele and 1 = derived allele.
#
# snvnames_D records the name of the SNVs.
#
# posns_D records the position of the SNVs on the geonme.
#
# hapMat_D is a hapMat object which will be passed to reconstructPPregion() from perfectphyloR package.
# 
sample_hap_data = t(sample_data$Haps$sample_haps)
snvnames        = colnames( sample_hap_data )
posns           = sample_data$Posn$SNV_Position
hapMat          = perfectphyloR::createHapMat(hapmat = sample_hap_data, 
                                              snvNames = snvnames, 
                                              hapNames = rownames(sample_hap_data),
                                              posns = posns)




#
# Step 2.
# Reconstruct the partitions with a minimum window of size 500 using perfectphyloR package in R.
# The output is a list. Each element of this list is a reconstructed partition along the genome.
# 
part              = perfectphyloR::reconstructPPregion(hapMat = hapMat, minWindow = 500)





#
# Step 3. 
# Calculate the distance matrices
# get_dist_matrices() function takes partitions as a list and return a list of distance matrices.
# 
dists = get_dist_matrices(partition_list = part)



# Check if the lengths match
length(dists) == length(part)



#
# Step 4.
# Save the reconstructed partitions and distance matrices
#
save( x = part,  file = "part.RData") 
save( x = dists, file = "dists.RData")



  