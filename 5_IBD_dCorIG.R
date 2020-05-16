# Load all the necessary functions for this script
source("5_IBD_Utillity_Functions.R")



# Load sample data, reconstructed partitions, and the distance matrices.
load("sample_data.RData")
load("part.RData")
load("dists.RData")



# Find the duplicated partitions and calculate the weight of each partition. 
#
# get_unique_part_and_weight() function take the partition_list and sample_data and returns:
#
# 1. Unique reconstructed partitions along the region
#
# 2. Weight of each reconstructed partition. The weight is the proportion of the genome 
#    that each reconstructed partition is spanned on.
# 
tt            = get_unique_part_and_weight(partition_list = part, sample_data = sample_data)
part          = tt$partition_list
weight.vector = tt$weights
rm(tt)




#
# dCorIG profile
#
# SeqID is a character vector and represents the label of the sequences in the sample
#
# Status is a vector, indicates the case/control status of the sample. 0 = Control , 1 = Case
#
# CClabels is a dataframe. 
# First column represents the label of the sequences.
# Second column is the case/control status of the sequences.
SeqID      = c( sample_data$CaseHapID, sample_data$ControlHapID )
Status     = sample_data$Haps$ccStatus
CClabels   = data.frame( SeqID = SeqID, Status = Status )





# dCorIG function takes the following inputs and returns the Naive dCor profile
# along the genome.
#
# 1. partition_list is a list of reconstructed partitions
# 2. w is a numeric vector representing the weight of each reconstrcuted partition.
# 3. distance_matrix_list is a list of distance matrices.
# 4. CClabels a dataframe with two columns. 
# First column is the sequence labels
# Second column is the case/control status of each sequence.
#
dCorIG_Profile = dCorIG(partition_list = part, w = weight.vector, distance_matrix_list = dists, CClabels = CClabels )
save( dCorIG_Profile, file = "dCorIG.RData")






#
# Permutation section
# Estimated time: 2 minutes for each iteration which is about 2000 minutes or 33.3 hours for 1000 permutations.
#
# We call the dCorIG_permute function to run the permutation for dCorIG
# 
number_of_permutations  = 1000
start.time = Sys.time()
dCorIG_permutation      = dCorIG_permute(nperm = number_of_permutations, sample_data, partition_list = part, w = weight.vector, distance_matrix_list = dists)
end.time = Sys.time()
print(end.time - start.time )



# Add the observed statistics to the last row of the matrix
dCorIG_permutation[number_of_permutations + 1,]  = dCorIG_Profile



# Save the results of permutation
save(x = dCorIG_permutation, file = "dCorIG_perm.RData")





# Optional: Saving profile plots automatically
pdf(file = "dCorIG.pdf", paper = "a4")
plot( x = sample_data$Posn$SNV_Position/1000, 
      y = dCorIG_Profile, 
      xlab = "SNV Positions (Kbp)", 
      ylab = "dCor", 
      main = "dCorIG"  )
abline(v = 900,    col = "red")
abline(v = 1100,   col = "red")
dev.off()


