# Load all necessary function for this script
source("4_IBD_Utilitty_Functions.R")



# Load population data, sample data, reconstructed partitions, and the distance matrices.
load("sample_data.RData")
load("part.RData")
load("dists.RData")



# Mantel profile
#
# SeqID is a character vector and represents the label of the sequences in the sample
#
# Status is a vector, indicates the case/control status of the sample. 0 = control , 1 = case
#
# CClabels is a dataframe. 
# First column represents the label of sequences
# Second column is the case/control status of the sequences
SeqID      = c( sample_data$CaseHapID, sample_data$ControlHapID )
Status     = sample_data$Haps$ccStatus
CClabels   = data.frame( SeqID = SeqID, Status = Status )



#
# Mantel_Profile function takes the following inputs and returns the Mantel profile
# along the genome.
#
# 1. cc_sample is a dataframe representing the case/control status of the sequences.
# 2. distance_matrix_list is a list containing the distance matrices of sequences
#    at each reconstructed partition.
#   
# 2. distance_matrix_list: a list of distance matrices calculated based on partitions
#
Mantel  = Mantel_Profile(cc_sample = CClabels, distance_matrix_list = dists )



# Save the Mantel profile
save(x = Mantel, file = "Mantel.RData")





#
# We call the Mantel_permute function to run the permutation for Mantel
# 
number_of_permutations = 1000
start.time = Sys.time()
Mantel_permutation      = Mantel_permute(nperm = number_of_permutations, sample_data, distance_matrix_list = dists)
end.time = Sys.time()
print(end.time - start.time )


# Add the observed statistics to the last row of the matrix
Mantel_permutation[ number_of_permutations + 1, ]  = Mantel



# Save the results of permutation
save(x = Mantel_permutation, file = "Mantel_perm.RData")




# Optional: Saving profile plots automatically
pdf(file = "Mantel.pdf", paper = "a4")
plot( x = sample_data$Posn$SNV_Position/1000, 
      y = Mantel, 
      xlab = "SNV Positions (Kbp)", 
      ylab = "Mantel", 
      main = "Mantel"  )
abline(v = 900,    col = "red")
abline(v = 1100,   col = "red")
dev.off()




