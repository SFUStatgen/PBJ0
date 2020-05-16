# Load all necessary function for this script
source("4_IBD_Utilitty_Functions.R")



# Load population data, sample data, reconstructed partitions, and the distance matrices.
load("sample_data.RData")
load("part.RData")
load("dists.RData")



# dCorN profile
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
# dCor_Profile function takes the following inputs and returns the Naive dCor profile
# along the genome.
#
# 1. cc_sample is a dataframe representing the case/control status of the sequences.
# 2. distance_matrix_list is a list containing the distance matrices of sequences
#    at each reconstructed partition.
#   
dCorN  = dCor_Profile(cc_sample = CClabels, distance_matrix_list = dists )



# Save the dCorN profile
save(x = dCorN, file = "dCorN.RData")




#
# Permutation section
# Estimated time
# 20 seconds for one permutation which is about 333 minutes or 5.5 hours for 1000 permutations.



#
# We call the dCorN_permute function to run the permutation for dCorN
# 
number_of_permutations = 1000
start.time = Sys.time()
dCorN_permutation      = dCorN_permute(nperm = number_of_permutations, sample_data, distance_matrix_list = dists)
end.time = Sys.time()
print(end.time - start.time )


# Add the observed statistics to the last row of the matrix
dCorN_permutation[ number_of_permutations + 1, ]  = dCorN



# Save the results of permutation
save(x = dCorN_permutation, file = "dCorN_perm.RData")




# Optional: Saving profile plots automatically
pdf(file = "dCorN.pdf", paper = "a4")
plot( x = sample_data$Posn$SNV_Position/1000, 
      y = dCorN, 
      xlab = "SNV Positions (Kbp)", 
      ylab = "dCor", 
      main = "dCorN"  )
abline(v = 900,    col = "red")
abline(v = 1100,   col = "red")
dev.off()




