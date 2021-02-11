# Load all the required functions for this script
source('3_NonIBD_Utillity_Functions.R')


# load the sample data
load('sample_data.RData')


# Run SKAT-O with new window size setting = 21 (Including target SNV)
SKATO_new = SKATO_TEST(sample_data = sample_data, CCLabel = sample_data$Genos$ccStatus, window.size = 21)


# Save new SKAT-O result
save( SKATO_new, file = "SKATO_new.RData")



# Define the number of permutation
number_of_permutation  = 1000


# Call the permute_SKATO function to run the permutataion.
# 
# Parameters of this function:
# 1. nperm is the desired number of permutations.
# 2. sample_data is the sample object generated in 1_SimulateData.R script.
start.time = Sys.time()
SKATO_permutations_new = permute_SKATO(nperm = number_of_permutation, sample_data = sample_data, window.size = 21)
end.time = Sys.time()
print(end.time - start.time )



# Add the observed statistics to the last row of the matrix
SKATO_permutations_new[number_of_permutation + 1, ] = SKATO_new$pvalue



# Save the permutation result
save( SKATO_permutations_new, file = "SKATO_perm_new.RData")








# Optional: Saving profile plots automatically
pdf(file = "SKATO_new.pdf", paper = "a4")
plot( x = SKATO_new$pos/1000, 
      y = -log(SKATO_new$pvalue, base = 10),
      xlab = "Genomic Position (Kbps)",
      main = "SKATO")
abline(v = 900,    col = "red")
abline(v = 1100,    col = "red")
dev.off()
