# Load all the required functions for this script
source('3_NonIBD_Utillity_Functions.R')


# load the sample data
load('sample_data.RData')





#
# The SKATO function takes the following inputs and run the SKATO test:
# 1. sample_data which is a list. We generated this in 1_SimulateData.R script.
# 2. A vector indicating the case/control status of individuals. 1 = Case , 0 = Control
#
SKATO     = SKATO_TEST(sample_data = sample_data, CCLabel = sample_data$Genos$ccStatus)
save(x = SKATO, file = "SKATO.RData" )







#
# Permutation section
#
# Estimated time on a laptop with intel core i7 2.56 GHZ , 8GB RAM:
# 15 seconds for each permutation which is about 250 minutes or 4 hours for 1000 permutations. 



# Define the number of permutation
number_of_permutation  = 1000


# Call the permute_SKATO function to run the permutataion.
# 
# Parameters of this function:
# 1. nperm is the desired number of permutations.
# 2. sample_data is the sample object generated in 1_SimulateData.R script.
start.time = Sys.time()
SKATO_permutations = permute_SKATO(nperm = number_of_permutation, sample_data = sample_data)
end.time = Sys.time()
print(end.time - start.time )



# Add the observed statistics to the last row of the matrix
SKATO_permutations[ number_of_permutation + 1, ] = SKATO$pvalue



# Save the permutation result
save( SKATO_permutations, file = "SKATO_perm.RData")








# Optional: Saving profile plots automatically
pdf(file = "SKATO.pdf", paper = "a4")
plot( x = sample_data$Posn$SNV_Position/1000, 
      y = -log(SKATO$pvalue, base = 10),
      xlab = "Genomic Position (Kbps)",
      main = "SKATO")
abline(v = 900,    col = "red")
abline(v = 1100,    col = "red")
dev.off()
