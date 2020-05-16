# Load all the required functions for this script
source('3_NonIBD_Utillity_Functions.R')



# load the sample data
load('sample_data.RData')


#
# The Fishers.Exact.Test() takes the following inputs and run the Fisher's exact test along the genome:
# 1. sample_data which is a list. We generated this in 1_SimulateData.R script.
# 2. A vector indicating the case/control status of individuals. 1 = Case , 0 = Control
#
# FET is a list of two:
# 
# 1. pvalue is a numeric vector and represents the calculated pvalue at each position
# 2. pos is a numeric vector showing the position of the SNVs on the genome. We use this for plotting purposes only.
# 
FET           = Fishers.Exact.Test(sample_data = sample_data, CCLabel = sample_data$Genos$ccStatus)
save(x = FET, file = "FET.RData" )







#
# Permutation section
# 
# Estimated time on a laptop with intel core i7 2.56 GHZ , 8GB RAM:
# 1 second for each permutation which is about 20 minutes for 1000 permutations. 

# Define the number of permutation
number_of_permutation = 1000


# Call the permute_FET function to run the permutataion.
# 
# Parameters of this function:
# 1. nperm is the desired number of permutations
# 2. sample_data is the sample object generated in 1_SimulateData.R script.
# 
start.time = Sys.time()
FET_permutations = permute_FET(nperm = number_of_permutation, sample_data = sample_data)
end.time = Sys.time()
print(end.time - start.time )


# Add the observed statistics to the last row of the matrix
FET_permutations[ number_of_permutation + 1, ] = FET$pvalue



# Save the results of permutation
save( FET_permutations, file = "FET_perm.RData")







# Optional: Saving profile plots automatically
pdf(file = "FET.pdf", paper = "a4")
plot( x = sample_data$Posn$SNV_Position/1000, 
      y = -log(FET$pvalue, base = 10), 
      xlab = "SNV Positions (Kbp)", 
      ylab = "-log10(p-value)", 
      main = "FET"  )
abline(v = 900,    col = "red")
abline(v = 1100,   col = "red")
dev.off()

