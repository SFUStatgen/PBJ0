#
# Load all the necessary functions to calculate the ECDF datasets both for localization and detection
#
source('ECDF_Utillity_Functions.R')




#
# Define the main directory where 500 datasets are located
#
main_path = "D:/Research/Project1/Simulations/msprime/WrightFisherModel/SearchForExampleDataset_July28/"




#####  Localization ####
#
# Define a matrix to record the average distance of the peak from the risk region for each method.
# Rows in this matrix represents each simulated dataset and columns represents each method.
# 
# We defined the order of columns as follows: 
# "FET", "SKATO", "dCorN", "Mantel" 
# 
Localization_results = matrix(NA, nrow = 500, ncol = 4) 
colnames(Localization_results) = c("FET", "SKATO", "dCorN", "Mantel")
for(i in 1:500){
  dir          = paste0( main_path, 'Dataset', i, '/')  
  Localization_results[i,] = get_localization_results(dir)
  print(dir)
}
save(x = Localization_results, file = 'Localization_results_1_500.RData')







#####  Detection ####
#
# Define a matrix to record the calculated pvalue based on permutation for:
# FET, SKATO, dCorN, Mantel
# 
# Rows in this matrix represents each simulated dataset and columns represents each method.
# 
# We defined the order of columns as follows: 
# "FET", "SKATO", "dCorN", "Mantel"
# 

Detection_results           = matrix(NA, nrow = 500, ncol = 4)
colnames(Detection_results) = c("FET", "SKATO", "dCorN", "Mantel")
for(i in 1:500){
  dir                   = paste0( main_path, 'Dataset', i, '/')  
  Detection_results[i,] = get_detection_results(dir)
  print(dir)
}
save(x = Detection_results, file = 'Detection_results_1_500.RData')


