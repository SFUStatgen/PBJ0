#
# Load all the necessary functions to calculate the ECDF datasets both for detection under H0
#
source('null_ECDF_Utillity_Functions.R')




#
# Define the main directory where 500 datasets are located
#
main_path = "D:/Research/Project1/Simulations/msprime/WrightFisherModel/SearchForExampleDataset_July28/"



#####  Detection ####
#
# Define a matrix to record the detection results:
# FET, SKATO, dCor, Mantel
# 
# Rows in this matrix represents each simulated dataset and columns represents each method.
# 
# We defined the order of columns as follows: 
# "FET", "SKATO", "dCor", "Mantel"
# 

null_Detection_results           = matrix(NA, nrow = 500, ncol = 4)
colnames(null_Detection_results) = c("FET", "SKATO", "dCor", "Mantel" )
for(i in 1:500){
  dir                   = paste0( main_path, 'Dataset', i, '/null/')  
  null_Detection_results[i,] = get_detection_results(dir)
  print(dir)
}
save(x = null_Detection_results, file = 'null_Detection_results_1_500.RData')


