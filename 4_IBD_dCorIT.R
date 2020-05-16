# Load all necessary functions for this script
source("4_IBD_Utilitty_Functions.R")



# Load population data, sample data, reconstructed partitions, and the distance matrices.
load("sample_data.RData")
load("dists.RData")




# dCorIT profile 
# Reclassifying case sequences based on their true carrier status
#
SeqID             = c( sample_data$CaseHapID, sample_data$ControlHapID )
Status            = Case_reclassify(sample_data)
CClabels          = data.frame( SeqID = SeqID, Status = Status )
dCorIT            = dCor_Profile(cc_sample = CClabels, distance_matrix_list = dists )



# Save the results
save(x = dCorIT, file = "dCorIT.RData")




# Optional: Saving profile plots automatically
pdf(file = "dCorIT.pdf", paper = "a4")
plot( x = sample_data$Posn$SNV_Position/1000, 
      y = dCorIT, 
      xlab = "SNV Positions (Kbp)", 
      ylab = "dCor", 
      main = "dCorIT"  )
abline(v = 900,    col = "red")
abline(v = 1100,   col = "red")
dev.off()


