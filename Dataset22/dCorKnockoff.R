start.time = Sys.time()
# dCor knockoff-based test, based on Nknockoff simulated knockoff datasets.
Nknockoff <- 1000 #Increase to 1000 to make comparable to nperm in permutation testing.
#------------------------------------------------------------------#
# Setup, part 1: read in data and functions needed to simulate a 
# knockoff null distribution for dCorIG (and dCorN, while we're at it).
#------------------------------------------------------------------#
# 1. As no msprime datasets are uploaded to GitHub, will simulate
# my own. To do so I ran the following scripts from the PBJ0 directory:
# (i) 1_SimulateData.R, (ii) 2_part_dist.R, (iii) the first 41 lines of
# 4_IBD_dCorN.R and (iv) the first 59 lines of 5_IBD_dCorIG.R.
# The results of these scripts are stored in .RData files which
# can be loaded into my workspace.
load("sample_data.RData") # from 1_SimulateData.R
load("part.RData") # from 2_part_dist.R
load("dCorN.RData") # from 4_IBD_dCorN.R
load("dCorIG.RData") # from 5_IBD_dCorIG.R
# Load functions that do dCor-related calculations...
source("2_part_dist_Utillity_Functions.R") # Authors PN/BK
source("4_IBD_Utilitty_Functions.R") # Authors PN/BK
source("5_IBD_Utillity_Functions.R") # Authors PN/BK
source("fitknock.R") # My fitknock() function to fit the knockoff model
#------------------------------------------------------------------#
# Setup, part 2: Load packages for parallel computing and set up a "cluster" 
# object to distribute computations across the cluster. Lines 27-37 below 
# are cut-and-pasted directly from BK/PN's 4_IBD_dCorN.R script.
#------------------------------------------------------------------#
library(foreach)
library(doParallel)
# To run parallel on compute canada:
# (Ref: https://docs.computecanada.ca/wiki/R)
# Create an array from the NODESLIST environnement variable
nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
# Create the cluster with the nodes name. One process per count of node name.
# nodeslist = node1 node1 node2 node2, means we are starting 2 processes on node1, likewise on node2.
# For testing on a laptop, change the above to:  
cl = makeCluster(nodeslist, type = "PSOCK") 
registerDoParallel(cl)
#------------------------------------------------------------------#
# Create the variables, such as the case/control labels, that are 
# common to all knockoff datasets. Use the same variable names as BK/PN.
#------------------------------------------------------------------#
control.inds <- which(sample_data$Haps$ccStatus == 0)
control.haps <- t(sample_data$Haps$sample_haps[,control.inds]) # cols are haplos
nhaps <- ncol(sample_data$Haps$sample_haps)
snvnames  <- row.names(sample_data$Haps$sample_haps)
hapnames <- colnames(sample_data$Haps$sample_haps)
posns <- sample_data$Posn$SNV_Position
CClabels <- sample_data$Haps$ccStatus
SeqID <- c(sample_data$CaseHapID, sample_data$ControlHapID) 
cc_sample <- data.frame(SeqID = SeqID, Status = CClabels) 
#------------------------------------------------------------------#
# Generate knockoff null distributions for both dCorN and dCorIG
#------------------------------------------------------------------#
hmm <- fitknock(control.haps) #fills in the hmm els pInit, Q, pEmit used below.
set.seed(123)
# Create a list of functions needed by each worker in the cluster.
export.list = c("get_GNN_statistic","get_dist_matrices","reclassify_cases_gnn",
                "dCor_Profile","GNN_Calculation","get_unique_part_and_weight",
                "dCorIG")
res = foreach(i = 1:Nknockoff, .export = export.list) %dopar%{
  ko_haps <- SNPknock::sampleHMM(hmm$pInit,hmm$Q,hmm$pEmit,nhaps)
  hapMat<-perfectphyloR::createHapMat(hapmat = ko_haps,snvNames = snvnames,
                                    hapNames = hapnames,posns = posns)
  #Get the full partition for all SNVs with many redundancies
  part = perfectphyloR::reconstructPPregion(hapMat = hapMat,minWindow = 500) 
  dists <- get_dist_matrices(part) 
  mN <- max(dCor_Profile(cc_sample = cc_sample, dists))
  tt <- get_unique_part_and_weight(part, sample_data) 
  part <- tt$partition_list #Need the more efficient minimal partition for GNN
  weight.vector <- tt$weights #Also for GNN-based predn of non-carrier case seqs.
  #Now get the max dCorIG statistic over the genomic region
  mIG <- max(dCorIG(part, weight.vector, dists, CClabels = cc_sample,q_cut_off = 0.25 ))
  c(mN, mIG) #Maximum dCorN and dCorIG statistics over the genomic region.
} #end %dopar%
# Results are in the list res with one element for each knockoff replicate. 
# Coerce to a data frame with rows corresp to ko reps and save to a .csv file.
ko<- t(data.frame(res))
row.names(ko) <- NULL
colnames(ko) <- c("dCorN","dCorIG")
write.csv(ko,file="dCorNIGko.csv",quote=FALSE,row.names=FALSE)
# Calculate p-value from the knockoff null distribution. NB: Unlike in
# permutation testing, we should *not* include the observed statistic in the 
# knockoff distribution
mean(ko[,"dCorN"] >= max(dCorN))
mean(ko[,"dCorIG"] >= max(dCorIG_Profile))

# Stop the cluster
stopCluster(cl)
end.time = Sys.time()
end.time - start.time