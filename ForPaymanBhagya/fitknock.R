# fitknock() uses tools from the SNPknock R package to fit a HMM that 
# can be used to simulate knockoff datasets. 
fitknock <- function(haps,fp="./fastPHASE",fpdir="FastPhaseDir"){
  #------------------------------------------------------------------
  # Inputs:
  #  * haps - matrix of haplotypes, with haplotypes as rows, SNVs as cols
  #  * fp - location of fastPHASE executable
  #  * fpdir - directory in which to save fastPHASE results
  # Output:
  #  * a list specifying the knockoff HMM that can be passed to SNPknock's
  #    sampleHMM() function to sample haplotypes.
  #------------------------------------------------------------------
  # Create a directory to hold the fastPHASE results, if necessary
  if(!dir.exists(fpdir)) dir.create(fpdir)
  # Create a fastPHASE input file
  hfile <- paste(fpdir,"haplotype.inp",sep="/")
  SNPknock::writeXtoInp(haps,phased=TRUE,out_file=hfile)
  # Run fastPHASE to fit the haplotype model and save output in the
  # fpdir directory. All output files are preceeded with "out_"
  outroot <- paste(fpdir,"out",sep="/")
  SNPknock::runFastPhase(fp_path=fp,X_file=hfile,
             out_path = outroot, phased=TRUE)
  # Create the data structure that SNPknock's sampleHMM requires
  r_file <- paste(outroot,"rhat.txt",sep="_")
  alpha_file <- paste(outroot,"alphahat.txt",sep="_")
  theta_file <- paste(outroot,"thetahat.txt",sep="_")
  char_file <- paste(outroot,"origchars",sep="_")
  hmm <- SNPknock::loadHMM(r_file, alpha_file, theta_file, char_file, 
                   compact=FALSE, phased=TRUE)
  return(hmm)
}