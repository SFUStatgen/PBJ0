# Step1: Load the msprime package in Python
import msprime



# Step2: Define object ts (tree sequence)
# Simulating the ancestory of 6200 sequences, with an diploid effective population size of 3100, 
# mutation and recombination rate of 2*e-8 and 1e-8 per base per generation, respectively. 
# We use discrete time Wright-Fisher (dtwf). 
# At 5000 generations back in the time, we switch to Hudson coalescent with recombination (CwR) to complete ARG. 
ts = msprime.simulate(sample_size=6200,Ne=3100,length=2e6,recombination_rate=1e-8,mutation_rate=2e-8,model="dtwf",demographic_events=[
        msprime.SimulationModelChange(time=5000, model="hudson")])
		
		
		
# Step3: Save the simulated tree sequence.
ts.dump("TreeSequence.trees")
