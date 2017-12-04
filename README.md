# abcsmc-batch
A batch implementation of abcsmc. If modifying for own use, your simulation and distance function should be contained in simulation.cpp
See parameter_distance.setup

# To compile

$ cd source

$ g++ main.cpp abcsmc.cpp simulation.cpp -std=gnu++11 -L/usr/local/lib/ -lgsl -lgslcblas -I/usr/local/include -o ABCSMC 

$ cp ABCSMC ../


# To run

The SMC algorithm proceeds with an automated schedule and the populations are submitted as independent batch files, with the particles distributed in parallel.

The running of the simulation is semi-automated and should be checked from time 
to time as populations can be re-run if a problem occurs halfway through the algorithm. 
Populations are submitted as a job chain with files contained in the folder JobSubmission. 
The submission is implemented in JobChain.sh using the qsub #$ -hold_jid $(JobName). 
Depending on how many populations you expect before the number of generations increases, 
you can set the number of consecutive jobs accordingly to ensure the random seed is varied between populations within the same generation.

There are seven input arguments, Seed: The gsl rng initializer. Particles: The number of particles accepted. Final: The final epsilon. Job: The job number. 
Alpha: The percentage of the distances that is used to compute the next epsilon. Nparameters: The number of parameters with prior distributions.
Time: This parameter is used for initialisation and modifying the seed between populations if naccepted<nparticles. If set to 0 the population is initialised, if > 0 then the population is resampled from the previous

# Example input arguments

./ABCSMC Seed=$SGE_TASK_ID Particles=100 Final=0.1 Job=$SGE_TASK_ID Alpha=0.95 Nparameters=9 Time=0
