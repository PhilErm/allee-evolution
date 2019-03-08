# Parameters
# The parameters for the model

# Invasion parameters
N <- 300 # The number of individuals in the initial population
nSpread <- 3 # The number of patches invaders the initial population is spread over
m <- 0.5 # The probability of an individual dispersing

# Growth function parameters
r <- 0.2 # The intrinsic growth rate of individuals
K <- 500 # The carrying capacity of each patch

# Trait parameters
A.mean <- 0 # The mean around which initial invaders' Allee traits are distributed
A.sd <- 10 # The SD of the distribution of initial invaders' Allee traits

# Evolution parameters
mut.prob <- 0.01 # Probability a given individual will mutate
mut.sd <- 10 # SD of the Allee trait's mutation

# Miscellaneous parameters
ex.prob <- 0 # The probability of any given patch being randomly wiped out

# Simulation parameters
nGens <- 500 # The number of generations of growth, dispersal, mutation, etc. per simulation
nSims <- 20  # The number of simulations that are run