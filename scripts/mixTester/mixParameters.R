# Mix parameters
# The parameters for the model when modified for examining mixing.

# Invasion parameters
nSpread <- 5 # The number of patches invaders the initial population is spread over
N <- 500*nSpread # The number of individuals per initially occupied patch
m <- 0.5 # The probability of an individual dispersing
dist <- 1 # The number of patches forwards of backwards possible to move when dispersing

# Growth function parameters
r <- 0.2 # The intrinsic growth rate of individuals
K <- 500 # The carrying capacity of each patch

# Trait parameters
A.mean <- 250 # The mean around which initial invaders' Allee traits are distributed
A.sd <- 0 # The SD of the distribution of initial invaders' Allee traits

# Mutation parameters
mut.prob <- 0 # Probability a given individual will mutate
mut.sd <- 0 # SD of the Allee trait's mutation

# Miscellaneous parameters
ex.prob <- 0 # The probability of any given patch being randomly wiped out

# Simulation parameters
nGens <- 250 # The number of generations of growth, dispersal, mutation, etc. per simulation
nSims <- 1 # The number of simulations that are run