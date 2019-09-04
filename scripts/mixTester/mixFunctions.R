# Mix functions
# The functions that are used in the simulation when modified for examining mixing.

# Modified Haond et al. (2018) function for growth subject to an Allee effect
alleeFunc <- function(u, t, K, r){
  exp(4*r*(K/(K-t)^2)*(1-(u/K))*(u-t))
}

# Creating a population matrix
initPop <- function(N){
  patch <- rep(0:(nSpread-1), length.out = N) # Spreads the starting population over a selected number of patches
  n <- ave(patch, patch, FUN = length) # The density of the patch each individual is in is calculated
  A <- rnorm(N, mean = A.mean, sd = A.sd) # All individuals are given a value for their Allee trait (A), which is drawn from a normal distribution. The value indicates how sensitive they are to the Allee effect
  NM <- patch # Individuals in each starting patch contain a discrete value for a neutral mutation, NM. For simplicity, that value is the same as the number of the patch they started in
  pop <- cbind(patch, n, A, NM) # Binds all the vectors into a matrix
  pop # The initial population
}

# Dispersal
ranWalk <- function(patch, m, dist){
  newPatch <- patch+sample(x = seq(from = 0-dist, to = 0+dist, by = 1), prob = c(rep(m/2/dist, dist),1-m,rep(m/2/dist, dist)), replace = TRUE, size = length(patch)) # Individuals migrate dist (or less than dist) patches backwards or forwards (or not at all) according to equation and value of m
  newPatch[newPatch<0] <- 0 # Individuals that enter patch -1 are bounced back to patch 0
  newPatch
}

# Mutating individuals with a mut.prob chance per individual
mutation <- function(pop){
  mutate <- rbinom(n = nrow(pop), size = 1, prob = mut.prob) # Flips a weighted coin to determine who will mutate
  mutate <- mutate==1
  pop[mutate,"A"] <- rnorm(n = sum(mutate), mean = pop[mutate,"A"], sd = mut.sd) # Individuals marked for mutation are mutated according to mut.sd
  pop # The final mutated population
}

# Reproduction and dispersal of individuals
pGrowth <- function(pop){
  indivRepro <- alleeFunc(pop[,"n"], pop[,"A"], K, r) # Calculates an expected reproductive output for each individual as a function of local density and A
  realN <- rpois(indivRepro, indivRepro) # Randomly generates number of offspring based off expected reproductive output
  off <- rep(1:nrow(pop), realN) # Prepares each individual to be replicated according to its realised reproductive output
  offMat <- pop[off,] # Replicates individuals with inheritance
  offMat <- mutation(offMat) # Mutates offspring with mut.prob probability
  offMat[,"patch"] <- ranWalk(offMat[,"patch"], m, dist) # Offspring disperse according to probability m and max dist patches
  offMat[,"n"] <- ave(offMat[,"patch"], offMat[,"patch"], FUN=length) # Calculates local density for each patch
  pop <- offMat # Replaces old population with new population consisting of offspring only (i.e. full adult mortality)
  pop # The population as it now stands
}

# Running the model
runModel <- function(N, nGens){
  genList <- vector(length = nGens, mode = 'list') # Creates an empty list that each generation's results will be saved into
  genList[[1]] <- initPop(N) # Initial population is established with N individuals and saved into the first spot on the list
  pb <- txtProgressBar(min = 0, max = nGens, style = 3)
  for (i in 2:nGens){ # Loop the below 2:nGens times
    setTxtProgressBar(pb, i)
    genList[[i]] <- pGrowth(genList[[i-1]]) # Takes the first member of the list and runs pGrowth on it. Puts the result into the second spot on list. Takes the second spot on list and runs pGrowth on it to nGens
  }
  close(pb)
  genList
}

# Running multiple realisations of the model
sim <- function(N, nGens, nSims) {
  simList <- vector(length = nSims, mode = 'list') # Creates an empty list that each invasion will be saved into
  for (i in 1:nSims) { # Loop the below 1:nSims times
    cat("Simulating invasion", i, "\n")
    simList[[i]] <- runModel(N, nGens) # Runs one realisation of the model and saves it into the list
  }
  simList
}
