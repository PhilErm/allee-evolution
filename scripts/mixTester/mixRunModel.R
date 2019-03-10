# Mix running the model
# Run this script to simulate invasions in which it is possible to exmaine mixing.

# Required packages
library(beepr)
library(tidyverse)
library(data.table)

# Required scripts
source("scripts/mixTester/mixParameters.R")
source("scripts/mixTester/mixFunctions.R")

# Running simulation
simResults <- sim(N, nGens, nSims)

# Saving results
#save(simResults, file="data/simResults.RData")
beep()

# Figures for output analysis ####

# The first invasion's mixing by the final generation

# Data manipulation
finalGen <- simResults[[1]][[nGens]]
finalGen <- as.data.frame(finalGen)

# Plot
ggplot(finalGen,aes(x=patch, fill = as.factor(NM))) +
  geom_bar(position = "stack") +
  labs(x="Patch",y="Density") +
  scale_fill_discrete(name="Neutral trait carried by individuals")

# Simpson's diversity index in the whole population by the end of one invasion

# Data manipulation
# Function for calculating Simpson's Diversity Index from a population matrix
simpson <- function(pop){
  clones <- count(pop, NM)
  v <- NULL
  for(i in 1:nrow(clones)){
    v[i] <- (clones[i,"nn"])*((clones[i,"nn"])-1)
  }
  v <- unlist(v)
  numerator <- sum(v)
  denominator <- sum(clones[,"nn"])*(sum(clones[,"nn"])-1)
  simpsons.D <- 1 - (numerator/denominator)
  simpsons.D
}

simpson(finalGen) # The diversity index

# Simpson's diversity index in the vanguard by the end of one invasion

# Data manipulation
van <- subset(finalGen, finalGen[,"patch"] > max(finalGen[,"patch"])-5) # Reduces population to vanguard (front 5 patches) alone
simpson(van)

# Simpson's Diversity Index in the vanguard over time in multiple invasions

# Data manipulation
# Function that finds Simpson's Diversity Index for each generation of each invasion
simpFinder <- function(data){
  vanSub <- lapply(data, function(x) lapply(x, function(x) subset(x, x[,"patch"] > max(x[,"patch"])-5)))
  vanSub.df <- lapply(vanSub, function(x) lapply(x, function(x) as.data.frame(x)))
  vanSimp <- lapply(vanSub.df, function(x) lapply(x, function(x) simpson(x)))
  df <- data.frame(invasion = rep(1:nSims, each = nGens), generation = rep(1:nGens, times = nSims))
  flat <- unlist(vanSimp)
  vanSimp.df <- cbind(df, flat)
  vanSimp.df
}

simpInvasions <- simpFinder(simResults)

# Plot

ggplot(simpInvasions, aes(x=generation, y=flat)) +
  stat_summary(fun.y=mean, geom="line")

# The raw number of species in the vanguard over time
# Function for calculating Simpson's Diversity Index from a population matrix
species <- function(pop){
  sp <- distinct(pop, NM)
  sp.n <- nrow(sp)
}

# Function that finds species number in vanguard for each generation of each invasion
spesFinder <- function(data){
  vanSub <- lapply(data, function(x) lapply(x, function(x) subset(x, x[,"patch"] > max(x[,"patch"])-5)))
  vanSub.df <- lapply(vanSub, function(x) lapply(x, function(x) as.data.frame(x)))
  vanSpes <- lapply(vanSub.df, function(x) lapply(x, function(x) species(x)))
  df <- data.frame(invasion = rep(1:nSims, each = nGens), generation = rep(1:nGens, times = nSims))
  flat <- unlist(vanSpes)
  vanSpes.df <- cbind(df, flat)
  vanSpes.df
}

spesInvasions <- spesFinder(simResults)

# Plot

ggplot(spesInvasions, aes(x=generation, y=flat)) +
  stat_summary(fun.y=mean, geom="line")