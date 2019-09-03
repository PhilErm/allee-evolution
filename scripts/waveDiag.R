# Wave diagnostics
# This is a multi-purpose script for performing a number of tests on waves to determine if they are likely pushed or pulled.

# Notes
# There are two basic criteria that we can use to see if a wave is pushed or pulled.
# The first is seeing if changes in K have an effect on speed. In a pulled wave, changes to K should have no effect
# on speed.
# The second is to see what level of "mixing" is occurring. This is achieved by tracking individual lineages in the vanguard over time,
# and seeing if they lose diversity over time. Mutation is paused for this. If there is a loss of diversity over time, then the wave is pulled. If there isn't a
# loss of diversity over time, then the wave is pushed.

# Required packages
library(beepr)
library(tidyverse)

# Required scripts
source("scripts/functions.R")

# Seeing if a change in K changes the speed of an invasion

# Loading data. This is the full set of 20 pushed-to-pulled invasions
load(file="data/pushedToPulled.RData")
pp <- simResults

# Extracting just one invasion to use as the basis for our tests
pp.inv.exp <- pp[[2]]

# Extracting the first generation, where the wave is likely pushed
pp.inv.gen1 <- pp.inv.exp[[1]]

# Extracting just the vanguard of the final generation, where the wave is likely pulled
pp.inv.gen500 <- pp.inv.exp[[500]]

# Function that extracts just the vanguard patches
distFinder <- function(data){
  subset(data, data[,1] > max(data[,1])-5)
}

# Extracting the vanguard patches
pp.inv.gen500.van <- distFinder(pp.inv.gen500)
pp.inv.gen500.van

# Taking a look at the average A value in the the vanguard
ggplot(as.data.frame(pp.inv.gen500.van), aes(x=patch,y=A)) +
  stat_summary(fun.y = mean, geom = "line")

# Taking a look at the average A value in the core
ggplot(as.data.frame(pp.inv.gen1), aes(x=patch,y=A)) +
  stat_summary(fun.y = mean, geom = "line")

# Checking invasion speed under different K levels for core invaders (probably pushed)

# Need a new function that allows us to plug an existing population into invasion and start simulation
# from there
runModelWaveDiag <- function(pop, nGens){
  genList <- vector(length = nGens, mode = 'list') # Creates an empty list that each generation's results will be saved into
  genList[[1]] <- pop # Initial population is established by importing an existing invasion in
  pb <- txtProgressBar(min = 0, max = nGens, style = 3)
  for (i in 2:nGens){ # Loop the below 2:nGens times
    setTxtProgressBar(pb, i)
    genList[[i]] <- pGrowth(genList[[i-1]]) # Takes the first member of the list and runs pGrowth on it. Puts the result into the second spot on list. Takes the second spot on list and runs pGrowth on it to nGens
  }
  close(pb)
  genList
}

# Testing the new invasion function
# Loading necessary parameters
m <- 0.5 # The probability of an individual dispersing
r <- 0.2 # The intrinsic growth rate of individuals
K <- 500 # The carrying capacity of each patch
mut.prob <- 0 # Probability a given individual will mutate
mut.sd <- 0 # SD of the Allee trait's mutation
ex.prob <- 0 # The probability of any given patch being randomly wiped out
nGens <- 100 # The number of generations of growth, dispersal, mutation, etc. per simulation

# Running the invasion
diagResCore <- runModelWaveDiag(pp.inv.gen1,nGens)

# Checking the results in terms of average A and speed
# Average A in final generation
ggplot(as.data.frame(diagResCore[[nGens]]), aes(x=patch,y=A)) +
  stat_summary(fun.y = mean, geom = "line")

# Final distance reached
max(diagResCore[[nGens]][,1])
# Increasing K slows things down. Wave is not pulled.

# Now testing the effect of different K on vanguard invasions
# First, we have to change the patch numbers of the old vanguard data so that they start at 0, rather
# than whatever distance that particular invasion reached
# Since we're just testing, I'm going to do it manually
pp.inv.gen500.van # Min distance is 197. So 197 needs to be subtracted from all patch counts
pp.inv.gen500.van[,"patch"] <- pp.inv.gen500.van[,"patch"]-197
pp.inv.gen500.van

# Loading necessary parameters
m <- 0.5 # The probability of an individual dispersing
r <- 0.2 # The intrinsic growth rate of individuals
K <- 100000 # The carrying capacity of each patch
mut.prob <- 0 # Probability a given individual will mutate
mut.sd <- 0 # SD of the Allee trait's mutation
ex.prob <- 0 # The probability of any given patch being randomly wiped out
nGens <- 25 # The number of generations of growth, dispersal, mutation, etc. per simulation

# Running the invasion
diagResVan <- runModelWaveDiag(pp.inv.gen500.van,nGens)

# Checking the results in terms of average A and speed
# Average A in final generation
ggplot(as.data.frame(diagResVan[[nGens]]), aes(x=patch,y=A)) +
  stat_summary(fun.y = mean, geom = "line")

# Final distance reached
max(diagResVan[[nGens]][,1])
# Seem to get some difference in average speed at extreme K values. But quite resistant to change.

# One way of doing things might be to specifically test a wave where everyone is at -500 and can't
# evolve, instead of these intermixed clone lines that average out to something like -500.
# That's something that can be done by the normal invasion code, so running off to do that now
# and see how speed changes over time

# Running 5 such invasions at K = 100 for 100 generations dist. reached was:
# 36, 37, 33, 32, 33
# Running 5 such invasions at K = 250 for 100 generations dist. reached was:
# 39, 41, 47, 47, 42
# Running 5 such invasions at K = 500 for 100 generations dist. reached was:
# 43, 44, 46, 43, 45
# Running 5 such invasions at K = 1000 for 100 generations dist. reached was:
# 44, 45, 45, 42, 47
# Running 5 such invasions at K = 5000 for 100 generations dist. reached was:
# 31, 32, 31, 33, 32

# Some variability. Will have to mull those results over.

# Let's do the same thing with invasions starting at A = 50.

# Running 5 such invasions at K = 100 for 100 generations dist. reached was:
# 6, 4, 8, 6, 7
# Running 5 such invasions at K = 250 for 100 generations dist. reached was:
# 28, 30, 28, 29, 29
# Running 5 such invasions at K = 500 for 100 generations dist. reached was:
# 30, 31, 30, 33, 29
# Running 5 such invasions at K = 1000 for 100 generations dist. reached was:
# 29, 13, 29, 28, 24
# Running 5 such invasions at K = 5000 for 100 generations dist. reached was:
# 23, 15, 17, 15, 12

# Would boosting r (current/default quite low at 0.2) exacerbate the effect seen? 
# Let's find out

# For -500 invasions with r = 0.5
#K 100
#49,46,54,47,51
#K 250
#61,65,64,63,62
#K 500
#68,67,69,67,65
#K 1000
#65,65,62,65,67
#K 5000 (takes a long time, only doing 3)
#54,53,54

# For A = 50 invasions with r = 0.5
#K 100
#EXT, EXT, EXT, 6, EXT
#K 250
#44,44,45,44,43
#K 500
#47,48,49,49,47
#K 1000
#48,48,48,50,50
#K 5000
#24,23,20,15,11

# A = 50 definitely more sensitive to changes in r. Both do exhibit some change however

# Running A/K inverse to see if speed changes. This maintains a monotonically decreasing function at N = 0
# All with -A +K
# +/- 100
#44,43,41,43,43
# +/- 250
#43,46,43,41,46
# +/- 500
#47,45,45,48,50
# +/- 1000
#47,46,47,45,45
# +/- 5000
#46,45,45,47,44