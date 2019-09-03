# Running the model
# Run this script to simulate invasions

# Required packages
library(beepr)

# Required scripts
source("scripts/parameters.R")
source("scripts/functions.R")

# Running simulation
simResults <- sim(N, nGens, nSims)

# Saving results
#save(simResults, file="data/simResults.RData")
beep()