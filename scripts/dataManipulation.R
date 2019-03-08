# Data manipulation

# Required packages
library(data.table)

# Required scripts
source("scripts/parameters.R")

# Importing data
#load(file="../Outputs/simResults.RData")

# Manipulating data into formats appropriate for figures and analysis

# Speed of single invasion over time
speed1 <- lapply(simResults[[1]], function(x) max(x[,"patch"])) # Determines how far each generation got by looking at the farthest patch
speed1 <- lapply(speed1, function(x) as.data.frame(x)) # Begins conversion of data to data frame
speed1 <- rbindlist(speed1, idcol = TRUE) # Continues conversion of data to data frame
speed1 <- as.data.frame(speed1) # Finishes conversion of data to data frame 
colnames(speed1) <- c("time", "distance") # Renames columns

# Change in A trait over invasion space for single invasion
alleeRes1 <- aggregate(simResults[[1]][[nGens]][,"A"]~simResults[[1]][[nGens]][,"patch"], FUN = mean) # Means of A as a function of patch

# Change in density over invasion space for single invasion
density1 <- table(simResults[[1]][[nGens]][,"patch"]) # Tabulates how many individuals are in each patch

# Number of individuals belonging to each clone line in the population for single invasion
cloneDensity1 <- table(simResults[[1]][[nGens]][,"A"]) # Tabulates how many individuals are in each clone line in the population
cloneDensity1 <- as.data.frame(cloneDensity1) # Converts table to data frame
colnames(cloneDensity1) <- c("clone", "frequency") # Renames columns

# Preparing data concerning final state of all invasions
genFinal <- sapply(simResults, '[', nGens) # Creates a list consisting of the final generation of each invasion

# Change in A trait over invasion space for all invasions
alleeResFinal <- lapply(genFinal, function(x) aggregate(x[,"A"]~x[,"patch"], FUN = mean)) # Calculates mean A per patch for each invasion
alleeResFinal <- rbindlist(alleeResFinal, idcol = TRUE) # Begins conversion of data to data frame
alleeResFinal <- as.data.frame(alleeResFinal) # Finishes conversion of data to data frame 
colnames(alleeResFinal) <- c("invasion", "patch", "A") # Renames columns

# Change in density over invasion space for all invasions
densityFinal <- lapply(genFinal, function(x) table(x[,"patch"])) # Calculates density per patch for each invasion
densityFinal <- lapply(densityFinal, function(x) as.data.frame(x)) # Begins conversion of data to data frame
densityFinal <- rbindlist(densityFinal, idcol = TRUE) # Continues conversion of data to data frame
densityFinal <- as.data.frame(densityFinal) # Finishes conversion of data to data frame
colnames(densityFinal) <- c("invasion", "patch", "density") # Renames columns
densityFinal <- transform(densityFinal, patch = as.numeric(patch)) # Transforms patch from a factor into a numeric

# Change in A trait over time for single invasion
alleeResTime <- simResults[[1]]
alleeResTime <- lapply(alleeResTime, function(x) aggregate(x[,"A"]~x[,"patch"], FUN = mean)) # Calculates mean A per patch for each generation
alleeResTime <- rbindlist(alleeResTime, idcol = TRUE) # Begins conversion of data to data frame
alleeResTime <- as.data.frame(alleeResTime) # Finishes conversion of data to data frame 
colnames(alleeResTime) <- c("generation", "patch", "A") # Renames columns

# Change in density over time for single invasion
densityTime <- simResults[[1]]
densityTime <- lapply(densityTime, function(x) table(x[,"patch"])) # Calculates density per patch for each generation
densityTime <- lapply(densityTime, function(x) as.data.frame(x)) # Begins conversion of data to data frame
densityTime <- rbindlist(densityTime, idcol = TRUE) # Continues conversion of data to data frame
densityTime <- as.data.frame(densityTime) # Finishes conversion of data to data frame
colnames(densityTime) <- c("generation", "patch", "density") # Renames columns
densityTime <- transform(densityTime, patch = as.numeric(patch)) # Transforms patch from a factor into a numeric

# Mean speed of all invasions
speed <- rapply(simResults, function(x) max(x[,"patch"]), how = "list") # Determines how far each invasion got at each generation
invasion <- rep(seq(1,nSims,1), each = nGens) # Creates vector to prepare speed data frame
generation <- rep(seq(1,nGens,1), times = nSims) # Creates vector to prepare speed data frame
distance <- unlist(speed) # All speeds moved to vector
speed <- cbind(invasion, generation, distance) # Binds vectors together into matrix
speed <- as.data.frame(speed) # Finishes conversion of data to data frame
