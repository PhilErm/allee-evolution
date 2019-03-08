# Analysis figures
# A range of figures exploring different characteristics of the simulated population

# Required packages
library(ggplot2)
library(latex2exp)

# Required scripts
source("scripts/parameters.R")
source("scripts/dataManipulation.R")

# Change in A trait over invasion space for single invasion
plot(alleeRes1, ylab = "Allee effect vulnerability", xlab = "Patch", type = "l", main = paste("Mean A trait after", nGens, "generations of single invasion", sep=" "))

# Change in density over invasion space for single invasion
plot(density1, type = "l", ylab = "Density", xlab = "Patch", main = paste("Density of individuals after", nGens, "generations of single invasion", sep = " "))

# Invasion speed for single invasion
plot(speed1, type = "l", ylab = "Distance", xlab = "Time", main = paste("Speed of single invasion", sep = " "))

# Number of individuals belonging to each clone line in the population for single invasion
#cloneDenisty1Title <- paste("Frequency of each clone line after", nGens, "generations of single invasion", sep = " ")
#ggplot(data = cloneDensity1) +
#  geom_col(mapping = aes(x = clone, y = frequency)) +
#  theme_classic() +
#  labs(title = cloneDenisty1Title)

# Change in A trait over time for single invasion
alleeResTimeTitle <- paste("Change in mean A trait over time for single invasion", sep = " ")
ggplot(data = alleeResTime) +
  geom_line(mapping = aes(x = patch, y = A, group = generation, colour = generation)) +
  theme_classic() +
  scale_colour_gradientn(colours = terrain.colors(nGens)) +
  labs(title = alleeResTimeTitle)

# Change in density over time for single invasion
densityTimeTitle <- paste("Change in density over time for single invasion", sep = " ")
ggplot(data = densityTime) +
  geom_line(mapping = aes(x = patch, y = density, group = generation, colour = generation)) +
  theme_classic() +
  scale_colour_gradientn(colours = terrain.colors(nGens)) +
  #scale_colour_gradient2(low="chartreuse4",high="indianred3",mid="goldenrod",midpoint=nGens/2) +
  labs(title = densityTimeTitle)

# Change in A trait over invasion space for all invasions
alleeResFinalTitle <- paste("Mean A trait for", nSims, "invasions after", nGens, "generations", sep = " ")
ggplot(data = alleeResFinal) +
  geom_line(mapping = aes(x = patch, y = A, group = invasion, colour = invasion)) +
  theme_classic() +
  scale_colour_gradientn(colours = terrain.colors(nSims)) +
  labs(title = alleeResFinalTitle)

# Change in density over invasion space for all invasions
densityFinalTitle <- paste("Density of individuals for", nSims, "invasions after", nGens, "generations", sep = " ")
ggplot(data = densityFinal) +
  geom_line(mapping = aes(x = patch, y = density, group = invasion, colour = invasion)) +
  theme_classic() +
  scale_colour_gradientn(colours = terrain.colors(nSims)) +
  labs(title = densityFinalTitle)

# Mean change in A trait over invasion space for all invasions
alleeResFinalColTitle <- paste("Mean A trait across", nSims, "invasions after", nGens, "generations", sep = " ")
ggplot(alleeResFinal, aes(patch, A)) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", fill = "grey") +
  stat_summary(fun.y = mean, geom = "line") + 
  theme_classic() +
  labs(title = alleeResFinalColTitle)

# Mean change in density over invasion space for all invasions
densityFinalColTitle <- paste("Mean density of individuals across", nSims, "invasions after", nGens, "generations", sep = " ")
ggplot(densityFinal, aes(patch, density)) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", fill = "grey") +
  stat_summary(fun.y = mean, geom = "line") + 
  theme_classic() +
  labs(title = densityFinalColTitle)

# Mean invasion speed for all invasions
speedColTitle <- paste("Mean speed across", nSims, "invasions", sep = " ")
ggplot(speed, aes(generation, distance)) +
  stat_summary(fun.data = mean_se, geom = "ribbon", fill = "grey") +
  stat_summary(fun.y = mean, geom = "line") +
  theme_classic() +
  labs(title = speedColTitle)