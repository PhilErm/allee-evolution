# Manuscript figures
# Code for generating figures used in manuscript from simulation output.
  
# Required packages
library(data.table)
library(colorspace)
library(latex2exp)
library(gridExtra)
library(grid)
library(directlabels)
library(scales)
library(tidyverse)
#library(plotly) # Uncomment if creating the animated figure in "The change in A and density over time and space for a single invasion"

# Required scripts
source("scripts/functions.R")
source("scripts/parameters.R") # Correctly generating some plots requires that the current values in parameters.R (particularly the value of nGens) are identical to those used in their particular simulation data

# Sensitivity analysis of growth function subject ot Allee effect ####

# Range of N used in plot
dens <- 0:550

# Values of A used in plot
A1 <- -500
A2 <- -250
A3 <- 0
A4 <- 125
A5 <- 250

plotA1 <- alleeFunc(dens, A1, K, r)
plotA2 <- alleeFunc(dens, A2, K, r)
plotA3 <- alleeFunc(dens, A3, K, r)
plotA4 <- alleeFunc(dens, A4, K, r)
plotA5 <- alleeFunc(dens, A5, K, r)
expectedRepro <- c(plotA1, plotA2, plotA3, plotA4, plotA5)
density <- rep(dens,5)
A.value <- rep(c(A1, A2, A3, A4, A5), each = 551)
sens <- cbind(expectedRepro, density, A.value)
sens <- as.data.frame(sens)

# Plot
p1 <- ggplot(data = sens) +
  geom_line(mapping = aes(x = density, y = expectedRepro, color = as.factor(A.value))) +
  geom_hline(yintercept = 1, linetype="dotted") +
  expand_limits(y = c(1,1.3)) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.9, 1.25)) +
  labs(x = TeX("Density ($\\textit{N_x})"), y = TeX("Expected reproductive output ($\\textit{E}(\\textit{W_{ix}})$)" )) +
  scale_color_hue(labels = c(substitute(paste(italic('A'[i]), " = -500")), substitute(paste(italic('A'[i]), " = -250")), substitute(paste(italic('A'[i]), " = 0")), substitute(paste(italic('A'[i]), " = 125")), substitute(paste(italic('A'[i]), " = 250")))) +
  theme_classic(base_size = 12) +
  theme(legend.text=element_text(size=12), legend.title=element_blank())
p1
ggsave(filename = "figures/figFuncSens.pdf", p1, width = 20, height = 10, units = "cm")

# Average mean Allee trait and density across space for typical suite of invasions ####

# Data manipulation functions
# Preparing data concerning final state of all invasions
invFinal <- function(data){
  genFinal <- sapply(data, '[', nGens) # Creates a list consisting of the final generation of each invasion
  genFinal
}

# Change in A trait over invasion space for all invasions
ADiffFinal <- function(data){
  alleeResFinal <- lapply(data, function(x) aggregate(x[,"A"]~x[,"patch"], FUN = mean)) # Calculates mean A per patch for each invasion
  alleeResFinal <- rbindlist(alleeResFinal, idcol = TRUE) # Begins conversion of data to data frame
  alleeResFinal <- as.data.frame(alleeResFinal) # Finishes conversion of data to data frame 
  colnames(alleeResFinal) <- c("invasion", "patch", "A") # Renames columns
  alleeResFinal
}

# Change in density over invasion space for all invasions
densFinal <- function(data){
  densityFinal <- lapply(data, function(x) table(x[,"patch"])) # Calculates density per patch for each invasion
  densityFinal <- lapply(densityFinal, function(x) as.data.frame(x)) # Begins conversion of data to data frame
  densityFinal <- rbindlist(densityFinal, idcol = TRUE) # Continues conversion of data to data frame
  densityFinal <- as.data.frame(densityFinal) # Finishes conversion of data to data frame
  colnames(densityFinal) <- c("invasion", "patch", "density") # Renames columns
  densityFinal <- transform(densityFinal, patch = as.numeric(patch)) # Transforms patch from a factor into a numeric
  densityFinal
}

# Loading data
load(file="data/defaultConditions.RData")
defCon <- simResults

# Manipulating data
defCon <- invFinal(defCon)
defCon.ADiff <- ADiffFinal(defCon)
densityFinal <- densFinal(defCon)
densADiff <- merge(densityFinal, defCon.ADiff)
densADiffmean<-aggregate(densADiff[,3:4], by=list(densADiff$patch), FUN=mean)
colnames(densADiffmean) <- c("patch", "density", "A") # Renames columns

# Plot
p1 <- ggplot(densADiffmean, aes(patch, density, fill = A, colour = A)) + 
  geom_col() +
  theme_classic(base_size = 12) +
  scale_fill_gradient(low="red3", high="royalblue2", breaks=c(0,-50,-100,-150,-200,-250), limits = c(-250, 0)) +
  scale_colour_gradient(low="red3", high="royalblue2", breaks=c(0,-50,-100,-150,-200,-250), limits = c(-250, 0)) +
  labs(x = "Patch", y = "Mean density") +
  guides(fill=guide_legend(title=expression("Mean "~bar(italic(A))[italic(x)]~"")),colour=guide_legend(title=expression("Mean "~bar(italic(A))[italic(x)]~"")))
p1
#ggsave(filename = "figures/figDefConFinal.pdf", p1, width = 20, height = 10, units = "cm")

# The evolution of A in the vanguard over time ####

# Data manipulation functions
# Finds the five (or failing that, below five) farmost occupied patches
distFinder <- function(data){
  subset(data, data[,1] > length(data[,1])-5) # Looks at a summarised population matrix and keeps only summaries that pertain to 5 foremost patches
}

# Finds the mean of the vanguard over time
meanFinder <- function(data){
  invMean <- lapply(data, function(x) lapply(x, function(x) aggregate(x[,"A"]~x[,"patch"], FUN = mean))) # Finds A-mean for each patch of each invasion at each generation
  vanMean <- lapply(invMean, function(x) lapply(x, function(x) distFinder(x))) # Keeps only the five farthest most patches in each invasion
  vanMean <- lapply(vanMean, function(x) rbindlist(x, idcol = "generation")) # Reduces each invasion to a single item in a list
  finMean <- rbindlist(vanMean, idcol = "invasion") # Puts all invasions into the same table
  df.Mean <- as.data.frame(finMean) # Converts the table into a data frame to facilitate plotting
  colnames(df.Mean) <- c("invasion", "generation", "patch", "A") # Renames columns
  df.Mean
}

# Loading data
load(file="data/pushedToPulledWithMixing.RData")
pp <- simResults

# Tracking vanguard A over time
# Manipulating data
pp.mean <- meanFinder(pp) # Applying meanFinder

# Plot
p1 <- ggplot(pp.mean, aes(x=generation, y=A)) +
  stat_summary(fun.y = mean, geom = "line") +
  labs(x = "Generation", y = expression("Mean "~bar(italic(A))[italic(van)]~""), tag = 'A') +
  theme_classic(base_size = 12) +
  theme(plot.tag = element_text(size=24)) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_hline(yintercept = -497, linetype="dotted") +
  scale_y_continuous(breaks=c(250,0,-250,-497,-750))
p1

# Observing growth curves for start and end vanguard A

# Create data frame to facilitate plot
dens.df <- data.frame(x = 0:550)
dens <- dens.df$x

# Calculating final mean of Vanguard
finMean <- aggregate(pp.mean[,"A"]~pp.mean[,"patch"], FUN = mean) # Calculating A mean of patches
finMean <- subset(finMean, finMean[,1] > length(finMean[,1])-6)
finMean <- mean(finMean[,2])

# Plot
p2 <- ggplot(dens.df, aes(x)) +
  stat_function(fun = alleeFunc, args = c(t=finMean, K=K, r=r), geom = "line", aes(colour = "finMean")) +
  stat_function(fun = alleeFunc, args = c(t=250, K=K, r=r), geom = "line", aes(colour = "initMean")) +
  theme_classic(base_size = 12) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.9, 1.25)) +
  geom_hline(yintercept = 1, linetype="dotted") +
  expand_limits(y = c(1,1.3)) +
  labs(tag = 'B', x = TeX("Density ($\\textit{N_x})"), y = TeX("Expected reproductive output ($\\textit{E}(\\textit{W_{ix}})$)" )) +
  scale_colour_manual(name = " ", labels = c(expression(~italic(A[i])~"= final mean "~bar(italic(A))[italic(van)]~""),expression(~italic(A[i])~"= initial mean "~bar(italic(A))[italic(van)]~"")), values = c("#F8766D", "#619CFF")) +
  theme(legend.text=element_text(size=12), legend.title=element_blank(), plot.tag = element_text(size=24)) +
  theme(legend.text.align = 0)
p2

# Multiplot
p1 <- arrangeGrob(p1)

p2 <- arrangeGrob(p2)

f1 <- grid.arrange(p1, p2, ncol = 1)
f1
ggsave(filename = "figures/figPP.pdf", f1, width = 20, height = 20, units = "cm")

# The change in A and density over time and space for a single invasion ####

# Data manipulation functions
# Change in A trait over invasion space for single invasion
singleASpace <- function(data){
  alleeRes1 <- aggregate(data[[1]][[nGens]][,"A"]~data[[1]][[nGens]][,"patch"], FUN = mean) # Means of A as a function of patch
}

# Change in A trait over time
singleATime <- function(data){
  alleeResTime <- data[[1]]
  alleeResTime <- lapply(alleeResTime, function(x) aggregate(x[,"A"]~x[,"patch"], FUN = mean)) # Calculates mean A per patch for each generation
  alleeResTime <- rbindlist(alleeResTime, idcol = TRUE) # Begins conversion of data to data frame
  alleeResTime <- as.data.frame(alleeResTime) # Finishes conversion of data to data frame 
  colnames(alleeResTime) <- c("Generation", "patch", "A") # Renames columns
  alleeResTime
}

# Change in density over invasion space for single invasion
singleDensSpace <- function(data){
  density1 <- table(data[[1]][[nGens]][,"patch"]) # Tabulates how many individuals are in each patch
  density1 <- as.data.frame(density1)
  density1[,1] <- as.numeric(density1[,1])
  density1
}

# Change in density over time for single invasion
singleDensTime <- function(data){
  densityTime <- data[[1]]
  densityTime <- lapply(densityTime, function(x) table(x[,"patch"])) # Calculates density per patch for each generation
  densityTime <- lapply(densityTime, function(x) as.data.frame(x)) # Begins conversion of data to data frame
  densityTime <- rbindlist(densityTime, idcol = TRUE) # Continues conversion of data to data frame
  densityTime <- as.data.frame(densityTime) # Finishes conversion of data to data frame
  colnames(densityTime) <- c("Generation", "patch", "density") # Renames columns
  densityTime <- transform(densityTime, patch = as.numeric(patch)) # Transforms patch from a factor into a numeric
  densityTime
}

# Loading data
load(file="data/defaultConditions.RData")
defCon <- simResults

# Manipulating data
sing.ATime <- singleATime(defCon)
sing.densTime <- singleDensTime(defCon)

# Plot
p1 <- ggplot(data = sing.ATime) +
  geom_line(mapping = aes(x = patch, y = A, group = Generation, colour = Generation)) +
  theme_classic(base_size = 12) +
  scale_colour_viridis_c() +
  labs(x = "Patch", y = expression(""~bar(italic(A))[italic(x)]~""), fill = "Generation", tag = "A") +
  theme(axis.title.x=element_blank(), plot.tag = element_text(size=24)) +
  theme(legend.position="none")
p1

p2 <- ggplot(data = sing.densTime) +
  geom_line(mapping = aes(x = patch, y = density, group = Generation, colour = Generation)) +
  theme_classic() +
  scale_colour_viridis_c() +
  labs(x = "Patch", y = "Density", fill = "Generation", tag = "B") +
  theme(legend.position='bottom', plot.tag = element_text(size=24))
p2

# # Animated plot
# # Load the package 'plotly' to run this
# sing.densTime$A <- sing.ATime$A
# p3 <- ggplot(data = sing.densTime) +
#   geom_point(mapping = aes(x = patch, y = density, colour = A, frame = Generation)) +
#   theme_classic() +
#   scale_colour_viridis_c() + 
#   labs(x = "Patch", y = "Density", colour = "Mean A trait value")
# p3 <- ggplotly(p3)
# p3

# Multiplot
p1 <- arrangeGrob(p1)

p2 <- arrangeGrob(p2)

f1 <- grid.arrange(p1, p2, ncol = 1)
f1
#ggsave(filename = "figures/figSingInvTimeSpace.pdf", f1, width = 20, height = 20, units = "cm")


# Sensitivity analysis of invasion parameters ####

# Data manipulation functions
ADiffFinder <- function(data){
  genFinal <- sapply(data, '[', nGens) # Creates a list consisting of the final generation of each invasion
  alleeResFinal <- lapply(genFinal, function(x) aggregate(x[,"A"]~x[,"patch"], FUN = mean)) # Calculates mean A per patch for each invasion
  alleeResFinal <- rbindlist(alleeResFinal, idcol = TRUE) # Begins conversion of data to data frame
  alleeResFinal <- as.data.frame(alleeResFinal) # Finishes conversion of data to data frame 
  colnames(alleeResFinal) <- c("invasion", "patch", "A") # Renames columns
  alleeResFinal # Mean of each patch in all invasions
  core5 <- alleeResFinal[alleeResFinal[,"patch"]<= 4,]
  maxDist <- split(alleeResFinal, alleeResFinal$invasion)
  maxDist <- as.data.frame(sapply(maxDist, function(x) max(x[,"patch"])))
  maxDist <- cbind(maxDist, rep(1:20, nSims, length = nSims))
  colnames(maxDist) <- c("speed", "invasion") # Renames columns
  alleeResFinal <- merge(alleeResFinal, maxDist)
  top5bottom5 <- subset(alleeResFinal, patch <= 4 | patch >= speed-4) # The key. Keeps first 5 patches and last 5 patches for each invasion. Might need to change it to get rid of very last patch.
  top5bottom5$speed <- NULL
  coreFront <- ifelse(top5bottom5$patch <= 4, "core", "front")
  top5bottom5 <- cbind(top5bottom5, coreFront)
  top5bottom5 # Average for each invasion's first 5 and last 5 patches.
}

ADiffByScen <- function(data){
  ADiff.data <- ADiffFinder(data)
  scenario <- rep(deparse(substitute(data)), length = nrow(ADiff.data))
  ADiff.data <- cbind(ADiff.data, scenario)
  ADiff.data
}

ADiffPlot <- function(data){
  ggplot(data = data) +
    geom_boxplot(mapping = aes(x = scenario, y = A, fill = coreFront)) +
    geom_hline(yintercept = 0, linetype="dotted") +
    theme_classic(base_size = 14)
}

speedByScen <- function(data){
  speed <- rapply(data, function(x) max(x[,"patch"]), how = "list") # Determines how far each invasion got at each generation
  invasion <- rep(seq(1,nSims,1), each = nGens) # Creates vector to prepare speed data frame
  generation <- rep(seq(1,nGens,1), times = nSims) # Creates vector to prepare speed data frame
  distance <- unlist(speed) # All speeds moved to vector
  speed <- cbind(invasion, generation, distance) # Binds vectors together into matrix
  speed <- as.data.frame(speed) # Finishes conversion of data to data frame
  scenario <- rep(deparse(substitute(data)), length = nrow(speed))
  speed <- cbind(speed, scenario)
  speed
}

speedPlot <- function(data){
  ggplot(data, aes(x = generation, y = distance, group = scenario, colour = scenario, fill = scenario)) +
    stat_summary(fun.y = mean, geom = "line") + 
    labs(x = "Generation", y = "Mean distance") +
    theme_classic(base_size = 14)
  #theme(legend.position="none")
}

# For extracting a legend as a plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# Function for sharing legend between plots
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}

# Response to r

# Loading data
load(file="data/defConR0_1.RData")
r0_1 <- simResults
load(file="data/defConR0_2.RData")
r0_2 <- simResults
load(file="data/defConR0_3.RData")
r0_3 <- simResults

# Manipulating data
r.ADiff <- rbind(ADiffByScen(r0_1),ADiffByScen(r0_2),ADiffByScen(r0_3))
r.Speed <- rbind(speedByScen(r0_1),speedByScen(r0_2),speedByScen(r0_3))
rm(r0_1, r0_2, r0_3)

# Manipulating figures
p1 <- ADiffPlot(r.ADiff)
p1 <- p1 + scale_x_discrete(labels=c("0.1","0.2","0.3"))
p1 <- p1 + labs(x = TeX("$\\textit{r}$"), y = expression(""~bar(italic(A))[italic(fin)]~""))
p1 <- p1 + scale_fill_discrete(labels = c("Core", "Vanguard"))
p1 <- p1 + theme(legend.title=element_blank()) + theme(legend.position="none")
p1

s1 <- speedPlot(r.Speed)
#s1 <- s1 + geom_dl(aes(label=scenario),method=list(box.color = NA, "angled.boxes"))
s1 <- s1 + theme(legend.title=element_blank())
s1 <- s1 + scale_color_hue(labels = c(substitute(paste(italic('r'), " = 0.1")), substitute(paste(italic('r'), " = 0.2")), substitute(paste(italic('r'), " = 0.3"))))
s1 <- s1 + theme(legend.text.align = 0) + theme(axis.title.x=element_blank())
s1

# Response to K*

# Loading data
load(file="data/defConK250.RData")
K250 <- simResults
load(file="data/defConK500.RData")
K500 <- simResults
load(file="data/defConK750.RData")
K750 <- simResults

# Manipulating data
K.ADiff <- rbind(ADiffByScen(K250),ADiffByScen(K500),ADiffByScen(K750))
K.Speed <- rbind(speedByScen(K250),speedByScen(K500),speedByScen(K750))
rm(K250, K500, K750)

# Manipulating figures
p2 <- ADiffPlot(K.ADiff)
p2 <- p2 + scale_x_discrete(labels=c("250","500","750"))
p2 <- p2 + labs(x = TeX("$\\textit{K}$"), y = expression(""~bar(italic(A))[italic(fin)]~""))
p2 <- p2 + scale_fill_discrete(labels = c("Core", "Vanguard"))
p2 <- p2 + theme(legend.title=element_blank()) + theme(legend.position="none") + theme(axis.title.y=element_blank())
p2

s2 <- speedPlot(K.Speed)
#s2 <- s2 + geom_dl(aes(label=scenario),method=list(box.color = NA, "angled.boxes"))
s2 <- s2 + theme(legend.title=element_blank())
s2 <- s2 + scale_color_hue(labels = c(substitute(paste(italic('K'), " = 250")), substitute(paste(italic('K'), " = 500")), substitute(paste(italic('K'), " = 750"))))
s2 <- s2 + theme(legend.text.align = 0) + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())
s2

# Response to p_mutation

# Loading data
load(file="data/defConMut0.RData")
mut0 <- simResults
load(file="data/defConMut0_001.RData")
mut0_001 <- simResults
load(file="data/defConMut0_01.RData")
mut0_01 <- simResults
load(file="data/defConMut0_1.RData")
mut0_1 <- simResults

# Manipulating data
mut.ADiff <- rbind(ADiffByScen(mut0),ADiffByScen(mut0_001),ADiffByScen(mut0_01),ADiffByScen(mut0_1))
mut.Speed <- rbind(speedByScen(mut0),speedByScen(mut0_001),speedByScen(mut0_01), speedByScen(mut0_1))
rm(mut0, mut0_001, mut0_01, mut0_1)

# Manipulating figures
p3 <- ADiffPlot(mut.ADiff)
p3 <- p3 + scale_x_discrete(labels=c("0","0.001","0.01","0.1"))
p3 <- p3 + labs(x = TeX("$\\textit{p_{mut}}$"), y = expression(""~bar(italic(A))[italic(fin)]~""))
p3 <- p3 + scale_fill_discrete(labels = c("Core", "Vanguard"))
p3 <- p3 + theme(legend.title=element_blank()) + theme(legend.position="none")
p3

s3 <- speedPlot(mut.Speed)
#s3 <- s3 + geom_dl(aes(label=scenario),method=list(box.color = NA, "angled.boxes"))
s3 <- s3 + theme(legend.title=element_blank())
s3 <- s3 + scale_color_hue(labels = c(expression(""~italic(p)[italic(mut)]~"= 0"), expression(""~italic(p)[italic(mut)]~"= 0.001"), expression(""~italic(p)[italic(mut)]~"= 0.01"), expression(""~italic(p)[italic(mut)]~"= 0.1")))
s3 <- s3 + theme(legend.text.align = 0) + theme(axis.title.x=element_blank())
s3

# Response to A

# Loading data
load(file="data/defConA-50.RData")
ANeg50 <- simResults
load(file="data/defConA-25.RData")
ANeg25 <- simResults
load(file="data/defConA0.RData")
A0 <- simResults
load(file="data/defConA25.RData")
A25 <- simResults
load(file="data/defConA50.RData")
A50 <- simResults

# Manipulating data
A.ADiff <- rbind(ADiffByScen(ANeg50),ADiffByScen(ANeg25),ADiffByScen(A0),ADiffByScen(A25),ADiffByScen(A50))
A.Speed <- rbind(speedByScen(ANeg50),speedByScen(ANeg25),speedByScen(A0),speedByScen(A25),speedByScen(A50))
rm(ANeg50, ANeg25, A0, A25, A50)

# Manipulating figures
p4 <- ADiffPlot(A.ADiff)
p4 <- p4 + scale_x_discrete(labels=c("-50","-25","0","25","50"))
p4 <- p4 + labs(x = expression(""~bar(italic(A))[italic(init)]~""), y = expression(""~bar(italic(A))[italic(fin)]~""))
p4 <- p4 + scale_fill_discrete(labels = c("Core", "Vanguard"))
p4 <- p4 + theme(legend.title=element_blank()) + theme(legend.position="none") + theme(axis.title.y=element_blank())
p4

s4 <- speedPlot(A.Speed)
#s4 <- s4 + geom_dl(aes(label=scenario),method=list(box.color = NA, "angled.boxes"))
s4 <- s4 + theme(legend.title=element_blank())
s4 <- s4 + scale_color_hue(labels = c(expression(""~bar(italic(A))[italic(init)]~"= -50"), expression(""~bar(italic(A))[italic(init)]~"= -25"), expression(""~bar(italic(A))[italic(init)]~"= 0"), expression(""~bar(italic(A))[italic(init)]~"= 25"), expression(""~bar(italic(A))[italic(init)]~"= 50")))
s4 <- s4 + theme(legend.text.align = 0) + theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank())
s4

# Response to m

# Loading data
load(file="data/defConM0_25.RData")
m0_25 <- simResults
load(file="data/defConM0_5.RData")
m0_5 <- simResults
load(file="data/defConM0_75.RData")
m0_75 <- simResults

# Manipulating data
m.ADiff <- rbind(ADiffByScen(m0_25),ADiffByScen(m0_5),ADiffByScen(m0_75))
m.Speed <- rbind(speedByScen(m0_25),speedByScen(m0_5),speedByScen(m0_75))
rm(m0_25, m0_5, m0_75)

# Manipulating figures
p5 <- ADiffPlot(m.ADiff)
p5 <- p5 + scale_x_discrete(labels=c("0.25","0.5","0.75"))
p5 <- p5 + labs(x = TeX("$\\textit{m}$"), y = expression(""~bar(italic(A))[italic(fin)]~""))
p5 <- p5 + scale_fill_discrete(labels = c("Core", "Vanguard"))
p5 <- p5 + theme(legend.title=element_blank()) + theme(legend.position="none")
p5

s5 <- speedPlot(m.Speed)
#s5 <- s5 + geom_dl(aes(label=scenario),method=list(box.color = NA, "angled.boxes"))
s5 <- s5 + theme(legend.title=element_blank())
s5 <- s5 + scale_color_hue(labels = c(substitute(paste(italic('m'), " = 0.25")), substitute(paste(italic('m'), " = 0.5")), substitute(paste(italic('m'), " = 0.75"))))
s5 <- s5 + theme(legend.text.align = 0)
s5

# Response to dist

# Loading data
load(file="data/defConDist1.RData")
Dist1 <- simResults
load(file="data/defConDist2.RData")
Dist2 <- simResults
load(file="data/defConDist3.RData")
Dist3 <- simResults

# Manipulating data
m.ADiff <- rbind(ADiffByScen(Dist1),ADiffByScen(Dist2),ADiffByScen(Dist3))
m.Speed <- rbind(speedByScen(Dist1),speedByScen(Dist2),speedByScen(Dist3))
rm(Dist1, Dist2, Dist3)

# Manipulating figures
p6 <- ADiffPlot(m.ADiff)
p6 <- p6 + scale_x_discrete(labels=c("1","2","3"))
p6 <- p6 + labs(x = TeX("$\\textit{dist}$"), y = expression(""~bar(italic(A))[italic(fin)]~""))
p6 <- p6 + scale_fill_discrete(labels = c("Core", "Vanguard"))
p6 <- p6 + theme(legend.title=element_blank()) + theme(legend.position="none") + theme(axis.title.y=element_blank())
p6

s6 <- speedPlot(m.Speed)
#s6 <- s6 + geom_dl(aes(label=scenario),method=list(box.color = NA, "angled.boxes"))
s6 <- s6 + theme(legend.title=element_blank())
s6 <- s6 + scale_color_hue(labels = c(substitute(paste(italic('dist'), " = 1")), substitute(paste(italic('dist'), " = 2")), substitute(paste(italic('dist'), " = 3"))))
s6 <- s6 + theme(legend.text.align = 0) + theme(axis.title.y=element_blank())
s6

# Labelling each plot

p1 <- arrangeGrob(p1, top = textGrob("A", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

p2 <- arrangeGrob(p2, top = textGrob("B", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

p3 <- arrangeGrob(p3, top = textGrob("C", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

p4 <- arrangeGrob(p4, top = textGrob("D", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

p5 <- arrangeGrob(p5, top = textGrob("E", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

p6 <- arrangeGrob(p6, top = textGrob("F", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

s1 <- arrangeGrob(s1, top = textGrob("A", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

s2 <- arrangeGrob(s2, top = textGrob("B", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

s3 <- arrangeGrob(s3, top = textGrob("C", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

s4 <- arrangeGrob(s4, top = textGrob("D", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

s5 <- arrangeGrob(s5, top = textGrob("E", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

s6 <- arrangeGrob(s6, top = textGrob("F", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))

# Combining all A-difference scenarios with seperate legend
# If moving to an even number of graphs, consider using grid_arrange_shared_legend code for printing one legend at bottom
pLegend <- ADiffPlot(r.ADiff)
pLegend <- pLegend + scale_x_discrete(labels=c("0.5","1","1.5"))
pLegend <- pLegend + labs(x = TeX("$\\textit{r}$"), y = expression(""~bar(italic(A))[italic(fin)]~""))
pLegend <- pLegend + scale_fill_discrete(labels = c("Core", "Vanguard"))
pLegend <- pLegend + theme(legend.title=element_blank(), legend.text=element_text(size=14))
pLegend <- pLegend + theme(legend.position="bottom")
pLegend <- g_legend(pLegend)

lay <- rbind(c(1,2),
             c(3,4),
             c(5,6),
             c(7,7))

f1 <- grid.arrange(p1, p2, p3, p4, p5, p6, pLegend, layout_matrix = lay, heights=c(10,10,10,2))
f1
ggsave(filename = "figures/figInvSensSpeedADiff.pdf", f1, width = 20, height = 20, units = "cm")


# Combining all speed scenarios with seperate legend
g1 <- grid.arrange(s1, s2, s3, s4, s5, s6, ncol = 2, nrow = 3)
g1
ggsave(filename = "figures/figInvSensSpeed.pdf", g1, width = 20, height = 20, units = "cm")

# The degree of mixing that occurs in a wave subject to a strong Allee effect versus a wave subject to no Allee effect, and for the pushed-pulled transition ####

# Data manipulation functions
# Calcultes Simpson's diversity index from a population matrix
simpson <- function(pop){
  clones <- count(pop, NM)
  v <- NULL
  for(i in 1:nrow(clones)){
    v[i] <- (clones[i,"n"])*((clones[i,"n"])-1)
  }
  v <- unlist(v)
  numerator <- sum(v)
  denominator <- sum(clones[,"n"])*(sum(clones[,"n"])-1)
  simpsons.D <- 1 - (numerator/denominator)
  simpsons.D
}

# Calculates Simpson's diversity index for each generation of each invasion
simpFinder <- function(data){
  vanSub <- lapply(data, function(x) lapply(x, function(x) subset(x, x[,"patch"] > max(x[,"patch"])-5)))
  vanSub.df <- lapply(vanSub, function(x) lapply(x, function(x) as.data.frame(x)))
  vanSimp <- lapply(vanSub.df, function(x) lapply(x, function(x) simpson(x)))
  df <- data.frame(invasion = rep(1:nSims, each = nGens), generation = rep(1:nGens, times = nSims))
  flat <- unlist(vanSimp)
  vanSimp.df <- cbind(df, flat)
  vanSimp.df
}

# Loading data
load(file="data/pushedToPulledWithMixing.RData")
pp <- simResults
load(file="data/mixingPushed.RData")
mixPushed <- simResults
load(file="data/mixingPulled.RData")
mixPulled <- simResults

# Manipulating data
simpPP <- simpFinder(pp)
simpPushed <- simpFinder(mixPushed)
simpPulled <- simpFinder(mixPulled)
simpMerged <- rbindlist(list(simpPushed, simpPP, simpPulled), idcol = "ID")

# Plots
# Simpson's diversity index in the vanguard over time in invasions subject to evolution
p1 <- ggplot(simpPP, aes(x = generation, y = flat)) +
  stat_summary(fun.y = mean, geom = "line") +
  labs(x = "Generation", y = expression(""~bar(italic(D))[italic(van)]~"")) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0.8, linetype="dotted") +
  scale_y_continuous(breaks=c(0.8,0.6,0.4,0.2,0))
p1
ggsave(filename = "figures/figMixPP.pdf", p1, width = 20, height = 20, units = "cm")

# Simpson's diversity index in the vanguard over time
p1 <- ggplot(simpMerged, aes(x = generation, y = flat, colour = ID)) +
  stat_summary(fun.y = mean, geom = "line", aes(colour = paste("mean", ID))) +
  labs(tag = "A", x = "Generation", y = expression(""~bar(italic(D))[italic(van)]~"")) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0.8, linetype="dotted") +
  scale_y_continuous(breaks=c(0.8,0.6,0.4,0.2,0)) +
  scale_colour_manual(name = " ", labels = c(expression(~bar(italic(A))[italic(init)]~"= 250"),expression(~bar(italic(A))[italic(init)]~"= 250, "~italic(A)[italic(i)]~"evolution possible"),expression(~bar(italic(A))[italic(init)]~"= -497")), values = c("#619CFF", "black", "#F8766D")) +
  theme(legend.text=element_text(size=12), legend.title=element_blank(), plot.tag = element_text(size=24)) +
  theme(legend.text.align = 0) +
  theme(legend.position = c(0.7, 0.5))
p1

# A single example of the final mixing state of a pushed wave
p2 <- ggplot(as.data.frame(mixPushed[[10]][[nGens]]), aes(x=patch, fill = as.factor(NM))) +
  geom_bar(position = "stack") +
  labs(tag = "B", x = "Patch", y = "Density") +
  scale_fill_discrete(name = "Neutral trait carried", labels = c("Trait 0 ","Trait 1 ","Trait 2 ","Trait 3 ","Trait 4 ")) +
  theme_classic(base_size = 12) +
  guides(fill=FALSE) +
  theme(legend.text=element_text(size=12), legend.title=element_blank(), plot.tag = element_text(size=24))
p2

# A single example of the final mixing state of a pulled wave
p3 <- ggplot(as.data.frame(mixPulled[[10]][[nGens]]), aes(x = patch, fill = as.factor(NM))) +
  geom_bar(position = "stack") +
  labs(tag = "C", x = "Patch", y = " ") +
  scale_fill_discrete(name = "Neutral trait carried", labels = c("Trait 0 ","Trait 1 ","Trait 2 ","Trait 3 ","Trait 4 ")) +
  theme_classic(base_size = 12) +
  guides(fill=FALSE) +
  theme(legend.text=element_text(size=12), legend.title=element_blank(), plot.tag = element_text(size=24))
p3

# Multiplot
# For extracting a legend as a plot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# Function for sharing legend between plots
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}

pz <- ggplot(as.data.frame(mixPushed[[10]][[nGens]]), aes(x=patch, fill = as.factor(NM))) +
  geom_bar(position = "stack") +
  labs(x = "Patch", y = "Density") +
  scale_fill_discrete(name = "Neutral trait carried", labels = c("Trait 0 ","Trait 1 ","Trait 2 ","Trait 3 ","Trait 4 ")) +
  theme_classic(base_size = 12) +
  theme(legend.text=element_text(size=12), plot.tag = element_text(size=24))
pLegend <- pz
pLegend <- pLegend + theme()
pLegend <- pLegend + theme(legend.position="bottom")
pLegend <- g_legend(pLegend)

a1 <- arrangeGrob(p1)

a2 <- arrangeGrob(p2)

a3 <- arrangeGrob(p3)

lay <- rbind(c(1,1),
             c(2,3),
             c(4,4))

f1 <- grid.arrange(a1, a2, a3, pLegend, layout_matrix = lay, heights=c(10,10,1))
f1
ggsave(filename = "figures/figMixing.pdf", f1, width = 20, height = 20, units = "cm")

# Animated plot
# Load the package 'plotly' to run this
sing.densTime$A <- sing.ATime$A
p3 <- ggplot(data = sing.densTime) +
  geom_point(mapping = aes(x = patch, y = density, colour = A, frame = Generation)) +
  theme_classic() +
  scale_colour_viridis_c() +
  labs(x = "Patch", y = "Density", colour = "Mean A trait value")
p3 <- ggplotly(p3)
p3

# The degree of mixing in a pushed wave for 5000 generations ####

# Data manipulation functions
# Calcultes Simpson's diversity index from a population matrix
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

# Calculates Simpson's diversity index for each generation of each invasion
simpFinder <- function(data){
  vanSub <- lapply(data, function(x) lapply(x, function(x) subset(x, x[,"patch"] > max(x[,"patch"])-5)))
  vanSub.df <- lapply(vanSub, function(x) lapply(x, function(x) as.data.frame(x)))
  vanSimp <- lapply(vanSub.df, function(x) lapply(x, function(x) simpson(x)))
  df <- data.frame(invasion = rep(1:nSims, each = nGens), generation = rep(1:nGens, times = nSims))
  flat <- unlist(vanSimp)
  vanSimp.df <- cbind(df, flat)
  vanSimp.df
}

# Loading data
load(file="data/mixingPushed5000Gens.RData")
pushed5000 <- simResults

# Manipulating data
simpP5000 <- simpFinder(pushed5000)

# Plots
# Simpson's diversity index in the vanguard over time

p1 <- ggplot(simpP5000, aes(x = generation, y = flat)) +
  stat_summary(fun.y = mean, geom = "line") +
  labs(tag = "A", x = "Generation", y = expression(""~bar(italic(D))[italic(van)]~"")) +
  theme_classic(base_size = 12) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0.8, linetype="dotted") +
  scale_y_continuous(breaks=c(0.8,0.6,0.4,0.2,0)) +
  theme(legend.text=element_text(size=12), legend.title=element_blank(), plot.tag = element_text(size=24)) +
  theme(legend.text.align = 0)
p1

# A single example of the final mixing state of a pushed wave

p2 <- ggplot(as.data.frame(pushed5000[[1]][[nGens]]), aes(x=patch, fill = as.factor(NM))) +
  geom_bar(position = "stack") +
  labs(tag = "B", x = "Patch", y = "Density") +
  scale_fill_discrete(name = "Neutral trait carried", labels = c("Trait 0 ","Trait 1 ","Trait 2 ","Trait 3 ","Trait 4 ")) +
  theme_classic(base_size = 12) +
  theme(legend.text=element_text(size=12), plot.tag = element_text(size=24),legend.position="bottom")
p2

a1 <- arrangeGrob(p1)

a2 <- arrangeGrob(p2)

lay <- rbind(c(1,1),
             c(2,2))

f1 <- grid.arrange(a1, a2, layout_matrix = lay, heights=c(10,10,1))
f1
ggsave(filename = "figures/figMixPushed5000.pdf", f1, width = 20, height = 20, units = "cm")

# The region of highest population growth in an evolving wave ####

# Data manipulation functions
# Preparing data concerning the gen state of all invasions
invState <- function(data, gen){
  genFinal <- sapply(data, '[', gen) # Creates a list consisting of the final generation of each invasion
  genFinal
}

# Change in A trait over invasion space for all invasions
ADiffFinal <- function(data){
  alleeResFinal <- lapply(data, function(x) aggregate(x[,"A"]~x[,"patch"], FUN = mean)) # Calculates mean A per patch for each invasion
  alleeResFinal <- rbindlist(alleeResFinal, idcol = TRUE) # Begins conversion of data to data frame
  alleeResFinal <- as.data.frame(alleeResFinal) # Finishes conversion of data to data frame 
  colnames(alleeResFinal) <- c("invasion", "patch", "A") # Renames columns
  alleeResFinal
}

# Change in density over invasion space for all invasions
densFinal <- function(data){
  densityFinal <- lapply(data, function(x) table(x[,"patch"])) # Calculates density per patch for each invasion
  densityFinal <- lapply(densityFinal, function(x) as.data.frame(x)) # Begins conversion of data to data frame
  densityFinal <- rbindlist(densityFinal, idcol = TRUE) # Continues conversion of data to data frame
  densityFinal <- as.data.frame(densityFinal) # Finishes conversion of data to data frame
  colnames(densityFinal) <- c("invasion", "patch", "density") # Renames columns
  densityFinal <- transform(densityFinal, patch = as.numeric(patch)) # Transforms patch from a factor into a numeric
  densityFinal
}

# Loading data
load(file="data/pushedToPulledWithMixing.RData")
pp <- simResults

# Manipulating data
# Generating early, pushed wave snapshot
earlyInv <- invState(pp,20)
earlyInv.ADiff <- ADiffFinal(earlyInv)
earlyInv.dens <- densFinal(earlyInv)
earlyInv.densADiff <- merge(earlyInv.dens, earlyInv.ADiff)
earlyInv.densADiffmean<-aggregate(earlyInv.densADiff[,3:4], by=list(earlyInv.densADiff$patch), FUN=mean)
earlyInv.densADiffmean$growth <- alleeFunc(u=earlyInv.densADiffmean$density, t=earlyInv.densADiffmean$A, K=K, r=r)
colnames(earlyInv.densADiffmean) <- c("patch", "density", "A", "growth") # Renames columns

# Generating late, pulled wave snapshot
lateInv <- invState(pp,220)
lateInv.ADiff <- ADiffFinal(lateInv)
lateInv.dens <- densFinal(lateInv)
lateInv.densADiff <- merge(lateInv.dens, lateInv.ADiff)
lateInv.densADiffmean<-aggregate(lateInv.densADiff[,3:4], by=list(lateInv.densADiff$patch), FUN=mean)
lateInv.densADiffmean$growth <- alleeFunc(u=lateInv.densADiffmean$density, t=lateInv.densADiffmean$A, K=K, r=r)
colnames(lateInv.densADiffmean) <- c("patch", "density", "A", "growth") # Renames columns

# Plots

# Pushed wave population growth
legendTitle <- guide_legend(title=TeX("Mean expected reproductive output ($\\textit{E}(\\textit{W_{ix}})$)"))
p1 <- ggplot(earlyInv.densADiffmean, aes(patch, density, fill = growth, colour = growth)) + 
  geom_col() +
  theme_classic(base_size = 12) +
  scale_fill_gradient(breaks=c(0.90,1.00,1.10,1.20,1.30), limits = c(0.9,1.3)) +
  scale_colour_gradient(breaks=c(0.90,1.00,1.10,1.20,1.30), limits = c(0.9,1.3)) +
  #scale_fill_gradient(low="red3", high="royalblue2", breaks=c(0,-50,-100,-150,-200,-250), limits = c(-250, 0)) +
  #scale_colour_gradient(low="red3", high="royalblue2", breaks=c(0,-50,-100,-150,-200,-250), limits = c(-250, 0)) +
  labs(x = " ", y = "Mean density", tag = "A") +
  guides(fill = legendTitle, colour = legendTitle) +
  theme(legend.position = c(0.6, 0.6)) +
  theme(plot.tag = element_text(size=24)) + 
  expand_limits(x = 95)
p1

# Pulled wave population growth
p2 <- ggplot(lateInv.densADiffmean, aes(patch, density, fill = growth, colour = growth)) + 
  geom_col() +
  theme_classic(base_size = 12) +
  scale_fill_gradient(breaks=c(0.90,1.00,1.10,1.20,1.30), limits = c(0.9,1.3)) +
  scale_colour_gradient(breaks=c(0.90,1.00,1.10,1.20,1.30), limits = c(0.9,1.3)) +
  #scale_fill_gradient(low="red3", high="royalblue2", breaks=c(0,-50,-100,-150,-200,-250), limits = c(-250, 0)) +
  #scale_colour_gradient(low="red3", high="royalblue2", breaks=c(0,-50,-100,-150,-200,-250), limits = c(-250, 0)) +
  labs(x = "Patch", y = "Mean density", tag = "B") +
  guides(fill = FALSE, colour = FALSE) +
  theme(legend.text=element_text(size=12), plot.tag = element_text(size=24)) + 
  expand_limits(x = 95)
p2

# Multiplot

a1 <- arrangeGrob(p1)

a2 <- arrangeGrob(p2)

lay <- rbind(c(1),
             c(2))

f1 <- grid.arrange(a1, a2, layout_matrix = lay, heights=c(10,10))
f1
ggsave(filename = "figures/figPPGrowth.pdf", f1, width = 20, height = 20, units = "cm")

