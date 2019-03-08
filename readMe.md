# Read me

## Description of model
This is an individual-based model of a spreading population with discrete generations in which a trait that governs resistance to the Allee effect may undergo selection on the invasion front.

An output object called 'simResults' is produced at the end of the simulation. 'simResults' is a list that keeps track of all generations of all invasions simulated.

Each generation of the invasion is stored as a matrix, in which each row is an individual. Attributes of each individual are tracked in the columns of the matrix, such as their location ('patch') and the local density of their location ('n').

The function and purpose of each script is described in its respective file.

## How to run the model
After opening '2018 Evolution transforms pushed waves into pulled waves.Rproj', in the scripts folder, run 'runModel.R' after setting your parameters in 'parameters.R'.

## How to explore simulation results
In the scripts folder, run 'analysisFigures.R' after simulating your invasions.