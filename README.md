# Reproducing "Food web structure alters ecological communities top-heaviness with little effect on the biodiversity-functioning relationship"

This project contains the scripts needed to reproduce the experiments and analyses for the following paper: "Food web structure alters ecological communities top-heaviness with little effect on the biodiversity-functioning relationship" by Eva Delmas, Daniel B. Stouffer and Timothée Poisot.

## Requirement

We used the `Julia` programming language (`v. 1.3`), so you will first need to install `Julia` (see [this page for platform specific instructions](https://julialang.org/downloads/platform/)). To install all the packages needed, follow these steps:
1. Clone or download this project
2. `cd` into the project
3. open `julia`
4. use this script:

```
import Pkg #import the package manager
Pkg.activate(pwd()) #tell the package manager where to find the Manifest and Project files that store the environment info
Pkg.instantiate() #install the packages
```

## Content

### Scripts

- `ADBM_model.jl`: functions used to simulate food webs using the Allometric Diet Breadth model and body mass data from the Benguela pelagic ecosystem (stored in the data folder)
- `utils.jl`: various functions used in the other scripts
- `foodwebs_generation.jl`: simulates and stores the food webs
- `insilico_experiment.jl`: perform and stores simulations
- `outputs_extraction.jl`: calculate various measures of structure and functioning from the raw simulation outputs
- `analysis_and_figures.jl`: uses these measures to make the figures for the paper
- `biomassVSintake_baseline.jl`: estimates the biomass to intake relationship for a food web with no consumers (data used in the `analysis_and_figures.jl` script)

### Data

This contains the body mass data and the simulations outputs.
