# bayesianproject
Code and stuff from our project for Bayesian functional ANOVA


## required packages
If you want to be sure of not missing anything, here there is a list of all the packages we have used
```r
needed_packages  <- c("fda.usc","fda","fields","LaplacesDemon","matlab","mvtnorm", "Rmpfr", "LearnBayes")
new_packages  <- needed_packages[!(needed_packages %in%installed.packages ()[, "Package"])] if (length(new_packages))
install.packages(new_packages)
lapply(needed_packages , require , character.only = TRUE)}
```

## structure: folders and files
We have tried to organise the code in a thoughtful way
- **data** contains:
  - **_dati.csv_** a table with useful health and rehab indicators
  - **_Codice_analisi_funzionale.RData_** the base of our work: the data!
- **code** contains:
  - **_core.R_** - the heart of our code, with the calls to all the external functions
  
  - **_Codice_analisi_funzionale.R_** we build up our dataset from given data (its output is stored in **data**)
  - **_Functional.R_** - where we worked on functional representation of the data
  - **_fPCA.R_** - some preprocessing work
  - **_graphics.R_** 
- **samplers** contains:
  - **_niw.RData_**
  - **_nor.RData_**
  - **_gib.RData_** 
- **functions** contains some functions we built ad hoc for the problem:
  - **_calcolo.loglike.RData_** to computhe the likelihood of the data
  - **_bieffe.RData_** to compute the Bayes Factor
  - **_bieffe.power.RData_** to compute the BF with power likelihood

## mathematical model

## structure: code.... SI METTE OPPURE BASTANO I FOLDER?????? MAGARI UNO PSEUDOCODE?

## results: an overview
Plots and spoiler BF
