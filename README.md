# Bayesian functional analysis
## Study of EEGs from brain responses in coma patients
The aim of this project is to study the electrical brain waves of patients in a state of coma. In particular, an external electrical stimulus is applied to the patient and then the electrical brain response is recorded. These responses are here studied and then compared to the actual recovery of the patients, in order to verify whether there are significant differences in the curves of patients who had a good recovery and in the curves of patients with a bad recovery. 
To do so, we followed a Baeysian Functional Anova approach.

## required packages
If you don't want to miss anything, here there is a list of all the packages we have used
```r
needed_packages  <- c("fda.usc","fda","fields","LaplacesDemon","matlab","mvtnorm", "Rmpfr", "LearnBayes")
new_packages  <- needed_packages[!(needed_packages %in%installed.packages ()[, "Package"])]
if (length(new_packages))install.packages(new_packages)
lapply(needed_packages , require , character.only = TRUE)
```

## structure: folders and files
We have tried to organise the code in a thoughtful way
- **data** contains:
  - **_Codice_analisi_funzionale.RData_** the base of our work: the data!
- **code** contains:
  - **_core.R_** - the heart of our code, with the calls to all the external functions
  - **_Functional.R_** - where we worked on functional representation of the data
  - **_fPCA.R_** - some preprocessing work
  - **_graphics.R_**

- **samplers** contains:
  - **_niw.RData_**
  - **_nor.RData_**
  - **_gib.RData_** 
  - **_update.RData_** additional function called by gibbs
- **functions** contains some functions we built ad hoc for the problem:
  - **_calcolo.loglike.RData_** to computhe the likelihood of the data
  - **_bieffe.RData_** to compute the Bayes Factor
  - **_bieffe.power.RData_** to compute the BF with power likelihood

## mathematical model
Having in mind our goal, namely to test if the curves of patients with a good recovery are statistically different from the others, we build the following hypothesis test.
\H_0
We assume our signals to be i.i.d samples from a 
## results: an overview
In this overview, we want to take you along the discoveries we made, to finally lead you to the result.
First of all, a glimpe over a couple of the signals we studied
* signal 5 (Central SX Medium Latency) for all patients
<img src="https://https://github.com/aliarca/bayesianfunctionalanalysis/edit/main/images/all_5.png" width="50%" height="50%">
* signal 4 (Fontal DX Short Latency) for all patients
<img src="https://https://github.com/aliarca/bayesianfunctionalanalysis/edit/main/images/all_4.png" width="50%" height="50%">
As you can see, some signals are more tidy than others, which insteas vary much more.



## Acknowledgments

* Bayesian Statistics course Professor Alessandra Guglielmi
* Project supervisors Dr.s Riccardo Corradin
* Project tutor Dr. Mario Beraha
