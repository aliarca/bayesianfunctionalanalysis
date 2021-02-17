# Bayesian functional analysis
## Study of EEGs from brain responses in coma patients
The aim of this project is to study the electrical brain waves of patients in a state of coma. In particular, an external electrical stimulus is applied to the patient and then the electrical brain response is recorded. These responses are here studied and then compared to the actual recovery of the patients, in order to verify whether there are significant differences in the curves of patients who had a good recovery (identified by the GOSE index equal to 2) and in the curves of patients with a bad recovery (GOSE equal to 1). 
To do so, we followed a Baeysian Functional Anova approach.

## required packages
If you don't want to miss anything, here there is a list of all the packages we have used
```r
needed_packages  <- c("fda.usc","fda","fields","LaplacesDemon","matlab","mvtnorm", "Rmpfr", "LearnBayes")
new_packages  <- needed_packages[!(needed_packages %in%installed.packages ()[, "Package"])]
if (length(new_packages))install.packages(new_packages)
lapply(needed_packages , require , character.only = TRUE)
```
## quick guide
_To directly run the code, just go to the core.R file, it contains all the code that will patientely guide you to the final result. But if you are curious about the preprocessing, and you want to understand better our approach, you can have a detailed perspective by running Functional.R and fPCA.R before the main one. In those files we work on several functional representations of the EEG curves, and then we reduced the optimal representation obtained before thanks to functional PCA. In the end, for the bravest, the file stan.R contains the code for an upgraded, more complex and computationally demanding model._

## structure: folders and files
- **data** contains:
  - **_Codice_analisi_funzionale.RData_** the base of our work: the data!
- **code** contains:
  - **_core.R_** - the heart of our code, with the calls to all the external functions
  - **_Functional.R_** - where we worked on functional representation of the data
  - **_fPCA.R_** - some preprocessing work
  - **_graphics.R_**
  - **_stan.R_** an additional file with some tentatives of implementing a more complicated model
- **samplers** contains:
  - **_niw.RData_**
  - **_nor.RData_**
  - **_gib.RData_** 
  - **_update.RData_** additional function called by gibbs
- **functions** contains some functions we built ad hoc for the problem:
  - **_calcolo.loglike.RData_** to computhe the likelihood of the data
  - **_bieffe.RData_** to compute the Bayes Factor
  - **_bieffe.power.RData_** to compute the BF with power likelihood
- **stan** contains:
  - **_precision.stan_**
  - **_niw.stan_**
  - **_cholesky.stan_**
  - **_gp.stan_**
  - **_gp_cholesky.stan_**

## mathematical model
The core of the model is the Bayes Factor computation. This because we are dealing with a hypothesis testing, with in the null hypothesis the fact that all the patients waves are equally distributes, and under the alternative they are separated according to the recovery

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/hyp.png" width="40%" height="40%">

To do so, we need to compute the joint likelihood, so here it is the model for the curves, assumed independent among patiens

First of all, we are dealing with dense data, that can be treated as functional, so instead of using a multivariate Gaussian model, we prefer a Gaussian Process.

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model1.png" width="20%" height="20%">

Secondly we divide the curve in two components, a random mean term, and a random (zero mean) error

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model2.png" width="40%" height="40%">

The first step is the sampling for the mean, according to the Bayesian approach. We represent it via suitable basis, to work with an equivalent but lighter model

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model3.png" width="15%" height="15%">

So that we study a small dimensional Normal Inverse Wishart distribution, from which we can direcly sample

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model4.png" width="40%" height="40%">

Then the error term, which we assume time independent

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model5.png" width="30%" height="30%">

is studied and sampled through its hyperparameters (in this case we simulate from the two posteriors via Gibbs Sampling)

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model6.png" width="50%" height="50%">


Finally we have to compute (via MC approximation) the likelihood

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model7.png" width="40%" height="40%">

which we translate in logarithmic form, and also thanks to independece, it assumes the followin form

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model8.png" width="40%" height="40%">

Ultimately, we computed the Bayes Factor, which will give us the final results

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/model9.png" width="30%" height="30%">

## results: an overview
In this overview, we want to take you along the discoveries we made, to finally lead you to the result.
First of all, a glimpe over a couple of the signals we studied
* signal 5 (Central SX Medium Latency) for all patients, GOSE1 patiens and GOSE2 patients
<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/all_5.png" width="50%" height="50%">
<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/gose1_5.png" width="50%" height="50%">
<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/gose2_5.png" width="50%" height="50%">

* signal 4 (Fontal DX Short Latency) for all patients, GOSE1 patiens and GOSE2 patients
<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/all_4.png" width="50%" height="50%">
<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/gose1_4.png" width="50%" height="50%">
<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/gose2_4.png" width="50%" height="50%">
As you can see, some signals are more tidy than others, which insteas vary much more.

As declared in the **mathematical model**, the BF is the tool that gives us the finale result, which is the statistical evidence in favour of the separated model vs the complete one. As you can see from this table, which actually reports a more robust version of the BF, obtained via a coarsened likelihood, **all the signals except the fourth** are evidently better explained by the "recovery model" (2*logBF >> 0). As regards the fourth, we have no evidence in favour of one model nor the other, as one could qualitatevely guess from the plots above.

<img src="https://github.com/aliarca/bayesianfunctionalanalysis/blob/main/images/result.png" width="50%" height="50%">



## Acknowledgments

* Bayesian Statistics course Professor Alessandra Guglielmi
* Project tutor Dr. Riccardo Corradin
