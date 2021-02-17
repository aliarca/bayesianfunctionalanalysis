// -------------------------------------------------------------------------------------------
// DATA
/* The data block is for the declaration of variables that are read in as data
   Data variables are not transformed in any way
   The data block does not allow statements
   Each variable's value is validated against its declaration as it is read */

// The input data is a matrix 'epsilon' of dimensions n (rows=patient), TT (columns=time).
data {
  int<lower=0> TT;      //n.time instants
  int<lower=0> n;       //n. patients
  vector[TT] epsilon[n];//x_it-mean
  vector[1] time[TT];
  real <lower=3> a0;
  real b0;
  real muL;
  real <lower=0> sigmaL;
  real a1;
}

// ------------------------------------------------------------------------------------------------
// TRANSFORMED DATA
/* The transformed data block is for declaring and defining variables that do not 
   need to be changed when running the program.
   Any variable that is defined wholly in terms of data or transformed data should 
   be declared and defined in the transformed data block. */

//I declare AND define epsilon mean vector (of zeros)   
transformed data{
  vector[TT] mu = rep_vector(0, TT); 
}


// ------------------------------------------------------------------------------------------------
// PARAMETERS
/* The variables declared in the parameters program block correspond directly to the variables 
   being sampled by Stan's samplers. Variables declared as parameters cannot be directly
   assigned values. So there is no block of statements in the parameters program block. 
   
   DISCLAIMER: There is a substantial amount of computation involved for parameter variables in a
   Stan program at each leapfrog step within the HMC or NUTS samplers, and a bit more computation 
   along with writes involved for saving the parameter values corresponding to a sample. */

//The parameters accepted by the model: sigma2[t] for t in 1:1600 and L 
parameters {
  real<lower=0> sigma2; // constraint on the instantaneous variance
  real L;                   // characteristic distance
  real<lower=0> nu;       // dof of the inverse wishart for Omega
}
// ------------------------------------------------------------------------------------------------
// TRANSFORMED PARAMETERS
/* The transformed parameters program block consists of optional variable declarations followed 
   by statements. After the statements are executed, the constraints on the transformed
   parameters are validated. Any variable declared as a transformed parameter is part of 
   the output produced for draws. */

// ------------------------------------------------------------------------------------------------
// MODEL
/* The model program block consists of optional variable declarations followed by statements. 
   The variables in the model block are local variables and are not written as part of the output. 
   
   Local variables may not be defined with constraints because there is no well-defined way to 
   have them be both flexible and easy to validate
   
   The statements in the model block typically define the model. This is the block in which
   probability (sampling notation) statements are allowed.*/
   
// I define two local variables: R covariance matrix, tridiagonal, and its Cholesky decomposition   
model{
  // first of all: all the local variables declared 
  matrix[TT,TT] R;
  matrix[TT,TT] Omega;
  // second of all: all the local variables defined
  R = cov_exp_quad(time, sqrt(sigma2), L);
  // for (i in 1:TT){
  //   R[i,i] = 1/sigma2[i];
  //   if(i>1 && i<=TT){
  //     R[i,i-1] = 1/sigma2[i]*exp(-0.5/L^2);
  //     R[i-1,i] = R[i,i-1];
  //   }
  // }
//   for (i in 1:TT){
// 	  sigma2[i] ~ inv_gamma(a0,b0); // iid sample
// 	}
 // CR = cholesky_decompose(R);
  sigma2 ~ gamma(a0,b0);

  // for(i in 1:TT){
  //   sigma2[i] ~ gamma(a0,b0);
  // }

	L ~ normal(muL, sigmaL);

	// nu ~ chi_square(a1);

 
  Omega ~ wishart(nu+TT+1, R);
	
	// Likelihood:
	for(i in 1:n){
		epsilon[i,] ~ multi_normal_prec(mu, R); 
	}
}


// PROVA MULTI_NORMAL_PRECISION (inverso di inv_whishart è whicshart) / CAR_PRIOR PER TRIDIAGONALI
