// -------------------------------------------------------------------------------------------
// DATA
/* The data block is for the declaration of variables that are read in as data
   Data variables are not transformed in any way
   The data block does not allow statements
   Each variable's value is validated against its declaration as it is read */

// The input data is a matrix 'epsilon' of dimensions n (rows=patient), TT (columns=time).data {
  int<lower=1> TT;      //n.time instants
  int<lower=1> n;       //n. patients
  matrix[n,TT] epsilon; //x_it-mean
  vector[1] time[TT];   
  real <lower=3> a0;
  real b0;
  real muL;
  real <lower=0> sigmaL;
}
// ------------------------------------------------------------------------------------------------
// TRANSFORMED DATA
/* The transformed data block is for declaring and defining variables that do not 
   need to be changed when running the program.
   Any variable that is defined wholly in terms of data or transformed data should 
   be declared and defined in the transformed data block. */

//I declare AND define the vector of weights of the multi GP 
transformed data{
  vector[n] w = rep_vector(1.0,n);
}

// The parameters accepted by the model: sigma2[t] for t in 1:1600 and L (characteristic distance)
parameters {
  real <lower=0> sigma2[TT];
  real <lower=0, upper=1> L;
}


model{ 
  
  matrix[TT,TT] R;
  matrix[TT,TT] Tri;
  Tri = rep_matrix(0,TT,TT);
  // second of all: all the local variables defined
  R = cov_exp_quad(time, sqrt(sigma2), L);
  
  for (i in 1:TT){
   Tri[i,i] = R[i,i];
   if(i>1 && i<=TT){
     Tri[i,i-1] = R[i,i-1];
     Tri[i-1,i] = Tri[i,i-1];
     }
  }
  // Prior:
 
	sigma2 ~ inv_gamma(a0,b0); // iid sample
	
	L ~ normal(muL, sigmaL);  
	
	// Likelihood:
	epsilon ~ multi_gp(Tri,w);
}

