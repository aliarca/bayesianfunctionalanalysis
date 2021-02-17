// The input data is a matrix 'epsilon' of dimensions n (rows=patient), TT (columns=time).
data {
  int<lower=1> TT;      //n.time instants
  int<lower=1> n;       //n. patients
  matrix[n,TT] epsilon; //x_it-mean
  vector[1] time[TT];      
  real <lower=3> a0;
  real b0;
  real muL;
  real <lower=0> sigma2L;
}
transformed data{
  vector[n] w = rep_vector(1.0,n);
}

// The parameters accepted by the model: sigma2[t] for t in 1:1600 and L (characteristic distance)
parameters {
  real <lower=0> sigma2[TT];
  real <lower=0, upper=1> L;
}


model{  
  matrix[TT,TT] R = rep_matrix(0,TT,TT);
  matrix[TT,TT] L_R;
  for (i in 1:TT){
    R[i,i] = sigma2[i];
    if(i>1 && i<=TT){
      R[i,i-1] = sigma2[i]*exp(-0.5/L^2);
      R[i-1,i] = R[i,i-1];
    }
  }
  L_R = cholesky_decompose(R);  //I put the matrix here because it uses L (unknown)
  
  // Prior:
  for (i in 1:TT){
	  sigma2[i] ~ inv_gamma(a0,b0); // iid sample
	}
	L ~ normal(muL, sigma2L);  
	
	// Likelihood:
	epsilon ~ multi_gp_cholesky(L_R,w);
}

