data {
  int <lower=0> n;
  real y[n];
  real W;
  real m0;
  real C0;
}

parameters {
  real mu; //Level
  // real rho;
  real theta_1;
  real epsilon[n];
}

transformed parameters{
  real<lower=0> sigma[n]; //Level
  real theta[n];
  
  theta[1]=theta_1;
  sigma[1]=exp(-(theta[1])/2);
  
  for (t in 2:n) {
    // theta[t]=tanh(rho)*theta[t-1]+epsilon[t];
    // theta[t]=rho*theta[t-1]+epsilon[t];
    theta[t]=m0*theta[t-1]+epsilon[t];
    sigma[t]=exp(-(theta[t])/2);
  }
  
  
}

model{
  mu ~ normal(0,1);
  // rho    ~ normal(m0,sqrt(C0));
  theta_1   ~ normal(0,1);
  epsilon   ~ normal(0,sqrt(W));
  y ~ normal(mu,sigma);
}