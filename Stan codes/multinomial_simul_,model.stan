data {
  int <lower=0> t;
  int <lower=0> n;
  int <lower=0> k;
  int y[t,k+1];
  matrix[k,n] F[t];
  matrix[n,n] G[t-1];
  vector[n] a1;
  cov_matrix[n] R1;
  cov_matrix[n] W;
}
parameters {
  vector[n] theta_1;
  vector[n] omega[t-1];
}
transformed parameters{
  vector[n] theta[t];
  vector[k+1] lambda[t];
  vector<lower=0>[k+1] eta[t];
  theta[1]=theta_1;
  lambda[1,1:k]=F[1]*theta[1];
  lambda[1,k+1]=0;
  eta[1]=exp(lambda[1])/sum(exp(lambda[1]));
  for (i in 2:t) {
    theta[i]=G[i-1]*theta[i-1]+omega[i-1];
    lambda[i,1:k]=F[i]*theta[i];
    lambda[i,k+1]=0;
    eta[i]=exp(lambda[i])/sum(exp(lambda[i]));
  }
}
model{
  theta_1 ~ multi_normal(a1,R1);
  for(i in 1:(t-1)){
    omega[i]   ~ multi_normal(rep_vector(0,n),W);
  }
  for(i in 1:t){
    y[i] ~ multinomial(eta[i]);
  }
}
