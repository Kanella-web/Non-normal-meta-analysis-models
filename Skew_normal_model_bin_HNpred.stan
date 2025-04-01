
data {
  int<lower=1> ns;  // Number of observations
  array[ns] int r1;
  array[ns] int r2;
  array[ns] int n1;
  array[ns] int n2;

}

parameters {
  array[ns] real m;
  
  real xi;
  real<lower=0> omega;
  real alpha;
  
  array[ns] real delta12;
}
transformed parameters {
  real tau_sqr;
  real mu;
  real skew;
  real beta;
  real delt;
  
  beta = sqrt(2/pi()) ; 
  
  delt = ( alpha / (sqrt(1 + (alpha^2))) ) ; 
  
  mu = (xi + (beta * omega * delt) ) ; 
  
  tau_sqr = ( (omega^2) * ( 1 - ((beta^2) * (delt^2)) ) ) ;
  
  skew = ( (4 - pi()) / 2 ) * ( (( beta*delt)^3) / (( 1 - ((beta*delt)^2) )^(1.5)) ) ; 
  
}

model {
  for (i in 1:ns) {
    m[i] ~ normal(0, 100);
    
    r1[i] ~ binomial_logit(n1[i], m[i]);
    r2[i] ~ binomial_logit(n2[i], m[i] + delta12[i]);
   
    delta12[i] ~ skew_normal(xi, omega, alpha);
  }
  
  xi ~ normal(0, 100);
  omega ~ normal(0, 1);
  alpha ~ normal(0, 5);
}

 generated quantities{
   real pred1;
   pred1 = skew_normal_rng(xi, omega, alpha);

 }

