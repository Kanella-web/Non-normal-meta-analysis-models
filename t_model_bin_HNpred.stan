 
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
  array[ns] real delta12;
  real<lower=2.5, upper=1000> nu;
}
transformed parameters {
  real tau_sqr;
  real mu;
  //real nu;

  //tau_sqr = (nu / (nu-1)) * omega ;
  //nu = (ns - 1) ; 
  tau_sqr = ((nu / (nu-2)) * (omega^2)) ;
  mu = xi ;

}
model{
  for (i in 1:ns) {
    m[i] ~ normal(0, 100);
    
    r1[i] ~ binomial_logit(n1[i], m[i]);
    r2[i] ~ binomial_logit(n2[i], m[i] + delta12[i]);
   
    delta12[i] ~ student_t(nu, xi, omega);


  }
  
  xi ~ normal(0, 100);
  omega ~ normal(0, 1);
  nu ~ exponential(0.10);
}

 generated quantities{
     real pred;
    pred = student_t_rng(nu, xi, omega);
  }
 
