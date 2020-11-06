//RSSp.stan
//10-20-2017
data {
int<lower=0> P; // number of SNPs
real q[P]; //transformed data
}
parameters {
  real<lower=0> sigma_u;
}
transformed parameters{
  real<lower=0> t_sigma_u;
  t_sigma_u = sqrt(square(sigma_u) +1);
}
model {
  q ~ normal(0,t_sigma_u);
}
