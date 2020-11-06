//RSSp.stan
//10-20-2017
data {
int<lower=0> P; // number of SNPs
real u_hat[P]; // data
}
parameters {
  real<lower=0> sigma_u;
  real u[P];
}
model {
    u ~ normal(0,sigma_u);
    u_hat ~ normal(u,1);
}
