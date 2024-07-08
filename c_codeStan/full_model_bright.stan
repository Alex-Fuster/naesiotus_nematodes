
data{
  int N;
  int N_spp;
  array[N] int<lower=1,upper=N_spp> spp_id;
  int N_islands;
  array[N] int<lower=1,upper=N_islands> island_id;
  vector[N] brightness;
}
parameters{
  vector[N_spp] spp_effects;
  vector[N_islands] island_effects;

  real mu;
  real<lower=0> sigma_spp;
  real<lower=0> sigma_islands;

}

model {
  brightness ~ poisson_log(mu + spp_effects[spp_id] + island_effects[island_id]);
  spp_effects ~ normal(0, sigma_spp);
  island_effects ~ normal(0, sigma_islands);
  mu ~ normal(5, 2);
  sigma_spp ~ exponential(1);
  sigma_islands ~ exponential(1);

}



