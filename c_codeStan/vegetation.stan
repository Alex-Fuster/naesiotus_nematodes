data {
  int<lower=0> N_nem;               // Number of observations for nematode load
  int<lower=0> N_veg;               // Number of observations for vegetation zone
  array[N_nem] int<lower=0, upper=100> load;  // Nematode load counts
  array[N_veg] int<lower=0, upper=1> vegetation; // Vegetation zone: 0 for humid, 1 for arid
  int<lower=0> N_sp_nem;             // Number of species in the nematode data
  array[N_nem] int<lower=1, upper=N_sp_nem> species_i;   // New species id for load dataset
  array[N_veg] int<lower=1, upper=N_sp_nem> species_v; // New species id for vegetation dataset
}

parameters {
  real intercept;
  real<upper=0> slope_vegetation;
  vector<lower=0, upper=1>[N_sp_nem] vegetation_prob;  // Vegetation probability for each species
}

model {
  vector[N_sp_nem] lambda;
  
  // Priors
  intercept ~ normal(0, 1);
  slope_vegetation ~ normal(0, 1);
  vegetation_prob ~ beta(2, 2);  // Prior for vegetation probabilities
  
  // Model vegetation zone
  vegetation ~ bernoulli(vegetation_prob[species_v]);
  
  // Model nematode load
  lambda = intercept + slope_vegetation * vegetation_prob;
  load ~ poisson_log(lambda[species_i]);
}

generated quantities {
  array[N_veg] real predicted_vegetation;  // Predicted vegetation probabilities
  array[N_nem] int predicted_load;  // Predicted nematode load
  
  for (n in 1:N_nem) {
    predicted_load[n] = poisson_log_rng(intercept + slope_vegetation * vegetation_prob[species_i[n]]);
  }
  
  for (v in 1:N_veg) {
    predicted_vegetation[v] = bernoulli_rng(vegetation_prob[species_v[v]]);
  }
}
