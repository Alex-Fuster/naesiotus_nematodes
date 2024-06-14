data {
  int<lower=0> N_nem;               // Number of observations for nematode load
  int<lower=0> N_bright;            // Number of observations for brightness
  array[N_nem] int<lower=0, upper=100> load;  // Nematode load counts
  vector[N_bright] brightness;       // Brightness scores
  int<lower=0> N_sp_nem;             // Number of species in the nematode data
  // int<lower=0> N_sp_bright;          // Number of species in the brightness data
  array[N_nem] int<lower=1, upper=N_sp_nem> species_i;   // New species id for load dataset
  array[N_bright] int<lower=1, upper=N_sp_nem> species_k; // New species id for brightness dataset
}

parameters {
  real intercept;
  real slope_bright;
  vector[N_sp_nem] mu_bright;       // Mean brightness for each species
  real<lower=0> sigma_bright; // Standard deviation of brightness for each species
}

model {
  vector[N_sp_nem] lambda;
  
  // Priors
  intercept ~ normal(0, 1);
  slope_bright ~ normal(0, 1);
  mu_bright ~ normal(0, 1);
  sigma_bright ~ exponential(1);
  
  // Model brightness
  brightness ~ normal(mu_bright[species_k], sigma_bright);
  
  // Model nematode load
  lambda = intercept + slope_bright * mu_bright;
  load ~ poisson_log(lambda[species_i]);
}

generated quantities {
  array[N_nem] int predicted_load;  // Predicted nematode load
  
  for (n in 1:N_nem) {
    predicted_load[n] = poisson_log_rng(intercept + slope_bright * mu_bright[species_i[n]]);
  }
}
