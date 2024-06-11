data {
  int<lower=0> N_nem;               // Number of observations for nematode load
  int<lower=0> N_bright;            // Number of observations for brightness
  array[N_nem] int<lower=0, upper=100> load;  // Nematode load counts
  vector[N_bright] brightness;       // Brightness scores
  int<lower=0> N_sp_nem;             // Number of species in the nematode data
  int<lower=0> N_sp_bright;          // Number of species in the brightness data
  array[N_nem] int<lower=1> species_i;   // New species id for load dataset
  array[N_bright] int<lower=1> species_k; // New species id for brightness dataset
  array[N_sp_nem] int<lower=1, upper=N_sp_bright> species_map; // Mapping from nematode species to brightness species
}

parameters {
  real intercept;
  real slope_bright;
  vector[N_sp_bright] mu_bright;       // Mean brightness for each species
  vector<lower=0>[N_sp_bright] sigma_bright; // Standard deviation of brightness for each species
}

model {
  vector[N_sp_nem] lambda;
  
  // Priors
  intercept ~ normal(0, 10);
  slope_bright ~ normal(0, 1);
  mu_bright ~ normal(0, 1);
  sigma_bright ~ exponential(1);
  
  // Model brightness
  brightness ~ normal(mu_bright[species_k], sigma_bright[species_k]);
  
  // Model nematode load
  lambda = intercept + slope_bright * mu_bright[species_map[species_i]];
  load ~ poisson(lambda[species_i]);
}
