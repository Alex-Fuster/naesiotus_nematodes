data {
  int<lower=0> N_nem;               // Number of observations for nematode load
  int<lower=0> N_bright;             // Number of observations for brightness
  
  array[N_nem] int<lower=0, upper=100>load;           // Nematode load counts
  vector[N_bright] brightness;        // Brightness scores
  
  int<lower=0> N_sp_nem;              // Number of species in the nematode data
  int<lower=0> N_sp_bright;           // Number of species in the brightness data
  
  array[N_nem] int<lower=1>species_i; // Species id for load dataset
  array[N_bright] int<lower=1>species_k; // Species id for brightness dataset
}


parameters {
  real intercept;
  real slope_bright;
  
  vector[N_sp_bright] mu_bright;      // Mean brightness for each species
  vector<lower=0>[N_sp_bright] sigma_bright; // Standard deviation of brightness for each species
}


model {
  
  vector[N_sp_nem] lamda;
  
 // Priors
  intercept ~ normal(0, 10);
  slope_bright ~ normal(0, 1);
  mu_bright ~ normal(0, 1);
  sigma_bright ~ exponential(1);
  
 // Model nematode load
  
  brightness[species_k] ~ normal(mu_bright[species_k], sigma_bright[species_k]);
  lamda[species_i] = intercept + slope_bright * mu_bright[species_i];
  load[species_i] ~ poisson(lamda[species_i]);
}

  
