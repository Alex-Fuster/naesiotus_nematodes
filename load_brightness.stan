data {
  int<lower=0> N_nem;                 // Number of observations for nematode load
  int<lower=0> N_bright;              // Number of observations for brightness
  
  int load[N_nem];           // Nematode load counts
  vector[N_bright] brightness;        // Brightness scores
  
  int<lower=0> N_sp_nem;              // Number of species in the nematode data
  int<lower=0> N_sp_bright;           // Number of species in the brightness data
  
  int<lower=1, upper=N_sp_nem> species_i[N_nem]; // Species id for load dataset
  int<lower=1, upper=N_sp_bright> species_k[N_bright]; // Species id for brightness dataset
}


parameters {
  
  real intercept;
  real slope_bright;
  
  vector[N_sp_bright] mu_bright;      // Mean brightness for each species
  vector<lower=0>[N_sp_bright] sigma_bright; // Standard deviation of brightness for each species
}


model {
  
  // Priors
  intercept ~ normal(0, 10);
  slope_bright ~ normal(0, 1);
  mu_bright ~ normal(0, 1);
  sigma_bright ~ exp(1);
  
  // Model nematode load
  
  brightness[species_k] ~ normal(mu_bright[species_k], sigma_bright[species_k]);
  lamda[species_i] ~ intercept + slope_bright * mu_bright[species_i];
  load[species_i] ~ poisson(lamda[species_i]);
}

