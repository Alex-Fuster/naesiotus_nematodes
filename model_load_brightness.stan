



data {
  int<lower=0> N_nem;            // Number of observations for nematode load
  int<lower=0> N_bright;         // Number of observations for habitat

  int<lower=0> load[N_nem];      // Nematode load counts;
  vector[N_bright] brightness;   // Brightness scores;

  int<lower=0> N_sp_nem;         // Number of species in the nematode data
  int<lower=0> N_sp_bright;      // Number of species in the habitat data

  // Species id for each dataset
  vector[N_sp_nem] species_i;
  vector[N_sp_bright] species_k;
}





parameters {
  
  // habitat is modeled for each species, so there is a slope, mu and sigma for each spp
  
  real intercept;
  real slope_hab[N_sp_bright]; // N_sp_hab or N_sp_nem?????
  real mu_hab[N_sp_bright]; // N_sp_hab or N_sp_nem?????
  real sigma_hab[N_sp_bright]; // N_sp_hab or N_sp_nem?????
  
}





model {
  
  
  // Priors
  
  intercept ~ normal(0,10);
  slope_hab ~ normal(0,1);
  mu_hab ~  normal(0,1);
  sigma_hab ~ exp(1);
  
  
  // Model nematode load
  
  hab[species_k] ~ normal(mu_hab[species_k], sigma_hab[species_k]);
  
  lamda[species_i] ~ intercept + slope_hab * mu_hab[species_i];
  
  load[species_i] ~ poisson(lamda[species_i]);
  
}

