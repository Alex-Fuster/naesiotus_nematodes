data {
  int<lower=0> N_nem;               // Number of observations for nematode load
  int<lower=0> N_habitat;           // Number of observations for habitat
  array[N_nem] int<lower=0, upper=100> load;  // Nematode load counts
  array[N_habitat] int<lower=0, upper=1> habitat; // Habitat: 0 for terrestrial, 1 for arboreal
  int<lower=0> N_sp_nem;             // Number of species in the nematode data
  array[N_nem] int<lower=1, upper=N_sp_nem> species_i;   // New species id for load dataset
  array[N_habitat] int<lower=1, upper=N_sp_nem> species_h; // New species id for habitat dataset
}

parameters {
  real intercept;
  real<upper=0> slope_habitat;
  vector<lower=0, upper=1>[N_sp_nem] habitat_prob;  // Habitat probability for each species
}

model {
  vector[N_sp_nem] lambda;
  
  // Priors
  intercept ~ normal(0, 1);
  slope_habitat ~ normal(0, 1);
  habitat_prob ~ beta(2, 2);  // Prior for habitat probabilities
  
  // Model habitat
  habitat ~ bernoulli(habitat_prob[species_h]);
  
  // Model nematode load
  lambda = intercept + slope_habitat * habitat_prob;
  load ~ poisson_log(lambda[species_i]);
}

generated quantities {
  array[N_habitat] real predicted_habitat;  // Predicted habitat probabilities
  array[N_nem] int predicted_load;  // Predicted nematode load
  
  for (n in 1:N_nem) {
    predicted_load[n] = poisson_log_rng(intercept + slope_habitat * habitat_prob[species_i[n]]);
  }
  
  for (h in 1:N_habitat) {
    predicted_habitat[h] = bernoulli_rng(habitat_prob[species_h[h]]);
  }
}
