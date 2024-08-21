data {
  int N_spp;
  array[N_spp] int<lower=1, upper=N_spp> sp_id;
  
  // habitat data
  int N_habitat;
  array[N_habitat] int habitat_arboreal;   // Total counts of arid habitat for each species
  array[N_habitat] int total_hab;  //total counts of habitat (arboreal or terrestrial)
  array[N_habitat] int sp_index_hab; 
}
parameters {
  real mu_arboreal;
  real<lower=0> sd_arboreal;
}
transformed parameters {
  vector[N_habitat] arboreal_prob = inv_logit(arboreal_prob_logit);
}
model {
  vector[N_habitat] arboreal_prob_logit;
  // Hierarchical priors for the probabilities
  arboreal_prob_logit ~ normal(mu_arboreal, sd_arboreal);
  
  habitat_arboreal ~ binomial_logit(total_hab, arboreal_prob_logit[sp_index_hab])
}
  