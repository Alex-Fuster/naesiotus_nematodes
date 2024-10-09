
library(BAMMtools)
library(coda)


# 1. Check Convergence with coda

# Read in the tree, event data and MCMC data
phy <- read.tree("C:/Users/alexf/Downloads/bamm-2.5.0-Windows/bamm-2.5.0-Windows/bamm-2.5.0-Windows/cleannames_tree_041024.txt")
event_data <- read.csv("C:/Users/alexf/Downloads/bamm-2.5.0-Windows/bamm-2.5.0-Windows/bamm-2.5.0-Windows/event_data.txt", sep = ",", header = TRUE)
mcmc_data <- read.table("C:/Users/alexf/Downloads/bamm-2.5.0-Windows/bamm-2.5.0-Windows/bamm-2.5.0-Windows/mcmc_out.txt", header = TRUE, sep = ",")

# Convert to mcmc object
mcmc_obj <- as.mcmc(mcmc_data)

# Check convergence with ESS and trace plots
effectiveSize(mcmc_obj)
plot(mcmc_obj)


# Summarize Posterior Distribution

summary(mcmc_obj)
# Plot histogram of the number of rate shifts
hist(mcmc_data$N_shifts, breaks = 30, main = "Posterior Distribution of Rate Shifts", xlab = "Number of Rate Shifts")


#####################################################


# Process and plot events
bamm_data <- getEventData(phy, eventdata = event_data, burnin = 0.25)
summary(bamm_data)


# Plot mean rate shifts for specific clades
plot.bammdata(bamm_data, lwd = 2, logcolor = TRUE)

# # Generate a rate-through-time plot
# # Extract rate-through-time data without additional arguments
# rtt <- getRateThroughTimeMatrix(bamm_data)
# 
# # Calculate mean net diversification rates over time
# mean_speciation <- apply(rtt$lambda, 2, mean)
# mean_extinction <- apply(rtt$mu, 2, mean)
# net_diversification <- mean_speciation - mean_extinction
# 
# # Plot net diversification rate over time with log scale on the y-axis
# plot(rtt$times, net_diversification, type = "l", log = "y",
#      xlab = "Time", ylab = "Net Diversification Rate",
#      main = "Net Diversification Rate Through Time (Log-Scale)")


# Plot the tree with labeled nodes
plot(phy, cex = 0.6)
nodelabels((Ntip(phy) + 1):(Ntip(phy) + Nnode(phy)), adj = c(1.2, -0.5), frame = "none", cex = 0.6)


# Compute rate for a specific clade

# List of nodes you want to process
nodes_of_interest <- c(118)  # You can add more nodes here

# Initialize an empty list to store results
results <- list()

# Loop over each node to compute metrics and store results
for (node in nodes_of_interest) {
  # Get clade rates for the node
  clade_rates <- getCladeRates(bamm_data, node = node, nodetype = "include", verbose = TRUE)
  
  # Extract lambda and mu
  lambda <- clade_rates$lambda
  mu <- clade_rates$mu
  
  # Calculate extinction fraction (?? = ?? / ??) and net diversification rate (r = ?? - ??)
  extinction_fraction <- mu / lambda
  net_diversification <- lambda - mu
  
  # Store results in a list
  results[[as.character(node)]] <- data.frame(
    Node = node,
    Lambda = lambda,
    Mu = mu,
    ExtinctionFraction = extinction_fraction,
    NetDiversification = net_diversification
  )
}

# Combine all results into a single data frame
summary_table <- do.call(rbind, results)
