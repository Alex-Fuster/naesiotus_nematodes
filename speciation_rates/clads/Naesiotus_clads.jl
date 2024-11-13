## Package loading

using PANDA

## Importing tree

my_tree = load_tree("naesiotus_rescaled.tre")

## Running ClaDS

output = infer_ClaDS(my_tree, n_trees = 100, print_state = 100, f = 0.8588)

## Saving to RData file

save_ClaDS_in_R(output, "naesiotus_clads_rates.RData")

