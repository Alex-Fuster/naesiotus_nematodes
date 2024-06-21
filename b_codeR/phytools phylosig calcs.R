#nematode phytools stuff
library(phytools)

#set working directory
setwd("~/R/nematodes")

#get the tree and load data from the 'phylog_distances_naesiotus' script from alex
tree<- as.phylo(tree)
df_load

#make the table we want to use for the phylosig arguments (will need the tidyverse packages from alex's script)
df_meanload <- df_load %>%
  group_by(spp) %>%
  summarise(n = n(),
            mean_load = mean(nematode_count, na.rm = TRUE),
            se_load = sd(nematode_count, na.rm = TRUE)/sqrt(length(nematode_count)), 
            sd_load = sd(nematode_count, na.rm = TRUE)) %>%
  arrange(desc(mean_load))

colnames(df_meanload) = c("species", "n" ,"mean_load", "se_load", "sd_load")

df_meanload$mean_load <- round(df_meanload$mean_load, 1)
df_meanload$sd_load <- round(df_meanload$sd_load, 1)
df_meanload$se_load <- round(df_meanload$se_load, 1)

meanload<-setNames(df_meanload$mean_load, df_meanload$species)
s.error<-setNames(df_meanload$se_load, df_meanload$species)



#prune the tree to match the data (71 tree tips and 47 species with data)
plot(tree)
pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, df_meanload$species));

#should be ready to use phylosig
?phylosig
K_nemaload <- phylosig(tree, meanload, method = 'K', test = TRUE, se = s.error)
K_nemaload
plot.phylosig(K_nemaload)

#why can't I use the standard error with lambda?
lambda_nemaload <- phylosig(tree, meanload, method = 'lambda', test = TRUE) #, se = s.error) - not sure why it doesn't like this
lambda_nemaload
plot.phylosig(lambda_nemaload)