#nematode phytools stuff
library(phytools)
library(tidyverse)

tree <- read.tree("a_data/phylogeny_2022/22Nov_Naesiotus.nwk")
df_rads <- read_excel("a_data/GPS RAD snails.xlsx")

#get the tree and load data from the 'phylog_distances_naesiotus' script from alex
tree<- as.phylo(tree)
#load
df_load <- read_excel("a_data/04april24_datalab.xlsx", sheet = 3)
df_load$nematode_count[which(df_load$nematode_count == ">100")] <- "100"
df_load$nematode_count <- as.numeric(df_load$nematode_count)

#make the table we want to use for the phylosig arguments (will need the tidyverse packages from alex's script)
df_meanload <- df_load %>%
  group_by(spp) %>%
  summarise(n = n(),
            mean_load = mean(nematode_count, na.rm = TRUE),
            se_load = sd(nematode_count, na.rm = TRUE)/sqrt(length(nematode_count)), 
            sd_load = sd(nematode_count, na.rm = TRUE)) %>%
  arrange(desc(mean_load)) %>% 
  rename(
    species = spp
  )



df_meanload$mean_load <- round(df_meanload$mean_load, 1)
df_meanload$sd_load <- round(df_meanload$sd_load, 1)
df_meanload$se_load <- round(df_meanload$se_load, 1)

meanload<-setNames(df_meanload$mean_load, df_meanload$species)
s.error<-setNames(df_meanload$se_load, df_meanload$species)



#prune the tree to match the data (71 tree tips and 47 species with data)
plot(tree)

# Names checks

#df_rads[which(df_rads$Taxon == "albermarlensis"),] <- "albemarlensis"
df_rads[which(df_rads$Taxon == "cf. albemarlensis"),"Taxon"] <- "albemarlensis"

#df_rads[which(df_rads$Taxon == "ustulatus pallescens"),] <- "ustulatus"
#df_rads[which(df_rads$Taxon == "ustulatus phlegonis"),] <- "ustulatus"
#df_rads[which(df_rads$Taxon == "ustulatus mahogany"),] <- "ustulatus"

df_rads[which(df_rads$Taxon == "invalidus 1"),"Taxon"] <- "invalidus"
df_rads[which(df_rads$Taxon == "invalidus 2"),"Taxon"] <- "invalidus"

df_rads[which(df_rads$Taxon == "sculpturatus 1"),"Taxon"] <- "sculpturatus"
df_rads[which(df_rads$Taxon == "sculpturatus 2"),"Taxon"] <- "sculpturatus"
df_rads[which(df_rads$Taxon == "sculpturatus 3"),"Taxon"] <- "sculpturatus"

df_rads[which(df_rads$Taxon == "wolfi 1"),"Taxon"] <- "wolfi"
df_rads[which(df_rads$Taxon == "wolfi 2"),"Taxon"] <- "wolfi"
df_rads[which(df_rads$Taxon == "wolfi 3"),"Taxon"] <- "wolfi"

df_rads[which(df_rads$Taxon == "cf. perspectivus"),"Taxon"] <- "perspectivus"
df_rads[which(df_rads$Taxon == "perspectivus 1"),"Taxon"] <- "perspectivus"
df_rads[which(df_rads$Taxon == "perspectivus 2"),"Taxon"] <- "perspectivus"

df_rads[which(df_rads$Taxon == "cf. amastroides"),"Taxon"] <- "amastroides"

df_rads[which(df_rads$Taxon == "canaliferus 1"),"Taxon"] <- "canaliferus"
df_rads[which(df_rads$Taxon == "canaliferus 2"),"Taxon"] <- "canaliferus"

df_rads[which(df_rads$Taxon == "simrothi 1"),"Taxon"] <- "simrothi"
df_rads[which(df_rads$Taxon == "simrothi 2"),"Taxon"] <- "simrothi"

df_rads[which(df_rads$Taxon == "cf. tortuganus"),"Taxon"] <- "tortuganus"

df_rads[which(df_rads$Taxon == "cf. nux"),"Taxon"] <- "nux"


# Create a mapping between Sample and Taxon
sample_to_taxon <- setNames(df_rads$Taxon, df_rads$Sample)



# Replace tip labels in the tree with corresponding taxa
tree$tip.label <- sample_to_taxon[tree$tip.label]

plot(tree)


pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, df_meanload$species))

plot(pruned.tree)

#should be ready to use phylosig

K_nemaload <- phylosig(pruned.tree, meanload, method = 'K', test = TRUE, se = s.error)
K_nemaload
plot.phylosig(K_nemaload)

#why can't I use the standard error with lambda?
lambda_nemaload <- phylosig(tree, meanload, method = 'lambda', test = TRUE) #, se = s.error) - not sure why it doesn't like this
lambda_nemaload
plot.phylosig(lambda_nemaload)