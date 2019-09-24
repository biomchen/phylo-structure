
#########################################################################
## Phylogenetic Community Structure of small-bodied mammal communities ##
#########################################################################

require(ape)
require(picante)

## Reading data
occur<-read.csv("PCS_Data_5Feb2016.csv",header=T, sep=",")

# create community by species occurrence matrix
spec_occur<-xtabs(~occur$Community+occur$Species)

# get rid of any values that might be greater than one
for(i in 1:nrow(spec_occur)){
  for(j in 1:ncol(spec_occur)){
    if(spec_occur[i,j]>1){
      spec_occur[i,j]<-1
    }
  }
}

write.csv(spec_occur,"Species_by_community_matrix.csv")

spec_occur<-read.csv("Species_by_community_matrix.csv",header=T, row.names=1)

# calculate phycomm
fritz_tree<-read.nexus("Fritz_tree.nex")
spec_Meng<-occur$Species

fritz_tree1<-fritz_tree[[1]]
taxa_in_both<-intersect(fritz_tree1$tip.label, spec_Meng) # 294 species
mismatches<-setdiff(spec_Meng,taxa_in_both)

drop_set<-setdiff(fritz_tree1$tip.label, taxa_in_both)
pruned_tree<-drop.tip(fritz_tree1, tip=drop_set)

spec_occur<-spec_occur[, pruned_tree$tip.label]

## Alternatives
occur<-read.csv("species_names.csv",header=T, sep=",")

# create community by species occurrence matrix
spec_occur<-xtabs(~occur$Community+occur$Species)

# get rid of any values that might be greater than one
for(i in 1:nrow(spec_occur)){
  for(j in 1:ncol(spec_occur)){
    if(spec_occur[i,j]>1){
      spec_occur[i,j]<-1
    }
  }
}

# calculate phycomm
fritz_tree<-read.nexus("Fritz_tree.nex")
spec_Meng<-occur$Species

fritz_tree1<-fritz_tree[[1]]
taxa_in_both<-intersect(fritz_tree1$tip.label, spec_Meng) # 294 species
mismatches<-setdiff(spec_Meng,taxa_in_both)

drop_set<-setdiff(fritz_tree1$tip.label, taxa_in_both)
pruned_tree<-drop.tip(fritz_tree1, tip=drop_set)

spec_occur<-spec_occur[,prunedtree$tip.label]

# Net relatedness Index
phy_dist<- cophenetic(pruned_tree)
system.time(ses_mpd_result <- ses.mpd(spec_occur, phy_dist, null.model = "taxa.labels",
                                      abundance.weighted = FALSE, runs = 100))
NRI<--1*(ses_mpd_result$mpd.obs.z)

community<-unique(occur$NO)

community_nri<-data.frame(cbind(community, NRI))

plot(community_nri$community,community_nri$NRI)

# By different climte
spec_occur<-xtabs(~occur$Climate.Code+occur$Species)

# get rid of any values that might be greater than one
for(i in 1:nrow(spec_occur)){
  for(j in 1:ncol(spec_occur)){
    if(spec_occur[i,j]>1){
      spec_occur[i,j]<-1
    }
  }
}

write.csv(spec_occur,"Species_by_climate_code_matrix.csv")

# calculate phycomm
fritz_tree<-read.nexus("Fritz tree.nex")

fritz_tree1<-fritz_tree[[1]]
taxa_in_both<-intersect(fritz_tree1$tip.label, colnames(spec_occur)) # 266 species
mismatches<-setdiff(colnames(spec_occur), taxa_in_both)

drop_set<-setdiff(fritz_tree1$tip.label, taxa_in_both)
pruned_tree<-drop.tip(fritz_tree1, tip=drop_set)

spec_occur<-spec_occur[,pruned_tree$tip.label]

# Net relatedness Index
phy_dist<- cophenetic(pruned_tree)
system.time(ses_mpd_result <- ses.mpd(spec_occur, phy_dist, null.model = "taxa.labels",
                                      abundance.weighted = FALSE, runs = 1000))

NRI<--1*(ses_mpd_result$mpd.obs.z)

community<-unique(occur$Climate.Code)

community_nri<-data.frame(cbind(community,NRI))

plot(community_nri$community,community_nri$NRI)
