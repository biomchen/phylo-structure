
#########################################################################
## Phylogenetic Community Structure of small-bodied mammal communities ##
#########################################################################

## Reading data
occur<-read.csv("PCS Data 5Feb2016.csv",header=T, sep=",")

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

write.csv(spec_occur,"Species by community matrix.csv")

spec_occur<-read.csv("Species by community matrix.csv",header=T,row.names=1)


# calculate phycomm
require(ape)

fritz_tree<-read.nexus("Fritz_tree.nex")
spec_Meng<-occur$Species

fritz_tree1<-fritz_tree[[1]]
taxa.in.both<-intersect(fritz_tree1$tip.label, spec_Meng) # 294 species
mismatches<-setdiff(spec_Meng,taxa.in.both)

dropset<-setdiff(fritz_tree1$tip.label,taxa.in.both)
prunedtree<-drop.tip(fritz_tree1,tip=dropset)

spec_occur<-spec_occur[,prunedtree$tip.label]

## Alternatives
occur<-read.csv("Chen_Wilson_Species_2015.csv",header=T, sep=",")

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
require(ape)

fritz_tree<-read.nexus("Fritz_tree.nex")
spec_Meng<-occur$Species

fritz_tree1<-fritz_tree[[1]]
taxa.in.both<-intersect(fritz_tree1$tip.label, spec_Meng) # 294 species
mismatches<-setdiff(spec_Meng,taxa.in.both)

dropset<-setdiff(fritz_tree1$tip.label,taxa.in.both)
prunedtree<-drop.tip(fritz_tree1,tip=dropset)

spec_occur<-spec_occur[,prunedtree$tip.label]

# Net relatedness Index
phydist<- cophenetic(prunedtree)

require(picante)

system.time(ses.mpd.result <- ses.mpd(spec_occur,phydist,null.model = "taxa.labels",
                                      abundance.weighted = FALSE, runs = 100))

NRI<--1*(ses.mpd.result$mpd.obs.z)

community<-unique(occur$NO)

community_nri<-data.frame(cbind(community,NRI))

plot(community_nri$community,community_nri$NRI)


# By climate code
spec_occur<-xtabs(~occur$Climate.Code+occur$Species)

# get rid of any values that might be greater than one

for(i in 1:nrow(spec_occur)){
  for(j in 1:ncol(spec_occur)){
    if(spec_occur[i,j]>1){
      spec_occur[i,j]<-1
    }
  }
}

write.csv(spec_occur,"Species by climate code matrix.csv")

# calculate phycomm

require(ape)

fritz_tree<-read.nexus("Fritz tree.nex")

fritz_tree1<-fritz_tree[[1]]
taxa.in.both<-intersect(fritz_tree1$tip.label, colnames(spec_occur)) # 266 species
mismatches<-setdiff(colnames(spec_occur),taxa.in.both)

dropset<-setdiff(fritz_tree1$tip.label,taxa.in.both)
prunedtree<-drop.tip(fritz_tree1,tip=dropset)

spec_occur<-spec_occur[,prunedtree$tip.label]

# Net relatedness Index
phydist<- cophenetic(prunedtree)

require(picante)

system.time(ses.mpd.result <- ses.mpd(spec_occur,phydist,null.model = "taxa.labels",
                                      abundance.weighted = FALSE, runs = 1000))

NRI<--1*(ses.mpd.result$mpd.obs.z)

community<-unique(occur$Climate.Code)

community_nri<-data.frame(cbind(community,NRI))

plot(community_nri$community,community_nri$NRI)

##
