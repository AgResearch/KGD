genofile <- "HapMap.hmc.txt.gz"

source("/Code/GBS-Chip-Gmatrix.R")
Gfull <- calcG()
GHWdgm.05 <- calcG(which(HWdis > -0.05),"HWdgm.05", npc=4)  # recalculate using Hardy-Weinberg disequilibrium cut-off at -0.05

pedfile <- "Ped-GBS.csv"
groupsfile <- "Ped-Groups.csv"

rel.thresh <- 0.2
GCheck <- "GHWdgm.05$G5"
source("/Code/GBSPedAssign.R")

#G5 <- GHWdgm.05$G5
#save(G5,file="G5.RData")
