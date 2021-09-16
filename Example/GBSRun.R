#!/usr/bin/env Rscript

genofile <- "HapMap.hmc.txt.gz"

#source("<source directory>/GBS-Chip-Gmatrix.R")
source("../../GBS-Chip-Gmatrix.R")
Gfull <- calcG()
GHWdgm.05 <- calcG(which(HWdis > -0.05),"HWdgm.05", npc=4)  # recalculate using Hardy-Weinberg disequilibrium cut-off at -0.05

pedfile <- "Ped-GBS.csv"
groupsfile <- "Ped-Groups.csv"

rel.thresh <- 0.2
emm.thresh <- 0.075  # to make results same as before emm used (to match original example)
GCheck <- "GHWdgm.05$G5"
#source("<source directory>/GBSPedAssign.R")
source("../../GBSPedAssign.R")

#G5 <- GHWdgm.05$G5
#save(G5,seqID,file="G5.RData")
