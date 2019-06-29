genofile <- "HapMap.hmc.txt.gz"
sampdepth.thresh <- 0.3
cex.pointsize <- 1.2
functions.only <- TRUE
sink("GBSParDeerOut.txt")
source("/home/doddsk/GBS/Code/GBS-Chip-Gmatrix.R")
readTassel()
outlevel <- 1  # reduced QC output
GBSsummary()

pedfile <- "../DeerPedGBS.csv"
groupsfile <- "../Ped-Groups.csv"

breed <- read.table(text=seqID,sep="_",stringsAsFactors=FALSE)[,1]
fcolo <- c("darkblue","darkred")[match(breed,c("W","R"))]

snpsubset <- which(HWdis > -0.05)
GHW <- calcG(npc=4,snpsubset=snpsubset,sfx="RWHW")
G5 <- GHW$G5
GCheck <- "G5"
set.seed(230985)  # to get same bootstrap results if rerun
source("/home/doddsk/GBS/Code/GBSPedAssign.R")
write.csv(BothMatches,"BothMatchesRW.csv",row.names=FALSE,quote=FALSE) # Combined breeds 

dir.create("W")
setwd("W")
indW <- which(breed=="W")
pW <- calcp(indsubset=indW)
snpsubset <- which(HWdis > -0.05 & pW > 0 & pW < 1)
 GHWW <- calcG(snpsubset,indsubset=indW,sfx="W",puse=pW,calclevel=1)  # using Hardy-Weinberg disequilibrium cut-off at -0.05
G5W <- GHWW$G5
seqIDW <- seqID[indW]; if(length(GHWW$samp.removed) > 0 ) seqIDW <-  seqIDW[-GHWW$samp.removed]
GCheck <- "G5W"
puse <- pW
indsubset <- indW
rm(minr4inb)
pedfile <- "../../DeerPedGBS.csv"
groupsfile <- "../../Ped-Groups.csv"
source("/home/doddsk/GBS/Code/GBSPedAssign.R")
MatchesW <- BothMatches
write.csv(MatchesW,"BothMatchesW.csv",row.names=FALSE,quote=FALSE) 

# Alt models
 uY <- which(MatchesW$BothAssign=="Y")
 bbopt <- optimize(ssbbmm,lower=0,upper=20, tol=0.001)
 depth2K <- depth2Kchoose (dmodel="bb", bbopt$minimum) # 4.609
 mmstatsW.bb <- mismatch.2par(BothMatches$IndivID,BothMatches$BestFatherMatch, BothMatches$BestMotherMatch) 
 names(mmstatsW.bb) <- paste0(names(mmstatsW.bb),".bb")
 mpopt <- optimize(ssmpmm,lower=0.5,upper=0.9, tol=0.001)
 depth2K <- depth2Kchoose (dmodel="modp", mpopt$minimum)  # 0.591
 mmstatsW.mp <- mismatch.2par(BothMatches$IndivID,BothMatches$BestFatherMatch, BothMatches$BestMotherMatch)
 names(mmstatsW.mp) <- paste0(names(mmstatsW.mp),".mp")
 MatchesW  <- cbind(MatchesW,mmstatsW.bb,mmstatsW.mp)
 write.csv(MatchesW,"BothMatchesW.csv",row.names=FALSE,quote=FALSE) 
 depth2K <- depth2Kchoose (dmodel="modp")  # back to default model
setwd("..")


dir.create("R")
setwd("R")
indR <- which(breed=="R")
pR <- calcp(indsubset=indR)
snpsubset <- which(HWdis > -0.05 & pR > 0 & pR < 1)
 GHWR <- calcG(snpsubset,indsubset=indR,sfx="R",puse=pR,calclevel=1)  # using Hardy-Weinberg disequilibrium cut-off at -0.05
G5R <- GHWR$G5
seqIDR <- seqID[indR]; if(length(GHWR$samp.removed) > 0 ) seqIDR <-  seqIDR[-GHWR$samp.removed]
GCheck <- "G5R"
puse <- pR
indsubset <- indR
rm(minr4inb)
source("/home/doddsk/GBS/Code/GBSPedAssign.R")
MatchesR <- BothMatches
write.csv(MatchesR,"BothMatchesR.csv",row.names=FALSE,quote=FALSE) 

# Alt models
 uY <- which(MatchesR$BothAssign=="Y")
 bbopt <- optimize(ssbbmm,lower=0,upper=20, tol=0.001)
 depth2K <- depth2Kchoose (dmodel="bb", bbopt$minimum) # 3.956
 mmstatsR.bb <- mismatch.2par(BothMatches$IndivID,BothMatches$BestFatherMatch, BothMatches$BestMotherMatch) 
 names(mmstatsR.bb) <- paste0(names(mmstatsR.bb),".bb")
 mpopt <- optimize(ssmpmm,lower=0.5,upper=0.8, tol=0.001)
 depth2K <- depth2Kchoose (dmodel="modp", mpopt$minimum)  # 0.604
 mmstatsR.mp <- mismatch.2par(BothMatches$IndivID,BothMatches$BestFatherMatch, BothMatches$BestMotherMatch)
 names(mmstatsR.mp) <- paste0(names(mmstatsR.mp),".mp")
 MatchesR  <- cbind(MatchesR,mmstatsR.bb,mmstatsR.mp)
 write.csv(MatchesR,"BothMatchesR.csv",row.names=FALSE,quote=FALSE) 
 depth2K <- depth2Kchoose (dmodel="modp")  # back to default model
setwd("..")

sink()
