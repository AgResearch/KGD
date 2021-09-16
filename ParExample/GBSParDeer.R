genofile <- "HapMap.hmc.txt.gz"
pedfile <- "DeerPedGBS.csv"
groupsfile <- "Ped-Groups.csv"
sampdepth.thresh <- 0.3
cex.pointsize <- 1.2
functions.only <- TRUE
sink("GBSParDeerOut.txt", split=TRUE)
source("GBS-Chip-Gmatrix.R")
readGBS()
outlevel <- 1  # reduced QC output
GBSsummary()

breed <- read.table(text=seqID,sep="_",stringsAsFactors=FALSE)[,1]
fcolo <- c("darkblue","darkred")[match(breed,c("W","R"))]

snpsubset <- which(HWdis > -0.05)
GHW <- calcG(npc=4,snpsubset=snpsubset,sfx="RWHW")
G5 <- GHW$G5
GCheck <- "G5"
set.seed(230985)  # to get same bootstrap results if rerun
source("GBSPedAssign.R")
RWResults <- GBSPed()
write.csv(RWResults$BothMatches,"BothMatchesRW.csv",row.names=FALSE,quote=FALSE) # Combined breeds 

dir.create("W")
setwd("W")
indW <- which(breed=="W")
pW <- calcp(indsubset=indW)
snpsubset <- which(HWdis > -0.05 & pW > 0 & pW < 1)
 GHWW <- calcG(snpsubset,indsubset=indW,sfx="W",puse=pW,calclevel=1, npc=-2)  # using Hardy-Weinberg disequilibrium cut-off at -0.05
G5W <- GHWW$G5
seqIDW <- seqID[indW]; if(length(GHWW$samp.removed) > 0 ) seqIDW <-  seqIDW[-GHWW$samp.removed]
GCheck <- "G5W"
puse <- pW
indsubset <- indW
rm(minr4inb)
pedfile <- "../../DeerPedGBS.csv"
groupsfile <- "../../Ped-Groups.csv"
WResults <- GBSPed()
MatchesW <- WResults$BothMatches
write.csv(MatchesW,"BothMatchesW.csv",row.names=FALSE,quote=FALSE) 
bestparPCA(GHWW, sfx="W",keypos="bottomright", pedinfo=WResults$pedinfo, BothMatches=MatchesW)

# Alt models
 uY <- which(MatchesW$BothAssign=="Y")
 bbopt <- optimize(ssbbmm,lower=0,upper=20, tol=0.001, uuse=uY, pedinfo=WResults$pedinfo, BothMatches=MatchesW)
 depth2K <- depth2Kchoose (dmodel="bb", bbopt$minimum) # 4.609
 mmstatsW.bb <- mismatch.2par(MatchesW$IndivID,MatchesW$BestFatherMatch, MatchesW$BestMotherMatch,pedinfo=WResults$pedinfo) 
 names(mmstatsW.bb) <- paste0(names(mmstatsW.bb),".bb")
 mpopt <- optimize(ssmpmm,lower=0.5,upper=0.9, tol=0.001, uuse=uY, pedinfo=WResults$pedinfo, BothMatches=MatchesW)
 depth2K <- depth2Kchoose (dmodel="modp", mpopt$minimum)  # 0.591
 mmstatsW.mp <- mismatch.2par(MatchesW$IndivID,MatchesW$BestFatherMatch, MatchesW$BestMotherMatch,pedinfo=WResults$pedinfo) 
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
 GHWR <- calcG(snpsubset,indsubset=indR,sfx="R",puse=pR,calclevel=1, npc=-2)  # using Hardy-Weinberg disequilibrium cut-off at -0.05
G5R <- GHWR$G5
seqIDR <- seqID[indR]; if(length(GHWR$samp.removed) > 0 ) seqIDR <-  seqIDR[-GHWR$samp.removed]
GCheck <- "G5R"
puse <- pR
indsubset <- indR
rm(minr4inb)
RResults <- GBSPed()
MatchesR <- RResults$BothMatches
write.csv(MatchesR, "BothMatchesR.csv", row.names=FALSE, quote=FALSE) 
bestparPCA(GHWR, sfx="R",keypos="bottomright", pedinfo=RResults$pedinfo, BothMatches=MatchesR)

# Alt models
 uY <- which(MatchesR$BothAssign=="Y")
 bbopt <- optimize(ssbbmm,lower=0,upper=20, tol=0.001, uuse=uY, pedinfo=RResults$pedinfo, BothMatches=MatchesR)
 depth2K <- depth2Kchoose (dmodel="bb", bbopt$minimum) # 3.956
 mmstatsR.bb <- mismatch.2par(MatchesR$IndivID,MatchesR$BestFatherMatch, MatchesR$BestMotherMatch,pedinfo=RResults$pedinfo)  
 names(mmstatsR.bb) <- paste0(names(mmstatsR.bb),".bb")
 mpopt <- optimize(ssmpmm,lower=0.5,upper=0.8, tol=0.001, uuse=uY, pedinfo=RResults$pedinfo, BothMatches=MatchesR)
 depth2K <- depth2Kchoose (dmodel="modp", mpopt$minimum)  # 0.604
 mmstatsR.mp <- mismatch.2par(MatchesR$IndivID,MatchesR$BestFatherMatch, MatchesR$BestMotherMatch,pedinfo=RResults$pedinfo) 
 names(mmstatsR.mp) <- paste0(names(mmstatsR.mp),".mp")
 MatchesR  <- cbind(MatchesR,mmstatsR.bb,mmstatsR.mp)
 write.csv(MatchesR, "BothMatchesR.csv", row.names=FALSE, quote=FALSE) 
 depth2K <- depth2Kchoose (dmodel="modp")  # back to default model

setwd("..")

sink()
