#!/bin/echo Source me don't execute me 

# assume all in pedfile are in the genotype results. To do: remove those that are not
OK4ped <- TRUE
if (!exists("rel.thresh"))   rel.thresh  <- 0.4
if (!exists("rel.threshF"))  rel.threshF <- rel.thresh
if (!exists("rel.threshM"))  rel.threshM <- rel.thresh
if (!exists("emm.thresh"))   emm.thresh  <- 0.01           # Excess for single parent match
if (!exists("emm.thresh2"))  emm.thresh2 <- 2*emm.thresh   # Excess for parent-pair match
if (!exists("emmdiff.thresh2"))  emmdiff.thresh2 <-0  # alternate parent-pair based on emm
if (!exists("inb.thresh"))   inb.thresh  <- 0.2       # par relatedness - inbreeding
if (!exists("boota.thresh")) boota.thresh <- 99       # assignment threshold
if (!exists("mindepth.mm")) mindepth.mm <- 1            # changed to 1 to coincide to change to using exp mm rate
if (!exists("indsubset")) indsubset <- seq_along(seqID)
if (!exists("snpsubset")) snpsubset <- seq(nsnps)
if (!exists("depth.min")) depth.min <- 0    # for bootstrapping
if (!exists("depth.max")) depth.max <- Inf  # for bootstrapping
if (!exists("puse")) puse <- p
if (!exists("nboot")) nboot <- 1000  # for bootstrapping
if (!exists("boot.thresh")) boot.thresh <- 0.05 # rel diff for invoking bootstrapping
cat("Parentage parameter settings\n----------------------------\n rel.threshF\t",rel.threshF,
    "\n rel.threshM\t",rel.threshM,
    "\n emm.thresh\t",emm.thresh,
    "\n emm.thresh2\t",emm.thresh2,
    "\n emmdiff.thresh2\t",emmdiff.thresh2,
    "\n inb.thresh\t",inb.thresh," (parent relatedness v inbreeding)",
    "\n boota.thresh\t",boota.thresh,
    "\n depth.min\t",depth.min," (for bootstrapping)",
    "\n depth.max\t",depth.max," (for bootstrapping)",
    "\n nboot\t\t",nboot,
    "\n boot.thresh\t",boot.thresh," (relatedness difference to invoke bootstrapping)",
    "\n")
if(length(indsubset) != nrow(eval(parse(text = GCheck)))) {
 OK4ped <- FALSE
 cat("Number of individuals",length(indsubset),"does not match G matrix",nrow(eval(parse(text = GCheck))),"\n")
 } 

panel.yeqx <- function(x,y,col.points="black",col.line="red",...){   #panel function for pairs with identity line added
    points(x,y,col=col.points,...)
    abline(a = 0,b = 1, col=col.line, ...)
}
coordprop <- function(propn,crange) crange[1]+propn*diff(crange)  # function to find proportional positions on plots

mismatch.par <- function(offspring.id, par.id) {
  # ids as in the pedigree file, if only 1 parent compare with all offspring
  noffspring <- length(offspring.id)
  if (length(par.id) == 1) par.id <- rep(par.id, noffspring)
  nmismatch <- ncompare <- exp.mmrate <- rep(NA, noffspring)
  opos <- match(pedinfo$seqID[match(offspring.id, pedinfo$IndivID)], seqID)
  ppos <- match(pedinfo$seqID[match(par.id, pedinfo$IndivID)], seqID)
  for (ioffspring in 1:noffspring) {
   depthi <- depth.orig[opos[ioffspring], ]
   depthj <- depth.orig[ppos[ioffspring], ]
   usnp <- intersect(snpsubset,which(depthi >= mindepth.mm & depthj >= mindepth.mm))
   pi <- genon[opos[ioffspring], usnp]/2
   pj <- genon[ppos[ioffspring], usnp]/2
   Ko  <-  depth2K(depthi[usnp])
   Kp  <-  depth2K(depthj[usnp])
   nmismatch[ioffspring] <- length(which(abs(pi - pj) == 1))
   ncompare[ioffspring] <- length(usnp)
   ptemp <- puse[usnp]
   P <- ptemp*(1-ptemp)
   expmm <- rep(NA,length(usnp))
   ug <- which(pi==1)
   if(length(ug)>0) expmm[ug] <- P[ug] *(ptemp[ug]*Kp[ug] + (1-ptemp[ug])*Ko[ug] + Kp[ug]*Ko[ug] ) / (ptemp[ug]^2+2*P[ug]*Ko[ug])
   ug <- which(pi==0.5)
   if(length(ug)>0) expmm[ug] <- 0
   ug <- which(pi==0)
   if(length(ug)>0) expmm[ug] <- P[ug] *(ptemp[ug]*Ko[ug] + (1-ptemp[ug])*Kp[ug] + Kp[ug]*Ko[ug] ) / ((1-ptemp[ug])^2+2*P[ug]*Ko[ug])
   exp.mmrate[ioffspring] <- mean(expmm,na.rm=TRUE)
  }
  mmrate <- nmismatch/ncompare
  list(mmrate=mmrate,ncompare=ncompare,exp.mmrate=exp.mmrate)
}

mismatch.2par <- function(offspring.id, par1.id, par2.id,alph=Inf) {
  noffspring <- length(offspring.id)
  nmismatch <- ncompare <- exp.mmrate <- rep(NA, noffspring)
  opos <- match(pedinfo$seqID[match(offspring.id, pedinfo$IndivID)], seqID)
  p1pos <- match(pedinfo$seqID[match(par1.id, pedinfo$IndivID)], seqID)
  p2pos <- match(pedinfo$seqID[match(par2.id, pedinfo$IndivID)], seqID)
  for (ioffspring in 1:noffspring) {
   depthi <- depth.orig[opos[ioffspring], ]
   depthj <- depth.orig[p1pos[ioffspring], ]
   depthk <- depth.orig[p2pos[ioffspring], ]
   usnp <- intersect(snpsubset,which(pmin(depthi,depthj,depthk) >= mindepth.mm))
   pi <- genon[opos[ioffspring], usnp]/2
   pj <- genon[p1pos[ioffspring], usnp]/2
   pk <- genon[p2pos[ioffspring], usnp]/2
   Ko  <-  depth2K(depthi[usnp])
   Kf  <-  depth2K(depthj[usnp])
   Km  <-  depth2K(depthk[usnp])
   ptemp <- puse[usnp]
   P <- ptemp*(1-ptemp)
   expmm <- rep(NA,length(usnp))
   ug <- which(pi==1)
   if(length(ug)>0) expmm[ug] <- 
        ( ptemp[ug]^2 * P[ug] *(Km[ug]+Kf[ug])*(1+ Ko[ug]) + P[ug]^2*(2*Ko[ug] + Km[ug] + Kf[ug] - Kf[ug] *Km[ug] +2*Km[ug] *Ko[ug] +2 *Kf[ug] *Ko[ug] - 2 *Kf[ug] *Km[ug] *Ko[ug]) +
                 2*P[ug] *(1-ptemp[ug])^2 * Ko[ug] ) / (ptemp[ug]^2+2*P[ug]*Ko[ug])
   ug <- which(pi==0.5)
   if(length(ug)>0) expmm[ug] <- ( (1-2*P[ug]) *(Km[ug]+Kf[ug]) +4*P[ug]*Kf[ug] *Km[ug] ) / 2
   ug <- which(pi==0)
   if(length(ug)>0) expmm[ug] <- 
       ( 2*ptemp[ug]^2 * P[ug] * Ko[ug]  + P[ug]^2*(2*Ko[ug] + Km[ug] + Kf[ug] - Kf[ug] *Km[ug] +2*Km[ug] *Ko[ug] +2 *Kf[ug] *Ko[ug] - 2 *Kf[ug] *Km[ug] *Ko[ug]) +
                 P[ug] *(1-ptemp[ug])^2 * (Km[ug]+Kf[ug])*(1+ Ko[ug])) / ((1-ptemp[ug])^2+2*P[ug]*Ko[ug])
   exp.mmrate[ioffspring] <- mean(expmm,na.rm=TRUE)
   nmismatch[ioffspring] <- length(which(abs(pi - pj) == 1 | abs(pi - pk) == 1 | (pj == pk & !pj==0.5 & pi == 0.5)))
   ncompare[ioffspring] <- length(usnp)
   }
  mmrate <- nmismatch/ncompare
  list(mmrate=mmrate,ncompare=ncompare,exp.mmrate=exp.mmrate)
}

parmatch <- function(partype, Gmatrix) {
  if (missing(partype)) partype <- "Father"
  ParseqID <- with(pedinfo, seqID[match(pedinfo[, paste0(partype, "ID")], IndivID)])
  offspringpos <- match(pedinfo$seqID, seqID[indsubset])  # all in pedigree file considered as offspring
  parpos <- match(ParseqID, seqID[indsubset])
  ParRel <- Gmatrix[cbind(offspringpos, parpos)]
  ParMatch <- (ParRel > ifelse(partype == "Father", rel.threshF,rel.threshM))
  png(paste0(partype, "Verify.png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
   pairs(cbind(ParRel, offspringpos, parpos), labels = c("Relatedness", "Offspring order", paste(partype, "order")))
   dev.off()
  ncompare <- sum(!is.na(ParMatch))
  nmatch <- sum(ParMatch, na.rm = TRUE)
  matchperc <- 100 * nmatch/ncompare
  cat(nmatch, "matches out of", ncompare, partype, "comparisons:", format(matchperc, digits = 3), "%\n")
  if (nmatch > 0) 
    cat("Mean relatedness for", partype, "matches", format(mean(ParRel[which(ParMatch)]), digits = 3), "\n")
  if (ncompare > nmatch) 
    cat("Mean relatedness for", partype, "non-matches", format(mean(ParRel[which(!ParMatch)]), digits = 3), "\n")
  matchinfo <- data.frame(ParRel, ParMatch)
  names(matchinfo) <- paste0(partype, c("Rel", "Match"))
  cbind(pedinfo, matchinfo)
}

bestmatch <- function(ospos, parpos, Guse, partype) {
  if (missing(partype)) 
    partype <- "Par"
  parchk <- Guse[ospos, parpos,drop=FALSE]
  maxpos <- apply(parchk, 1, which.max)
  parchktemp <- parchk
  parchktemp[cbind(1:nrow(parchk), maxpos)] <- -1
  maxpos.2 <- apply(parchktemp, 1, which.max)
  rm(parchktemp)
  maxrel <- cbind(parchk[cbind(1:nrow(parchk), maxpos)], parchk[cbind(1:nrow(parchk), maxpos.2)])
  rel12 <- Guse[cbind(parpos[maxpos],parpos[maxpos.2])]
  out.df <- data.frame(seqID[indsubset][ospos], seqID[indsubset][parpos[maxpos]], seqID[indsubset][parpos[maxpos.2]], 
            relatedness = maxrel[, 1], rel2nd = maxrel[, 2], rel12=rel12, stringsAsFactors = FALSE)
  names(out.df) <- c("seqID", paste0("Best", partype, "Match"), paste0(partype, "Match2nd"), paste0(partype, "rel"), paste0(partype, "rel2nd"), paste0(partype, "12rel") )
  out.df
}

groupmatch <- function(Guse, partype) {
  rel.thresh <-  ifelse(partype == "Father", rel.threshF,rel.threshM)
  groupIDs <- unique(pedinfo[, paste0(partype, "Group")])
  groupIDs <- na.omit(groupIDs)
  groupIDs <- groupIDs[!groupIDs == ""]
  ngroups <- length(groupIDs)
  if (ngroups > 0) {
    for (g in 1:ngroups) {
      group <- groupIDs[g]
      offspringID <- pedinfo$IndivID[which(pedinfo[, paste0(partype, "Group")] == group)]
      ParGroupID <- groupsinfo$IndivID[which(groupsinfo$ParGroup == group)]
      offspringseqID <- with(pedinfo, seqID[match(offspringID, IndivID)])
      ParGroupseqID <- with(pedinfo, seqID[match(ParGroupID, IndivID)])
      offspringpos <- match(offspringseqID, seqID[indsubset])
      parpos <- match(ParGroupseqID, seqID[indsubset])
      gmatch <- bestmatch(offspringpos, parpos, Guse, partype)
      if (g == 1) 
        allgmatch <- gmatch else allgmatch <- rbind(allgmatch, gmatch)
    }
    allgmatch$IndivID <- pedinfo$IndivID[match(allgmatch$seqID, pedinfo$seqID)]
    ncolallg <- ncol(allgmatch)
    allgmatch <- allgmatch[, c(ncolallg, 1:(ncolallg-1))]
    allgmatch[, paste0("Best", partype, "Match")] <- pedinfo$IndivID[match(allgmatch[, paste0("Best", partype, "Match")], pedinfo$seqID)]
    allgmatch[, paste0(partype, "Match2nd")] <- pedinfo$IndivID[match(allgmatch[, paste0(partype, "Match2nd")], pedinfo$seqID)]
    mmstats <- mismatch.par(allgmatch$IndivID, allgmatch[, paste0("Best", partype, "Match")])
    allgmatch[, paste0("mmrate", partype)] <- mmstats$mmrate
    allgmatch[, paste0("mmnum", partype)] <- mmstats$ncompare
    allgmatch[, paste0("exp.mmrate", partype)] <- mmstats$exp.mmrate
    mmstats <- mismatch.par(allgmatch$IndivID, allgmatch[, paste0(partype, "Match2nd")])
    allgmatch[, paste0("mmrate", partype,"2")] <- mmstats$mmrate
    allgmatch[, paste0("exp.mmrate", partype,"2")] <- mmstats$exp.mmrate
    EMMrate <- allgmatch[, paste0("mmrate", partype)] - allgmatch[, paste0("exp.mmrate", partype)]
    EMMrate2 <- allgmatch[, paste0("mmrate", partype, "2")] - allgmatch[, paste0("exp.mmrate", partype, "2")] 
    ### bootstrap section
    bootpos <- which(allgmatch[,paste0(partype, "rel")] > rel.thresh & allgmatch[,paste0(partype, "rel")] - allgmatch[,paste0(partype, "rel2nd")] < boot.thresh)
    nsnpsub <- length(snpsubset)
    if(length(bootpos) > 0) {
     for(bcase in seq_along(bootpos)) {
      offspringseqID <- allgmatch$seqID[bootpos[bcase]]
      ParseqID <- with(pedinfo, seqID[match(allgmatch[bootpos[bcase], paste0("Best", partype, "Match")], IndivID)])
      Par2seqID <- with(pedinfo, seqID[match(allgmatch[bootpos[bcase], paste0(partype, "Match2nd")], IndivID)])
      offspringpos <- match(offspringseqID, seqID[indsubset])
      parpos <- match(ParseqID, seqID[indsubset])
      par2pos <- match(Par2seqID, seqID[indsubset])
      offs.depth <- depth.orig[indsubset[offspringpos], snpsubset]
      par.depth <-  depth.orig[indsubset[parpos], snpsubset]
      par2.depth <-  depth.orig[indsubset[par2pos], snpsubset]
      offs.genon0 <- genon[indsubset[offspringpos], snpsubset]
      par.genon0 <- genon[indsubset[parpos], snpsubset]
      par2.genon0 <- genon[indsubset[par2pos], snpsubset]
      if (depth.min > 1 | depth.max < Inf) {
       offs.genon0[offs.depth < depth.min] <- NA
       offs.genon0[offs.depth > depth.max] <- NA
       par.genon0[par.depth < depth.min] <- NA
       par.genon0[par.depth > depth.max] <- NA
       offs.depth[is.na(offs.genon0)] <- 0
       par.depth[is.na(par.genon0)] <- 0
       }
      offs.usegeno <- !is.na(offs.genon0)
      par.usegeno <- !is.na(par.genon0)
      par2.usegeno <- !is.na(par2.genon0)
      offs.genon0 <- offs.genon0 - rep.int(2 * puse[snpsubset], rep(1, nsnpsub))
      offs.genon0[is.na(offs.genon0)] <- 0     # equivalent to using 2p for missing genos
      par.genon0 <- par.genon0 - rep.int(2 * puse[snpsubset], rep(1, nsnpsub))
      par.genon0[is.na(par.genon0)] <- 0     # equivalent to using 2p for missing genos
      par2.genon0 <- par2.genon0 - rep.int(2 * puse[snpsubset], rep(1, nsnpsub))
      par2.genon0[is.na(par2.genon0)] <- 0     # equivalent to using 2p for missing genos
#      relcheck <- tcrossprod(matrix(offs.genon0,nrow=1),matrix(par.genon0,nrow=1)) / 
#                  sum((2*puse[snpsubset]*(1-puse[snpsubset]))[offs.usegeno & par.usegeno])
      bootrels <-     bootrels2 <- double(nboot)
      for (b in seq(nboot)) {
       bootsnps <- sample.int(nsnpsub,replace=TRUE)
       bootrels[b] <- tcrossprod(matrix(offs.genon0[bootsnps],nrow=1),matrix(par.genon0[bootsnps],nrow=1)) / 
                  sum((2*puse[snpsubset[bootsnps]]*(1-puse[snpsubset[bootsnps]]))[offs.usegeno[bootsnps] & par.usegeno[bootsnps]])
       bootrels2[b] <- tcrossprod(matrix(offs.genon0[bootsnps],nrow=1),matrix(par2.genon0[bootsnps],nrow=1)) / 
                  sum((2*puse[snpsubset[bootsnps]]*(1-puse[snpsubset[bootsnps]]))[offs.usegeno[bootsnps] & par2.usegeno[bootsnps]])
       }
      allgmatch[bootpos[bcase],paste0(partype,"sd")] <- sd(bootrels) 
      allgmatch[bootpos[bcase],paste0(partype,"Reliability")] <- 100*sum(bootrels > bootrels2) /nboot
      }
     } #bootpos
    ### end of bootstrap section
    tempAssign <- rep("Y",nrow(allgmatch))
    if(length(bootpos) > 0) tempAssign[which(allgmatch[,paste0(partype,"Reliability")] < boota.thresh )] <- "B"
    tempAssign[which(EMMrate > emm.thresh)] <- "E"
    tempAssign[which(allgmatch[, paste0(partype, "rel")] < rel.thresh)] <- "N"
     # for E assigns, check if the 2nd parent is possible.  
    tempAssign[which(allgmatch[, paste0(partype, "rel2nd")] >= rel.thresh & EMMrate2 < emm.thresh & tempAssign == "E")] <- "A"
    allgmatch[, paste0(partype, "Assign")] <- tempAssign
    png(paste0("Best", partype, "Matches.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
     plot(allgmatch[, paste0("mmrate", partype)] ~ allgmatch[, paste0(partype, "rel")], main = paste("Best", partype, "Matches"), xlab = "Estimated Relatedness", 
         ylab = "Raw mismatch rate",col=fcolo[match(allgmatch$seqID,seqID)], cex=0.8)
     abline(v=rel.thresh, col="grey")
     dev.off()
    png(paste0("Best", partype, "MatchesE.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
     plot(EMMrate ~ allgmatch[, paste0(partype, "rel")], main = paste("Best", partype, "Matches"), xlab = "Estimated Relatedness", 
         ylab = "Excess mismatch rate",col=fcolo[match(allgmatch$seqID,seqID)], cex=0.8)
     abline(v=rel.thresh, col="grey")
     abline(h=emm.thresh, col="grey")
     dev.off()
    tempch <- assign.pch[match(tempAssign,assign.rank)]
    png(paste0("ExpMM-", partype, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
     plot(allgmatch[, paste0("mmrate", partype)] ~ allgmatch[, paste0("exp.mmrate", partype)] , main = paste("Best", partype, "Matches"), xlab = "Expected mismatch rate", 
         ylab = "Raw mismatch rate",col=fcolo[match(allgmatch$seqID,seqID)], cex=0.8, pch=tempch)
     pch.used <- sort(match(unique(tempAssign),assign.rank))
     legend("bottomright",title="Assign",cex=0.75,pch=assign.pch[pch.used],legend=assign.rank[pch.used])
     abline(a=0,b=1,col="red")
     abline(a=emm.thresh,b=1,col="grey")
     dev.off()
    mmpalette <- colorRampPalette(c("blue","red"))(50)
    mmcol <- mmpalette[trunc(1+50*(EMMrate-min(EMMrate,na.rm=TRUE))/(diff(range(EMMrate,na.rm=TRUE))+1E-6))]
    legend_image <- as.raster(matrix(rev(mmpalette), ncol = 1))
    xyrange <- range(c(allgmatch[, paste0(partype, "rel")],allgmatch[, paste0(partype,"rel2nd")]))
    png(paste0("Best2", partype, "Matches.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
     plot(allgmatch[, paste0(partype,"rel2nd")] ~ allgmatch[, paste0(partype, "rel")], main = paste("Best", partype, "Matches"), xlab = "Estimated Relatedness", 
         ylab = "Relatedness to 2nd best",col=mmcol,xlim=xyrange,ylim=xyrange, cex=0.8)
     abline(a=0,b=1)
     abline(v=rel.thresh, col="grey")
     abline(h=rel.thresh, col="grey")
     rasterImage(legend_image, coordprop(0.05,xyrange), coordprop(0.7,xyrange), coordprop(0.1,xyrange), coordprop(0.9,xyrange))
     text(x=coordprop(0.11,xyrange),y=coordprop(0.7,xyrange),signif(min(EMMrate,na.rm=TRUE),2),pos=4,cex=0.8)
     text(x=coordprop(0.11,xyrange),y=coordprop(0.9,xyrange),signif(max(EMMrate,na.rm=TRUE),2),pos=4,cex=0.8)
     text(x=coordprop(0,xyrange),y=coordprop(0.95,xyrange),"Excess MM rate best",pos=4,cex=0.8)
     dev.off()
    write.csv(allgmatch, paste0(partype, "Matches.csv"), row.names = FALSE)
    noffspringpar <- data.frame(table(allgmatch[, paste0("Best", partype, "Match")]))
    colnames(noffspringpar)[2] <- paste0(partype, "Freq")
    groupsinfo <<- merge(groupsinfo, noffspringpar, by.x = "IndivID", by.y = "Var1", all = TRUE)
    groupsinfo[is.na(groupsinfo[, paste0(partype, "Freq")]), paste0(partype, "Freq")] <<- 0
    allgmatch
  } else {
    NULL
  }
}

if (OK4ped & exists("pedfile") & exists("GCheck")) {
  pedinfo <- read.csv(pedfile, stringsAsFactors = FALSE, colClasses=c(FatherGroup="character", MotherGroup="character"))
  pedinfo <- pedinfo[!is.na(pedinfo$seqID), ]
  pedinfo <- pedinfo[!is.na(match(pedinfo$seqID, seqID[indsubset])), ]
  dupids <- which(duplicated(pedinfo$IndivID))
  if(length(dupids)>0) cat("Warning: dupicates in IndivID",pedinfo$IndivID[dupids],"\n")
  dupids <- which(duplicated(pedinfo$seqID))
  if(length(dupids)>0) cat("Warning: dupicates in seqID",pedinfo$seqID[dupids],"\n")
  if ("FatherID" %in% colnames(pedinfo)) 
    pedinfo <- parmatch("Father", eval(parse(text = GCheck)))
  if ("MotherID" %in% colnames(pedinfo)) 
    pedinfo <- parmatch("Mother", eval(parse(text = GCheck)))
  if ("FatherID" %in% colnames(pedinfo) & "MotherID" %in% colnames(pedinfo)) {
    if (is.character(pedinfo$IndivID)) {
      umiss <- which(pedinfo$FatherID == "")
      if (length(umiss) > 0) 
        pedinfo$FatherID[umiss] <- NA
      umiss <- which(pedinfo$MotherID == "")
      if (length(umiss) > 0) 
        pedinfo$MotherID[umiss] <- NA
    }
    famtable <- with(pedinfo, table(FatherID, MotherID))
    fampos <- which(famtable > 1, arr.ind = TRUE)
    famfathers <- dimnames(famtable)$FatherID[fampos[, 1]]
    fammothers <- dimnames(famtable)$MotherID[fampos[, 2]]
    if (is.numeric(pedinfo$FatherID)) 
      famfathers <- as.numeric(famfathers)
    if (is.numeric(pedinfo$FatherID)) 
      famfathers <- as.numeric(famfathers)
    noffspring <- famtable[fampos]
    nfamilies <- length(noffspring)
    famnumber <- rep(NA, nrow(pedinfo))
    famresults <- rep(NA, nfamilies)
    for (ifam in 1:nfamilies) {
      famnumber[which(pedinfo$FatherID == famfathers[ifam] & pedinfo$MotherID == fammothers[ifam])] <- ifam
      uoffspring <- match(pedinfo$seqID[which(famnumber == ifam)], seqID[indsubset])
      famresults[ifam] <- mean(eval(parse(text = GCheck))[uoffspring, uoffspring][upper.tri(diag(nrow = length(uoffspring)))])
    }
    cat("Mean relatedness for full-sib families (as given)\n")
    print(data.frame(famfathers, fammothers, noffspring, meanrel = famresults))
    cat("Mean relatedness within all full-sib families", weighted.mean(famresults, noffspring), "\n")
    
    uoffspring <- which(!is.na(famnumber))
    Fatherset <- unique(na.omit(pedinfo$FatherID))
    Motherset <- unique(na.omit(pedinfo$MotherID))
    udiff.fathers <- which(!match(pedinfo$FatherID[uoffspring], Fatherset) %*% t(rep(1, length(uoffspring))) == rep(1, length(uoffspring)) %*% 
                             t(match(pedinfo$FatherID[uoffspring], Fatherset)), arr.ind = T)
    udiff.mothers <- which(!match(pedinfo$MotherID[uoffspring], Motherset) %*% t(rep(1, length(uoffspring))) == rep(1, length(uoffspring)) %*% 
                             t(match(pedinfo$MotherID[uoffspring], Motherset)), arr.ind = T)
    udiff <- as.matrix(merge(udiff.fathers, udiff.mothers))
    cat("Mean relatedness between individuals in full-sib families with different parents", mean(eval(parse(text = GCheck))[uoffspring, 
                                                                                                                            uoffspring][udiff]), "\n")
  }
  write.csv(pedinfo, "PedVerify.csv", row.names = FALSE)
  if (exists("groupsfile")) {
    assign.rank <- c("Y","I","B","A","E","F","M","N")
    assign.pch <- c(16,2,6,1,15,13,13,4)
    groupsinfo <- read.csv(groupsfile, stringsAsFactors = FALSE, colClasses=(ParGroup="character"))
    groupsinfo <- groupsinfo[!duplicated(groupsinfo),]  # remove row duplicates
    groupsinfo$Genotyped <- ifelse(is.na(match(pedinfo$seqID[match(groupsinfo$IndivID,pedinfo$IndivID)],seqID)),"N","Y")
    if ("FatherGroup" %in% colnames(pedinfo)) 
      FatherMatches <- groupmatch(eval(parse(text = GCheck)), "Father")
    if ("MotherGroup" %in% colnames(pedinfo)) 
      MotherMatches <- groupmatch(eval(parse(text = GCheck)), "Mother")
    write.csv(groupsinfo, "GroupsParentCounts.csv", row.names = FALSE)
    if ("FatherGroup" %in% colnames(pedinfo) & "MotherGroup" %in% colnames(pedinfo)) {
     BothMatches <- merge(FatherMatches,MotherMatches)
     if(nrow(BothMatches) > 0) {
      BothMatches <- BothMatches[order(BothMatches$IndivID),]
      mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$BestFatherMatch, BothMatches$BestMotherMatch)
      BothMatches$mmrateF1M1 <- mmstats$mmrate
      BothMatches$mmnumF1M1 <- mmstats$ncompare
      BothMatches$exp.mmrateF1M1 <- mmstats$exp.mmrate
      mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$FatherMatch2nd, BothMatches$BestMotherMatch)
      BothMatches$mmrateF2M1 <- mmstats$mmrate
      BothMatches$mmnumF2M1 <- mmstats$ncompare
      BothMatches$exp.mmrateF2M1 <- mmstats$exp.mmrate
      mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$BestFatherMatch, BothMatches$MotherMatch2nd)
      BothMatches$mmrateF1M2 <- mmstats$mmrate
      BothMatches$mmnumF1M2 <- mmstats$ncompare
      BothMatches$exp.mmrateF1M2 <- mmstats$exp.mmrate
      mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$FatherMatch2nd, BothMatches$MotherMatch2nd)
      BothMatches$mmrateF2M2 <- mmstats$mmrate
      BothMatches$mmnumF2M2 <- mmstats$ncompare
      BothMatches$exp.mmrateF2M2 <- mmstats$exp.mmrate
      # parent relatedness & Inbreeding
      uf <- match(pedinfo$seqID[match(BothMatches$BestFather,pedinfo$IndivID)],seqID[indsubset])
      um <- match(pedinfo$seqID[match(BothMatches$BestMother,pedinfo$IndivID)],seqID[indsubset])
      BothMatches$relF1M1 <- eval(parse(text = GCheck))[cbind(uf,um)]
      uf <- match(pedinfo$seqID[match(BothMatches$BestFather,pedinfo$IndivID)],seqID[indsubset])
      um <- match(pedinfo$seqID[match(BothMatches$MotherMatch2nd,pedinfo$IndivID)],seqID[indsubset])
      BothMatches$relF1M2 <- eval(parse(text = GCheck))[cbind(uf,um)]
      uf <- match(pedinfo$seqID[match(BothMatches$FatherMatch2nd,pedinfo$IndivID)],seqID[indsubset])
      um <- match(pedinfo$seqID[match(BothMatches$BestMother,pedinfo$IndivID)],seqID[indsubset])
      BothMatches$relF2M1 <- eval(parse(text = GCheck))[cbind(uf,um)]
      uf <- match(pedinfo$seqID[match(BothMatches$FatherMatch2nd,pedinfo$IndivID)],seqID[indsubset])
      um <- match(pedinfo$seqID[match(BothMatches$MotherMatch2nd,pedinfo$IndivID)],seqID[indsubset])
      BothMatches$relF2M2 <- eval(parse(text = GCheck))[cbind(uf,um)]
      uo <- match(BothMatches$seqID,seqID[indsubset])
      BothMatches$Inb <- diag(eval(parse(text = GCheck)))[uo] - 1
      BothMatches$BothAssign <- assign.rank[pmax(match(BothMatches$FatherAssign,assign.rank),match(BothMatches$MotherAssign,assign.rank))]
      EMMrates <- with(BothMatches,cbind(mmrateF2M2-exp.mmrateF2M2,mmrateF1M2-exp.mmrateF1M2,mmrateF2M1-exp.mmrateF2M1,mmrateF1M1-exp.mmrateF1M1))
      EMMrate.min <- apply(EMMrates, MARGIN=1, min)
      BothMatches$BothAssign[which(EMMrates[,4] > emm.thresh2 & BothMatches$BothAssign %in% assign.rank[1:4])] <- "E"
      BothMatches$BothAssign[EMMrates[,4]-EMMrate.min > emmdiff.thresh2 & EMMrate.min > emm.thresh2 & BothMatches$BothAssign %in% assign.rank[1:5]] <- "A"
      BothMatches$BothAssign[which(BothMatches$relF1M1 - BothMatches$Inb > inb.thresh & BothMatches$BothAssign == "Y")] <- "I"
      BothMatches$BothAssign[which(BothMatches$FatherAssign == "Y" & BothMatches$BothAssign == "N")] <- "F"
      BothMatches$BothAssign[which(BothMatches$MotherAssign == "Y" & BothMatches$BothAssign == "N")] <- "M"
      BothMatches$BothAssign[which(BothMatches$FatherAssign == "Y" & BothMatches$BothAssign == "E" & BothMatches$MotherAssign == "E")] <- "F"
      BothMatches$BothAssign[which(BothMatches$MotherAssign == "Y" & BothMatches$BothAssign == "E" & BothMatches$FatherAssign == "E")] <- "M"

      BothMatches$Alternate <- ""
      Apos <- which(BothMatches$BothAssign=="A")
      for (ipos in Apos) {
       altpar = c("F2M2","F1M2","F2M1")[which(EMMrates[ipos,] == EMMrate.min[ipos])]
       altOK <- TRUE
       if (BothMatches[ipos,paste0("rel",altpar)] - BothMatches$Inb[ipos] > inb.thresh | EMMrate.min[ipos] > emm.thresh2) altOK <- FALSE
       if (grepl("F2",altpar)) {
        if(BothMatches[ipos, "Fatherrel2nd"] < rel.threshF | BothMatches[ipos, "mmrateFather2"] - BothMatches[ipos, "exp.mmrateFather2"] > emm.thresh ) altOK <- FALSE
        }
       if (grepl("M2",altpar)) {
        if(BothMatches[ipos, "Motherrel2nd"] < rel.threshM | BothMatches[ipos, "mmrateMother2"] - BothMatches[ipos, "exp.mmrateMother2"] > emm.thresh ) altOK <- FALSE
        }
       if(altOK) BothMatches$Alternate[ipos] <- altpar
       }
      write.csv(BothMatches,"BothMatches.csv", row.names = FALSE)
      uo <- match(BothMatches$seqID,seqID)
      plotch <- assign.pch[match(BothMatches$BothAssign,assign.rank)]
      plotcol <- rep("mediumblue",length(uo))
      plotcol[which(sampdepth[uo] < 1 )] <- "skyblue2"
      plotcol[which(sampdepth[uo] < 0.5 ) ] <- "grey75"
      png("ParRel-Inb.png", width = 960, height = 960, pointsize = cex.pointsize *  21)
       plot(BothMatches$relF1M1 ~ BothMatches$Inb,pch=plotch,col=plotcol,sub="0.5 <= mean depth < 1",col.sub="skyblue2",
            xlab="Estimated Inbreeding",ylab="Estimated (best match) parent relatedness")
       title(sub="X: unassigned parent(s)",adj=0)
       title(sub="mean depth < 0.5",col.sub="grey75",adj=0.95)
       dev.off()
      png(paste0("ExpMM-Both.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
       plot(mmrateF1M1 ~ exp.mmrateF1M1, data=BothMatches, main = paste("Best Parent Matches"), xlab = "Expected mismatch rate", 
           ylab = "Raw mismatch rate",col=fcolo[uo], pch=plotch, cex=0.8)
       abline(a=0,b=1,col="red")
       abline(a=emm.thresh2,b=1,col="grey")
       pch.used <- sort(match(unique(BothMatches$BothAssign),assign.rank))
       legend("bottomright",title="Assign",cex=0.75,pch=assign.pch[pch.used],legend=assign.rank[pch.used])
       dev.off()
      png("MMrateBoth.png", width = 640, height = 640, pointsize = cex.pointsize *  15)
       pairs(with(BothMatches,cbind(mmrateF2M2,mmrateF1M2,mmrateF2M1,mmrateF1M1)),upper.panel=panel.yeqx,lower.panel=NULL,
                  main="Raw Mismatch Rates", labels=c("Father2,\nMother2","Father1,\nMother2","Father2,\nMother1","Father1,\nMother1"),
                  col.points=fcolo[uo],pch=plotch)
       dev.off()
      png("MMrateBothE.png", width = 640, height = 640, pointsize = cex.pointsize *  15)
       pairs(EMMrates, main="Excess Mismatch Rates", labels=c("Father2,\nMother2","Father1,\nMother2","Father2,\nMother1","Father1,\nMother1"),
                  upper.panel=panel.yeqx,lower.panel=NULL,col.points=fcolo[uo],pch=plotch)
       dev.off()
      print( table(BothMatches$BothAssign) )
      }
     }
  }
}


if (FALSE) { # temporary working code
  # 
 }

