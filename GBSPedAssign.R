#!/bin/echo Source me don't execute me 
pedver <- "1.4.0"
cat("GBS-PedAssign for KGD version:",pedver,"\n")

  verif.ch <- c(".","Y","N")  # NA, Y, N
   assign.rank <<- c("Y","I","B","A","E","F","M","N")
#   assign.pch <- c(16,2,6,1,15,13,13,4)
   assign.pch <<- c(16,2,6,1,15,70,77,4)

# assume all in pedfile are in the genotype results. To do: remove those that are not
pedsetup <- function() {
 OK4ped <<- TRUE
 if (!exists("developer"))    developer   <<- FALSE
 if (!exists("rel.thresh"))   rel.thresh  <<- 0.4
 if (!exists("rel.threshF"))  rel.threshF <<- rel.thresh
 if (!exists("rel.threshM"))  rel.threshM <<- rel.thresh
 if (!exists("emm.thresh"))   emm.thresh  <<- 0.01           # Excess for single parent match
 if (!exists("emm.thresh2"))  emm.thresh2 <<- 2*emm.thresh   # Excess for parent-pair match
 if (!exists("doublemm"))     doublemm <<- FALSE   # count 2 for AA x AA = BB types?
 if (!exists("emmdiff.thresh2"))  emmdiff.thresh2 <<- 0  # alternate parent-pair based on emm
 if (!exists("inb.thresh"))   inb.thresh  <<- 0.2       # par relatedness - 2 * inbreeding
 if (!exists("minr4inb"))     minr4inb    <<- NULL      # par relatedness - inbreeding
 if (!exists("boota.thresh")) boota.thresh <<- 99       # assignment threshold
 if (!exists("mindepth.mm")) mindepth.mm <<- 1          # changed to 1 to coincide to change to using exp mm rate
 if (!exists("matchmethod")) matchmethod <<- "rel"      # choose best 2 parents based on "rel" or "EMM"
 if (!exists("indsubset")) indsubset <<- seq_along(seqID)
 if (!exists("snpsubset")) snpsubset <<- seq(nsnps)
 if (!exists("depth.min")) depth.min <<- 0    # for bootstrapping
 if (!exists("depth.max")) depth.max <<- Inf  # for bootstrapping
 if (!exists("puse")) puse <<- p
 if (!exists("nboot")) nboot <<- 1000  # for bootstrapping
 if (!exists("boot.thresh")) boot.thresh <<- 0.05 # rel diff for invoking bootstrapping
 if (!exists("allow.selfing")) allow.selfing <<- FALSE # if FALSE find best parent pair without selfing
 cat("Parentage parameter settings\n----------------------------\n rel.threshF\t",rel.threshF,
     "\n rel.threshM\t",rel.threshM,
     "\n emm.thresh\t",emm.thresh,
     "\n doublemm\t",doublemm,
     "\n emm.thresh2\t",emm.thresh2,
     "\n emmdiff.thresh2",emmdiff.thresh2,
     "\n mindepth.mm\t",mindepth.mm,
     "\n inb.thresh\t",inb.thresh,"(parent relatedness - 2 * inbreeding)",
     "\n minr4inb\t",minr4inb,
     "\n boota.thresh\t",boota.thresh,
     "\n depth.min\t",depth.min,"(for bootstrapping)",
     "\n depth.max\t",depth.max,"(for bootstrapping)",
     "\n nboot\t\t",nboot,
     "\n boot.thresh\t",boot.thresh,"(relatedness difference to invoke bootstrapping)",
     "\n matchmethod\t",matchmethod,
     "\n allow.selfing\t",allow.selfing,
     "\n")
 if(length(indsubset) != nrow(eval(parse(text = GCheck)))) {
  OK4ped <<- FALSE
  cat("Number of individuals",length(indsubset),"does not match G matrix",nrow(eval(parse(text = GCheck))),"\n")
  } 
 if (exists("pedfile")) if(!file.exists(pedfile)) {
  OK4ped <<- FALSE
  cat("Warning: Pedigree file", pedfile, "not found\n")
  }
 } #pedsetup

panel.yeqx <- function(x,y,col.points="black",col.line="red",...){   #panel function for pairs with identity line added
    points(x,y,col=col.points,...)
    abline(a = 0,b = 1, col=col.line, ...)
}
coordprop <- function(propn,crange) crange[1]+propn*diff(crange)  # function to find proportional positions on plots

mismatch.par <- function(offspring.id, par.id, pedinfo) {
  # ids as in the pedigree file, if only 1 parent compare with all offspring
  noffspring <- length(offspring.id)
  if (length(par.id) == 1) par.id <- rep(par.id, noffspring)
  nmismatch <- ncompare <- exp.mmrate <- rep(NA, noffspring)
  opos <- match(pedinfo$seqID[match(offspring.id, pedinfo$IndivID)], seqID)
  ppos <- match(pedinfo$seqID[match(par.id, pedinfo$IndivID)], seqID)
  for (ioffspring in 1:noffspring) {
   depthi <- depth[opos[ioffspring], ]
   depthj <- depth[ppos[ioffspring], ]
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

cemult <- function(x,A) t(x*t(A))  # columnwise elementwise mult of A by x (usually length(x) = nrow(A))
ceadd  <- function(x,A) t(x+t(A))  # columnwise elementwise add  of A by x (usually length(x) = nrow(A))

mismatch.par.comb <- function(offspring.id, par.id, pedinfo) {
  # ids as in the pedigree file, all combs of offspring and parent (for given parent type)
  noffspring <- length(offspring.id)
  nmismatch <- ncompare <- exp.mmrate <- matrix(NA, nrow=noffspring,ncol=length(par.id))
  opos <- match(pedinfo$seqID[match(offspring.id, pedinfo$IndivID)], seqID)
  ppos <- match(pedinfo$seqID[match(par.id, pedinfo$IndivID)], seqID)
  depthj <- depth[ppos, ]
  Kpall  <-  depth2K(depthj)
  for (ioffspring in 1:noffspring) {
   depthi <- depth[opos[ioffspring], ]
   usnp <- intersect(snpsubset,which(depthi >= mindepth.mm)) 
   pi <- genon[opos[ioffspring], usnp]/2
   pj <- genon[ppos, usnp]/2
   Ko  <-  depth2K(depthi[usnp])
   Kp  <-  Kpall[,usnp]
   nmismatch[ioffspring,] <- rowSums((abs(matrix(pi,nrow=nrow(pj),ncol=ncol(pj),byrow=TRUE) - pj) == 1),na.rm=TRUE)
   ncompare[ioffspring,] <- rowSums(!is.na(pj)) 
   ptemp <- puse[usnp]
   P <- ptemp*(1-ptemp)
   expmm <- matrix(NA,ncol=length(usnp),nrow=length(par.id))
   ug <- which(pi==1)
   if(length(ug)>0) expmm[,ug] <- cemult(P[ug]/ (ptemp[ug]^2+2*P[ug]*Ko[ug]), ceadd((1-ptemp[ug])*Ko[ug],cemult(ptemp[ug],Kp[,ug])  + cemult(Ko[ug],Kp[,ug])) ) 
   ug <- which(pi==0.5)
   if(length(ug)>0) expmm[,ug] <- 0
   ug <- which(pi==0)
   if(length(ug)>0) expmm[,ug] <- cemult(P[ug]/ ((1-ptemp[ug])^2+2*P[ug]*Ko[ug]) , ceadd(ptemp[ug]*Ko[ug], cemult((1-ptemp[ug]),Kp[,ug]) + cemult(Ko[ug],Kp[,ug]) )) 
   expmm[which(is.na(pj) | (depthj[,usnp] < mindepth.mm)) ] <- NA
   exp.mmrate[ioffspring,] <- rowMeans(expmm,na.rm=TRUE)
  }
  mmrate <- nmismatch/ncompare
  list(mmrate=mmrate,ncompare=ncompare,exp.mmrate=exp.mmrate)
}

mismatch.2par <- function(offspring.id, par1.id, par2.id,alph=Inf, pedinfo) {
  noffspring <- length(offspring.id)
  nmismatch <- ncompare <- exp.mmrate <- rep(NA, noffspring)
  opos <- match(pedinfo$seqID[match(offspring.id, pedinfo$IndivID)], seqID)
  p1pos <- match(pedinfo$seqID[match(par1.id, pedinfo$IndivID)], seqID)
  p2pos <- match(pedinfo$seqID[match(par2.id, pedinfo$IndivID)], seqID)
  for (ioffspring in 1:noffspring) {
   depthi <- depth[opos[ioffspring], ]
   depthj <- depth[p1pos[ioffspring], ]
   depthk <- depth[p2pos[ioffspring], ]
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
   if(!doublemm) {
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
    } else { # doublemm
    ug <- which(pi==1)
    if(length(ug)>0) expmm[ug] <- 
        ( ptemp[ug]^2 * P[ug] *(Km[ug]+Kf[ug])*(1+ Ko[ug]) + P[ug]^2*(2*Ko[ug] + Km[ug] + Kf[ug] +2*Km[ug] *Ko[ug] +2 *Kf[ug] *Ko[ug]) +
                 2*P[ug] *(1-ptemp[ug])^2 * Ko[ug] * (1+ Km[ug] + Kf[ug]) ) / (ptemp[ug]^2+2*P[ug]*Ko[ug])
    ug <- which(pi==0.5)
    if(length(ug)>0) expmm[ug] <- ( (1-2*P[ug]) *(Km[ug]+Kf[ug]) +4*P[ug]*Kf[ug] *Km[ug] ) / 2
    ug <- which(pi==0)
    if(length(ug)>0) expmm[ug] <- 
       ( 2*ptemp[ug]^2 * P[ug] * Ko[ug] * (1+ Km[ug] + Kf[ug])  + P[ug]^2*(2*Ko[ug] + Km[ug] + Kf[ug] +2*Km[ug] *Ko[ug] +2 *Kf[ug] *Ko[ug]) +
                 P[ug] *(1-ptemp[ug])^2 * (Km[ug]+Kf[ug])*(1+ Ko[ug])) / ((1-ptemp[ug])^2+2*P[ug]*Ko[ug])
    }
   exp.mmrate[ioffspring] <- mean(expmm,na.rm=TRUE)
   nmismatch[ioffspring] <- length(which(abs(pi - pj) == 1 | abs(pi - pk) == 1 | (pj == pk & !pj==0.5 & pi == 0.5)))
   if(doublemm) nmismatch[ioffspring] <-  nmismatch[ioffspring] + length(which(abs(pi - pj) == 1 & pj==pk ))
   ncompare[ioffspring] <- length(usnp)
   }
  mmrate <- nmismatch/ncompare
  list(mmrate=mmrate,ncompare=ncompare,exp.mmrate=exp.mmrate)
}

parmatch <- function(partype, Gmatrix, pedinfo) {
  if (missing(partype)) partype <- "Father"
  ParseqID <- with(pedinfo, seqID[match(pedinfo[, paste0(partype, "ID")], IndivID)])
  offspringpos <- match(pedinfo$seqID, seqID[indsubset])  # all in pedigree file considered as offspring
  parpos <- match(ParseqID, seqID[indsubset])
  ParRel <- Gmatrix[cbind(offspringpos, parpos)]
  Parmm <- mismatch.par(pedinfo$IndivID,pedinfo[, paste0(partype, "ID")], pedinfo)
  ParEMM <- Parmm$mmrate - Parmm$exp.mmrate
  rel.threshpar <- ifelse(partype == "Father", rel.threshF,rel.threshM)
  ParMatch <- (ParRel > rel.threshpar & ParEMM < emm.thresh)
  noffsp <- length(offspringpos)
  tempch <- verif.ch[match(ParMatch, c(NA,TRUE,FALSE))]
  if(any(!is.na(ParRel))) {
   relrange <- seq(min(ParRel,na.rm=TRUE),max(ParRel,na.rm=TRUE),diff(range(ParRel,na.rm=TRUE))/20)
   emmrange <- seq(min(ParEMM,na.rm=TRUE),max(ParEMM,na.rm=TRUE),diff(range(ParEMM,na.rm=TRUE))/20)
   png(paste0(partype, "Verify.png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
    nthresh <- length(c(relrange,emmrange))
    plotch <- c(rep(1,noffsp),rep(46,nthresh)) # circles, dots for thresholds
    pairs(cbind(c(ParRel, relrange,rep(rel.threshpar, length(emmrange))),c(ParEMM, rep(emm.thresh,length(relrange)),emmrange), c(offspringpos, rep(NA,nthresh)), c(parpos, rep(NA,nthresh))), 
          labels = c("Relatedness", "EMM", "Offspring order", paste(partype, "order")), gap=0, pch=plotch)
    dev.off()
   }
  ncompare <- sum(!is.na(ParMatch))
  nmatch <- sum(ParMatch, na.rm = TRUE)
  matchperc <- 100 * nmatch/ncompare
  cat(nmatch, "matches out of", ncompare, partype, "comparisons:", format(matchperc, digits = 3), "%\n")
  if (nmatch > 0) 
    cat("Mean relatedness for", partype, "matches", format(mean(ParRel[which(ParMatch)]), digits = 3), "\n")
  if (ncompare > nmatch) 
    cat("Mean relatedness for", partype, "non-matches", format(mean(ParRel[which(!ParMatch)]), digits = 3), "\n")
  ParInb <- diag(Gmatrix)[parpos] - 1
  matchinfo <- data.frame(ParRel, ParEMM, ParMatch,ParInb)
  names(matchinfo) <- paste0(partype, c("Rel", "EMM", "Match","Inb"))
  pedinfo <- cbind(pedinfo, matchinfo)
  tempinfo <- pedinfo; tempinfo[,paste0(partype, "rel")] <- tempinfo[,paste0(partype, "Rel")]
  parEplot(partype,ParEMM,tempinfo[,paste0(partype, "Rel")],matchtype="Rec",
     plotcol=fcolo[match(pedinfo$seqID,seqID)],relthresh=rel.threshpar) 

  png(paste0("ExpMM-Rec", partype, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
   plot(Parmm$mmrate ~ Parmm$exp.mmrate , main = paste("Rec", partype, "Matches"), xlab = "Expected mismatch rate", 
       ylab = "Raw mismatch rate",col=fcolo[match(pedinfo$seqID,seqID)], cex=0.8, pch=tempch)
#   legend("bottomright",title="Verify",cex=0.75,pch=verif.pch,legend=c("","Y","N")
   abline(a=0,b=1,col="red")
   abline(a=emm.thresh,b=1,col="grey")
   edges <- par("usr") # xl,xr,yb,yt
   poly1 <- data.frame(x1 = c(edges[1], edges[1], edges[2],edges[2],edges[1]), 
                       y1 = c(emm.thresh+edges[1],edges[4],edges[4],emm.thresh+edges[2],emm.thresh+edges[1]))
   polygon(poly1,col = rgb(0,0,0,alpha=0.1),border=NA) 
   dev.off()

  if(developer) {
   Exprel <- 0.5 + pedinfo$Inb + ParInb/2   # par-offspr rel from theory
   reldevn <- ParRel - Exprel
   lmrel <- lm(ParRel ~ pedinfo$Inb + ParInb )
   print(summary(lmrel))
   reldevplots(partype,ParEMM,reldevn,matchtype="Rec",gmatch=tempinfo,Exprelvar=Exprel,plotcol=fcolo[match(pedinfo$seqID,seqID)],plotch=1) 
   }
  pedinfo
}

bestmatch <- function(ospos, parpos, Guse, partype, matchcriterion = "rel", groupname=group, pedinfo) {
  #matchcriterion == "EMM" added later, but made to work exactly the same (will redo EMM for best 2 later)
  if(!matchcriterion == "EMM") matchcriterion <- "rel"
  if (missing(partype)) partype <- "Par"
  groupsize <- length(na.omit(parpos))
  if(groupsize > 0 ) {
   diag(Guse) <- -1     # prevent self-parenting
   parchk <- Guse[ospos, parpos,drop=FALSE]
   maxpos.2 <- rep(NA,nrow(parchk))
   if(matchcriterion == "rel") {
    maxpos <- apply(parchk, 1, which.max)
    if(groupsize > 1) {
     parchktemp <- parchk
     parchktemp[cbind(1:nrow(parchk), maxpos)] <- -1
     maxpos.2 <- apply(parchktemp, 1, which.max)
     rm(parchktemp)
     }
    }
   if(matchcriterion == "EMM") {
    offspringID.bm <- pedinfo$IndivID[match(seqID[indsubset][ospos],pedinfo$seqID)]
    parGroupID.bm <- pedinfo$IndivID[match(seqID[indsubset][parpos],pedinfo$seqID)]
    mm.bm <- mismatch.par.comb(offspringID.bm,parGroupID.bm, pedinfo)
    EMMchk <- with(mm.bm,mmrate-exp.mmrate)
    maxpos <- apply(EMMchk, 1, which.min)
    EMMchk[cbind(1:nrow(EMMchk), maxpos)] <- 1
    if(groupsize > 1) maxpos.2 <- apply(EMMchk, 1, which.min)  
    }
   maxrel <- cbind(parchk[cbind(1:nrow(parchk), maxpos)], parchk[cbind(1:nrow(parchk), maxpos.2)])
   rel12 <- Guse[cbind(parpos[maxpos],parpos[maxpos.2])]
   out.df <- data.frame(seqID[indsubset][ospos], seqID[indsubset][parpos[maxpos]], seqID[indsubset][parpos[maxpos.2]], 
             relatedness = maxrel[, 1], rel2nd = maxrel[, 2], rel12=rel12, stringsAsFactors = FALSE)
   } else {
#   cat("Warning: a",partype,"group has no genotyped individuals, includes offspring seqID:",seqID[indsubset][ospos][1],"\n")
   cat("Warning:",partype,"group",groupname, "has no genotyped individuals\n")
   nprog <- length(ospos)
   out.df <- data.frame(seqID[indsubset][ospos],character(nprog),character(nprog),rep(NA,nprog),rep(NA,nprog),rep(NA,nprog),stringsAsFactors = FALSE)
   }
  names(out.df) <- c("seqID", paste0("Best", partype, "Match"), paste0(partype, "Match2nd"), paste0(partype, "rel"), paste0(partype, "rel2nd"), paste0(partype, "12rel") )
  out.df
}

bestmatesmatch <- function(ospos, matespos, Guse, matchcriterion = "rel", pedinfo) {
  #matespos has a column for each mate (male, female)
  if(!matchcriterion == "EMM") matchcriterion <- "rel"
  if(nrow(matespos) > 0 ) {
   selfrel <- diag(Guse)
   diag(Guse) <- -1     # prevent self-parenting
   par1chk <- Guse[ospos, matespos[,1],drop=FALSE]
   par2chk <- Guse[ospos, matespos[,2],drop=FALSE]
   if(matchcriterion == "rel") {
    maxpos <- apply(pmin(par1chk,par2chk), 1, which.max)  # find the best min parent rel (for each offspring)
    parchktemp <- par1chk
    parchktemp[cbind(1:nrow(parchktemp), maxpos)] <- -1  # only need to do this for one mate
    maxpos.2 <- apply(pmin(parchktemp,par2chk), 1, which.max)
    rm(parchktemp)
    }
   if(matchcriterion == "EMM") {
    offspringID.bm <- pedinfo$IndivID[match(seqID[indsubset][ospos],pedinfo$seqID)]
    parGroupID1.bm <- pedinfo$IndivID[match(seqID[indsubset][matespos[,1]],pedinfo$seqID)]
    parGroupID2.bm <- pedinfo$IndivID[match(seqID[indsubset][matespos[,2]],pedinfo$seqID)]
    # in a loop for now - is there a better way (similar to mismatch.par.comb ?)
    maxpos <- maxpos.2 <- integer(length(ospos))
    for(ipos in 1:length(ospos)) {
     mm.bm <- mismatch.2par(rep(offspringID.bm[ipos],length(matespos)),parGroupID1.bm,parGroupID2.bm,pedinfo=pedinfo)   # check offspring against all possible mate pairs
     EMMchk <- with(mm.bm,mmrate-exp.mmrate)
     maxpos[ipos] <- which.min(EMMchk)
     EMMchk[maxpos[ipos]] <- 1
     maxpos.2[ipos] <- which.min(EMMchk)
     }
    }
   diag(Guse) <- selfrel
   maxrel <- cbind(par1chk[cbind(1:nrow(par1chk), maxpos)],  par2chk[cbind(1:nrow(par2chk), maxpos)],
                   par1chk[cbind(1:nrow(par1chk), maxpos.2)],par2chk[cbind(1:nrow(par2chk), maxpos.2)]) #F1,M1, F2, M2
   rel12 <- cbind(Guse[cbind(matespos[maxpos,1],matespos[maxpos.2,1])],Guse[cbind(matespos[maxpos,2],matespos[maxpos.2,2])]) #F1F2, M1M2
   out.df <- data.frame(seqID[indsubset][ospos], seqID[indsubset][matespos[maxpos,1]], seqID[indsubset][matespos[maxpos,2]], 
             seqID[indsubset][matespos[maxpos.2,1]], seqID[indsubset][matespos[maxpos.2,2]], 
             maxrel, rel12, stringsAsFactors = FALSE)
   } else {
   out.df <- data.frame(matrix(character(0),ncol=5),matrix(numeric(0),ncol=6), stringsAsFactors = FALSE)
   }
  names(out.df) <- c("seqID", "BestFatherMatch", "BestMotherMatch","FatherMatch2nd","MotherMatch2nd",
                     "Fatherrel","Motherrel","Fatherrel2nd","Motherrel2nd","Father12rel","Mother12rel")
  out.df
}

parEplot <- function(partype,EMMvar,relvar,matchtype="Best",plotcol="black",relthresh=rel.thresh,emmthresh=emm.thresh) {
      png(paste0(matchtype, partype, "MatchesE.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
      plot(EMMvar ~ relvar, main = paste(matchtype, partype, "Matches"), xlab = "Estimated Relatedness", 
          ylab = "Excess mismatch rate",col=plotcol, cex=0.8)
      abline(v=relthresh, col="grey")
      abline(h=emmthresh, col="grey")
      edges <- par("usr") # xl,xr,yb,yt
      poly1 <- data.frame(x1 = c(edges[1],relthresh,relthresh,edges[2],edges[2],edges[1],edges[1]), 
                          y1 = c(edges[3],edges[3],emm.thresh,emm.thresh,edges[4],edges[4],edges[3]))
      polygon(poly1,col = rgb(0,0,0,alpha=0.1),border=NA) 
      dev.off()
 }

reldevplots <- function(partype,EMMvar,relvar,matchtype="Best",gmatch,Exprelvar,plotcol="black",plotch=1) { 
      png(paste0(matchtype, partype, "MatchesDE.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
       plot(EMMvar ~ relvar, main = paste(matchtype, partype, "Matches"), xlab = "Relatedness Deviation", 
           ylab = "Excess mismatch rate",col=plotcol, cex=0.8)
       abline(v=rel.thresh-0.5, col="grey")
       abline(h=emm.thresh, col="grey")
       dev.off()
      png(paste0(matchtype, partype, "MatchesExprel.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
       plot(gmatch[, paste0(partype, "rel")] ~ Exprelvar, pch=plotch, main = paste(matchtype, partype, "Matches"),
          xlab="Expected relatedness using estimated Inbreeding", ylab="Estimated relatedness")
       abline(a=0,b=1,col="red")
       abline(h=0.5,col="red",lty=2)
       dev.off()
 }

trioplots <- function(BothMatches) {
    uo <- match(BothMatches$seqID,seqID)
    plotch <- assign.pch[match(BothMatches$BothAssign,assign.rank)]
    plotcol <- rep("mediumblue",length(uo))
    plotcol[which(sampdepth[uo] < 1 )] <- "skyblue2"
    plotcol[which(sampdepth[uo] < 0.5 ) ] <- "grey75"
    uplot <- which(!is.na(BothMatches$relF1M1))
    if(length(uplot) > 0) {
      if(sum(!is.na(BothMatches$Inb[uplot])) > 0) {
      png("ParRel-Inb.png", width = 960, height = 960, pointsize = cex.pointsize *  21)
#       plot(BothMatches$relF1M1[uplot] ~ BothMatches$Inb[uplot],pch=plotch[uplot],col=plotcol[uplot],sub="0.5 <= mean depth < 1",col.sub="skyblue2",
#            xlab="Estimated Inbreeding",ylab="Estimated (best match) parent relatedness", cex.sub=0.9)
#       abline(a=0, b=2, col="red")
#       lines(x=c(min(min(BothMatches$Inb[uplot]), (minr4inb-inb.thresh)/2), (minr4inb-inb.thresh)/2,max(BothMatches$Inb[uplot])),
#             y=c(minr4inb,minr4inb,2*max(BothMatches$Inb[uplot])+inb.thresh),col="grey")
       plot(BothMatches$Inb[uplot] ~ BothMatches$relF1M1[uplot] ,pch=plotch[uplot],col=plotcol[uplot],sub="0.5 <= mean depth < 1",col.sub="skyblue2",
            ylab="Estimated Inbreeding",xlab="Estimated (best match) parent relatedness", cex.sub=0.9)
       abline(a=0, b=1/2, col="red")
       lines(y=c(min(min(BothMatches$Inb[uplot]), (minr4inb-inb.thresh)/2), (minr4inb-inb.thresh)/2,max(BothMatches$Inb[uplot])),
             x=c(minr4inb,minr4inb,2*max(BothMatches$Inb[uplot])+inb.thresh),col="grey")
       edges <- par("usr") # xl,xr,yb,yt
       poly1 <- data.frame(x1 = c(minr4inb, minr4inb, edges[2],edges[2],minr4inb), 
                           y1 = c(edges[3],(minr4inb-inb.thresh)/2,(edges[2]-inb.thresh)/2,edges[3],edges[3]))
       polygon(poly1,col = rgb(0,0,0,alpha=0.1),border=NA) 
       title(sub="X: unassigned parent(s)",adj=0,cex.sub=0.9)
       title(sub="mean depth < 0.5",col.sub="grey75",adj=0.95,cex.sub=0.8)
       dev.off()
       }
      png(paste0("ExpMM-Both.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
       plot(mmrateF1M1 ~ exp.mmrateF1M1, data=BothMatches[uplot,,drop=FALSE], main = paste("Best Parent Matches"), xlab = "Expected mismatch rate", 
           ylab = "Raw mismatch rate",col=fcolo[uo][uplot], pch=plotch[uplot], cex=0.8)
       abline(a=0,b=1,col="red")
       abline(a=emm.thresh2,b=1,col="grey")
       pch.used <- sort(match(unique(BothMatches$BothAssign[uplot]),assign.rank))
       legend("bottomright",title="Assign",cex=0.75,pch=assign.pch[pch.used],legend=assign.rank[pch.used])
       edges <- par("usr") # xl,xr,yb,yt
       poly1 <- data.frame(x1 = c(edges[1], edges[1], edges[2],edges[2],edges[1]), 
                           y1 = c(emm.thresh2+edges[1],edges[4],edges[4],emm.thresh2+edges[2],emm.thresh2+edges[1]))
       polygon(poly1,col = rgb(0,0,0,alpha=0.1),border=NA) 
       dev.off()
      }
 }

groupmatch <- function(Guse, partype, pedinfo, groupsinfo) {
  rel.thresh <-  ifelse(partype == "Father", rel.threshF,rel.threshM)
  groupIDs <- unique(pedinfo[, paste0(partype, "Group")])
  groupIDs <- na.omit(groupIDs)
  groupIDs <- groupIDs[!groupIDs == ""]
  ngroups <- length(groupIDs)
  nsnpsub <- length(snpsubset)
  if (ngroups > 0) {
    for (g in 1:ngroups) {
      group <- groupIDs[g]
      offspringID <- pedinfo$IndivID[which(pedinfo[, paste0(partype, "Group")] == group)]
      ParGroupID <- groupsinfo$IndivID[which(groupsinfo$ParGroup == group)]
      offspringseqID <- with(pedinfo, seqID[match(offspringID, IndivID)])
      ParGroupseqID <- with(pedinfo, seqID[match(ParGroupID, IndivID)])
      offspringpos <- match(offspringseqID, seqID[indsubset])
      parpos <- match(ParGroupseqID, seqID[indsubset])
      gmatch <- bestmatch(offspringpos, parpos, Guse, partype, matchcriterion = matchmethod, groupname=group,pedinfo=pedinfo)
      if (g == 1) allgmatch <- gmatch else allgmatch <- rbind(allgmatch, gmatch)
      }
    allgmatch$IndivID <- pedinfo$IndivID[match(allgmatch$seqID, pedinfo$seqID)]
    ncolallg <- ncol(allgmatch)
    allgmatch <- allgmatch[, c(ncolallg, 1:(ncolallg-1))]
    if (nrow(allgmatch) > 0) {
     allgmatch[, paste0("Best", partype, "Match")] <- pedinfo$IndivID[match(allgmatch[, paste0("Best", partype, "Match")], pedinfo$seqID)]
     allgmatch[, paste0(partype, "Match2nd")] <- pedinfo$IndivID[match(allgmatch[, paste0(partype, "Match2nd")], pedinfo$seqID)]
     mmstats <- mismatch.par(allgmatch$IndivID, allgmatch[, paste0("Best", partype, "Match")], pedinfo)
     allgmatch[, paste0("mmrate", partype)] <- mmstats$mmrate
     allgmatch[, paste0("mmnum", partype)] <- mmstats$ncompare
     allgmatch[, paste0("exp.mmrate", partype)] <- mmstats$exp.mmrate
     mmstats <- mismatch.par(allgmatch$IndivID, allgmatch[, paste0(partype, "Match2nd")], pedinfo)
     allgmatch[, paste0("mmrate", partype,"2")] <- mmstats$mmrate
     allgmatch[, paste0("exp.mmrate", partype,"2")] <- mmstats$exp.mmrate
     EMMrate <- allgmatch[, paste0("mmrate", partype)] - allgmatch[, paste0("exp.mmrate", partype)]
     EMMrate2 <- allgmatch[, paste0("mmrate", partype, "2")] - allgmatch[, paste0("exp.mmrate", partype, "2")] 
     ### bootstrap section
     bootpos <- which(allgmatch[,paste0(partype, "rel")] > rel.thresh & EMMrate < emm.thresh & allgmatch[,paste0(partype, "rel")] - allgmatch[,paste0(partype, "rel2nd")] < boot.thresh)
     bootpos <- intersect(bootpos,which(!is.na(allgmatch[, paste0(partype, "Match2nd")])))
     if(length(bootpos) > 0) {
      for(bcase in seq_along(bootpos)) {
       offspringseqID <- allgmatch$seqID[bootpos[bcase]]
       ParseqID <- with(pedinfo, seqID[match(allgmatch[bootpos[bcase], paste0("Best", partype, "Match")], IndivID)])
       Par2seqID <- with(pedinfo, seqID[match(allgmatch[bootpos[bcase], paste0(partype, "Match2nd")], IndivID)])
       offspringpos <- match(offspringseqID, seqID[indsubset])
       parpos <- match(ParseqID, seqID[indsubset])
       par2pos <- match(Par2seqID, seqID[indsubset])
       offs.depth <- depth[indsubset[offspringpos], snpsubset]
       par.depth <-  depth[indsubset[parpos], snpsubset]
       par2.depth <-  depth[indsubset[par2pos], snpsubset]
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
#       relcheck <- tcrossprod(matrix(offs.genon0,nrow=1),matrix(par.genon0,nrow=1)) / 
#                   sum((2*puse[snpsubset]*(1-puse[snpsubset]))[offs.usegeno & par.usegeno])
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
#     tempAssign[which(allgmatch[, paste0(partype, "rel")] < rel.thresh | is.na(allgmatch[, paste0(partype, "rel")]))] <- "N"  # replace with next 2 for NA
     tempAssign[which(allgmatch[, paste0(partype, "rel")] < rel.thresh )] <- "N" 
     tempAssign[is.na(allgmatch[, paste0(partype, "rel")])] <- NA
      # for E assigns, check if the 2nd parent is possible.  
     tempAssign[which(allgmatch[, paste0(partype, "rel2nd")] >= rel.thresh & EMMrate2 < emm.thresh & tempAssign == "E")] <- "A"
     allgmatch[, paste0(partype, "Assign")] <- tempAssign
     cat("\nSummary of",partype,"Assignments\n")
     print(addmargins(table(allgmatch[, paste0(partype, "Assign")],useNA="ifany")))
     png(paste0("Best", partype, "Matches.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
      plot(allgmatch[, paste0("mmrate", partype)] ~ allgmatch[, paste0(partype, "rel")], main = paste("Best", partype, "Matches"), xlab = "Estimated Relatedness", 
          ylab = "Raw mismatch rate",col=fcolo[match(allgmatch$seqID,seqID)], cex=0.8)
      abline(v=rel.thresh, col="grey")
      dev.off()
     parEplot(partype,EMMrate,allgmatch[, paste0(partype, "rel")],plotcol=fcolo[match(allgmatch$seqID,seqID)],relthresh=rel.thresh)
     tempch <- assign.pch[match(tempAssign,assign.rank)]
     png(paste0("ExpMM-", partype, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
      plot(allgmatch[, paste0("mmrate", partype)] ~ allgmatch[, paste0("exp.mmrate", partype)] , main = paste("Best", partype, "Matches"), xlab = "Expected mismatch rate", 
          ylab = "Raw mismatch rate",col=fcolo[match(allgmatch$seqID,seqID)], cex=0.8, pch=tempch)
      pch.used <- sort(match(unique(tempAssign),assign.rank))
      legend("bottomright",title="Assign",cex=0.75,pch=assign.pch[pch.used],legend=assign.rank[pch.used])
      abline(a=0,b=1,col="red")
      abline(a=emm.thresh,b=1,col="grey")
      edges <- par("usr") # xl,xr,yb,yt
      poly1 <- data.frame(x1 = c(edges[1], edges[1], edges[2],edges[2],edges[1]), 
                          y1 = c(emm.thresh+edges[1],edges[4],edges[4],emm.thresh+edges[2],emm.thresh+edges[1]))
      polygon(poly1,col = rgb(0,0,0,alpha=0.1),border=NA) 
      dev.off()
     mmpalette <- colorRampPalette(c("blue","red"))(50)
     mmcol <- mmpalette[trunc(1+50*(EMMrate-min(EMMrate,na.rm=TRUE))/(diff(range(EMMrate,na.rm=TRUE))+1E-6))]
     legend_image <- as.raster(matrix(rev(mmpalette), ncol = 1))
     xyrange <- range(c(allgmatch[, paste0(partype, "rel")],allgmatch[, paste0(partype,"rel2nd")]),na.rm=TRUE)
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
     noffspringpar <- data.frame(table(allgmatch[, paste0("Best", partype, "Match")]))
     colnames(noffspringpar)[2] <- paste0(partype, "Freq")
     groupsinfo <- merge(groupsinfo, noffspringpar, by.x = "IndivID", by.y = "Var1", all = TRUE)
     groupsinfo[is.na(groupsinfo[, paste0(partype, "Freq")]), paste0(partype, "Freq")] <- 0
     uo <- match(allgmatch$seqID,seqID[indsubset])
     allgmatch$Inb <- diag(Guse)[uo] - 1
     Parposped <- match(allgmatch[,paste0("Best", partype, "Match")],pedinfo$IndivID)
     Parposg <- match(pedinfo$seqID[Parposped],seqID[indsubset])
     allgmatch[,paste0(partype,"Inb")] <- diag(Guse)[Parposg] - 1
     if(developer) {
      Exprel <- 0.5 + allgmatch$Inb+allgmatch[,paste0(partype,"Inb")]/2   # par-offspr rel from theory
      reldevn <- allgmatch[,paste0(partype,"rel")] - Exprel
      uY <- which(tempAssign == "Y")
      if(length(uY) > 0) {
       lmrel <- lm(allgmatch[uY,paste0(partype,"rel")] ~ allgmatch$Inb[uY] + allgmatch[uY,paste0(partype,"Inb")] )
       print(summary(lmrel))
       }
      reldevplots(partype,EMMrate,reldevn,gmatch=allgmatch,Exprelvar=Exprel,plotcol=fcolo[match(allgmatch$seqID,seqID)],plotch=tempch) 
      }
     }
    write.csv(allgmatch, paste0(partype, "Matches.csv"), row.names = FALSE)
    list(groupmatch=allgmatch, groupsinfo=groupsinfo)
  } else {
    NULL
  }
}

matesmatch <- function(Guse, pedinfo, matesinfo) {
  groupIDs <- na.omit(unique(pedinfo$MatesGroup))
  groupIDs <- groupIDs[!groupIDs == ""]
  ngroups <- length(groupIDs)
  if (ngroups > 0) {
    for (g in 1:ngroups) {
      group <- groupIDs[g]
      offspringID <- pedinfo$IndivID[which(pedinfo$MatesGroup == group)]
      GroupID <- matesinfo[which(matesinfo$MatesGroup == group),c("MaleID","FemaleID")]
      offspringseqID <- with(pedinfo, seqID[match(offspringID, IndivID)])
      GroupseqID <- with(pedinfo, cbind(seqID[match(GroupID[,1], IndivID)],seqID[match(GroupID[,2], IndivID)]))
      offspringpos <- match(offspringseqID, seqID[indsubset])
      matespos <- cbind(match(GroupseqID[,1], seqID[indsubset]),match(GroupseqID[,2], seqID[indsubset]))
      gmatch <- bestmatesmatch(offspringpos, matespos, Guse, matchcriterion = matchmethod,pedinfo=pedinfo)
      if (g == 1) allgmatch <- gmatch else allgmatch <- rbind(allgmatch, gmatch)
      }
    allgmatch$IndivID <- pedinfo$IndivID[match(allgmatch$seqID, pedinfo$seqID)]
    ncolallg <- ncol(allgmatch)
    allgmatch <- allgmatch[, c(ncolallg, 1:(ncolallg-1))]
    if (nrow(allgmatch) > 0) {
     allgmatch$BestFatherMatch <- pedinfo$IndivID[match(allgmatch$BestFatherMatch, pedinfo$seqID)]
     allgmatch$BestMotherMatch <- pedinfo$IndivID[match(allgmatch$BestMotherMatch, pedinfo$seqID)]
     allgmatch$FatherMatch2nd <- pedinfo$IndivID[match(allgmatch$FatherMatch2nd, pedinfo$seqID)]
     allgmatch$MotherMatch2nd <- pedinfo$IndivID[match(allgmatch$MotherMatch2nd, pedinfo$seqID)]
     mmstats <- mismatch.par(allgmatch$IndivID, allgmatch$BestFatherMatch, pedinfo)
     allgmatch$mmrateFather <- mmstats$mmrate
     allgmatch$mmnumFather <- mmstats$ncompare
     allgmatch$exp.mmrateFather <- mmstats$exp.mmrate
     EMMrateFather <- allgmatch$mmrateFather  - allgmatch$exp.mmrateFather 
     #not doing 2nds - its the combo that is important here (so no single parent A's either)
     mmstats <- mismatch.par(allgmatch$IndivID, allgmatch$BestMotherMatch, pedinfo)
     allgmatch$mmrateMother <- mmstats$mmrate
     allgmatch$mmnumMother <- mmstats$ncompare
     allgmatch$exp.mmrateMother <- mmstats$exp.mmrate
     EMMrateMother <- allgmatch$mmrateMother  - allgmatch$exp.mmrateMother 
      # closerel / boostrapping implemented for mating pairs analysis
     tempAssign <- rep("Y",nrow(allgmatch))
     tempAssign[which(pmax(EMMrateFather,EMMrateMother) > emm.thresh)] <- "E"
     tempAssign[which(allgmatch$Fatherrel < rel.threshF)] <- "N"
     tempAssign[which(allgmatch$Motherrel < rel.threshM)] <- "N"
     allgmatch$BothAssign <- tempAssign 
     tempch <- assign.pch[match(tempAssign,assign.rank)]
      #havent done joint EMM yet ...
     parEplot("Father",EMMrateFather,allgmatch$Fatherrel,plotcol=fcolo[match(allgmatch$seqID,seqID)],relthresh=rel.threshF)
     parEplot("Mother",EMMrateMother,allgmatch$Motherrel,plotcol=fcolo[match(allgmatch$seqID,seqID)],relthresh=rel.threshM)
     noffspringmates <- data.frame(table(with(allgmatch,paste(BestFatherMatch,BestMotherMatch,sep=":"))))
     colnames(noffspringmates)[2] <- "MatesPairFreq"
     matesinfo$pairID <- with(matesinfo,paste(MaleID,FemaleID,sep=":"))
     matesinfo <- merge(matesinfo, noffspringmates, by.x = "pairID", by.y = "Var1", all = TRUE)
     matesinfo$MatesPairFreq[is.na(matesinfo$MatesPairFreq)] <- 0
     matesinfo$pairID <- NULL
     uo <- match(allgmatch$seqID,seqID[indsubset])
     allgmatch$Inb <- diag(Guse)[uo] - 1
     Parposped <- match(allgmatch$BestFatherMatch,pedinfo$IndivID)
     Parposg <- match(pedinfo$seqID[Parposped],seqID[indsubset])
     allgmatch$FatherInb <- diag(Guse)[Parposg] - 1
     Parposped <- match(allgmatch$BestMotherMatch,pedinfo$IndivID)
     Parposg <- match(pedinfo$seqID[Parposped],seqID[indsubset])
     allgmatch$MotherInb <- diag(Guse)[Parposg] - 1
     if(developer) {
      uY <- which(tempAssign == "Y")
      Exprel <- 0.5 + allgmatch$Inb+allgmatch$FatherInb/2   # par-offspr rel from theory
      reldevn <- allgmatch$Fatherrel - Exprel
      if(length(uY) > 0) {
       lmrel <- lm(allgmatch$Fatherrel[uY] ~ allgmatch$Inb[uY] + allgmatch$FatherInb[uY] )
       print(summary(lmrel))
       }
      reldevplots("Father",EMMrateFather,reldevn,gmatch=allgmatch,Exprelvar=Exprel,plotcol=fcolo[match(allgmatch$seqID,seqID)],plotch=tempch) 
      Exprel <- 0.5 + allgmatch$Inb+allgmatch$MotherInb/2   # par-offspr rel from theory
      reldevn <- allgmatch$Motherrel - Exprel
      if(length(uY) > 0) {
       lmrel <- lm(allgmatch$Motherrel[uY] ~ allgmatch$Inb[uY] + allgmatch$MotherInb[uY] )
       print(summary(lmrel))
       }
      reldevplots("Mother",EMMrateMother,reldevn,gmatch=allgmatch,Exprelvar=Exprel,plotcol=fcolo[match(allgmatch$seqID,seqID)],plotch=tempch) 
      }
     }
    write.csv(allgmatch, "MatePairMatches.csv", row.names = FALSE)
    list(matesmatch=allgmatch, matesinfo=matesinfo)
  } else {
    NULL
  }
}

ssbbmm <- function(bbpar,uuse, pedinfo, BothMatches, quiet=FALSE) {
 depth2K <<- depth2Kchoose (dmodel="bb", bbpar)
 mmstatsbb <- mismatch.2par(BothMatches$IndivID,BothMatches$BestFatherMatch, BothMatches$BestMotherMatch,pedinfo=pedinfo)
 mmssbb <- sum((mmstatsbb$mmrate-mmstatsbb$exp.mmrate)[uuse]^2)
 if(!quiet) cat("bb param = ",bbpar,"ss = ",mmssbb,"\n")
 mmssbb
 }

ssmpmm <- function(mppar,uuse, pedinfo, BothMatches, quiet=FALSE) {
 depth2K <<- depth2Kchoose (dmodel="modp", mppar)
 mmstatsmp <- mismatch.2par(BothMatches$IndivID,BothMatches$BestFatherMatch, BothMatches$BestMotherMatch,pedinfo=pedinfo)
 mmssmp <- sum((mmstatsmp$mmrate-mmstatsmp$exp.mmrate)[uuse]^2)
 if(!quiet) cat("mp param = ",mppar,"ss = ",mmssmp,"\n")
 mmssmp
 }

addtagIDs <- function(sampinfo,indvar,tagvar,matchtype="both", pedresults) {
  matchtype <- tolower(matchtype)
#  if(matchtype=="both") pedresults <- BothMatches
#  if(matchtype=="father") pedresults <- FatherMatches
#  if(matchtype=="mother") pedresults <- MotherMatches
  progpos <- match(pedresults$IndivID,sampinfo[,indvar])
  pedresults$IndivTag <- sampinfo[progpos,tagvar]
  if(matchtype=="both" | matchtype=="father") {
    fpos <-  match(pedresults$BestFatherMatch,sampinfo[,indvar])
    pedresults$FatherTag <- sampinfo[fpos,tagvar]
    }
  if(matchtype=="both" | matchtype=="mother") {
    mpos <-  match(pedresults$BestMotherMatch,sampinfo[,indvar])
    pedresults$MotherTag <- sampinfo[mpos,tagvar]
  }
  pedresults
 }

bestparPCA <- function(Gobj,sfx="",keypos=NULL, pedinfo, BothMatches) {
 plotch <- assign.pch[match(BothMatches$BothAssign,assign.rank)]
 ukeep <- which(seqID[Gobj$indsubset] %in% pedinfo$seqID)
 uf <- match(seqID[Gobj$indsubset][ukeep],BothMatches$seqID)
 uo <- match(BothMatches$seqID,seqID[Gobj$indsubset][ukeep])
 nprog=length(uo)
 pchuse <- plotch[uf]
 pchuse[is.na(pchuse)] <- 16
 png(paste0("PC-BestParents",sfx,".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
  with(Gobj$PC,  plot(x[ukeep, 2] ~ x[ukeep, 1], cex = 1, col = fcolo[Gobj$indsubset][ukeep], pch=pchuse, xlab = "Principal component 1", ylab = "Principal component 2") )
  ParentLines <- data.frame(x=rep(NA,3*nprog),y=rep(NA,3*nprog))
  ParentLines[seq(1,(3*nprog-2),3),] <- Gobj$PC$x[ukeep[uo],1:2]
  ParentLines[seq(2,(3*nprog-1),3),] <- Gobj$PC$x[match(pedinfo$seqID[match(BothMatches$BestFatherMatch,pedinfo$IndivID)],seqID[Gobj$indsubset]),1:2]
  lines(ParentLines,col="blue")
  ParentLines[seq(2,(3*nprog-1),3),] <- Gobj$PC$x[match(pedinfo$seqID[match(BothMatches$BestMotherMatch,pedinfo$IndivID)],seqID[Gobj$indsubset]),1:2]
  lines(ParentLines,col="deeppink")
  pch.used <- sort(match(unique(BothMatches$BothAssign),assign.rank))
  if(!is.null(keypos)) legend(keypos,title="Assign",cex=0.75,pch=assign.pch[pch.used],legend=assign.rank[pch.used])
  dev.off()
 invisible(NULL)
 }

GBSPed <- function () {
pedsetup()
outobj <- list()
if (OK4ped & exists("pedfile") & exists("GCheck")) {
  pedinfo <- suppressWarnings(read.csv(pedfile, stringsAsFactors = FALSE, colClasses=c(FatherGroup="character", MotherGroup="character")))
  pedinfo <- pedinfo[!is.na(pedinfo$seqID), ]
  pedinfo <- pedinfo[!is.na(match(pedinfo$seqID, seqID[indsubset])), ]
  uo <- match(pedinfo$seqID,seqID)
  pedinfo$Inb <- diag(eval(parse(text = GCheck)))[uo] - 1
  dupids <- which(duplicated(pedinfo$IndivID))
  if(length(dupids)>0) cat("Warning: dupicates in IndivID",pedinfo$IndivID[dupids],"\n")
  dupids <- which(duplicated(pedinfo$seqID))
  if(length(dupids)>0) cat("Warning: dupicates in seqID",pedinfo$seqID[dupids],"\n")
  if ("FatherID" %in% colnames(pedinfo)) 
    pedinfo <- parmatch("Father", eval(parse(text = GCheck)), pedinfo=pedinfo)
  if ("MotherID" %in% colnames(pedinfo)) 
    pedinfo <- parmatch("Mother", eval(parse(text = GCheck)), pedinfo=pedinfo)
  if ("FatherID" %in% colnames(pedinfo) & "MotherID" %in% colnames(pedinfo)) {
    Par2mm <- mismatch.2par(pedinfo$IndivID,pedinfo$FatherID,pedinfo$MotherID, pedinfo=pedinfo)
    pedinfo$FandMEMM <- with(Par2mm,mmrate-exp.mmrate)
    pedinfo$FandMmatch <- with(pedinfo,FatherMatch & MotherMatch & FandMEMM < emm.thresh2)
    if (is.character(pedinfo$IndivID)) {
      umiss <- which(pedinfo$FatherID == "")
      if (length(umiss) > 0) pedinfo$FatherID[umiss] <- NA
      umiss <- which(pedinfo$MotherID == "")
      if (length(umiss) > 0) pedinfo$MotherID[umiss] <- NA
    }
    tempch <- verif.ch[match(pedinfo$FandMmatch, c(NA,TRUE,FALSE))]
    png(paste0("ExpMM-RecBoth.png"), width = 640, height = 640, pointsize = cex.pointsize *  18)
     plot(Par2mm$mmrate ~ Par2mm$exp.mmrate, main = paste("Rec Parent Matches"), xlab = "Expected mismatch rate", 
         ylab = "Raw mismatch rate",col=fcolo[match(pedinfo$seqID,seqID)], pch=tempch, cex=0.8)
     abline(a=0,b=1,col="red")
     abline(a=emm.thresh2,b=1,col="grey")
#     legend("bottomright",title="Assign",cex=0.75,pch=assign.pch[pch.used],legend=assign.rank[pch.used])
     edges <- par("usr") # xl,xr,yb,yt
     poly1 <- data.frame(x1 = c(edges[1], edges[1], edges[2],edges[2],edges[1]), 
                         y1 = c(emm.thresh2+edges[1],edges[4],edges[4],emm.thresh2+edges[2],emm.thresh2+edges[1]))
     polygon(poly1,col = rgb(0,0,0,alpha=0.1),border=NA) 
     dev.off()

    famnumber <- rep(NA, nrow(pedinfo))
    if(is.character(pedinfo$FatherID)) pedinfo$FatherID[pedinfo$FatherID==""] <- NA
    if(is.character(pedinfo$MotherID)) pedinfo$MotherID[pedinfo$MotherID==""] <- NA
    famtable <- with(pedinfo, table(FatherID, MotherID))
    fampos <- which(famtable > 1, arr.ind = TRUE)
    if(nrow(fampos)>0) {
     famfathers <- dimnames(famtable)$FatherID[fampos[, 1]]
     fammothers <- dimnames(famtable)$MotherID[fampos[, 2]]
     if (is.numeric(pedinfo$FatherID)) 
       famfathers <- as.numeric(famfathers)
     if (is.numeric(pedinfo$FatherID)) 
       famfathers <- as.numeric(famfathers)
     noffspring <- famtable[fampos]
     nfamilies <- length(noffspring)
     famresults <- rep(NA, nfamilies)
     for (ifam in 1:nfamilies) {
       famnumber[which(pedinfo$FatherID == famfathers[ifam] & pedinfo$MotherID == fammothers[ifam])] <- ifam
       uoffspring <- match(pedinfo$seqID[which(famnumber == ifam)], seqID[indsubset])
       famresults[ifam] <- mean(eval(parse(text = GCheck))[uoffspring, uoffspring][upper.tri(diag(nrow = length(uoffspring)))])
       }
     cat("Mean relatedness for full-sib families (as given)\n")
     print(data.frame(famfathers, fammothers, noffspring, meanrel = famresults))
     cat("Mean relatedness within all full-sib families", weighted.mean(famresults, noffspring), "\n")
     }
    
    uoffspring <- which(!is.na(famnumber))
    Fatherset <- unique(na.omit(pedinfo$FatherID))
    Motherset <- unique(na.omit(pedinfo$MotherID))
    udiff.fathers <- which(!match(pedinfo$FatherID[uoffspring], Fatherset) %*% t(rep(1, length(uoffspring))) == rep(1, length(uoffspring)) %*% 
                             t(match(pedinfo$FatherID[uoffspring], Fatherset)), arr.ind = T)
    udiff.mothers <- which(!match(pedinfo$MotherID[uoffspring], Motherset) %*% t(rep(1, length(uoffspring))) == rep(1, length(uoffspring)) %*% 
                             t(match(pedinfo$MotherID[uoffspring], Motherset)), arr.ind = T)
    udiff <- as.matrix(merge(udiff.fathers, udiff.mothers))
    opos <- match(pedinfo$seqID[uoffspring], seqID[indsubset])
    cat("Mean relatedness between individuals in full-sib families with different parents", mean(eval(parse(text = GCheck))[opos,opos][udiff]), "\n")
  }
  outobj$pedinfo <- pedinfo
  write.csv(pedinfo, "PedVerify.csv", row.names = FALSE)

  if (exists("matesfile")) if(!file.exists(matesfile)) {
   cat("Warning: Mates file", matesfile, "not found\n")
   rm(matesfile)
   }
  if (exists("matesfile")) {  ######### check for matching mating pair ###########
    suppressWarnings(rm(groupsfile)) # dont use groups
    matesinfo <- read.csv(matesfile, stringsAsFactors = FALSE, colClasses=(MatesGroup="character"))
    matesinfo <- matesinfo[!duplicated(matesinfo),]  # remove row duplicates
    matesinfo$Genotyped <- ifelse(is.na(match(pedinfo$seqID[match(matesinfo$MaleID,pedinfo$IndivID)],seqID))
                                | is.na(match(pedinfo$seqID[match(matesinfo$FemaleID,pedinfo$IndivID)],seqID)),"N","Y")
    if ("MatesGroup" %in% colnames(pedinfo)) {
     mresults <- matesmatch(eval(parse(text = GCheck)),pedinfo=pedinfo, matesinfo=matesinfo)
     if(is.null(mresults)) BothMatches <- NULL else {
      matesinfo <- mresults$matesinfo
      BothMatches <- mresults$matesmatch
      }
     }
    if(nrow(BothMatches) > 0) {
     BothMatches <- BothMatches[order(BothMatches$IndivID),,drop=FALSE]
     mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$BestFatherMatch, BothMatches$BestMotherMatch,pedinfo=pedinfo)
     BothMatches$mmrateF1M1 <- mmstats$mmrate
     BothMatches$mmnumF1M1 <- mmstats$ncompare
     BothMatches$exp.mmrateF1M1 <- mmstats$exp.mmrate
     mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$FatherMatch2nd, BothMatches$MotherMatch2nd,pedinfo=pedinfo)
     BothMatches$mmrateF2M2 <- mmstats$mmrate
     BothMatches$mmnumF2M2 <- mmstats$ncompare
     BothMatches$exp.mmrateF2M2 <- mmstats$exp.mmrate
     uf <- match(pedinfo$seqID[match(BothMatches$BestFather,pedinfo$IndivID)],seqID[indsubset])
     um <- match(pedinfo$seqID[match(BothMatches$BestMother,pedinfo$IndivID)],seqID[indsubset])
     BothMatches$relF1M1 <- eval(parse(text = GCheck))[cbind(uf,um)]
     uf <- match(pedinfo$seqID[match(BothMatches$FatherMatch2nd,pedinfo$IndivID)],seqID[indsubset])
     um <- match(pedinfo$seqID[match(BothMatches$MotherMatch2nd,pedinfo$IndivID)],seqID[indsubset])
     BothMatches$relF2M2 <- eval(parse(text = GCheck))[cbind(uf,um)]
     uo <- match(BothMatches$seqID,seqID[indsubset])
     BothMatches$Inb <- diag(eval(parse(text = GCheck)))[uo] - 1
     tempInb <- BothMatches$Inb; tempInb[is.na(tempInb)] <- 100 # arbitrary high so not a fail in Inb tests
     EMMrates <- with(BothMatches,cbind(mmrateF2M2-exp.mmrateF2M2,mmrateF1M1-exp.mmrateF1M1))
     EMMrate.min <- apply(EMMrates, MARGIN=1, min)
     if (is.null(minr4inb)) minr4inb <<- min(BothMatches$relF1M1) - 0.001
     BothMatches$BothAssign[which(EMMrates[,2] > emm.thresh2 & BothMatches$BothAssign %in% assign.rank[1:4])] <- "E"
     BothMatches$BothAssign[EMMrates[,2]-EMMrate.min > emmdiff.thresh2 & EMMrate.min < emm.thresh2 & BothMatches$BothAssign %in% assign.rank[1:5]] <- "A"
     BothMatches$BothAssign[which(BothMatches$relF1M1 - 2 * tempInb > inb.thresh & BothMatches$relF1M1 > minr4inb & BothMatches$BothAssign == "Y")] <- "I"
     BothMatches$Alternate <- ""
     Apos <- which(BothMatches$BothAssign %in% c("A","I"))
     for (ipos in Apos) {
      altOK <- TRUE
      if (BothMatches$relF2M2[ipos] - 2 * tempInb[ipos] > inb.thresh | EMMrate.min[ipos] > emm.thresh2) altOK <- FALSE
      if(BothMatches$Fatherrel2nd[ipos] < rel.threshF) altOK <- FALSE
      if(BothMatches$Motherrel2nd[ipos] < rel.threshM) altOK <- FALSE
      if(altOK) BothMatches$Alternate[ipos] <- "F2M2"
      }
     outobj$BothMatches <- BothMatches
     write.csv(BothMatches,"MatePairMatches.csv", row.names = FALSE)
     uo <- match(BothMatches$seqID,seqID)
     trioplots(BothMatches=BothMatches)
     plotch <- assign.pch[match(BothMatches$BothAssign,assign.rank)]
     png("MMrateBothE.png", width = 640, height = 640, pointsize = cex.pointsize *  15)
      plot(EMMrates[,1] ~ EMMrates[,2], main="Excess Mismatch Rates", ylab="Father2, Mother2",xlab="Father1, Mother1", col=fcolo[uo],pch=plotch)
      abline(a = 0,b = 1, col="red")
      dev.off()
     cat("\nSummary of joint Assignments\n")
     print( addmargins(table(BothMatches$BothAssign, useNA="ifany")) )
     }
    }
  if (exists("groupsfile")) if(!file.exists(groupsfile)) {
   cat("Warning: Groups file", groupsfile, "not found\n")
   rm(groupsfile)
   }
  if (exists("groupsfile")) {  ######### find fathers and mothers from possibles ###########
    groupsinfo <- read.csv(groupsfile, stringsAsFactors = FALSE, colClasses=(ParGroup="character"))
    groupsinfo <- groupsinfo[!duplicated(groupsinfo),]  # remove row duplicates
    groupsinfo$Genotyped <- ifelse(is.na(match(pedinfo$seqID[match(groupsinfo$IndivID,pedinfo$IndivID)],seqID)),"N","Y")
    if ("FatherGroup" %in% colnames(pedinfo)) {
     mresults <- groupmatch(eval(parse(text = GCheck)), "Father",pedinfo=pedinfo, groupsinfo=groupsinfo)
     if(is.null(mresults)) FatherMatches <- NULL else {
      groupsinfo <- mresults$groupsinfo
      FatherMatches <- mresults$groupmatch
      }
     outobj$FatherMatches <- FatherMatches
     }
    if ("MotherGroup" %in% colnames(pedinfo)) {
     mresults <- groupmatch(eval(parse(text = GCheck)), "Mother",pedinfo=pedinfo, groupsinfo=groupsinfo)
     if(is.null(mresults)) MotherMatches <- NULL else {
      groupsinfo <- mresults$groupsinfo
      MotherMatches <- mresults$groupmatch
      }
     outobj$MotherMatches <- MotherMatches
     }
    if ("FatherGroup" %in% colnames(pedinfo) & "MotherGroup" %in% colnames(pedinfo)) {
     BothMatches <- merge(FatherMatches,MotherMatches,all=TRUE)
     if(nrow(BothMatches) > 0) {
      BothMatches <- BothMatches[order(BothMatches$IndivID),,drop=FALSE]
      if(!allow.selfing) {
       uselfed <- with(BothMatches,which(BestFatherMatch==BestMotherMatch))
       if(matchmethod=="EMM") { #whichbetter: 1 retain F1, 2 retain M1
        whichbetter <- apply(with(BothMatches[uselfed,],cbind(mmrateMother2-exp.mmrateMother2,mmrateFather2-exp.mmrateFather2)),1,which.min)
        } else { #rel method
        whichbetter <- apply(with(BothMatches[uselfed,],cbind(Motherrel2nd,Fatherrel2nd)),1,which.max)
        }
       BothMatches$tempna <- NA
       uchange <- which(whichbetter==2)
       if(length(uchange)>0) { # promote 2nd father to 1st
        new0 <- BothMatches[uselfed[uchange],c("FatherMatch2nd","tempna","Fatherrel2nd","tempna","tempna","mmrateFather2",
                            "tempna","exp.mmrateFather2","tempna","tempna","tempna","tempna","tempna","tempna")]
        changecols <- c("BestFatherMatch","FatherMatch2nd","Fatherrel","Fatherrel2nd","Father12rel","mmrateFather",
                            "mmnumFather","exp.mmrateFather","mmrateFather2","exp.mmrateFather2","Fathersd","FatherReliability","FatherAssign","FatherInb")
        colnames(new0) <- changecols
        new0$FatherAssign <- "Y"
        new0$FatherAssign[which(new0$mmrateFather - new0$exp.mmrateFather  > emm.thresh)] <- "E"
        new0$FatherAssign[which(new0$Fatherrel  < rel.threshF)] <- "N"
        new0$FatherAssign[which(is.na(new0$Fatherrel))] <- "N"
        new0$FatherInb <- diag(eval(parse(text = GCheck)))[match(pedinfo$seqID[match(new0$BestFatherMatch, pedinfo$IndivID)], seqID[indsubset])] - 1
        BothMatches[uselfed[uchange],changecols] <- new0
        # if 2nd mother is the promoted father, remove that and change B results to Y
        u2same <- intersect(which(BothMatches$BestFatherMatch ==  BothMatches$MotherMatch2nd), uselfed[uchange])
        if(length(u2same)>0) {
         BothMatches[u2same,c("MotherMatch2nd","Mother12rel","mmrateMother2","exp.mmrateMother2")] <- NA
         uB <- which(BothMatches$MotherAssign[u2same]=="B")
         if(length(uB)>0) {
          BothMatches$MotherAssign[u2same[uB]] <- "Y"  # restore
          BothMatches[u2same[uB],c("Mothersd","MotherReliability")] <- NA
          }
         }
        } 
       uchange <- which(whichbetter==1)
       if(length(uchange)>0) {
        new0 <- BothMatches[uselfed[uchange],c("MotherMatch2nd","tempna","Motherrel2nd","tempna","tempna","mmrateMother2",
                            "tempna","exp.mmrateMother2","tempna","tempna","tempna","tempna","tempna","tempna")]
        changecols <- c("BestMotherMatch","MotherMatch2nd","Motherrel","Motherrel2nd","Mother12rel","mmrateMother",
                            "mmnumMother","exp.mmrateMother","mmrateMother2","exp.mmrateMother2","Mothersd","MotherReliability","MotherAssign","MotherInb")
        colnames(new0) <- changecols
        new0$MotherAssign <- "Y"
        new0$MotherAssign[which(new0$mmrateMother - new0$exp.mmrateMother  > emm.thresh)] <- "E"
        new0$MotherAssign[which(new0$Motherrel  < rel.threshM)] <- "N"
        new0$MotherAssign[which(is.na(new0$Motherrel))] <- "N"
        new0$MotherInb <- diag(eval(parse(text = GCheck)))[match(pedinfo$seqID[match(new0$BestMotherMatch, pedinfo$IndivID)], seqID[indsubset])] - 1
        BothMatches[uselfed[uchange],changecols] <- new0
        # if 2nd father is the promoted mother, remove that and change B results to Y
        u2same <- intersect(which(BothMatches$BestMotherMatch ==  BothMatches$FatherMatch2nd), uselfed[uchange])
        if(length(u2same)>0) {
         BothMatches[u2same,c("FatherMatch2nd","Father12rel","mmrateFather2","exp.mmrateFather2")] <- NA
         uB <- which(BothMatches$FatherAssign[u2same]=="B")
         if(length(uB)>0) {
          BothMatches$FatherAssign[u2same[uB]] <- "Y"  # restore
          BothMatches[u2same[uB],c("Fathersd","FatherReliability")] <- NA
          }
         }
        } 
       BothMatches$tempna <- NULL
       }
      mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$BestFatherMatch, BothMatches$BestMotherMatch,pedinfo=pedinfo)
      BothMatches$mmrateF1M1 <- mmstats$mmrate
      BothMatches$mmnumF1M1 <- mmstats$ncompare
      BothMatches$exp.mmrateF1M1 <- mmstats$exp.mmrate
      mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$FatherMatch2nd, BothMatches$BestMotherMatch,pedinfo=pedinfo)
      BothMatches$mmrateF2M1 <- mmstats$mmrate
      BothMatches$mmnumF2M1 <- mmstats$ncompare
      BothMatches$exp.mmrateF2M1 <- mmstats$exp.mmrate
      mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$BestFatherMatch, BothMatches$MotherMatch2nd,pedinfo=pedinfo)
      BothMatches$mmrateF1M2 <- mmstats$mmrate
      BothMatches$mmnumF1M2 <- mmstats$ncompare
      BothMatches$exp.mmrateF1M2 <- mmstats$exp.mmrate
      mmstats <- mismatch.2par(BothMatches$IndivID, BothMatches$FatherMatch2nd, BothMatches$MotherMatch2nd,pedinfo=pedinfo)
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
      tempInb <- BothMatches$Inb; tempInb[is.na(tempInb)] <- 100 # arbitrary high so not a fail in Inb tests
      BothMatches$BothAssign <- assign.rank[pmax(match(BothMatches$FatherAssign,assign.rank),match(BothMatches$MotherAssign,assign.rank),na.rm=TRUE)]
      BothMatches$BothAssign[BothMatches$FatherAssign=="Y" & is.na(BothMatches$MotherAssign)] <- "F"
      BothMatches$BothAssign[BothMatches$MotherAssign=="Y" & is.na(BothMatches$FatherAssign)] <- "M"
      EMMrates <- with(BothMatches,cbind(mmrateF2M2-exp.mmrateF2M2,mmrateF1M2-exp.mmrateF1M2,mmrateF2M1-exp.mmrateF2M1,mmrateF1M1-exp.mmrateF1M1))
      EMMrate.min <- apply(EMMrates, MARGIN=1, min, na.rm=TRUE)
      if (is.null(minr4inb)) minr4inb <<- min(BothMatches$relF1M1) - 0.001
      BothMatches$BothAssign[which(EMMrates[,4] > emm.thresh2 & BothMatches$BothAssign %in% assign.rank[1:4])] <- "E"
      BothMatches$BothAssign[EMMrates[,4]-EMMrate.min > emmdiff.thresh2 & EMMrate.min < emm.thresh2 & BothMatches$BothAssign %in% assign.rank[1:5]] <- "A"
      BothMatches$BothAssign[which(BothMatches$relF1M1 - 2 * tempInb > inb.thresh & BothMatches$relF1M1 > minr4inb & BothMatches$BothAssign == "Y")] <- "I"
      BothMatches$BothAssign[which(BothMatches$FatherAssign == "Y" & BothMatches$BothAssign == "N")] <- "F"
      BothMatches$BothAssign[which(BothMatches$MotherAssign == "Y" & BothMatches$BothAssign == "N")] <- "M"
      BothMatches$BothAssign[which(BothMatches$FatherAssign == "Y" & BothMatches$BothAssign %in% c("E") & BothMatches$MotherAssign %in% c("E"))] <- "F"
      BothMatches$BothAssign[which(BothMatches$MotherAssign == "Y" & BothMatches$BothAssign %in% c("E") & BothMatches$FatherAssign %in% c("E"))] <- "M"

      BothMatches$Alternate <- ""
      Apos <- which(BothMatches$BothAssign %in% c("A","I"))
      for (ipos in Apos) {
       altpar <- c("F2M2","F1M2","F2M1")[which.min(EMMrates[ipos,1:3])]
       if(length(altpar) > 0) {
        altOK <- TRUE
        if (BothMatches[ipos,paste0("rel",altpar)] - 2 * tempInb[ipos] > inb.thresh | EMMrate.min[ipos] > emm.thresh2) altOK <- FALSE
        if (grepl("F2",altpar)) {
         if(BothMatches[ipos, "Fatherrel2nd"] < rel.threshF | BothMatches[ipos, "mmrateFather2"] - BothMatches[ipos, "exp.mmrateFather2"] > emm.thresh ) altOK <- FALSE
         }
        if (grepl("M2",altpar)) {
         if(BothMatches[ipos, "Motherrel2nd"] < rel.threshM | BothMatches[ipos, "mmrateMother2"] - BothMatches[ipos, "exp.mmrateMother2"] > emm.thresh ) altOK <- FALSE
         }
        if(altOK) BothMatches$Alternate[ipos] <- altpar
        }
       }
      outobj$BothMatches <- BothMatches
      write.csv(BothMatches,"BothMatches.csv", row.names = FALSE)
      uo <- match(BothMatches$seqID,seqID)
      trioplots(BothMatches=BothMatches)
      plotch <- assign.pch[match(BothMatches$BothAssign,assign.rank)]
      upairs <- which(!is.na(BothMatches$mmrateF1M1))
      if(length(upairs)>0 & any(!is.na(rowSums(EMMrates[upairs,,drop=FALSE])))) {
       png("MMrateBoth.png", width = 640, height = 640, pointsize = cex.pointsize *  15)
        pairs(with(BothMatches[upairs,,drop=FALSE],cbind(mmrateF2M2,mmrateF1M2,mmrateF2M1,mmrateF1M1)),upper.panel=panel.yeqx,lower.panel=NULL,
                  main="Raw Mismatch Rates", labels=c("Father2,\nMother2","Father1,\nMother2","Father2,\nMother1","Father1,\nMother1"),
                  col.points=fcolo[uo][upairs],pch=plotch[upairs])
        dev.off()
       png("MMrateBothE.png", width = 640, height = 640, pointsize = cex.pointsize *  15)
        pairs(EMMrates[upairs,,drop=FALSE], main="Excess Mismatch Rates", labels=c("Father2,\nMother2","Father1,\nMother2","Father2,\nMother1","Father1,\nMother1"),
                  upper.panel=panel.yeqx,lower.panel=NULL,col.points=fcolo[uo][upairs],pch=plotch[upairs])
        dev.off()
      }
      cat("\nSummary of joint Assignments\n")
      print( addmargins(table(BothMatches$BothAssign, useNA="ifany")) )
      }
     }
    write.csv(groupsinfo, "GroupsParentCounts.csv", row.names = FALSE)
  }
 }
 invisible(outobj)
} # GBSPed

if (!functions.only) PedResults <- GBSPed()


if (FALSE) { # temporary working code
  # 
 }

