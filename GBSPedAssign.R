# assume all in pedfile are in the genotype results. To do: remove those that are not

if (!exists("rel.thresh")) rel.thresh <- 0.4
if (!exists(" mindepth.mm")) mindepth.mm <- 5

mismatch.par <- function(offspring.id,par.id) { #ids as in the pedigree file, if only 1 parent compare with all offspring
 noffspring <- length(offspring.id)
 if (length(par.id)==1) par.id <- rep(par.id,noffspring)
 nmismatch <- ncompare <- rep(NA,noffspring)
 opos <- match(pedinfo$seqID[match(offspring.id,pedinfo$IndivID)],seqID)
 ppos <- match(pedinfo$seqID[match(par.id,      pedinfo$IndivID)],seqID)
 for (ioffspring in 1:noffspring) {
  depthi <- depth[opos[ioffspring],]
  depthj <- depth[ppos[ioffspring],]
  usnp <- which(depthi >= mindepth.mm & depthj>= mindepth.mm)
  pi <- genon[opos[ioffspring],usnp]/2
  pj <- genon[ppos[ioffspring],usnp]/2
  nmismatch[ioffspring] <- length(which(abs(pi-pj)==1))
  ncompare[ioffspring] <- length(usnp)    
 }
 mmrate <- nmismatch/ncompare
}

parmatch <- function(partype,Gmatrix) {
  if(missing(partype)) partype <- "Father"
  ParseqID <- with(pedinfo,seqID[match(pedinfo[,paste0(partype,"ID")],IndivID)])
  offspringpos <- match(pedinfo$seqID,seqID) # all in pedigree file considered as offspring
  parpos <- match(ParseqID, seqID)
  ParRel <- Gmatrix[cbind(offspringpos,parpos)]
  ParMatch <- (ParRel > rel.thresh)
  png(paste0(partype,"Verify.png"),width = 640, height = 640, pointsize = 15)
   pairs(cbind(ParRel,offspringpos,parpos),labels=c("Relatedness","Offspring order",paste(partype,"order")))
   dev.off()
  ncompare <- sum(!is.na(ParMatch))
  nmatch <- sum(ParMatch,na.rm=TRUE)
  matchperc <- 100*nmatch/ncompare
  cat(nmatch,"matches out of",ncompare,partype,"comparisons:",format(matchperc,digits=3),"%\n")
  if (nmatch>0) cat("Mean relatedness for",partype,"matches",format(mean(ParRel[which(ParMatch)]),digits=3),"\n")
  if (ncompare>nmatch) cat("Mean relatedness for",partype,"non-matches",format(mean(ParRel[which(!ParMatch)]),digits=3),"\n")
  matchinfo <- data.frame(ParRel,ParMatch)
  names(matchinfo) <- paste0(partype,c("Rel","Match"))
  cbind(pedinfo,matchinfo)
 }

bestmatch <- function(ospos,parpos,Guse,partype) {
 if (missing(partype)) partype <- "Par"
 parchk <- Guse[ospos,parpos]
 maxpos <- apply(parchk,1,which.max)
 parchktemp <- parchk
 parchktemp[cbind(1:nrow(parchk),maxpos)] <- -1
 maxpos.2 <- apply(parchktemp,1,which.max)
 rm(parchktemp)
 maxrel <- cbind(parchk[cbind(1:nrow(parchk),maxpos)],parchk[cbind(1:nrow(parchk),maxpos.2)])
 out.df <- data.frame(seqID[ospos], seqID[parpos[maxpos]], seqID[parpos[maxpos.2]],relatedness=maxrel[,1],rel2nd=maxrel[,2],stringsAsFactors =FALSE)
 names(out.df) <- c("seqID",paste0("Best",partype,"Match"),paste0(partype,"Match2nd"),paste0(partype,"rel"),paste0(partype,"rel2nd"))
 out.df
 }

groupmatch <- function(Guse,partype) {
  groupIDs <- unique(pedinfo[,paste0(partype,"Group")])
  groupIDs <- na.omit(groupIDs)
  groupIDs <- groupIDs[!groupIDs ==""]
  ngroups <- length(groupIDs)
  if(ngroups > 0) {
   for (g in 1:ngroups) {
    group <- groupIDs[g]
    offspringID <- pedinfo$IndivID[which(pedinfo[,paste0(partype,"Group")]==group)]
    ParGroupID <- groupsinfo$IndivID[which(groupsinfo$ParGroup==group)]
    offspringseqID <- with(pedinfo,seqID[match(offspringID,IndivID)])
    ParGroupseqID <- with(pedinfo,seqID[match(ParGroupID,IndivID)])
    offspringpos <- match(offspringseqID,seqID)
    parpos <- match(ParGroupseqID, seqID)
    gmatch <- bestmatch(offspringpos,parpos,Guse,partype) 
    if(g==1) allgmatch <- gmatch else allgmatch <- rbind(allgmatch,gmatch)
    }
   allgmatch$IndivID <- pedinfo$IndivID[match(allgmatch$seqID,pedinfo$seqID)]
   allgmatch <- allgmatch[,c(6,1:5)]
   allgmatch[,paste0("Best",partype,"Match")] <- pedinfo$IndivID[match(allgmatch[,paste0("Best",partype,"Match")],pedinfo$seqID)]
   allgmatch[,paste0(partype,"Match2nd")] <- pedinfo$IndivID[match(allgmatch[,paste0(partype,"Match2nd")],pedinfo$seqID)]
   allgmatch[,paste0("mmrate",partype)] <- mismatch.par(allgmatch$IndivID,allgmatch[,paste0("Best",partype,"Match")])
   png(paste0("Best",partype,"Matches.png"),width = 640, height = 640, pointsize = 15)
    plot(allgmatch[,paste0("mmrate",partype)] ~ allgmatch[,paste0(partype,"rel")],
         main=paste("Best",partype,"Matches"),xlab="Estimated Relatedness",ylab="Raw mismatch rate")
    dev.off()
   write.csv(allgmatch,paste0(partype,"Matches.csv"),row.names=FALSE)
   noffspringpar <- data.frame(table(allgmatch[,paste0("Best",partype,"Match")]))
   colnames(noffspringpar)[2] <- paste0(partype,"Freq")
   groupsinfo <<- merge(groupsinfo,noffspringpar,by.x="IndivID",by.y="Var1",all=TRUE)
   groupsinfo[is.na(groupsinfo[,paste0(partype,"Freq")]),paste0(partype,"Freq")] <- 0
   allgmatch
   } else {
   NULL
   }
}

if(exists("pedfile") & exists("GCheck")) {
 pedinfo <- read.csv(pedfile,stringsAsFactors=FALSE)
 pedinfo <- pedinfo[!is.na(pedinfo$seqID),]
 pedinfo <- pedinfo[!is.na(match(pedinfo$seqID,seqID)),]
 if("FatherID" %in% colnames(pedinfo)) pedinfo <- parmatch("Father",eval(parse(text=GCheck)))
 if("MotherID" %in% colnames(pedinfo)) pedinfo <- parmatch("Mother",eval(parse(text=GCheck)))
 if("FatherID" %in% colnames(pedinfo) & "MotherID" %in% colnames(pedinfo)) {
  if (is.character(pedinfo$IndivID)) {
   umiss <- which(pedinfo$FatherID=="")
   if(length(umiss)>0) pedinfo$FatherID[umiss] <- NA
   umiss <- which(pedinfo$MotherID=="")
   if(length(umiss)>0) pedinfo$MotherID[umiss] <- NA
   }
  famtable <-  with(pedinfo,table(FatherID,MotherID))
  fampos <- which(famtable>1,arr.ind=TRUE)
  famfathers <- dimnames(famtable)$FatherID[fampos[,1]]
  fammothers <- dimnames(famtable)$MotherID[fampos[,2]]
  if (is.numeric(pedinfo$FatherID)) famfathers <- as.numeric(famfathers)
  if (is.numeric(pedinfo$FatherID)) famfathers <- as.numeric(famfathers)
  noffspring <- famtable[fampos]
  nfamilies <- length(noffspring)
  famnumber <- rep(NA,nrow(pedinfo))
  famresults <- rep(NA,nfamilies)
  for (ifam in 1:nfamilies) {
   famnumber[which(pedinfo$FatherID==famfathers[ifam] & pedinfo$MotherID==fammothers[ifam])] <- ifam
   uoffspring <- match(pedinfo$seqID[which(famnumber==ifam)],seqID)
   famresults[ifam] <- mean(eval(parse(text=GCheck))[uoffspring,uoffspring][upper.tri(diag(nrow=length(uoffspring)))])
   }
  cat("Mean relatedness for full-sib families (as given)\n")
  print(data.frame(famfathers,fammothers,noffspring,meanrel=famresults))
  cat("Mean relatedness within all full-sib families",weighted.mean(famresults,noffspring),"\n")

  uoffspring <- which(!is.na(famnumber))
  Fatherset <- unique(na.omit(pedinfo$FatherID))
  Motherset <- unique(na.omit(pedinfo$MotherID))
  udiff.fathers <- which(!match(pedinfo$FatherID[uoffspring],Fatherset) %*% t(rep(1,length(uoffspring))) == rep(1,length(uoffspring))  %*% t(match(pedinfo$FatherID[uoffspring],Fatherset)),arr.ind=T)
  udiff.mothers <- which(!match(pedinfo$MotherID[uoffspring],Motherset) %*% t(rep(1,length(uoffspring))) == rep(1,length(uoffspring))  %*% t(match(pedinfo$MotherID[uoffspring],Motherset)),arr.ind=T)
  udiff <- as.matrix(merge(udiff.fathers,udiff.mothers))
  cat("Mean relatedness between individuals in full-sib families with different parents",mean(eval(parse(text=GCheck))[uoffspring,uoffspring][udiff]),"\n")
  }
 write.csv(pedinfo,"PedVerify.csv",row.names=FALSE)
 if(exists("groupsfile")) {
  groupsinfo <- read.csv(groupsfile,stringsAsFactors=FALSE)
  if("FatherGroup" %in% colnames(pedinfo)) FatherMatches <- groupmatch(eval(parse(text=GCheck)),"Father")
  if("MotherGroup" %in% colnames(pedinfo)) MotherMatches <- groupmatch(eval(parse(text=GCheck)),"Mother")
  write.csv(groupsinfo,"GroupsParentCounts.csv",row.names=FALSE)
  }
 }



