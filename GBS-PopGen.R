PopGenver <- "0.9.7"
cat("GBS-PopGen for KGD version:",PopGenver,"\n")

heterozygosity <- function(indsubsetgf=1:nind,snpsubsetgf=1:nsnps,maxiter=100,convtol=0.001){
 nsnpsgf <- length(snpsubsetgf)
 nindgf <- length(indsubsetgf)
 genongf <- genon[indsubsetgf,snpsubsetgf,drop=FALSE]
 depth.use <- depth.orig[indsubsetgf,snpsubsetgf,drop=FALSE]
 depth.use[depth.use==0] <- NA # dont want to average over obs with depth 0
 obsgfreq <- cbind(colSums(genongf==0,na.rm=TRUE),colSums(genongf==1,na.rm=TRUE),colSums(genongf==2,na.rm=TRUE))/colSums(!is.na(genongf))
 obshet <- mean(obsgfreq[,2],na.rm=TRUE)
 psub <- calcp(indsubset=indsubsetgf)[snpsubsetgf]
 mafsub <- pmin(psub,1-psub)
 ehettrue <- mean(2*mafsub*(1-mafsub),na.rm=TRUE)  # exp het if true genos observed
 ehetstar <- mean(2*mafsub*(1-mafsub)*(1-2*colMeans(depth2K(depth.use),na.rm=TRUE)),na.rm=TRUE)
 ohet2 <- mean(obsgfreq[,2],na.rm=TRUE)/mean(1-2*depth2K(depth.use),na.rm=TRUE)
 ohet <- rep(NA,nsnpsgf)
 for(isnp in 1:nsnpsgf) {
  genodepth <- depth.use[,isnp]
  pnew <- obsgfreq[isnp,]
  ng <- sum(!is.na(genodepth))
  convtest<- 1; itcount <- 0
  while(convtest>convtol & itcount<maxiter) {
   itcount <- itcount+1
   pcurrent <- pnew
   paanew <- sum(1/(1+depth2K(genodepth[which(genongf[,isnp]==2)])*pcurrent[2]/pcurrent[3]))/ng
   pbbnew <- sum(1/(1+depth2K(genodepth[which(genongf[,isnp]==0)])*pcurrent[2]/pcurrent[1]))/ng
   pnew <- c(pbbnew,  1-paanew-pbbnew, paanew)
   convtest <- sum(abs(pnew-pcurrent))
   if(is.na(convtest)) convtest <- 0
   }
  ohet[isnp] <- pnew[2]
  }
 data.frame(ohetstar=obshet, ehetstar=ehetstar, ohet=mean(ohet, na.rm=TRUE),ohet2=ohet2, ehet=ehettrue)
}

Fst.GBS0 <- function(snpsubset, indsubset, populations) {
 if (missing(snpsubset))   snpsubset <- 1:nsnps
 if (missing(indsubset))   indsubset <- 1:nind
 npops <- length(unique(populations[indsubset]))
 Fst.results <- numeric(length(snpsubset))
 snppopdepth <- rowsum(depth.orig[indsubset, snpsubset],populations[indsubset])
 usnp <- which(!apply(snppopdepth,MARGIN=2,min) == 0)
 snpsubset <- snpsubset[usnp]
 aX2 <-  rbind(round(genon[indsubset,snpsubset]/2),ceiling(genon[indsubset,snpsubset]/2))
 X2results <- apply(aX2,MARGIN=2,function(x,y) chisq.test(table(x,y))$statistic,y=rep(populations[indsubset],2))
# effnuma <- colSums(2*(1-depth2K(depthsub <- depth.orig[indsubset, snpsubset])))  # 2(1-K)
 effnuma <- colSums(2*(1-depth2K(depth.orig[indsubset, snpsubset])))  # 2(1-K)
 Fst.results[usnp] <- npops * X2results / (effnuma * (npops - 1))
 Fst.results
 }

 chisq.adj <- function(x,y) {
    # x has alleles in first 2* num ind rows and popn eff num in last npops rows, 1 col per SNP (but passed by column?)
   chiresult <- NA
   if(is.factor(y)) y <- factor(y)  # resets levels if some are unused
   npopstemp <-  length(unique(y))
   temptable <- table(x[1:(length(x)-npopstemp)],y)
   if(nrow(temptable) == 2) {
    temptable2 <- temptable * matrix(x[length(x)+((1-npopstemp):0)]/colSums(temptable),nrow=2,ncol=npopstemp,byrow=TRUE)
    temptable2[is.na(temptable2)] <- 0  # note: if pop has 0 reads will result in NA X2 value
    chiresult <- suppressWarnings(chisq.test(temptable2)$statistic)
    }
   chiresult
   }

Fst.GBS <- function(snpsubset, indsubset, populations, varadj=0, SNPtest=FALSE) {
  #use varadj=1 to get Fst as Weir p166, varadj=0 for usual Fst
  if (missing(snpsubset))   snpsubset <- 1:nsnps
  if (missing(indsubset))   indsubset <- 1:nind
  npops <- length(unique(populations[indsubset]))
  Fst.results <- NA*numeric(length(snpsubset))
  snppopdepth <- rowsum(depth.orig[indsubset, snpsubset,drop=FALSE],populations[indsubset])
  pgsub <- colMeans(genon[indsubset,snpsubset,drop=FALSE],na.rm=TRUE) /2
  usnp <- which(!apply(snppopdepth,MARGIN=2,min) == 0 & pgsub > 0 & pgsub < 1 )
  snpsubset <- snpsubset[usnp]
  aX2 <-  rbind(round(genon[indsubset,snpsubset,drop=FALSE]/2),ceiling(genon[indsubset,snpsubset,drop=FALSE]/2))
  effnuma1 <- 2*(1-depth2K(depth.orig[indsubset, snpsubset,drop=FALSE]))  # 2(1-K)
  snppopeffn <- rowsum(effnuma1,populations[indsubset])  # effective number reads popn x snps
  X2results <- apply(rbind(aX2,snppopeffn),MARGIN=2,chisq.adj,y=rep(populations[indsubset],2))
  effnuma <- colSums(2*(1-depth2K(depth.orig[indsubset, snpsubset,drop=FALSE])))  # 2(1-K)
  Fst.results[usnp] <- npops * X2results / (effnuma * (npops - varadj))
  pvalue <- pchisq(X2results,df=npops-1,lower.tail=FALSE) 
  cat("Fst Mean:",mean(Fst.results,na.rm=TRUE),"Median:",median(Fst.results,na.rm=TRUE),"p-value:",mean(pvalue,na.rm=TRUE),"\n")
  if(SNPtest) Fst.results <- list(Fst=Fst.results, pvalue=pvalue)
  Fst.results
}


Fst.GBS.pairwise <- function(snpsubset, indsubset, populations,sortlevels=TRUE, SNPtest=FALSE, ...) {
 if (missing(snpsubset))   snpsubset <- 1:nsnps
 if (missing(indsubset))   indsubset <- 1:nind
 popnames <- unique(populations[indsubset])
 if(sortlevels) popnames <- sort(popnames)
 npops <- length(popnames)
 Fst.results <- array(dim=c(npops,npops,length(snpsubset)))
 if(SNPtest) pvalue <- Fst.results
 Fst.means <- Fst.medians <- array(dim=c(npops,npops),dimnames=list(popnames,popnames))
  for (ipop in 1:(npops-1)) {
  indsubseti <- indsubset[which(populations[indsubset] == popnames[ipop])]
  for (jpop in (ipop+1):npops) {
   indsubsetj <- indsubset[which(populations[indsubset] == popnames[jpop])]
   Fsttemp <- Fst.GBS(snpsubset,c(indsubseti,indsubsetj),populations, SNPtest = SNPtest, ...)
   if(SNPtest) Fst.results[ipop,jpop,] <- Fsttemp$Fst else Fst.results[ipop,jpop,] <- Fsttemp
   if(SNPtest) pvalue[ipop,jpop,] <- Fsttemp$pvalue
   Fst.means[ipop,jpop] <- mean(Fst.results[ipop,jpop,],na.rm=TRUE)
   Fst.medians[ipop,jpop] <- median(Fst.results[ipop,jpop,],na.rm=TRUE)
   }
  }
 cat("Pairwise Fst Means\n"); print(Fst.means[-npops,-1])
 cat("Pairwise Fst Medians\n");print(Fst.medians[-npops,-1])
 if(SNPtest) 
 Fst.results <- list(Fst=Fst.results, pvalue=pvalue)
 }




popmaf <- function(snpsubset, indsubset, populations=NULL, subpopulations=NULL, indcol, colobj, minsamps=10, mafmin=0, sortlevels=TRUE, unif = FALSE) {
 #populations defined relative to full set (length nind)
 if (missing(snpsubset))   snpsubset <- 1:nsnps
 if (missing(indsubset))   indsubset <- 1:nind
 if (!missing(colobj)) {populations=colobj$collabels[match(colobj$sampcol,colobj$collist)]; indcol <- colobj$sampcol }
 if (missing(indcol)) indcol <- rep("black",nind)
 if(is.null(populations)) populations <- rep("",nind)
 if(is.null(subpopulations)) subpopulations <- rep("",nind)
 popnames <- unique(populations[indsubset])
 sublevs <- unique(subpopulations[indsubset])
 if(sortlevels) {popnames <- sort(popnames); sublevs <- sort(sublevs) }
 nsub <- length(sublevs)
 if(nsub > 1) histdensity=30 else histdensity=NULL
 anglespan <- 180*(1-1/nsub)
 maxfreq <- 0
 if(unif) {
  for (i in seq_along(popnames)) {
   thigroup <- popnames[i]
   for (j in 1:nsub){
    thissub <- sublevs[j]
    indgroup <- intersect(indsubset,which(populations==thigroup & subpopulations==thissub))
    if(length(indgroup) >= minsamps) {
     plev <- calcp(indsubset=indgroup,pmethod="G")
     mafgroup <- pmin(plev,1-plev)
     snpsgroup <- intersect(snpsubset,which(mafgroup >= mafmin))
     histinfo <-  suppressWarnings(mafplot(mafgroup[snpsgroup], doplot=FALSE ) )
     maxfreq <- max(maxfreq,histinfo$counts)
     }
    } 
   }
  }
 for (i in seq_along(popnames)) {
  thigroup <- popnames[i]
  for (j in 1:nsub){
   thissub <- sublevs[j]
   indgroup <- intersect(indsubset,which(populations==thigroup & subpopulations==thissub))
   if(length(indgroup) >= minsamps) {
    plev <- calcp(indsubset=indgroup,pmethod="G")
    mafgroup <- pmin(plev,1-plev)
    snpsgroup <- intersect(snpsubset,which(mafgroup >= mafmin))
    groupcol <- unique(indcol[indgroup]); if (length(groupcol) > 1) groupcol <- "black"
    if (!unif) mafplot(mafgroup[snpsgroup],plotname=paste0("MAF-",thigroup,thissub),barcol=groupcol,main=paste("MAF for",thigroup,thissub),
            density=histdensity, angle=anglespan*((j-1)/(max(2,nsub)-1) - 0.5) )  # angle is irrelevant when nsub=1
    if (unif) mafplot(mafgroup[snpsgroup],plotname=paste0("MAF-",thigroup,thissub),barcol=groupcol,main=paste("MAF for",thigroup,thissub),
            density=histdensity, angle=anglespan*((j-1)/(max(2,nsub)-1) - 0.5), ylim=c(0,maxfreq) )  # angle is irrelevant when nsub=1
    nnz <- sum(mafgroup[snpsubset]>0,na.rm=TRUE)
    ng.2 <- sum(mafgroup[snpsubset]>0.2,na.rm=TRUE)
    cat(thigroup,thissub,"n=",length(indgroup),"# SNPs with MAF>0",nnz,"# SNPs with MAF>0.2",ng.2,"Proportion",ng.2/nnz,"\n")
    }
   }
  }
 }


manhatplot <- function(value, chrom, pos, plotname, qdistn=qunif, keyrot=0, symsize=0.8, legendm = NULL, ...) {
 chromcol <- colourby(chrom)
 colkey(chromcol,"chrom",srt=keyrot)
 plotord <- order(chrom,pos)
 valuetext <- gsub("[^[:alnum:]]", ".", deparse(substitute(value)))
 symsize0 <- symsize
 if(length(symsize)==length(value)) symsize <- symsize[plotord]
 png(paste0(plotname,"-Manhat.png"),width=1200)
  plot(value[plotord], col=chromcol$sampcol[plotord],xlab="Position",ylab=valuetext,cex=symsize, xaxt="n")
  if(!is.null(legendm)) legendm()
  dev.off()
 png(paste0(plotname,"-QQ.png"))
  if(length(symsize)==length(value)) symsize <- symsize0[order(value)]
  qqplot(qdistn(ppoints(length(value)),...), y=value, xlab="Theoretical quantiles", ylab=paste(valuetext,"quantiles"), cex=symsize,
         sub="Line for mid 98% of values", col=chromcol$sampcol[order(value)])
  qqline(value,col=2, distribution = function(p) qdistn(p, ...),prob=c(0.01,0.99))
  if(!is.null(legendm)) legendm()
  dev.off()
 }

#
#For the beta distribution see dbeta.
#For the binomial (including Bernoulli) distribution see dbinom.
#For the Cauchy distribution see dcauchy.
#For the chi-squared distribution see dchisq.
#For the exponential distribution see dexp.
#For the F distribution see df.
#For the gamma distribution see dgamma.
#For the geometric distribution see dgeom. (This is also a special case of the negative binomial.)
#For the hypergeometric distribution see dhyper.
#For the log-normal distribution see dlnorm.
#For the multinomial distribution see dmultinom.
#For the negative binomial distribution see dnbinom.
#For the normal distribution see dnorm.
#For the Poisson distribution see dpois.
#For the Student's t distribution see dt.
#For the uniform distribution see dunif.
#For the Weibull distribution see dweibull.


# select & pair (from different chromosomes) snps
# based on T Bilton code
snpselection <- function(chromosome,position,nsnpperchrom=100,seltype="centre",randseed=NULL, snpsubset,chromuse) {
 #seltype is centre, even or random
 if (missing(snpsubset))   snpsubset <- 1:length(chromosome)
 if (missing(chromuse)) chromuse <- unique(chromosome)
 if(seltype=="center") seltype <- "centre"
 usnp <- intersect(snpsubset, which(chromosome %in% chromuse))
 chromlist <- unique(chromosome[usnp])
 chromnone <- setdiff(chromuse,chromlist)
 if(length(chromnone) > 0 ) cat("Warning: no SNPs requested on chromosome(s)", chromnone, "\n")
 if(seltype=="random") {
  set.seed(randseed)
  snpchoose <- function(x) {
   indx <- which(chromosome[usnp] == x)
   nsnp <- min(c(length(indx),nsnpperchrom)) ## if there are fewer than nsnpperchrom SNPs on the chromosome
   return(sample(indx,nsnp))
   }
  }
 if(seltype=="even") {
  snpchoose <- function(x) {
   indx <- which(chromosome[usnp] == x)
   nsnp <- min(c(length(indx),nsnpperchrom)) ## if there are fewer than nsnpperchrom SNPs on the chromosome
   snpsep <- length(indx)/nsnpperchrom
   temp <- sort(position[usnp][indx], index.return=TRUE)
   return(indx[temp$ix[round(seq(snpsep/2,by=snpsep,length.out=nsnpperchrom))]])
   }
  }
 if(seltype=="centre") {
  meanpos <- aggregate(position[usnp] ~ chromosome[usnp],FUN=mean)
  colnames(meanpos) <- c("chr","pos")
  chromlist <- unique(chromosome[usnp])
  snpchoose <- function(x) {
   indx <- which(chromosome[usnp] == x)
   nsnp <- min(c(length(indx),nsnpperchrom)) ## if there are fewer than nsnpperchrom SNPs on the chromosome
   temp <- sort(abs(position[usnp][indx] - meanpos$pos[which(meanpos$chr==x)]), index.return=TRUE)
   return(indx[sort(temp$ix[1:nsnp])])
   }
  }
 snplist <- sapply(chromlist, snpchoose, simplify = FALSE)
 nchrom <- length(snplist)
 pairs <- do.call("rbind", sapply(1:(nchrom-1), function(x) as.matrix(expand.grid(snplist[[x]], unlist(snplist[(x+1):nchrom]), KEEP.OUT.ATTRS = FALSE)) ))
 pairs[,1] <- usnp[pairs[,1]]
 pairs[,2] <- usnp[pairs[,2]]
 return(pairs)
 }

### T Bilton, select SNPs (for LD analysis) from a UR object (modified)
snpselectionUR <- function(URobj, nsnpperchrom=100, nchrom, ...){
  chromlist <- unique(URobj$.__enclos_env__$private$chrom)
  if (!missing(nchrom)) chromlist <- chromlist[1:nchrom]
  URpairs <- snpselection (chromosome=URobj$.__enclos_env__$private$chrom,position=URobj$.__enclos_env__$private$pos,nsnpperchrom=nsnpperchrom,chromuse=chromlist, ...)
  return(URpairs)
}

Nefromr2 <- function(r2auto,nLD, alpha=1, weighted=FALSE,minN=1) {
 #r2auto = set of r2 values across different autosomes
 #nLD = # individuals for r2 calcs
 #alpha = mutation parameter
 #beta=1 (2) for phase unknown (known)
 if(length(nLD) ==1) nLD <- rep(nLD,length(r2auto))
 uN <- which(nLD >= minN)
 r2auto <- r2auto[uN]
 nLD <- nLD[uN]
 wt <- rep(1,length(r2auto))
 if(weighted) wt <- nLD
 meanN <- mean(nLD)
 Neauto <- (1/weighted.mean(r2auto,wt,na.rm=TRUE) - 1) /2   
 beta <- 1; Neauto.adj.b1 <- (1/(weighted.mean(r2auto,wt,na.rm=TRUE) -1/(beta*meanN)) - alpha) /2 
 beta <- 2; Neauto.adj.b2 <- (1/(weighted.mean(r2auto,wt,na.rm=TRUE) -1/(beta*meanN)) - alpha) /2 
 Neauto.med <- (1/median(r2auto,na.rm=TRUE) - 1) /2
 beta <- 1; Neauto.med.adj.b1 <- (1/(median(r2auto,na.rm=TRUE) -1/(beta*meanN)) - alpha) /2 
 beta <- 2; Neauto.med.adj.b2 <- (1/(median(r2auto,na.rm=TRUE) -1/(beta*meanN)) - alpha) /2 
 data.frame(n=meanN,Neauto,Neauto.adj.b1,Neauto.adj.b2,Neauto.med,Neauto.med.adj.b1,Neauto.med.adj.b2)
 }
