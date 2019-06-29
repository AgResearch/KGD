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
 effnuma <- colSums(2*(1-depth2K(depthsub <- depth.orig[indsubset, snpsubset])))  # 2(1-K)
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

Fst.GBS <- function(snpsubset, indsubset, populations, varadj=0) {
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
  cat("Fst Mean:",mean(Fst.results,na.rm=TRUE),"Median:",median(Fst.results,na.rm=TRUE),"\n")
  Fst.results
}


Fst.GBS.pairwise <- function(snpsubset, indsubset, populations,sortlevels=TRUE, ...) {
 if (missing(snpsubset))   snpsubset <- 1:nsnps
 if (missing(indsubset))   indsubset <- 1:nind
 popnames <- unique(populations[indsubset])
 if(sortlevels) popnames <- sort(popnames)
 npops <- length(popnames)
 Fst.results <- array(dim=c(npops,npops,length(snpsubset)))
 Fst.means <- Fst.medians <- array(dim=c(npops,npops),dimnames=list(popnames,popnames))
  for (ipop in 1:(npops-1)) {
  indsubseti <- indsubset[which(populations[indsubset] == popnames[ipop])]
  for (jpop in (ipop+1):npops) {
   indsubsetj <- indsubset[which(populations[indsubset] == popnames[jpop])]
   Fst.results[ipop,jpop,] <- Fst.GBS(snpsubset,c(indsubseti,indsubsetj),populations, ...)
   Fst.means[ipop,jpop] <- mean(Fst.results[ipop,jpop,],na.rm=TRUE)
   Fst.medians[ipop,jpop] <- median(Fst.results[ipop,jpop,],na.rm=TRUE)
   }
  }
 cat("Pairwise Fst Means\n"); print(Fst.means[-npops,-1])
 cat("Pairwise Fst Medians\n");print(Fst.medians[-npops,-1])
 Fst.results
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


manhatplot <- function(value, chrom, pos, plotname, qdistn=qunif, keyrot=0, ...) {
 chromcol <- colourby(chrom)
 colkey(chromcol,"chrom",srt=keyrot)
 plotord <- order(chrom,pos)
 png(paste0(plotname,"-Manhat.png"),width=1200)
  plot(value[plotord], col=chromcol$sampcol[plotord],xlab="Position",ylab=substitute(value),cex=0.8, xaxt="n")
  dev.off()
 png(paste0(plotname,"-QQ.png"))
  qqplot(qdistn(ppoints(length(value)),...), y=value, xlab="Theoretical quantiles", ylab=paste0(substitute(value)," quantiles"), 
         sub="Line for mid 98% of values", col=chromcol$sampcol[order(value)])
  qqline(value,col=2, distribution = function(p) qdistn(p, ...),prob=c(0.01,0.99))
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

