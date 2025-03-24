### Simulation code for Population differentiation (Fst) paper (in review).

# fade (Lumley)
fade <- function(colors,alpha) {
  rgbcols <- col2rgb(colors)
  rgb(rgbcols[1,],rgbcols[2,],rgbcols[3,],alpha,max=255)
}

# jpeg with new defaults
jpegj <- function(fname,width=15, height=13, units="cm", res=300, quality=95,...) jpeg(filename=fname, width=width, height=height, units=units, res=res, quality=quality, ...)

codedir <- "<KGD folder>"
functions.only <- TRUE
source(paste0(codedir,"GBS-PopGen.R"))
source(paste0(codedir,"GBS-Chip-Gmatrix.R"))


#-----------------
### Simulation 1
set.seed(2018)
popsize <- 100   # currently set to use one value only
popp <- c(0.2, 0.4, 0.6, 0.8)
#popp <- c(0.5,0.5,0.5,0.5)
npop <- length(popp)
meandepth <- c(2,4,10)
nrep <- 1000
nsnps <- nrep  # nsnps used by Fst.GBS
dcollist <- c("grey","steelblue","blue")

results <- data.frame(popsize=integer(0),rep=integer(0),FstAct=double(0), depth=double(0), 
                      Fststar=double(0), Fstunadj=double(0),pstar=double(0),punadj=double(0))
FstTrue = (npop-1)*var(popp)/(npop*mean(popp)*(1-mean(popp)))
for (ipopsize in seq_along(popsize)) {
  popsize0 <- popsize[ipopsize]
  nind <- popsize0*npop
  poplabels <- rep(1:npop,each=popsize0)
  for (ip in seq(npop)) {
   p0 <- popp[ip]
   genospop <- matrix(sample(0:2,size=popsize0*nrep,prob=c((1-p0)^2,2*p0*(1-p0),p0^2),replace=TRUE),ncol=nrep)
   if(ip==1) genos <- genospop else genos <- rbind(genos,genospop)
  } 
  
  aX2 <-  rbind(round(genos/2),ceiling(genos/2))
  X2results <- apply(rbind(aX2,matrix(2*popsize0,nrow=npop,ncol=nrep)),MARGIN=2,chisq.adj,y=rep(poplabels,2))
  FstAct <- X2results/(2*popsize0*npop)
  cat("mean FstAct = ",mean(FstAct),"\n")
  
  for (idepth in seq_along(meandepth)) {
    depth0 <- meandepth[idepth]
    cat("mean depth=",depth0,"\n")
    depth = matrix(rpois(popsize0*nrep*npop,depth0),ncol = nrep)
    Acounts <- matrix(rbinom(popsize0*nrep*npop,depth,genos/2),ncol=nrep)  # A allele
    genon <- trunc(2*Acounts/depth-1)+1
    Fststar <- Fst.GBS(populations=poplabels,SNPtest = TRUE)
    cat("propn with p<0.05 (Adj)  ",mean(Fststar$pvalue<0.05,na.rm=TRUE),"\n")
    depth[depth>0] <- Inf
    Fstunadj <- Fst.GBS(populations=poplabels,SNPtest = TRUE)
    cat("propn with p<0.05 (Unadj)",mean(Fstunadj$pvalue<0.05,na.rm=TRUE),"\n")
    cat("Adj.   sddiff, MAE",sd(Fststar$Fst-FstAct,na.rm=TRUE),
        mean(abs(Fststar$Fst-FstAct),na.rm=TRUE),"\n")
    cat("Unadj. sddiff, MAE",sd(Fstunadj$Fst-FstAct,na.rm=TRUE),
        mean(abs(Fstunadj$Fst-FstAct),na.rm=TRUE),"\n")
    results <- rbind(results,cbind.data.frame(popsize0,1:nrep,FstAct,depth0,
                     Fststar=Fststar$Fst,Fstunadj=Fstunadj$Fst,pstar=Fststar$pvalue,punadj=Fstunadj$pvalue))
    #print(summary(lm(Fststar$Fst~FstAct)))
    #print(summary(lm(Fststar$Fst~Fstunadj$Fst)))
  }
}

dcolour <- dcollist[match(results$depth0,meandepth)]

jpegj("Figure1.jpg")
 par(mar = par("mar") + c(0,1,-2,0))
 boxplot(c(results$Fststar,FstAct) ~ c(results$depth0,rep(Inf,length(FstAct))),col=c(dcollist,"darkblue"),
        xlab=expression(mean~depth~(italic(d))) ,ylab = expression(F[ST]^"*"), cex.lab=1.2,,xaxt="n" )
 axis(side=1,at=1:4,labels=c(meandepth,expression(infinity)))
 abline(h=FstTrue,col="red",lwd=2)
 abline(h=FstTrue+1/(2*popsize),lty = 2, col="red",lwd=2)
 dev.off()

jpegj("Figure2.jpg")  
 par(mar = par("mar") + c(0,1,-2,0))
 plot(Fststar ~ FstAct, data=results, col=fade(dcolour,127),pch=16, xlab=expression(paste(F[ST]," with true genotypes")),
      ylab=expression(paste(F[ST]^"*"," with sequencing genotypes")),cex.lab=1.2)
 abline(a=0,b=1,col="red",lwd=3)
 par(family="serif") # set and save orig setting
 legend("topleft",legend=meandepth,title="d",title.font=3,col=dcollist,pch=16)
 dev.off()

jpegj("Figure3.jpg", width=22)
 layout(matrix(1:2,nrow=1))
 par(mar = par("mar") + c(0,1,-2,0))
 plot(Fststar ~ Fstunadj, data=results, col=fade(dcolour,127),pch=16, 
     xlab=expression(paste(F[ST]," (without depth adjustment)")),
     ylab=expression(paste(F[ST]^"*"," (with depth adjustment)")),cex.lab=1.2)
 abline(a=0,b=1,col="red",lwd=3)
 text(x=max(results$Fstunadj), y=max(results$Fststar),"A", cex=1.5, adj=c(-1,-1),xpd=TRUE)
 op <- par(family="serif") # set and save orig setting
 legend("topleft",legend=meandepth,title="d",title.font=3,col=dcollist,pch=16)
 par(op) #reset
 plot(pstar ~ punadj, data=results, col=fade(dcolour,127),pch=16, 
     xlab="p-value without depth adjustment",
     ylab="p-value with depth adjustment",cex.lab=1.2, log="xy")
 abline(a=0,b=1,col="red",lwd=3)
 text(x=max(results$punadj), y=max(results$pstar),"B", cex=1.5, adj=c(-1,-1),xpd=TRUE)
 op <- par(family="serif") # set and save orig setting
 legend("topleft",legend=meandepth,title="d",title.font=3,col=dcollist,pch=16)
 dev.off()

write.csv(results,"FstResultsSim1.csv",row.names=FALSE,quote=FALSE) 


#-----------------
# find a situation where not adjusting is a problem  # Sim 2 
set.seed(20246)
popsize <- 10   # currently set to use one value only
popp <- c(0.3, 0.7)
npop <- length(popp)
meandepth <- c(0.5,1,10)
nrep <- 1000
nsnps <- nrep
testlevel <- 0.05
dcollist <- c("yellow3","grey","blue")

results <- data.frame(popsize=integer(0),rep=integer(0),FstAct=double(0), depth=double(0), 
                      Fststar=double(0), Fstunadj=double(0),pstar=double(0),punadj=double(0))
(FstTrue = (npop-1)*var(popp)/(npop*mean(popp)*(1-mean(popp))))
 
 for (ipopsize in seq_along(popsize)) {
   popsize0 <- popsize[ipopsize]
   nind <- popsize0*npop
   poplabels <- rep(1:npop,each=popsize0)
   for (ip in seq(npop)) {
     p0 <- popp[ip]
     genospop <- matrix(sample(0:2,size=popsize0*nrep,prob=c((1-p0)^2,2*p0*(1-p0),p0^2),replace=TRUE),ncol=nrep)
     if(ip==1) genos <- genospop else genos <- rbind(genos,genospop)
   } 
   
   aX2 <-  rbind(round(genos/2),ceiling(genos/2))
   X2results <- apply(rbind(aX2,matrix(2*popsize0,nrow=npop,ncol=nrep)),MARGIN=2,chisq.adj,y=rep(poplabels,2))
   FstAct <- X2results/(2*popsize0*npop)
   cat("mean FstAct = ",mean(FstAct),"\n")
   FststarAll <- list()
   for (idepth in seq_along(meandepth)) {
     depth0 <- meandepth[idepth]
     cat("mean depth=",depth0,"\n")
     depth = matrix(rpois(popsize0*nrep*npop,depth0),ncol = nrep)
     Acounts <- matrix(rbinom(popsize0*nrep*npop,depth,genos/2),ncol=nrep)
     genon <- trunc(2*Acounts/depth-1)+1
     Fststar <- Fst.GBS(populations=poplabels, SNPtest = TRUE)
     FststarAll[[idepth]] <- Fststar
     depth[depth>0] <- Inf
     Fstunadj <- Fst.GBS(populations=poplabels, SNPtest = TRUE)
     cat("Adj.   sddiff, MAE, power at",testlevel,":",sd(Fststar$Fst-FstAct,na.rm=TRUE),
         mean(abs(Fststar$Fst-FstAct),na.rm=TRUE),
         mean(Fststar$pvalue<testlevel,na.rm=TRUE),"\n")
     cat("Unadj. sddiff, MAE, power at",testlevel,":",sd(Fstunadj$Fst-FstAct,na.rm=TRUE),
         mean(abs(Fstunadj$Fst-FstAct),na.rm=TRUE),
         mean(Fstunadj$pvalue<testlevel,na.rm=TRUE),"\n")
     results <- rbind(results,cbind.data.frame(popsize0,1:nrep,FstAct,depth0,
                                               Fststar=Fststar$Fst,Fstunadj=Fstunadj$Fst,pstar=Fststar$pvalue,punadj=Fstunadj$pvalue))
   }
 }
 
dcolour <- dcollist[match(results$depth0,meandepth)]
jpegj("Figure4.jpg",width= 22)
 layout(matrix(1:2,nrow=1))
 par(mar = par("mar") + c(0,1,-2,0))
 plot(Fststar ~ Fstunadj, data=results, col=fade(dcolour,127),pch=16, xlab=expression(paste(F[ST]," (without depth adjustment)")),
      ylab=expression(paste(F[ST]^"*"," (with depth adjustment)")),cex.lab=1.2)
 abline(a=0,b=1,col="red",lwd=3)
 text(x=max(results$Fstunadj,na.rm=TRUE), y=max(results$Fststar,na.rm=TRUE),"A", cex=1.5, adj=c(-1,-1),xpd=TRUE)
 op <- par(family="serif") # set and save orig setting
 legend("topleft",legend=meandepth,title="d",title.font=3,col=dcollist,pch=16)
 par(op) #reset
 plot(pstar ~ punadj, data=results, col=fade(dcolour,127),pch=16, 
      xlab="p-value without depth adjustment",
      ylab="p-value with depth adjustment",cex.lab=1.2, log="xy")
 abline(a=0,b=1,col="red",lwd=3)
 text(x=max(results$punadj,na.rm=TRUE), y=max(results$pstar,na.rm=TRUE),"B", cex=1.5, adj=c(-1,-1),xpd=TRUE)
 op <- par(family="serif") # set and save orig setting
 legend("topleft",legend=meandepth,title="d",title.font=3,col=dcollist,pch=16)
 dev.off()


#-----------------
#Look at p values under null   # Sim 3
dcollist <- c("green") # not used
set.seed(2000)
popsize <- 100   # currently set to use one value only
popp <- c(0.5,0.5,0.5,0.5)  
npop <- length(popp)
meandepth <- c(0.5,1,10)
nrep <- 10000
nsnps <- nrep

results <- data.frame(popsize=integer(0),rep=integer(0),FstAct=double(0), depth=double(0), Fststar=double(0), Fstunadj=double(0))
FstTrue = (npop-1)*var(popp)/(npop*mean(popp)*(1-mean(popp)))
for (ipopsize in seq_along(popsize)) {
  popsize0 <- popsize[ipopsize]
  nind <- popsize0*npop
  poplabels <- rep(1:npop,each=popsize0)
  for (ip in seq(npop)) {
    p0 <- popp[ip]
    genospop <- matrix(sample(0:2,size=popsize0*nrep,prob=c((1-p0)^2,2*p0*(1-p0),p0^2),replace=TRUE),ncol=nrep)
    if(ip==1) genos <- genospop else genos <- rbind(genos,genospop)
  } 
  
  aX2 <-  rbind(round(genos/2),ceiling(genos/2))
  X2results <- apply(rbind(aX2,matrix(2*popsize0,nrow=npop,ncol=nrep)),MARGIN=2,chisq.adj,y=rep(poplabels,2))
  pvalresults <- pchisq(X2results,df=npop-1,lower.tail = FALSE)
  FstAct <- X2results/(2*popsize0*npop)
  cat("mean FstAct = ",mean(FstAct),"\n")
  cat("propn with p<0.05",mean(pvalresults<0.05,na.rm=TRUE),"\n")
  
  for (idepth in seq_along(meandepth)) {
    depth0 <- meandepth[idepth]
    cat("mean depth=",depth0,"\n")
    depth = matrix(rpois(popsize0*nrep*npop,depth0),ncol = nrep)
    Acounts <- matrix(rbinom(popsize0*nrep*npop,depth,genos/2),ncol=nrep)  # A allele
    genon <- trunc(2*Acounts/depth-1)+1
    Fststar <- Fst.GBS(populations=poplabels,SNPtest = TRUE)
    cat("propn with p<0.05",mean(Fststar$pvalue<0.05,na.rm=TRUE),"\n")
    depth[depth>0] <- Inf
    Fstunadj <- Fst.GBS(populations=poplabels, SNPtest = TRUE)
    cat("propn with p<0.05",mean(Fstunadj$pvalue<0.05,na.rm=TRUE),"\n")
    results <- rbind(results,cbind.data.frame(popsize0,1:nrep,FstAct,depth0,Fststar$Fst,Fstunadj$Fst))
  }
}


#-----------------
# sim for optimal depth Sim4
set.seed(20241002) #date redone
popp <- c(0.2, 0.4, 0.6, 0.8)
npop <- length(popp)
nrep <- 1000
nsnps <- nrep
readsperpop <- 200  # 200 reads per population
meandepth <- c(0.1,0.25, 0.5, 1,2,4,10)
resultsopt <- data.frame(popsize=integer(0),rep=integer(0),FstAct=double(0), depth=double(0), Fststar=double(0),
                         Fstunadj=double(0),pstar=double(0),punadj=double(0))
for (idepth in seq_along(meandepth)) {
  depth0 <- meandepth[idepth]
  popsize0 <- readsperpop / depth0
  nind <- popsize0*npop
  poplabels <- rep(1:npop,each=popsize0)
  for (ip in seq(npop)) {
    p0 <- popp[ip]
    genospop <- matrix(sample(0:2,size=popsize0*nrep,prob=c((1-p0)^2,2*p0*(1-p0),p0^2),replace=TRUE),ncol=nrep)
    if(ip==1) genos <- genospop else genos <- rbind(genos,genospop)
    } 
  aX2 <-  rbind(round(genos/2),ceiling(genos/2))
  X2results <- apply(rbind(aX2,matrix(2*popsize0,nrow=npop,ncol=nrep)),MARGIN=2,chisq.adj,y=rep(poplabels,2))
  FstAct <- X2results/(2*popsize0*npop)
  depth = matrix(rpois(popsize0*nrep*npop,depth0),ncol = nrep)
  Acounts <- matrix(rbinom(popsize0*nrep*npop,depth,genos/2),ncol=nrep)
  genon <- trunc(2*Acounts/depth-1)+1
  Fststar <- Fst.GBS(populations=poplabels,SNPtest = TRUE)
  depth[depth>0] <- Inf
  Fstunadj <- Fst.GBS(populations=poplabels,SNPtest = TRUE)
  resultsopt <- rbind(resultsopt,cbind.data.frame(popsize0,1:nrep,FstAct,depth0,
                                            Fststar=Fststar$Fst,Fstunadj=Fstunadj$Fst,pstar=Fststar$pvalue,punadj=Fstunadj$pvalue))
  }

depthsdopt <- aggregate(resultsopt$Fststar,list(resultsopt$depth0),sd)
#depthsdoptunadj <- aggregate(resultsopt$Fstunadj,list(resultsopt$depth0),sd)
aggregate(resultsopt$Fststar,list(resultsopt$depth0),mean)
powerstar <- resultsopt$pstar < 1e-20
depthpowopt <- aggregate(powerstar,list(resultsopt$depth0),mean, na.rm=TRUE)

jpegj("Figure5.jpg",height=14)
 par(mar = par("mar") + c(0,1,0,2))
 plot(depthsdopt,pch=16,xlab="Depth",ylab=expression(sd(F[ST]^'*')),cex.lab=1.2)
 lines(depthsdopt,lwd=2)
 par(new = TRUE)
 plot(depthpowopt,pch=16,col="blue",axes=FALSE,xlab="",ylab="",bty="n")
 lines(depthpowopt,lwd=2,col="blue")
 axis(side=4,at=seq(0,1,0.2), col.axis="blue")
 topticks <- c(0.1,2,4,10)
 axis(side=3,at=topticks, labels = readsperpop /topticks )
 mtext(expression(paste("Power (",alpha,"=",10^-20,")")),line=3, side=4, col="blue",cex=1.2)
 mtext("Number of individuals / subpopulation",line=2, side=3,cex=1.2)
 dev.off()
