#!/bin/echo Source me don't execute me 

KGDver <- "1.1.0"
cat("KGD version:",KGDver,"\n")
if (!exists("nogenos"))          nogenos          <- FALSE
if (!exists("gform"))            gform            <- "uneak"
if (!exists("genofile"))         genofile         <- "HapMap.hmc.txt"
if (!exists("sampdepth.thresh")) sampdepth.thresh <- 0.01
if (!exists("snpdepth.thresh"))  snpdepth.thresh  <- 0.01
if (!exists("hirel.thresh"))     hirel.thresh     <- 0.9
if (!exists("iemm.thresh"))      iemm.thresh      <- 0.05
if (!exists("triallelic.thresh")) triallelic.thresh <- 0.005  # if third allele present in higher than this proportion, then discard SNP, otherwise discard 3rd & 4th allele calls (currently only for ANGSD data)
if (!exists("cex.pointsize"))    cex.pointsize    <- 1
if (!exists("functions.only"))   functions.only   <- FALSE
if (!exists("alleles.keep"))     alleles.keep     <- FALSE
if (!exists("outlevel"))         outlevel         <- 9
if (!exists("use.Rcpp"))         use.Rcpp         <- TRUE
if (!exists("nThreads"))         nThreads         <- 4  # 0 means use all available
if (!exists("negC"))             negC             <- ""   # empty string will bypass negC checks
if (!exists("negCsettings"))     negCsettings     <- list()
if (!exists("QQprobpts"))        QQprobpts        <- c(0.5,0.8,0.9,0.95,0.99)

# function to locate Rcpp file (assume it is in the same directory as this file and this file was 'sourced')
pathToCppFile = function() {
    cpp.name <- "GBS-Rcpp-functions.cpp"
    this.file <- NULL
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) {
            f <- sys.frame(i)
            if (!f$chdir)
                this.file <- normalizePath(f$ofile)
            break
        }
    }
    if (!is.null(this.file)) {
        source.dir <- dirname(this.file)
        return(file.path(source.dir, cpp.name))
    }
    else {
        # assume it is in the current working directory
        return(cpp.name)
    }
}

# load C++ functions
# compiling can take ~10 seconds, so we cache the file (under ~/R/RcppCache)
# it will only be recompiled when the C++ file changes or is moved
have_rcpp <- FALSE
if (require(Rcpp) & use.Rcpp) {
    # RcppArmadillo is also required, but we don't need to load it here
    if (is.element("RcppArmadillo", installed.packages()[,"Package"])) {
        have_rcpp <- TRUE
        cpp.path <- pathToCppFile()
        cat("Loading C++ functions:", cpp.path, "\n")
        sourceCpp(file = cpp.path, showOutput = TRUE,
                  cacheDir = file.path(path.expand("~"), "R", "RcppCache"))
    }
}

### depth2K functions
 r_depth2K <- function(depthvals)  1/2^depthvals   # convert depth to K value assuming binomial 
 # select R or Rcpp version of depth2K depending on whether Rcpp is installed
 if (have_rcpp) {
     depth2K <- function(depthvals) {
         # Rcpp version only works with matrix as input, so fallback to R version otherwise
         if (is.matrix(depthvals)) {
             result <- rcpp_depth2K(depthvals, nThreads)
         } else {
             result <- r_depth2K(depthvals)
         }
         return(result)
     }
 } else {
     depth2K <- r_depth2K
 }

# convert depth to K value assuming beta-binomial with parameters alpha=beta=alph. Inf gives binomial
 r_depth2Kbb <- function(depthvals, alph=Inf) {
  if (alph==Inf) 1/2^depthvals else beta(alph,depthvals+alph)/beta(alph,alph)
  }
 # select R or Rcpp version depending on whether Rcpp is installed
 if (have_rcpp) {
    depth2Kbb <- function(depthvals, alph=Inf) {
        # Rcpp version only works with matrix as input, so fallback to R version otherwise
        if (is.matrix(depthvals) & alph < Inf) {
           if(alph < Inf) {
            result <- rcpp_depth2Kbb(depthvals, nThreads, alph)
            } else {
              result <- depth2K(depthvals)
            }
          } else {    
            result <- r_depth2Kbb(depthvals, alph)
        }
        return(result)
    }
 } else {
     depth2Kbb <- r_depth2Kbb
 }

# convert depth to K value modp model. prob of seeing same allele as last time is modp (usually >= 0.5)
 r_depth2Kmodp <- function(depthvals, modp=0.5 ) {
  Kvals <- 0.5*modp^(depthvals-1)
  Kvals[depthvals==0] <- 1
  Kvals
  }
 # select R or Rcpp version depending on whether Rcpp is installed
 if (have_rcpp) {
    depth2Kmodp <- function(depthvals, modp=0.5) {
        # Rcpp version only works with matrix as input, so fallback to R version otherwise
        if (is.matrix(depthvals)) {
            result <- rcpp_depth2Kmodp(depthvals, modp, nThreads)
        } else {
            result <- r_depth2Kmodp(depthvals, modp)
        }
        return(result)
    }
 } else {
     depth2Kmodp <- r_depth2Kmodp
 }

depth2Kchoose <- function(dmodel="bb",param) {  # function to choose redefine depth2K
 if (!dmodel=="modp") dmodel <- "bb"
 if (missing(param) & dmodel=="bb") param <- Inf
 if (missing(param) & dmodel=="modp") param <- 0.5
 if (dmodel=="bb") depth2K <- function(depthvals) depth2Kbb(depthvals,alph=param)
 if (dmodel=="modp") depth2K <- function(depthvals) depth2Kmodp(depthvals,modp=param)
 depth2K
 }


readGBS <- function(genofilefn = genofile, usedt="recommended") {
 if (gform == "chip") readChip(genofilefn)
 if (gform == "ANGSDcounts") readANGSD(genofilefn)
 if (gform == "TagDigger") readTD(genofilefn)
 if (toupper(gform) == "VCF") read.vcf(genofilefn)
 if (gform %in% c("uneak","Tassel")) readTassel(genofilefn, usedt=usedt)
 ndups <- sum(duplicated(seqID))
 if(ndups>0) cat("Warning: ",ndups,"duplicated seqIDs, 1st one is:",seqID[duplicated(seqID)][1],"\n")
 }

readTD <- function(genofilefn0 = genofile, skipcols=0) {
  havedt <- require("data.table")
  ghead <- scan(genofilefn0, what = "", nlines = 1, sep = ",", quote="", quiet=TRUE)
  if(skipcols > 0) ghead <- ghead[-(1:skipcols)]
  nsnps <<- (length(ghead) - 1)/2
  SNPinfo <- read.table(text=ghead[-1],sep="_",fill=TRUE,stringsAsFactors=FALSE)
  SNP_Names <<- SNPinfo[2*seq(nsnps),1]
  refalleles <<- SNPinfo[2*seq(nsnps)-1,2]
  altalleles <<- SNPinfo[2*seq(nsnps),2]
  cat("Data file has", nsnps, "SNPs \n")
  print(table(refalleles,altalleles))
  isgzfile <- grepl(".gz",genofilefn0) #gz unzipping will only work on linux systems
  if(nogenos) {
   if(.Platform$OS.type == "unix") {  # faster to cut in unix than to scan whole file for first field
    if(isgzfile) seqID <<- read.table(text=system(paste("zcat",genofilefn0, "| cut -f 1 -d,"), intern = TRUE, ignore.stderr = TRUE),header=FALSE, stringsAsFactors=FALSE)[,1]
    if(!isgzfile) seqID <<- read.table(text=system(paste("cut",genofilefn0, "-f 1 -d,"), intern = TRUE, ignore.stderr = TRUE),header=FALSE, stringsAsFactors=FALSE)[,1]
    } else {
    seqID <<-  scan(genofilefn0, what="",skip=1,quote="",flush=TRUE, quiet=TRUE) # 1st column only
    }
   nind <<- length(seqID)
   }
  if(!nogenos) {
   if (havedt) {
    if ( packageVersion("data.table") < "1.12" & isgzfile) {
     genosin <- fread(paste("gunzip -c",genofilefn0),sep=",",header=TRUE,showProgress=FALSE)
     } else {    
     genosin <- fread(genofilefn0,sep=",",header=TRUE)
     }
    if(skipcols > 0) genosin <- genosin[,-(1:skipcols)]
    seqID <<- as.character(as.matrix(genosin[,1])[,1])  # (Allows any name for first col), also convert to vector from data.table
    nind <<- length(seqID)
    alleles <<- as.matrix(genosin[,-1,with=FALSE])
    } else {
    genosin <- scan(genofilefn0, skip = 1, sep = ",", what = c(list(seqID = ""), rep(list(0), 2*nsnps)), quote="", quiet=TRUE) 
    seqID <<- as.character(genosin[[1]])
    nind <<- length(seqID)
    alleles <<- matrix(0, nrow = nind, ncol = 2 * nsnps)
    for (isnp in seq(2*nsnps)) alleles[, isnp] <<- genosin[[isnp+1]] 
    }
   }
  cat("Data file has", nind, "samples \n")
  invisible(NULL)
}

readANGSD <- function(genofilefn0 = genofile) {
  ghead <- scan(genofilefn0, what = "", nlines = 1, sep = "\t", quote="")
  nind <<- (length(ghead) - 1) / 4
  if (nind != floor(nind)) nind  <<- length(ghead) / 4
  if (nind != floor(nind)) print("Incorrect number of columns")
  seqID <<- substr(ghead,1,nchar(ghead)-2)[seq(4,4*nind,4)]
  if(nogenos) cat("Warning: No SNP info yet with gform ANGSDcounts and nogenos\n")
  if(!nogenos) {
   genosin <- scan(genofilefn0, skip = 1, sep = "\t", flush=TRUE, what = rep(list(integer(0)), 4*nind), quote="")
   nsnps <<- length(genosin[[1]])
   SNP_Names <<- paste0("SNP",formatC(seq(nsnps),width=nchar(nsnps),flag="0"))  # with leading zeros
   alleles4 <- matrix(aperm(array(unlist(genosin,use.names=FALSE),dim =c(nsnps,4,nind)), c(3,2,1)), nrow=nind)
    # cols in sets of 4 alleles, rows are samples  (temporary, until relevant alleles extracted)
   acounts <- matrix(colSums(alleles4),ncol=4,byrow=TRUE)
   acountord <- t(apply(acounts,1,order,decreasing=TRUE))
   SNP.discard <- which(acounts[cbind(1:nsnps,acountord[,3])] / rowSums(acounts) > triallelic.thresh)
   snpcols <- sort(c(seq(0,4*(nsnps-1),4)+acountord[,1],seq(0,4*(nsnps-1),4)+acountord[,2]))
   alleles <<- alleles4[,snpcols]
   if(length(SNP.discard) > 0) {
    alleles <<- alleles[,-c(2*SNP.discard-1,2*SNP.discard)]
    nsnps <<- nsnps - length(SNP.discard)
    SNP_Names <<- SNP_Names[-SNP.discard]
    cat(length(SNP.discard),"SNP(s) removed with 3rd allele frequency >",triallelic.thresh,"\n")
    }
  }
 invisible(NULL)
}

 read.vcf <- function(vcffile=genofile) {  # uses AD then GT if either/both available
  # Needs data.table and (for .gz files) R.utils
  if(!require(data.table)) stop("data.table package required for read.vcf")
  if(nogenos) cat("Warning: No SNP or sample info yet with gform VCF and nogenos\n")
  if(!nogenos) {
   vcfin <- fread(vcffile,skip="#CHROM")   # skip headers until #CHROM found on a row
   if(colnames(vcfin)[1]=="V1") { 
    cat("Correcting added first column name, remove trailing (blank?) column (ignore fread warning below)\n")
    colnames(vcfin)[1:ncol(vcfin)-1] <-  colnames(vcfin)[2:ncol(vcfin)] 
    vcfin[,ncol(vcfin)] <- NULL
    }
   SNP_Names <<- vcfin$ID
   refalleles <<- vcfin$REF
   altalleles <<- vcfin$ALT
   nsnps <<- length(SNP_Names)
   formatcol <- which(colnames(vcfin)=="FORMAT")
   seqID <<-  colnames(vcfin)[-(1:formatcol)]
   nind <<- length(seqID)
   chrom <<- vcfin$`#CHROM`
   pos <<- vcfin$POS
   if(all(SNP_Names==".")) SNP_Names <<- paste(chrom,pos,sep="_")
   nfields <- 1+ lengths(regmatches(vcfin$FORMAT,gregexpr(":",vcfin$FORMAT)))
   tempformats <- read.table(text=vcfin$FORMAT,sep=":",fill=TRUE,stringsAsFactors=FALSE,col.names=paste0("V",1:max(nfields)))
   gthave = apply(tempformats=="GT",1,any)
   gtpos = unlist(apply(tempformats=="GT",1,which))
   adhave = apply(tempformats=="AD",1,any)
   adpos = unlist(apply(tempformats=="AD",1,which))
   genon <<- matrix(NA,nrow=nind,ncol=nsnps)
   if(any(gthave)) {
    tempgt <- read.table(text=as.matrix(vcfin[gthave,-(1:formatcol)]),sep=":",fill=TRUE,stringsAsFactors=FALSE)[matrix(c(1:(nind*sum(gthave)),rep(gtpos,nind)),ncol=2,dimnames=list(NULL,c("row","col")))]
    tempgt[tempgt=="."] <- "./."
    genongt <- t(matrix(rowSums(read.table(text=gsub("/","|",tempgt),sep="|",na.strings=".")),nrow=sum(gthave)))
    genon[,which(gthave)] <<- genongt
    }
   depth <<- matrix(Inf, nrow = nind, ncol = nsnps)
   if(any(adhave)) {
    tempad <- read.table(text=as.matrix(vcfin[adhave,-(1:formatcol)]),sep=":",fill=TRUE,stringsAsFactors=FALSE)[matrix(c(1:(nind*sum(adhave)),rep(adpos,nind)),ncol=2,dimnames=list(NULL,c("row","col")))]
    tempad2 <- read.table(text=sub(".","0,0",sub(".,.","0,0",tempad,fixed=TRUE),fixed=TRUE),sep=",")
    ref <- t(matrix(tempad2[,1],nrow=sum(adhave)))
    alt <- t(matrix(tempad2[,2],nrow=sum(adhave)))
    genonad <- trunc((2*ref/(ref+alt))-1)+1
    genon[,which(adhave)] <<- genonad
    alleles <<- matrix(Inf, nrow = nind, ncol = 2 * nsnps)
    alleles[,seq(1, 2 * nsnps - 1, 2)[which(adhave)]] <<- ref
    alleles[,seq(2, 2 * nsnps, 2)[which(adhave)]] <<- alt
    depth[,which(adhave)] <<- ref + alt
    }
   depth[is.na(genon)] <<- 0
   p <<- colMeans(genon, na.rm = TRUE)/2 # same as pg further down
   if(!any(adhave)) gform <- "chip"
   }
  invisible(NULL)
 }

readChip <- function(genofilefn0 = genofile) {
  ghead <- scan(genofilefn0,what="",nlines=1,sep=",", quote="", quiet=TRUE)
  SNP_Names <<- ghead[-1]
  nsnps <<- length(SNP_Names)
  cat("Data file has", nsnps, "SNPs \n")
  if(nogenos) cat("Warning: No sample info yet with gform chip data and nogenos\n")
  if(!nogenos) {
   genost <- scan(genofilefn0,what="",skip=1,sep=",", quote="", quiet=TRUE) # read as text ... easier to pull out elements than list of nsnps+1
   snpnums <- ((1:length(genost))-1) %% (nsnps+1)
   genon <<- matrix(as.numeric(genost[which(snpnums !=0)]) ,ncol=nsnps,byrow=TRUE)
   seqID <<-  genost[which(snpnums ==0)]
   nind <<- length(seqID)
   cat("Data file has", nind, "samples \n")
   rm(genost)
   depth <<- matrix(Inf, nrow = nind, ncol = nsnps)
   depth[is.na(genon)] <<- 0
   p <<- colMeans(genon, na.rm = TRUE)/2 # same as pg further down
   }
  invisible(NULL)
}

readTassel <- function(genofilefn0 = genofile, usedt="recommended") {
  havedt <- FALSE       # problem with read.table(text= when > 2bill
  if(usedt=="always")   havedt <- require("data.table")
  gsep <- switch(gform, uneak = "|", Tassel = ",")
  ghead <- scan(genofilefn0, what = "", nlines = 1, sep = "\t", quote="", quiet=TRUE)
  nind <<- length(ghead) - switch(gform, uneak = 6, Tassel = 2)
  cat("Data file has", nind, "samples \n")
  seqID <<- switch(gform, uneak = ghead[2:(nind + 1)], Tassel = ghead[-(1:2)])
   # find SNPs without reading all genos (perhaps do this only if nogenos is TRUE if its taking too long)
   if(gform=="uneak") SNP_Names <<- scan(genofilefn0, what="",skip=1,quote="",flush=TRUE, quiet=TRUE) # 1st column only
   if(gform=="Tassel") {
    tempin <- scan(genofilefn0, what=list(CHROM = "", POS = 0),skip=1,quote="",flush=TRUE, quiet=TRUE) # 1st 2 columns only
    chrom <<- tempin$CHROM
    pos <<- tempin$POS
    SNP_Names <<- paste(chrom,pos,sep="_")
    }
   nsnps <<- length(SNP_Names)
   cat("Data file has", nsnps, "SNPs \n")
   if(usedt=="recommended" & nsnps * nind > 2^30) havedt <- FALSE    # problem with read.table(text= when > 2bill
  if(!nogenos) {
  if (havedt & gform=="Tassel") {
    # placeholder for code to use data.table
    isgzfile <- grepl(".gz",genofilefn0) #gz unzipping will only work on linux systems
    if ( packageVersion("data.table") < "1.12" & isgzfile) {
     genosin <- fread(paste("gunzip -c",genofilefn0),sep="\t",header=TRUE,showProgress=FALSE)
     } else {    
     genosin <- fread(genofilefn0,sep="\t",header=TRUE)
     }
   # alleles <<- as.matrix(reshape(data.frame(ids=rep(1:nind,nsnps),snps=rep(1:nsnps,each=nind),
   #                               read.table(text= t(as.matrix(genosin[,-(1:2),with=FALSE])),sep=",")),
   #                               direction="wide",idvar="ids",timevar="snps"))[,-1]
    # dcast syntax dcast(data.table(ids ...), ids ~ snps,value.var = c("V1","V2" )) puts all V1 then all V2 so need to rearrange ...
    tempalleles <- as.matrix(dcast(data.table(ids=rep(1:nind,nsnps),snps=rep(1:nsnps,each=nind),
                                 read.table(text= t(as.matrix(genosin[,-(1:2),with=FALSE])),sep=",")),
                                 ids ~ snps,value.var = c("V1","V2" ))[,-1])
    alleles <<- matrix(0,nrow=nind, ncol=2*nsnps)
    alleles[,seq(1, 2 * nsnps - 1, 2)] <<- tempalleles[,1:nsnps]
    alleles[,seq(2, 2 * nsnps, 2)] <<- tempalleles[,nsnps+(1:nsnps)]
    } else {  # use scan to read input
    if (gform == "Tassel") genosin <- scan(genofilefn0, skip = 1, sep = "\t", what = c(list(chrom = "", coord = 0), rep(list(""), nind)), quote="", quiet=TRUE)
    if (gform == "uneak") genosin <- scan(genofilefn0, skip = 1, sep = "\t", what = c(list(chrom = ""), rep(list(""), nind), list(hetc1 = 0, hetc2 = 0, acount1 = 0, acount2 = 0, p = 0)), quote="", quiet=TRUE)
    alleles <<- matrix(0, nrow = nind, ncol = 2 * nsnps)
    for (iind in 1:nind) alleles[iind, ] <<- matrix(as.numeric(unlist(strsplit(genosin[[iind + switch(gform, uneak = 1, Tassel = 2)]], split = gsep, 
                                                   fixed = TRUE))), nrow = 1)
    if (gform == "uneak") AFrq <<- genosin[[length(genosin)]]
    }
   }
  unadepth <- which(is.na(alleles))
  if(length(unadepth)>0) alleles[unadepth] <<- 0
  invisible(NULL)
}

parkGBS <- function() {
 #refalleles altalleles
 if(!exists("alleles")) {
    alleles <- matrix(0,nrow=nind, ncol=2*nsnps)
    tempalleles <- matrix(0,nrow=nind,ncol=nsnps)
    tempalleles[which(genon == 0)] <- depth[which(genon==0)]
    tempalleles[which(genon == 1)] <- depth[which(genon==1)]/2
    alleles[,seq(1, 2 * nsnps - 1, 2)] <- tempalleles
    tempalleles[which(genon == 2)] <- depth[which(genon==2)]
    tempalleles[which(genon == 0)] <- 0
    alleles[,seq(2, 2 * nsnps, 2)] <- tempalleles
  }
 parkeddata <- list(nsnps=nsnps,SNP_Names=SNP_Names, seqID = seqID, nind=nind, alleles = alleles) 
 }

activateGBS <- function(GBSobj) {
 nsnps <<- GBSobj$nsnps; SNP_Names <<- GBSobj$SNP_Names; seqID <<- GBSobj$seqID; nind <<- GBSobj$nind; alleles <<- GBSobj$alleles
 if(exists("depth")) rm(depth, p, genon, inherits=TRUE)
 invisible(NULL)
 }

joinGBS <- function(join1, join2=NULL, replace=TRUE, uniqueSNPs=FALSE) {
 if(is.null(join2)) join2 <- parkGBS()
 cat("Joining set 1",join1$nind,"ind x",join1$nsnps,"SNPs, set 2",join2$nind,"ind x",join2$nsnps,"SNPs\n")
 if(uniqueSNPs) join1$SNP_Names <- make.names(c(join2$SNP_Names,join1$SNP_Names),unique=TRUE)[join2$nsnps+(1:join1$nsnps)]
 SNP_Names <- unique( c(join1$SNP_Names,join2$SNP_Names))
 nsnps <- length(SNP_Names)
 nind <- join1$nind+join2$nind
 seqID <- c(join1$seqID,join2$seqID)
 snppos1 <- match(join1$SNP_Names,SNP_Names)
 snppos2 <- match(join2$SNP_Names,SNP_Names)
 alleles <- matrix(0,nrow=nind,ncol=2*nsnps)
 alleles[1:join1$nind,as.vector(rbind(snppos1*2-1,snppos1*2))] <- join1$alleles
 alleles[join1$nind+(1:join2$nind),as.vector(rbind(snppos2*2-1,snppos2*2))] <- join2$alleles
 if(any(duplicated(seqID))) cat(sum(duplicated(seqID)),"samples with seqID already present (kept separate). Consider using mergeSamples(seqID)\n")
 # scoping problems if mergeSamples used inside this function
 outobj <- NULL
 if(replace) {
  nsnps <<- nsnps; SNP_Names <<- SNP_Names; seqID <<- seqID; nind <<- nind; alleles <<- alleles
  if(exists("depth")) rm(depth,genon, pos=1)
  } else {
  outobj <- list(nsnps=nsnps,SNP_Names=SNP_Names, seqID = seqID, nind=nind, alleles = alleles)
  }
 invisible(outobj)
 }

samp.remove <- function (samppos=NULL, keep=FALSE) {
 if(keep) samppos <- setdiff(1:nind,samppos)
 if(length(samppos)>0) {
    if(gform != "chip") alleles <<- alleles[-samppos, ]
    if(gform == "chip" & exists("alleles") & alleles.keep) alleles <- alleles[-samppos,]
    if(exists("depth")) depth <<- depth[-samppos, ]
    if(exists("genon")) genon <<- genon[-samppos,]
    if(exists("sampdepth")) sampdepth <<- sampdepth[-samppos]
    seqID <<- seqID[-samppos]
    nind <<- nind - length(samppos)
  }
 }

snp.remove <- function(snppos=NULL, keep=FALSE) {
 if(keep) snppos <- setdiff(1:nsnps,snppos)
 if (length(snppos) > 0) {
   SNP_Names <<- SNP_Names[-snppos]
   nsnps <<- length(SNP_Names)
   if(exists("p")) p <<- p[-snppos]
   if(exists("depth")) depth <<- depth[, -snppos]
   if(exists("genon")) genon <<- genon[, -snppos]
   if(exists("chrom")) chrom <<- chrom[-snppos]
   if(exists("pos")) pos <<- pos[-snppos]
   if(exists("refalleles")) refalleles <<- refalleles[-snppos]
   if(exists("altalleles")) altalleles <<- altalleles[-snppos]
   if (exists("alleles")) {
     uremovea <- sort(c(2 * snppos, 2 * snppos - 1))  # allele positions
     if(exists("RAcounts")) RAcounts <<- RAcounts[-snppos, ]
     alleles <<- alleles[, -uremovea]
     if(exists("allelecounts")) allelecounts <<- allelecounts[uremovea]
     if (gform == "uneak") AFrq <<- AFrq[-snppos]
   }
  }
 }

HWpops <- function(snpsubset, indsubset, populations=NULL, depthmat = depth) {
 if (missing(snpsubset))   snpsubset <- 1:nsnps
 if (missing(indsubset))   indsubset <- 1:nind
 if (is.null(populations)) populations <- rep("A",nind)
 popnames <- unique(populations[indsubset])
 npops <- length(popnames)
 HWdis <- matrix(NA,nrow=npops,ncol=length(snpsubset)); rownames(HWdis) <- popnames  #initialise
 l10LRT <- x2star <- l10pstar <- maf <- HWdis
 for(ipop in 1:npops) {
  thigroup <- popnames[ipop]
  indgroup <- intersect(indsubset,which(populations==thigroup))
  naa <- colSums(genon[indgroup,snpsubset,drop=FALSE] == 2, na.rm = TRUE)
  nab <- colSums(genon[indgroup,snpsubset,drop=FALSE] == 1, na.rm = TRUE)
  nbb <- colSums(genon[indgroup,snpsubset,drop=FALSE] == 0, na.rm = TRUE)
  n1 <- 2 * naa + nab
  n2 <- nab + 2 * nbb
  n <- n1 + n2  #n alleles
  p1 <- n1/n
  p2 <- 1 - p1
  HWdis[ipop,] <- naa/(naa + nab + nbb) - p1 * p1
  x2 <- (naa + nab + nbb) * HWdis[ipop,]^2/(p1^2 * p2^2)
  LRT <- 2 * (n * log(n) + naa * log(pmax(1, naa)) + nab * log(pmax(1, nab)) + nbb * log(pmax(1, nbb)) - (n/2) * log(n/2) - n1 * log(n1) - n2 * 
               log(n2) - nab * log(2))  # n is # alleles = 2* n obs
#  l10p <- -log10(exp(1)) * pchisq(x2, 1, lower.tail = FALSE, log.p = TRUE)
  l10LRT[ipop,] <- -log10(exp(1)) * pchisq(LRT, 1, lower.tail = FALSE, log.p = TRUE)
  Kdepth <- depth2K(depthmat[indgroup,snpsubset,drop=FALSE])
  Kdepth[depthmat[indgroup,snpsubset,drop=FALSE]==0] <- NA
  esnphetstar <- 2*p1*p2*(1-2*colMeans(Kdepth,na.rm=TRUE))
  osnphetstar <- nab/(naa + nab + nbb)
  x2star[ipop,] <- colSums(1-2*Kdepth,na.rm=TRUE)*(1-osnphetstar/esnphetstar)^2 # corrected Nov 2018
  l10pstar[ipop,] <- -log10(exp(1)) * pchisq(x2star[ipop,], 1, lower.tail = FALSE, log.p = TRUE)
  maf[ipop,] <- ifelse(p1 > 0.5, p2, p1)
  }
  outobj <-  list( HWdis=HWdis, l10LRT=l10LRT, x2star=x2star, l10pstar=l10pstar, maf=maf)
  if(npops>1) outobj <- c(outobj, list(l10pstar.pop =  -log10(exp(1)) * pchisq(colSums(x2star,na.rm=TRUE), colSums(!is.na(x2star)), lower.tail = FALSE, log.p = TRUE) ))
  outobj
  }

finplot <- function(HWdiseq=HWdis, MAF=maf,  plotname="finplot", finpalette=palette.aquatic, finxlim=c(0,0.5), finylim=c(-0.25, 0.25)) {
 depthtrans <- function(x) round(20 * log(-log(1/(x + 0.9)) + 1.05))  # to compress colour scale at higher depths
 depthpoints <- c(0.5, 5, 50, 250)  # legend points
 transpoints <- depthtrans(depthpoints)
 mindepthplot <- 0.1
 maxdepthplot <- 256
 maxtrans <- depthtrans(maxdepthplot)
 legend_image <- as.raster(matrix(rev(finpalette[1:maxtrans]), ncol = 1))
 png(paste0(plotname,".png"), width = 960, height = 960, pointsize = cex.pointsize *  18)
  if(whitedist(finpalette) < 25) par(bg="grey")
  plot(HWdiseq ~ MAF, col = finpalette[depthtrans(pmax(mindepthplot, pmin(snpdepth, maxdepthplot)))], cex = 0.8, 
       xlab = "Minor allele frequency", ylab = "Hardy-Weinberg disequilibrium", cex.lab = 1.5, xlim=finxlim, ylim=finylim)
  rasterImage(legend_image, 0.05, -0.2, 0.07, -0.1)
  text(x = 0.1, y = -0.2 + 0.1 * transpoints/maxtrans, labels = format(depthpoints))
  text(x = 0.075, y = -0.075, labels = "SNP Depth", cex = 1.2)
  dev.off()
 }

HWsigplot <- function(HWdiseq=HWdis, MAF=maf, ll=l10LRT, plotname="HWdisMAFsig", finpalette=palette.aquatic, finxlim=c(0,0.5), finylim=c(-0.25, 0.25), 
                      llname=expression('-log'[10]*' LRT'), sortord=ll) {
 sigtrans <- function(x) round(sqrt(x) * 40/max(sqrt(ll),na.rm=TRUE)) + 1  # to compress colour scale at higher LRT
 sigpoints <- c(0.5, 2, 5, 20, 50, 200, 500)  # legend points
 sigpoints <- sigpoints[union(1:2,which(sigpoints < max(ll,na.rm=TRUE)))]
 if(length(sigpoints) > 5) sigpoints <- sigpoints[seq(1,length(sigpoints),2)]
 transpoints <- sigtrans(sigpoints)
 maxtrans <- sigtrans(max(ll,na.rm=TRUE))
 legend_image <- as.raster(matrix(rev(finpalette[1:maxtrans]), ncol = 1))
 plotord <- 1:length(MAF)
 if(!is.null(sortord)) plotord <- order(sortord)
 png(paste0(plotname,".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
  if(whitedist(finpalette) < 25) par(bg="grey")
  plot(HWdiseq[plotord] ~ MAF[plotord], col = finpalette[sigtrans(ll)][plotord], cex = 0.8, xlab = "Minor Allele Frequency", 
       ylab = "Hardy-Weinberg disequilibrium", cex.lab = 1.5, xlim=finxlim, ylim=finylim)
  rasterImage(legend_image, 0.05, -0.2, 0.07, -0.1)
  text(x = 0.1, y = -0.2 + 0.1 * transpoints/maxtrans, labels = format(sigpoints))
  text(x = 0.075, y = -0.075, labels = llname, cex = 1.2)
  dev.off()
 }

finclass <- function(HWdiseq=HWdis, MAF=maf, colobj, classname=NULL, plotname="finclass", finxlim=c(0,0.5), finylim=c(-0.25, 0.25)) {
 if(missing(colobj)) colobj <-  list(collabels="",collist="black",sampcol=rep("black",length(MAF)))
 plotord <- 1:length(MAF)  # fixed for now, keep in case want to allow reordering
 png(paste0(plotname,".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
  plot(HWdiseq[plotord] ~ MAF[plotord], col = colobj$sampcol[plotord], cex = 0.8, xlab = "Minor Allele Frequency", 
       ylab = "Hardy-Weinberg disequilibrium", cex.lab = 1.5, xlim=finxlim, ylim=finylim)
  legend(0.05, -0.075, legend=colobj$collabels,col=colobj$collist,pch=1,title=classname)
  dev.off()
 }

x2starplots <- function () {
 HWsigplot(ll=l10pstar,llname=expression('-log'[10]*'p X'^2*'*'), finpalette=colorRampPalette(c("deepskyblue2","red"))(50))
 yaxpts <- quantile(x2star,QQprobpts, na.rm=TRUE)  # NA when max depth is 1
 png("X2star-QQ.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
  par(mar = c(5.1, 4.1, 4.1, 4.1))
  qqplot(qchisq(ppoints(nsnps), df = 1), x2star, cex=0.75, main = parse(text = "Hardy-Weinberg ~~ X^2~~ '* Q-Q Plot'"), 
         xlab = parse(text = "Theoretical ~~ (chi[1]^2) ~~  Quantiles"), ylab = "Sample Quantiles")
  ytck <- axis(side=4, tick=FALSE, labels=FALSE)
  axis(side=4, at=ytck, labels = signif(-log10(exp(1)) * pchisq(ytck, 1, lower.tail = FALSE, log.p = TRUE),3))
  mtext(text= expression('-log'[10]*'p X'^2*'*'), side=4, line=3)
  if(length(QQprobpts) > 0) {
   for(ipt in 1:length(QQprobpts)) {
    lines(x=c(0,xaxpts[ipt],xaxpts[ipt]),y=c(yaxpts[ipt],yaxpts[ipt],0), col="grey")
    }
   text(x=xaxpts/4,y=yaxpts,labels=QQprobpts, pos=3, col="grey")
   }
  dev.off()
 }


HDplot <- function(plotname="HDplot", colourtype = "depth", finpalette=palette.aquatic, HDxlim=c(0,1), HDylim=c(-Inf, Inf), HDcol=NULL,
          sortcol = "asc") {
 # from Mckinney Mol Ecol Res 2017
 # colourtype can be "depth", "HW" (use l10LRT), "HW*" (use l10pstar)
 # sortcol can be "asc", "desc" or ""
 if(!exists("alleles")) cat("Cannot produce HD plots without alleles object\n")
 if(exists("alleles")) {
  if(!colourtype %in% c("depth", "HW", "HW*")) colourtype <- "depth"
  if(!sortcol %in% c("asc","desc")) sortcol <- ""
  if(!exists("HD.saved.KGD")) {
   HD.saved.KGD <<- TRUE
   Het <<- colSums(genon == 1, na.rm=TRUE) / (nind * SNPcallrate)
   SNPN <<- colSums(depth * (genon ==1), na.rm=TRUE)
   SNPNA <<- colSums(alleles[, seq(1, 2 * nsnps - 1, 2)] * (genon ==1), na.rm=TRUE)
   storage.mode(SNPN) <- storage.mode(SNPNA) <- "integer"
   }
  SNPz <- (SNPN/2 - SNPNA)/sqrt(SNPN * 0.5 * 0.5)  # denom is SNPsd
  if(colourtype == "depth") {
   colourvar <- snpdepth
   legendlab <- "SNP Depth"
   colvtrans <- function(x) round(20 * log(-log(1/(x + 0.9)) + 1.05))  # to compress colour scale at higher depths
   colvpoints <- c(0.5, 5, 50, 250)  # legend points
   mincolvplot <- 0.1
   maxcolvplot <- 256
   }
  if(colourtype=="HW") {colourvar <- l10LRT; legendlab <- expression('-log'[10]*' LRT')}
  if(colourtype=="HW*") {colourvar <- l10pstar; legendlab <- expression('-log'[10]*'p X'^2*'*')}
  if(grepl("^HW", colourtype)) {
   mincolvplot <- 0
   maxcolvplot <- max(colourvar,na.rm=TRUE)
   colvtrans <- function(x) round(sqrt(x) * 40/max(sqrt(colourvar),na.rm=TRUE)) + 1  # to compress colour scale at higher LRT
   colvpoints <- c(0.5, 2, 5, 20, 50, 200, 500)  # legend points
   colvpoints <- colvpoints[union(1:2,which(colvpoints < maxcolvplot))]
   if(length(colvpoints) > 5) colvpoints <- colvpoints[seq(1,length(colvpoints),2)]
   }
  maxtrans <- colvtrans(maxcolvplot)
  mintrans <- colvtrans(mincolvplot)
  transpoints <- colvtrans(colvpoints)

  HDxlim[2] <- HDxlim[2] + 0.13  # make some room for legend
  legend_image <- as.raster(matrix(rev(finpalette[1:maxtrans]), ncol = 1))
  HDc <- HDcol
  if(is.null(HDcol)) HDc <- finpalette[colvtrans(pmax(mincolvplot, pmin(colourvar, maxcolvplot)))]
  uplot <- which(SNPN > 0)
  if( HDylim[1] == -Inf) HDylim[1] <- min(SNPz[uplot])
  if( HDylim[2] == Inf) HDylim[2] <- max(SNPz[uplot])
  plotch <- rep(1,nsnps)
  plotch[which(SNPz > HDylim[2])] <- 94   # ^
  plotch[which(SNPz < HDylim[1])] <- 118  # v
  plotord <- 1:nsnps
  if(is.null(HDcol)) {
   if(sortcol == "asc") plotord <- order(colourvar)
   if(sortcol == "desc") plotord <- order(colourvar, decreasing = TRUE)
   }
  png(paste0(plotname,".png"), width = 960, height = 960, pointsize = cex.pointsize *  18)
   if(whitedist(finpalette) < 25) par(bg="grey")
   plot(pmax(pmin(SNPz,HDylim[2]),HDylim[1])[plotord][uplot] ~ Het[plotord][uplot], col = HDc[plotord][uplot], cex = 0.8,  pch = plotch[plotord][uplot],
        xlab = "Proportion of heterozygotes (H)", ylab = "Read ratio deviation (D)", cex.lab = 1.5, xlim=HDxlim, ylim=HDylim)
      # frame.plot = FALSE,
   if(is.null(HDcol)) {
    x0 <- 0.96  # lhs of key
    rect((x0-.05)*HDxlim[2], 1.1*HDylim[1]+0.05*HDylim[2], (x0+.13)*HDxlim[2], 0.85*HDylim[1]+0.15*HDylim[2], col="white", border="grey", xpd=NA)
    rasterImage(legend_image, x0*HDxlim[2], HDylim[1], (x0+0.04)*HDxlim[2], 0.9*HDylim[1]+0.1*HDylim[2])
    text(x = (x0+0.07)*HDxlim[2], y = HDylim[1] + 0.1*(HDylim[2]-HDylim[1]) * (transpoints-mintrans)/(maxtrans-mintrans), labels = format(colvpoints), xpd=NA, cex=0.8)
    text(x = (x0+0.03)*HDxlim[2], y =0.88*HDylim[1]+0.12*HDylim[2], labels = legendlab, cex = 1, xpd=NA)
    }
   dev.off()
  }
 }

mafplot <- function(MAF=maf,plotname="MAF", barcol="grey", doplot=TRUE, ...) {
 if(doplot) png(paste0(plotname,".png"), pointsize = cex.pointsize * 12)
  histinfo <- hist(MAF, breaks = 50, xlab = "Minor Allele Frequency", col = barcol, plot= doplot, ...)
  if(doplot) dev.off()
  if(doplot) histinfo <- paste(plotname," png output created")
  histinfo
 }

na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
}

calcp <- function(indsubset, pmethod="A") {
 if(!pmethod == "G") pmethod <- "A"
 if (missing(indsubset))   indsubset <- 1:nind
 if(exists("depth")) if(any(depth==Inf)) pmethod <- "G"
 if (pmethod == "A") {
  if (nrow(alleles) == nind) {
   RAcountstemp <- matrix(colSums(alleles[indsubset,,drop=FALSE]), ncol = 2, byrow = TRUE)  # 1 row per SNP, ref and alt allele counts
   afreqs <- RAcountstemp[, 1]/rowSums(RAcountstemp)  # p for ref allele - based on # reads, not on inferred # alleles
   } else {
   afreqs <- NULL
   print("Error: alleles is wrong size for pmethod A, using pmethod G instead")
   pmethod <- "G"
   }
  }
 if (pmethod == "G") afreqs <- colMeans(genon[indsubset,,drop=FALSE], na.rm = TRUE)/2  # allele freq assuming genotype calls
 afreqs
 }

GBSsummary <- function() {
 havedepth <- exists("depth")  # if depth present, assume it and genon are correct & shouldn't be recalculated (as alleles may be the wrong one)
 if(havedepth & !gform %in% c("VCF","chip")) cat("Warning: depth object already exists - reusing\n")
 if(gform != "chip") {
  if (!havedepth) depth <<- alleles[, seq(1, 2 * nsnps - 1, 2)] + alleles[, seq(2, 2 * nsnps, 2)]
#  storage.mode(depth) <- "integer"
  if (nchar(negC) > 0) { # check and report negative controls
   uneg <- do.call(grep,c(list(negC,seqID),negCsettings))
   if(length(uneg) > 0 ) {
    negCstats <<- data.frame(seqID=seqID[uneg], callrate= 1 - rowSums(depth[uneg,,drop=FALSE] == 0)/nsnps, sampdepth=rowMeans(depth[uneg,,drop=FALSE]), stringsAsFactors = FALSE)
    write.csv(negCstats, "negCStats.csv", row.names = FALSE)
    if(length(uneg)>0) cat(length(uneg),"Negative controls removed\n  mean call rate = ",mean(negCstats$callrate), 
                       "\n  max  call rate = ",max(negCstats$callrate),
                       "\n  mean depth = ",mean(negCstats$sampdepth),
                       "\n  max  depth = ",max(negCstats$sampdepth), "\n")
   samp.remove(uneg)
    }
   }
  if (have_rcpp) {
   sampdepth.max <- rcpp_rowMaximums(depth, nThreads)
  }
  else {
   sampdepth.max <- apply(depth, 1, max)
  }
  sampdepth <<- rowMeans(depth)
  u0 <- which(sampdepth.max == 0)
  u1 <- setdiff(which(sampdepth.max == 1 | sampdepth < sampdepth.thresh), u0)
  nmax0 <- length(u0)
  nmax1 <- length(u1)
  seqID.removed <<- character(0)
  if (nmax0 > 0) {
   cat(nmax0, "samples with no calls (maximum depth = 0) removed:\n")
   seqID.removed <<- seqID[u0]
   print(data.frame(indnum = u0, seqID = seqID.removed, sampdepth = sampdepth[u0]))
   }
  if (nmax1 > 0) {
   cat(nmax1, "samples with maximum depth of 1 and/or mean depth <", sampdepth.thresh, "removed:\n")
   print(data.frame(indnum = u1, seqID = seqID[u1], sampdepth = sampdepth[u1]))
   seqID.removed <<- c(seqID.removed, seqID[u1])
   }
  samp.remove(union(u0, u1))

 if (!gform == "chip") {
  if (!havedepth) {
   genon <<- alleles[, seq(1, 2 * nsnps - 1, 2)]/depth
   if(max(alleles)==Inf) {
    atemp <- alleles
    atemp[alleles==Inf] <- 1e6
    dtemp <- atemp[, seq(1, 2 * nsnps - 1, 2)] + atemp[, seq(2, 2 * nsnps, 2)]
    genon <<- atemp[, seq(1, 2 * nsnps - 1, 2)]/dtemp
    }
   genon <<-  trunc(2*genon-1)+1
   }
#  uhet <- which(!genon^2 == genon)
#  genon <<- 2*genon
#  genon[uhet] <<- 1
  if(outlevel > 7) {
   samples <<- genon
   samples[na.zero(genon) == 1] <<- 2* (sample.int(2, sum(genon == 1, na.rm=TRUE), replace = TRUE) - 1)  # allows genon to have > .Machine$integer.max elements
   }
#  rm(uhet)
  }
 gc()

  docalcp <- FALSE
  if (!exists("p") | length(seqID.removed) > 0) docalcp <- TRUE
  if (exists("p")) if (length(p) != nsnps) docalcp <- TRUE
  if (docalcp) { # not redone e.g. after merge
   p <<- calcp()
   if (is.null(p)) {
    cat("alleles not available, using genotype method for p\n")
    p <<- calcp(pmethod="G")
    }
   }
  if(exists("genosin")) rm(genosin)
  } #end GBS-specific
 write.csv(data.frame(seqID = seqID), "seqID.csv", row.names = FALSE)
 snpdepth <<- colMeans(depth)
 uremove <- which(p == 0 | p == 1 | is.nan(p) | snpdepth < snpdepth.thresh)
 nmaf0 <- length(which(p == 0 | p == 1))
 cat(nmaf0, "SNPs with MAF=0 and", length(uremove)-nmaf0, "SNPs with depth <", snpdepth.thresh, "removed\n")
 snp.remove(uremove)
cat("Analysing", nind, "individuals and", nsnps, "SNPs\n")

 ###### compare allele frequency estimates from allele counts and from genotype calls (& from input file, if uneak format)
 pg <<- colMeans(genon, na.rm = TRUE)/2  # allele freq assuming genotype calls
 if(outlevel > 4) {
  png("AlleleFreq.png", width = 960, height = 960, pointsize = cex.pointsize *  18)
   p.lab <- "Allele frequency from allele counts"
   pg.lab <- "Allele frequency from genotype calls"
   AF.lab <- "Allele frequency given"
   if (gform == "uneak" & outlevel > 6) pairs(cbind(pg, p, AFrq), col = "#80808020", pch = 16, cex = 0.8, labels = c(pg.lab, p.lab, AF.lab))
   if (gform != "uneak" | outlevel < 7) plot(pg ~ p, col="#80808020", pch=16, cex=0.8, xlab=p.lab, ylab=pg.lab)
   dev.off()
  }

 # calc some overall snp stats
 HWstats <- HWpops(depthmat = depth)
 HWdis <<- HWstats$HWdis[1,];  l10LRT <<- HWstats$l10LRT[1,]; x2star <<- HWstats$x2star[1,]; l10pstar <<- HWstats$l10pstar[1,]; maf <<- HWstats$maf[1,]
 LRT <- qchisq(-log(10)*l10LRT, 1, lower.tail = FALSE, log.p = TRUE)

 sampdepth <<- rowMeans(depth)  # recalc after removing SNPs and samples
 #if(outlevel > 4) sampdepth.med <<- apply(depth, 1, median)
 if(outlevel > 4) {
   if (have_rcpp) {
     sampdepth.med <<- rcpp_rowMedians(depth, nThreads)
   }
   else {
     sampdepth.med <<- apply(depth, 1, median)
   }
 }
 depth0 <- rowSums(depth == 0)
 snpdepth <<- colMeans(depth)
 missrate <- sum(as.numeric(depth == 0))/nrow(depth)/ncol(depth)
 cat("Proportion of missing genotypes: ", missrate, "Callrate:", 1-missrate,"\n")

 callrate <<- 1 - rowSums(depth == 0)/nsnps  # sample callrate, after removing SNPs, samples 
 SNPcallrate <<- 1 - colSums(depth == 0)/nind  
 png("CallRate.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
  hist(callrate, seq(0,1,0.02), col = "cornflowerblue", border = "cornflowerblue", main = "Histogram of sample call rates", xlab = "Call rate (proportion of SNPs scored)")
  dev.off()
 png("SNPCallRate.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
  # suggested by Jaroslav Klapste (Scion) 
  hist(SNPcallrate, seq(0,1,0.02), col = "cornflowerblue", border = "cornflowerblue", main = "Histogram of SNP call rates", xlab = "Call rate (proportion of samples scored)")
  dev.off()
 depthinfo <- TRUE; if(min(depth)==Inf) depthinfo <- FALSE
 if (!depthinfo) write.csv(data.frame(seqID, callrate), "SampleStats.csv", row.names = FALSE)
 if (depthinfo) {
  write.csv(data.frame(seqID, callrate, sampdepth), "SampleStats.csv", row.names = FALSE)
  sampdepth.scored <- sampdepth * nsnps/(nsnps - depth0)
  ufinite <- which(sampdepth < Inf)
  haschip <- (length(ufinite) < nind)
  noninf <- ""; if (haschip) noninf <- "(non-chip)"
  cat("Mean sample depth:", mean(sampdepth), "\n")
  if(haschip)   cat("Mean sample depth (non-chip):", mean(sampdepth[ufinite]), "\n")
  if(length(ufinite) > 0) {
   if(outlevel > 4) {
    png("SampDepth.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
     plot(sampdepth[ufinite] ~ sampdepth.med[ufinite], col = "#80808080", pch = 16, cex = 1.2, main = paste("Sample Depth",noninf), xlab = "Median", ylab = "Mean")
     dev.off()
    }
   png("SampDepth-scored.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
    plot(sampdepth.scored[ufinite] ~ sampdepth[ufinite], col = "#80808080", pch = 16, cex = 1.2, main = paste("Sample Depth",noninf), xlab = "Mean", ylab = "Mean with depth>0")
    dev.off()
   png("SampDepthHist.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
    hist(sampdepth[ufinite], 100, col = "cornflowerblue", border = "cornflowerblue", main = paste("Histogram of mean sample depth",noninf), xlab = "Mean sample depth")
    dev.off()
   png("SampDepthCR.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
    plot(sampdepth[ufinite] ~ callrate[ufinite], col = "#80808080", pch = 16, cex = 1.2, main = paste("Sample Depth v Call rate",noninf), xlab = "Sample call rate", ylab = "Mean sample depth")
    dev.off()
   }
  ufinite <- which(snpdepth < Inf)
  haschip <- (length(ufinite) < nsnps)
  noninf <- ""; if (haschip) noninf <- "(non-chip)"
  if(length(ufinite) >0 ) {
   png("SNPDepthHist.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
    hist(snpdepth[ufinite], 100, col = "cornflowerblue", border = "cornflowerblue", main = paste("Histogram of mean SNP depth",noninf), xlab = "Mean SNP depth")
    dev.off()
   png("SNPDepth.png", width = 640, height = 640, pointsize = cex.pointsize * 15)
    plot(SNPcallrate[ufinite] ~ snpdepth[ufinite], log="x", col = "#80808080", pch = 16, cex = 1, main = paste("SNP Depth",noninf), ylab = "SNP Call rate (proportion of samples scored)", 
         xlab = "Mean SNP depth (log scale)")
    dev.off()
   }
  }
 xaxpts <<- qchisq(QQprobpts,df=1)
 finplot()
 x2starplots()
 if(outlevel > 9) HWsigplot(plotname="HWdisMAFsig-raw")
 if(outlevel > 4) {
  yaxpts <- quantile(LRT,QQprobpts, na.rm=TRUE)
  png("LRT-QQ.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
   qqplot(qchisq(ppoints(nsnps), df = 1), LRT, cex=0.75, main = "Hardy-Weinberg LRT Q-Q Plot", 
         xlab = parse(text = "Theoretical ~~ (chi[1]^2) ~~  Quantiles"), ylab = "Sample Quantiles")
   if(length(QQprobpts) > 0) {
    for(ipt in 1:length(QQprobpts)) {
     lines(x=c(0,xaxpts[ipt],xaxpts[ipt]),y=c(yaxpts[ipt],yaxpts[ipt],0), col="grey")
     }
    text(x=xaxpts/4,y=yaxpts,labels=QQprobpts, pos=3, col="grey")
    }

   dev.off()
  png("LRT-hist.png", width = 480, height = 480, pointsize = cex.pointsize * 12)
   hist(LRT, breaks = 50, col = "grey", xlab = "Hardy Weinberg likelihood ratio test statistic")
   dev.off()
  }
 mafplot()
######### depth[depth < 2] <<- 1.1  # not using depth values <2 after this so set to >1 to avoid 0 divisor note do not use depth.max <2 though
 fcolo <<- rep("black", nind)  # modify this to specify colours for individuals in the plots
 }


palette.aquatic <- colorRampPalette(c(rgb(200, 200, 200, max = 255), "blue"))(50)  # grey to blue 
palette.terrain <- terrain.colors(50)
palette.temperature <- colorRampPalette(c("blue","white","red"))(50)
whitedist <- function(pal) min(colSums(abs(col2rgb(pal)-matrix(255,nrow=3,ncol=length(pal)))))
if(!functions.only) {
 readGBS()
 GBSsummary()
} # !functions.only


################## functions
upper.vec <- function(sqMatrix,diag=FALSE) as.vector(sqMatrix[upper.tri(sqMatrix,diag=diag)])
corner <- function(mtx, size=6L) print(mtx[1:size,1:size])
#seq2samp <- function(seqIDvec=seqID) read.table(text=seqIDvec,sep="_",stringsAsFactors=FALSE,fill=TRUE)[,1] # might not get number of cols right
seq2samp1 <- function(seqIDvec=seqID, splitby="_", ...) sapply(strsplit(seqIDvec,split=splitby, ...),"[",1)
seq2samp <- function(seqIDvec=seqID, splitby="_",nparts=NULL,dfout=FALSE, ...){
  if(is.null(nparts)) {
    sampout <- sapply(strsplit(seqIDvec,split=splitby, ...),"[",1)
  } else { 
    nparts <- as.integer(nparts)
    if(nparts < 2) nparts <- 2
    nsep <- nchar(gsub(paste0("[^",splitby,"]"),"",seqIDvec))
    seqIDvec1 <-  paste0(seqIDvec, strrep(splitby,pmax(0,nparts-1-nsep)))
    if(splitby %in% c(".","$","^","|","?","*","+","(",")")) splitby <- paste0("\\",splitby)
    stringpart <- paste(rep("(\\S*)",nparts),collapse=splitby)  # * means 0 or more times
    stringpart <- paste(rep("(.*)",nparts),collapse=splitby)
    seqIDexpr <- paste0("^",stringpart,"$")
    sampout <- gsub(seqIDexpr,"\\1",seqIDvec1)
    if(dfout) {
      sampout <- as.data.frame(sampout,stringsAsFactors=FALSE)
      for(ipart in 2:nparts) sampout <- cbind(sampout, gsub(seqIDexpr,paste0("\\",ipart),seqIDvec1),stringsAsFactors=FALSE)
      colnames(sampout) <- paste0("V",1:nparts)
    }
  }
  sampout
}


colourby <- function(colgroup, nbreaks=0, col.name=NULL, symbgroup=NULL, symb.name=NULL,groupsort=FALSE,maxlight=1,alpha=1,reverse=FALSE,symbset=NULL,
            hclpals=character(0),pal.upper=1, nacolour="black") {
 #maxlight is maximum lightness between 0 and 1
 #nbreaks is suggested # breaks
 isgroups <- (nbreaks==0)
 npals <- length(hclpals)
 if(isgroups) {
  collabels <- unique(colgroup)
  if(groupsort) collabels <- sort(collabels,na.last=TRUE)
  }
 if(!isgroups) {
  histinfo <- hist(colgroup,nbreaks,plot=FALSE)
  collabels <- histinfo$mids
  colgroup <- collabels[cut(colgroup,breaks = histinfo$breaks, labels=FALSE)]  # redefine input
  }
 ncol <- length(collabels)
 collist <- rainbow(ncol)
 if(ncol > 8) collist[seq(2,ncol,2)] <-  rgb((t(col2rgb(collist[seq(2,ncol,2)]))+matrix(127,ncol=3,nrow=floor(ncol/2)))/2,maxColorValue = 255)  # darken every 2nd one
 if(npals > 0 & exists("hcl.colors")) {
   for(ipal in 1:npals) {
    addcols <- hcl.colors(ceiling(ncol/(npals*pal.upper)),pal=hclpals[ipal])[1:ceiling(ncol/npals)]
    if(ipal==1) collist <- addcols else collist <- c(collist,addcols)
    }
   collist <- collist[1:ncol] # trim if needed
   }
 lightness <- colSums( col2rgb(collist))/(3*255)
 collist <- rgb(t(col2rgb(collist) %*% diag(pmin(lightness,rep(maxlight,length(lightness))) / lightness)), maxColorValue = 255)
 if(alpha < 1) {
  rgbcols <- col2rgb(collist)
  collist <- rgb(rgbcols[1,],rgbcols[2,],rgbcols[3,],alpha*255,maxColorValue=255)
  }
 if(reverse) collist <- rev(collist)
 sampcol <- collist[match(colgroup,collabels)]
 outlist <- list(collabels=collabels,collist=collist,sampcol=sampcol)
 if(!is.null(col.name)) outlist <- c(outlist,col.name=col.name)
 if (all(as.integer(symbset)==symbset)) {
  if(is.null(symbgroup)) {
   symblist <- suppressWarnings(symbset + rep(0,ncol))  # uses R recycle 
   sampsymb <- symblist[match(colgroup,collabels)]
   outlist <- c(outlist,list(symblist=symblist,sampsymb=sampsymb))
   } else { 
   symblabels <- unique(symbgroup)
   if(groupsort) symblabels <- sort(symblabels,na.last=TRUE)
   nsymb <- length(symblabels)
   if(length(symbset) < nsymb) symbset <- union(symbset,1:nsymb)[1:nsymb] # augment with unused from 1,2, ... 
   sampsymb <- symbset[match(symbgroup,symblabels)]
   outlist <- c(outlist,list(symblabels = symblabels, symblist=symbset,sampsymb=sampsymb))
   if(!is.null(symb.name)) outlist <- c(outlist,symb.name=symb.name)
   }
  }
 outlist$sampcol[is.na(outlist$sampcol)] <- nacolour
 outlist
 }

changecol <- function(colobject,colposition,newcolour) {# provide new colours to colourby object
 if(length(colposition)>0) {
  if(length(colposition) != length(newcolour))  stop("Error: colposition and newcolour are different lengths") 
  oldcolour <- colobject$collist
  colobject$collist[colposition] <- newcolour
  colobject$sampcol <-  colobject$collist[match(colobject$sampcol,oldcolour)] 
  }
 colobject
 }

coloursub <- function(colobj, indsubset) {
 #subset a colour object
 if(missing(indsubset)) indsubset <- 1:length(colobj$sampcol)
 colobj$sampcol <-  colobj$sampcol[indsubset]
 colobj$sampsymb <-  colobj$sampsymb[indsubset]
 ncol <- length(colobj$collist)
 uused <- which(colobj$collist %in% colobj$sampcol)
 colobj$collabels <-  colobj$collabels[uused]
 colobj$collist <-  colobj$collist[uused]
 if(is.null(colobj$symblabels)) {
  if(length(colobj$symblist)==ncol) colobj$symblist <- colobj$symblist[uused] 
  } else {
  uused <- which(colobj$symblist %in% colobj$sampsymb)
  colobj$symblabels <-  colobj$symblabels[uused]
  colobj$symblist <-  colobj$symblist[uused]
  }
 colobj
 }

colkey <- function(colobj, sfx="", srt.lab=0, plotch=16, horiz = TRUE, freq=FALSE) {  # plot a key for colours
 nlevels <- length(colobj$collabels)
 legname1 <- ""; if(!is.null(colobj$col.name)) legname1 <- colobj$col.name
 legname2 <- ""; if(!is.null(colobj$symb.name)) legname2 <- colobj$symb.name
 longdim <- 480 + max(0,nlevels-20) * 10
 labtext <- colobj$collabels
 temptab <- table(colobj$sampcol)
 dosymb <- !is.null(colobj$symblabels)  # if there are symblabels then plot separately
 if( length(colobj$symblist) == length(colobj$collist) & !dosymb) plotch <- colobj$symblist
 # possibly: check that colours match symbols 1-1: sum(with(colobj,table(sampcol,sampsymb))!=0) == length(colobj$collist) 
 if(freq) labtext <- paste(labtext, as.character(temptab[match(colobj$collist,names(temptab))]), sep="\t")
 if(dosymb) {
  symbtext=colobj$symblabels
  temptab <- table(colobj$sampsymb)
  if(freq) symbtext <- paste(symbtext, as.character(temptab[match(colobj$symblist,names(temptab))]), sep="\t")
  }
 if(horiz) {
  png(paste0("ColourKey",sfx,".png"),height=240, width = longdim)
  par(mar=0.1+c(1, 1, 1, 1))
  plot( 1:nlevels, rep(1,nlevels),col=colobj$collist,pch=plotch, cex=1.5,ylab="", axes=FALSE,xlab="", xpd=NA, ylim=c(0,1),xlim=c(0,nlevels))
  text(labels=c(legname1,labtext),x=(0:nlevels)+0.4,y=0.9,srt=srt.lab,pos=2,xpd=TRUE)
  if(dosymb) legend(x=nlevels/2,y=0.85,legend=symbtext,col="black",pch=colobj$symblist,ncol=length(colobj$symblist),
        bty="n", xjust=0.5, pt.cex=1.5, xpd=NA, title=legname2)
  } else {
  png(paste0("ColourKey",sfx,".png"),width=240, height = longdim)
  par(mar=0.1+c(1, 1, 1, 1))
  plot( rep(0,nlevels),nlevels:1, col=colobj$collist,pch=plotch, cex=1.5,ylab="", axes=FALSE,xlab="", xpd=NA, xlim=c(0,1),ylim=c(1,nlevels+1))
  text(labels=c(labtext[nlevels:1],legname1),y=1:(nlevels+1),x=0.1,pos=4,srt=srt.lab,xpd=TRUE)
  if(dosymb) legend(y=nlevels/2, x=0.2 + max(strwidth(labtext)),legend=symbtext,col="black",pch=colobj$symblist,
        bty="n", yjust=0.5, pt.cex=1.5, xpd=NA, title=legname2)
  }
 dev.off()
 }

collegend <- function(colobj,legpos="topleft", plotx=NULL, ploty=NULL, cex.leg=0.9) {
     # being tested
     legendsym <- 1  #default value
     if (length(unique(colobj$symblist))==1) legendsym <- colobj$symblist[1] else
     if(length(colobj$symblist) == length(colobj$collist)) if(
         length(colobj$collist) == qr(model.matrix(~ colobj$sampcol +  as.factor(colobj$sampsymb)))$rank) legendsym <- colobj$symblist
     legname1 <- ""; if(!is.null(colobj$col.name)) legname1 <- colobj$col.name
     legname2 <- ""; if(!is.null(colobj$symb.name)) legname2 <- colobj$symb.name
     ncolumns <- 1; if(length(colobj$collist) > 6 & !is.numeric(colobj$collabels)) ncolumns <- 2
     ncolumns2 <- 1; if(length(colobj$symblist) > 6) ncolumns2 <- 2
     if(!is.numeric(colobj$collabels)) leginfo <- legend(legpos, legend=colobj$collabels, col=colobj$collist, pch=legendsym, ncol=ncolumns, cex=cex.leg, title=legname1)$rect
     if(is.numeric(colobj$collabels)) { #raster instead of legend
      leginfo <- legend(legpos,legend=rep("test",3),plot=FALSE,title=legname1)$rect # find where a 3 element legend would go
      topadj <- 1.5 * strheight(legname1)
      leginfo$top <- leginfo$top - topadj; leginfo$h <- leginfo$h - topadj  #find legend placement adj for title
      legend_image <- as.raster(matrix(rev(colobj$collist), ncol = 1))
      rasterImage(legend_image, leginfo$left+0.05*leginfo$w, leginfo$top-0.95*leginfo$h, 
                                 leginfo$left+0.55*leginfo$w, leginfo$top-0.05*leginfo$h)
      rticks <- pretty(colobj$collabels)
      rticks <- rticks[2:(length(rticks)-1)]
      nbars <- length(colobj$collist)
      rasth <- 0.9*leginfo$h
      rast0 <- leginfo$top-0.95*leginfo$h+0.5*rasth/nbars
      rasth <- rasth*(nbars-1)/nbars  # height between end mid-pints
      lines(x=rep(leginfo$left+0.6*leginfo$w,2),y=c(rast0,rast0+rasth))
      rtickss <- (rticks - min(colobj$collabels))/diff(range(colobj$collabels))
      for(itick in seq_along(rticks)){
       lines(x=leginfo$left+c(0.6,0.63)*leginfo$w,y=rep(rast0+rtickss[itick]*rasth,2))
       text(x=leginfo$left+0.65*leginfo$w, y= rast0+rtickss[itick]*rasth, labels=format(rticks)[itick],adj=c(0,0.5), cex=0.7)
       }
      text(x=leginfo$left+0.5*leginfo$w,y=leginfo$top, labels=legname1, xpd=NA, adj=c(0.5,0), cex=cex.leg)
      } #raster
    if(!is.null(plotx) & !is.null(ploty)) { #check that plot points not covered
     ncovered <- length(which(plotx > leginfo$left & plotx < leginfo$left + leginfo$w & ploty < leginfo$top & ploty > leginfo$top - leginfo$h))
     if(ncovered > 0) cat("Warning: ",ncovered,"data points covered by colour legend\n")
     }
    if(!is.null(colobj$symblabels)) {
     leginfo2 <- legend(legpos, legend=colobj$symblabels, ncol=ncolumns2, cex=cex.leg, title=legname2, plot=FALSE)$rect
     x2 <- leginfo$left - leginfo2$w *1.05 
     y2 <- min(leginfo$top,max(ploty,na.rm=TRUE))
     if(!is.null(plotx) & !is.null(ploty)) { #direct 2nd legend to sparse location beside 1st legend & check that plot points not covered
      if(leginfo$left < mean(plotx,na.rm=TRUE)) x2 <- leginfo$left + leginfo$w *1.05
      if(leginfo$top < mean(ploty,na.rm=TRUE)) y2 <- min(ploty,na.rm=TRUE) + leginfo2$h
      ncovered <- length(which(plotx > leginfo2$left & plotx < leginfo2$left + leginfo2$w & ploty < leginfo2$top & ploty > leginfo2$top - leginfo2$h))
      if(ncovered > 0) cat("Warning: ",ncovered,"data points covered by symbol legend\n")
      }
     leginfo2 <- legend(x=x2, y=y2, legend=colobj$symblabels, col="black", pch=colobj$symblist, ncol=ncolumns2, cex=cex.leg, title=legname2)$rect
     }
 } #collegend


### calculate expected identity mismatch rate given observed geno & depths of indiv1 and depths of indiv2
mismatch.ident <- function(seqid1, seqid2,snpsubset=1:nsnps, puse=p, mindepth.mm=1) {
  if(mindepth.mm < 1) mindepth.mm <- 1
  pos1 <- match(seqid1, seqID)
  pos2 <- match(seqid2, seqID)
  depthi <- depth[pos1, ]
  depthj <- depth[pos2, ]
  usnp <- intersect(snpsubset,which(depthi >= mindepth.mm & depthj >= mindepth.mm))
  pi <- genon[pos1, usnp]/2
  pj <- genon[pos2, usnp]/2
  K1  <-  depth2K(depthi[usnp])
  K2  <-  depth2K(depthj[usnp])
  mismatched <- (pi != pj)
  #nmismatch <- length(which(pi != pj))
  nmismatch <- sum(mismatched)
  ncompare <- length(usnp)
  ptemp <- puse[usnp]
  P <- ptemp*(1-ptemp)
  expmm <- rep(NA,length(usnp))
  ug <- which(pi==1)
  if(length(ug)>0) expmm[ug] <- 2 * P[ug] * K1[ug] * (1 - K2[ug]) / (ptemp[ug]^2+2*P[ug]*K1[ug])
  ug <- which(pi==0.5)
  if(length(ug)>0) expmm[ug] <- 2 * K2[ug] 
  ug <- which(pi==0)
  if(length(ug)>0) expmm[ug] <- 2 * P[ug] * (1 - K1[ug]) * K2[ug] / ((1-ptemp[ug])^2+2*P[ug]*K1[ug])
  exp.mmrate <- mean(expmm,na.rm=TRUE)
  mmrate <- nmismatch/ncompare
  snpmmstats <- data.frame(mm=logical(nsnps),expmm=numeric(nsnps)); snpmmstats[,] <- NA
  snpmmstats[usnp,1] <- mismatched; snpmmstats[usnp,2] <- pmin(expmm,rep(1,ncompare))  # why can expmm be >1 ?
  list(mmrate=mmrate,ncompare=ncompare,exp.mmrate=exp.mmrate,snpmmstats=snpmmstats)
}


posCreport <- function(mergeIDs,Guse,sfx = "",indsubset,Gindsubset,snpsubset=1:nsnps,puse=p, snpreport=FALSE, quiet=FALSE) {
 if (missing(Gindsubset)) Gindsubset <- 1:nind
 if (missing(indsubset)) indsubset <- 1:nrow(Guse)
 havedepth <- exists("sampdepth")
 mergeIDs <- mergeIDs[indsubset]
 csvout <- paste0("posCreport", sfx, ".csv")
 cat("Positive Control Checks (also see", csvout, ")\n")
 seqIDtemp <- seqID[Gindsubset][indsubset]
 multiIDs <- unique(mergeIDs[which(duplicated(mergeIDs))])
 if(havedepth) depthsub <- depth[Gindsubset,,drop=FALSE]
 posCstats <- data.frame(mergeID=character(0),nresults=integer(0),selfrel=numeric(0),meanrel=numeric(0),minrel=numeric(0),
              meandepth=numeric(0),mindepth=numeric(0),meanCR=numeric(0),mmrate=numeric(0),exp.mmrate=numeric(0),EMM=numeric(0),stringsAsFactors=FALSE)
 meandepth <- mindepth <- meanCR <- mmrate <- exp.mmrate <- IEMM <- NA
 sink(paste0("posCchecks",sfx,".txt"),split=!quiet)
 snpstats <- data.frame(snpcount = integer(nsnps), snpmm=integer(nsnps), expmm=numeric(nsnps))
 for (i in seq_along(multiIDs)) {
  thisID <- multiIDs[i]
  thispos <- which(mergeIDs==thisID)
  nresults <- length(thispos)
  thisdepth <- rep(NA,nresults)
  thisG <- Guse[indsubset,indsubset][thispos,thispos]
  selfrel <- mean(diag(thisG))
  meanrel <- mean(upper.vec(thisG))
  minrel <- min(upper.vec(thisG))
  if(havedepth) {
   thisdepth <- sampdepth[Gindsubset][indsubset][thispos]
   meandepth <- mean(thisdepth)
   mindepth <- min(thisdepth)
   meanCR <- mean( 1 - rowSums(depthsub[thispos,,drop=FALSE] == 0)/nsnps )
   }
  #emmstats <- data.frame(mmrate=numeric(0),exp.mmrate=numeric(0))
  mmmat <- expmmmat <- matrix(NA,nrow=nresults,ncol=nresults)
  if(length(snpsubset) > 0) {
   for(ipos in 1:nresults) for (jpos in setdiff(1:nresults,ipos)) {
    emmstatstemp <- mismatch.ident(seqIDtemp[thispos[ipos]],seqIDtemp[thispos[jpos]],snpsubset=snpsubset,puse=puse)  # do it both ways
    mmmat[ipos,jpos] <-emmstatstemp$mmrate
    expmmmat[ipos,jpos] <- emmstatstemp$exp.mmrate
    ucompare <- which(!is.na(emmstatstemp$snpmmstats$mm))
    snpstats$snpcount[ucompare] <- snpstats$snpcount[ucompare] + 1
    snpstats$snpmm[ucompare] <- snpstats$snpmm[ucompare] + emmstatstemp$snpmmstats$mm[ucompare]
    snpstats$expmm[ucompare] <- snpstats$expmm[ucompare] + emmstatstemp$snpmmstats$expmm[ucompare]
    }
   expmmmat <- (expmmmat+t(expmmmat))/2
   mmrate <- mean(mmmat,na.rm=TRUE)
   exp.mmrate <- mean(expmmmat,na.rm=TRUE)
   IEMM <- mmrate-exp.mmrate
   }
  posCstats <- rbind(posCstats, data.frame(mergeID=thisID,nresults=nresults,selfrel=selfrel,meanrel=meanrel,minrel=minrel,
                     meandepth=meandepth,mindepth=mindepth,meanCR=meanCR,mmrate=mmrate,exp.mmrate =exp.mmrate, IEMM = IEMM, stringsAsFactors=FALSE))
  ulorel <- which( ((thisG < 1 & thisG - selfrel <= hirel.thresh - 1) | mmmat-expmmmat>iemm.thresh) & upper.tri(thisG), arr.ind = TRUE)
  if (nrow(ulorel) > 0) print(data.frame(Indiv1 = seqIDtemp[thispos[ulorel[, 1]]], Indiv2 = seqIDtemp[thispos[ulorel[, 2]]], rel = thisG[ulorel],
                        depth1 = thisdepth[ulorel[, 1]],depth2 = thisdepth[ulorel[, 2]], IEMM = (mmmat-expmmmat)[ulorel]))
  }
 sink()
 snpstats$mmrate <- snpstats$snpmm / snpstats$snpcount 
 snpstats$expmm <- snpstats$expmm / snpstats$snpcount 
 write.csv(posCstats,file=csvout,row.names=FALSE,quote=FALSE)
 if(nrow(posCstats) > 1) {
  png(paste0("SelfRel",sfx,".png"), width = 960, height = 960, pointsize = cex.pointsize *  18)
   with(posCstats,plot(meanrel~selfrel,xlab="Mean within run",ylab="Mean between run",main="Self-relatedness"))
   abline(a=0,b=1,col="red",lwd=2)
   lines(x=c(min(posCstats$selfrel), 2-hirel.thresh,max(max(posCstats$selfrel),2-hirel.thresh)),y=c(min(posCstats$selfrel)+hirel.thresh-1,1,1),col="grey",lwd=2)
   dev.off()
  if(length(snpsubset) > 0) {
   png(paste0("posC-MM", sfx, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
    with(posCstats,plot(mmrate ~ exp.mmrate, xlab = "Expected mismatch rate", ylab = "Raw mismatch rate", cex=0.8))
    abline(a=0,b=1,col="red",lwd=2)
    abline(a=iemm.thresh,b=1,col="grey",lwd=2)
    dev.off()
   png(paste0("posC-EMM", sfx, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
    with(posCstats,pairs(cbind(meanrel,selfrel,IEMM)))
    dev.off()
   if(snpreport) {
    depthtrans <- function(x) round(20 * log(-log(1/(x + 0.9)) + 1.05))  # to compress colour scale at higher depths
    depthpoints <- c(0.5, 5, 50, 250)  # legend points
    transpoints <- depthtrans(depthpoints)
    mindepthplot <- 0.1
    maxdepthplot <- 256
    maxtrans <- depthtrans(maxdepthplot)
    png(paste0("posC-SNPMM", sfx, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
     with(snpstats,plot(mmrate ~ expmm, col = palette.aquatic[depthtrans(pmax(mindepthplot, pmin(snpdepth, maxdepthplot)))],
                        xlab = "Expected SNP mismatch rate", ylab = "Raw SNP mismatch rate", cex=0.8, 
                        main="SNP mismatch rate vs expected mismatch rate"))
     title(main="using default finplot colours", line=1, cex.main=0.75)
     abline(a=0,b=1,col="red",lwd=2)
     dev.off()
    }
   }
  }
 if(snpreport) list(posCstats=posCstats, snpstats=snpstats) else posCstats
 }

mergeSamples <- function(mergeIDs, indsubset, replace=FALSE) {  
 # doesn't do samples0, so cant do G3 
 if (missing(indsubset)) indsubset <- 1:nind
 mergeIDs <- mergeIDs[indsubset]
 hasg <- exists("genon")
 if(hasg) {
  aggr.msum <- rowsum(genon[indsubset,,drop=FALSE],mergeIDs,na.rm=TRUE)   # rowsum very fast
  temp <- 1 * !is.na(genon[indsubset,,drop=FALSE])
  aggr.mn <- rowsum(temp,mergeIDs) 
# aggr.mn <- rowsum(1 * !is.na(genon[indsubset,,drop=FALSE]),mergeIDs) 
  rm(temp)
  genon.m <- aggr.msum/aggr.mn
  genon.m <- trunc(genon.m-1)+1
  depth.m <- rowsum(depth[indsubset,,drop=FALSE],mergeIDs,na.rm=TRUE)
  ID.m <- rownames(aggr.msum)
  sampdepth.m <- rowMeans(depth.m)
  snpdepth.m <- colMeans(depth.m)
  pg.m <- colMeans(genon.m, na.rm = TRUE)/2  # allele freq assuming genotype calls
  }
 if(alleles.keep | !hasg) alleles.m <- rowsum(alleles[indsubset,,drop=FALSE],mergeIDs,na.rm=TRUE) else alleles.m <- NULL
 if(!hasg) ID.m <- rownames(alleles.m)
 nind.m <- length(ID.m)
 if(exists("SNPcallrate") & hasg) { #GBSsummary has been run
  callrate.m <- 1 - rowSums(depth.m == 0)/nsnps  # sample callrate, after removing SNPs, samples 
  SNPcallrate.m <- 1 - colSums(depth.m == 0)/nind.m  
  }
 nseq <- rowsum(rep(1,length(indsubset)),mergeIDs)[,1] # results being merged
 seqID.m <- seqID[indsubset][match(ID.m,mergeIDs)]
 seqinfo <- read.table(text=seqID[indsubset][match(ID.m,mergeIDs)],sep="_",fill=TRUE,stringsAsFactors=FALSE)
 if(ncol(seqinfo)==5) { #Assume formated as ID_Flowcell_Lane_plate_X and return ID_merged_nsamples_0_X
  umerged <- which(nseq>1)
  seqinfo[umerged,2] <- "merged"
  seqinfo[umerged,3] <- nseq[umerged]
  seqinfo[umerged,4] <- 0
  seqID.m <- paste(seqinfo[,1],seqinfo[,2],seqinfo[,3],seqinfo[,4],seqinfo[,5],sep="_")
  }
 if(hasg) mergelist <- list(mergeIDs=ID.m, nind=nind.m, seqID=seqID.m, genon=genon.m, depth = depth.m, sampdepth=sampdepth.m, snpdepth=snpdepth.m, pg=pg.m, nmerged=nseq)
 if(alleles.keep & hasg) mergelist <- list(mergeIDs=ID.m, nind=nind.m, seqID=seqID.m, genon=genon.m, depth = depth.m, alleles=alleles.m, sampdepth=sampdepth.m, snpdepth=snpdepth.m, pg=pg.m, nmerged=nseq)
 if(!hasg) mergelist <- list(mergeIDs=ID.m, nind=nind.m, seqID=seqID.m, alleles=alleles.m, nmerged=nseq)
 if(exists("SNPcallrate") & hasg) mergelist <- c(mergelist,list(callrate=callrate.m, SNPcallrate=SNPcallrate.m))
 if(replace) {
  nind <<- nind.m; seqID <<- seqID.m; mergelist$nind <- NULL; mergelist$seqID <- NULL
  if(hasg) {
    genon <<- genon.m; depth <<- depth.m; sampdepth <<- sampdepth.m; snpdepth <<- snpdepth.m; pg <<- pg.m; 
    mergelist$genon <- mergelist$depth <- mergelist$sampdepth <- mergelist$snpdepth <- mergelist$pg <- NULL
    }
  if(exists("SNPcallrate") & hasg) { callrate <<- callrate.m; SNPcallrate <<- SNPcallrate.m; mergelist$callrate <- mergelist$SNPcallrate <- NULL }
  if(alleles.keep | !hasg) {alleles <<- alleles.m; mergelist$alleles <- NULL}
  }
 mergelist
 }

mergeSamples2 <- function(mergeIDs, indsubset) {  
 # this version to merge duplicates only ...
 # doesn't do samples0, so cant do G3 
 if (missing(indsubset)) indsubset <- 1:nind
 mergeIDs <- mergeIDs[indsubset]
 nseq <- rowsum(rep(1,length(indsubset)),mergeIDs) # results being merged
 singleIDs <- row.names(nseq)[nseq==1]
 usingle <- which(mergeIDs %in% singleIDs)
 umultiple <- setdiff(1:length(indsubset),usingle)
 aggr.msum <- rowsum(genon[indsubset[umultiple],,drop=FALSE],mergeIDs[umultiple],na.rm=TRUE)   # rowsum very fast
 temp <- 1 * !is.na(genon[indsubset[umultiple],,drop=FALSE])
 aggr.mn <- rowsum(temp,mergeIDs[umultiple]) 
# aggr.mn <- rowsum(1 * !is.na(genon[indsubset,,drop=FALSE]),mergeIDs) 
 rm(temp)
 genon.m <- aggr.msum/aggr.mn
 genon.m <- trunc(genon.m-1)+1
 ID.m <- rownames(aggr.msum)
 depth.m <- rowsum(depth[indsubset[umultiple],,drop=FALSE],mergeIDs[umultiple],na.rm=TRUE)
 nind.m <- nrow(genon.m)
 seqID.m <- seqID[indsubset[umultiple]][match(ID.m,mergeIDs[umultiple])]
 seqinfo <- read.table(text=seqID.m ,sep="_",fill=TRUE,stringsAsFactors=FALSE)
 if(ncol(seqinfo)==5) { #Assume formated as ID_Flowcell_Lane_plate_X and return ID_merged_nsamples_0_X
  seqinfo[,2] <- "merged"
  seqinfo[,3] <- nseq[nseq>1]
  seqinfo[,4] <- 0
  seqID.m <- paste(seqinfo[,1],seqinfo[,2],seqinfo[,3],seqinfo[,4],seqinfo[,5],sep="_")
  }
 if(length(usingle) > 0) {
  genon.m <- rbind(genon.m,genon[indsubset[usingle],,drop=FALSE])
  ID.m <- c(ID.m,mergeIDs[usingle])
  depth.m <- rbind(depth.m,depth[indsubset[usingle],,drop=FALSE])
  nind.m <- nind.m + length(usingle)
  seqID.m <- c(seqID.m,seqID[indsubset[usingle]])
  }
 sampdepth.m <- rowMeans(depth.m)
 snpdepth.m <- colMeans(depth.m)
 pg.m <- colMeans(genon.m, na.rm = TRUE)/2  # allele freq assuming genotype calls
 list(mergeIDs=ID.m, nind=nind.m, seqID=seqID.m, genon=genon.m, depth = depth.m, sampdepth=sampdepth.m, snpdepth=snpdepth.m, pg=pg.m, nmerged=nseq)
 }

plateplot <- function(plateinfo,plotvar=sampdepth,vardesc="", sfx="", neginfo,negcol="grey", colpal = hcl.colors(80,palette="YlGnBu", rev=TRUE)) {
 if(vardesc=="") vardesc <- as.character(substitute(plotvar))
 infonames <- names(plateinfo)
 if(!"row"  %in% infonames | ! "column" %in% infonames) stop("row and/or column not in plateinfo")
 if(nrow(plateinfo) != nind) stop (paste("Number of rows in plateinfo",nrow(plateinfo),"does not match # individuals",nind))
 if(missing(colpal) & !exists("hcl.colors")) colpal <- colorRampPalette(c("#FCFFDD","#BCE9C5","#18BDB0","#007EB3","#26185F"))(80)  # similar to YlGnBu
 colpalgray <- rev(gray.colors(60, end=0.975))  # keep colours away from white,

 markneg <- function(negcol=negcol) {
  if(nrow(neginfo) > 0) {
   for (ineg in 1:nrow(neginfo)) {
    thisrow <- which(rownames(zmat)==neginfo$row[ineg])
    thiscol <- which(colnames(zmat)==neginfo$column[ineg])
    if(length(thisrow) == 1 & length(thiscol)==1) { # only for rows and cols in the main plate plot.
     lines(y=c(thisrow-3/2, thisrow-1/2)/(maxrow-1),x=c(thiscol-3/2,thiscol-1/2)/(maxcol-1),col=negcol)
     lines(y=c(thisrow-1/2, thisrow-3/2)/(maxrow-1),x=c(thiscol-3/2,thiscol-1/2)/(maxcol-1),col=negcol)
     }
    }
   }
  }

 addlegend <- function(barmax=0.7,plotvar,main="",colpalette=colpal, subid=NULL, subcol=NULL) {
  ymax <- max(1,(maxrow-0.5)/(maxrow-1))  # 1 to cope with nrow=1
  parplt <- par("plt")
  paryratio <- (parplt[4]-parplt[3])/(parplt[2]-parplt[1])
  paryratio <- paryratio* max(2,maxcol)*max(1,(maxrow-1))/(max(2,maxrow)*max(1,(maxcol-1)))   # adjust for image x & y dimensions
  rasterImage(as.raster(matrix(colpalette, ncol = 1)), barmax, ymax+ 0.01, barmax+0.05, ymax+ 0.41, angle=90, xpd=TRUE)
  varstats <- summary(plotvar)
  varmin <- varstats[1]; varmax = varstats[6]; varmean <- varstats[4]; varrange <- varmax-varmin
  xlow <- barmax - 0.4*paryratio
  xmean <- xlow +  0.4*paryratio* (varmean-varmin)/varrange
  barheight <- 0.05/paryratio
  bartop <- ymax + 0.01 + barheight
  lines(x=c(barmax,barmax),y=c(bartop,ymax + 0.08),xpd=TRUE)
  lines(x=c(xlow,xlow),y=c(bartop,ymax + 0.08),xpd=TRUE)
  lines(x=c(xmean,xmean),y=c(bartop,ymax + 0.08),xpd=TRUE)
  text(x = xlow+.02, y = ymax + 0.12, pos=2, labels = paste("Min =",signif(varmin,3)),xpd=TRUE, cex=0.8)
  text(x = barmax-.05, y = ymax + 0.12, pos=4, labels = paste("Max =",signif(varmax,3)),xpd=TRUE, cex=0.8)
  text(x = max(xlow+0.07,xmean), y = ymax + 0.12  , pos=NULL, labels = paste("Mean =",signif(varmean,3)), xpd=TRUE, cex=0.8)
  text(x = -0.5/(max(1.5,maxcol)-1), y = ymax + 0.03, pos=4, labels = main, cex = 1, xpd=TRUE)  # place from left margin
  nsub <- length(subid)
  if(!is.null(subid)) for(isub in seq_along(subid))  {
   rect(xlow,  ymax + 0.01 + barheight*(isub-1)/nsub, barmax, ymax + 0.01 + barheight*isub/nsub, col = subcol[isub], border=NA, xpd=TRUE)
   text(x=barmax,y=ymax + 0.01 + barheight*(isub-1)/nsub, labels=subid[isub],cex=0.5, adj=c(-0.2,0), xpd=TRUE)
   if(isub==1) text(x=barmax,y=bartop, labels="Subplate",cex=0.7, adj=c(-0.2,0), xpd=TRUE)
   }
  }

 zmat <- xtabs( plotvar ~ plateinfo$row + plateinfo$column)  # sums (if needed) plotvar
 zmat[which(zmat==0)] <- NA
 maxrow <- nrow(zmat); maxcol <- ncol(zmat)

 png(paste0("Plate",sfx,".png"), pointsize = cex.pointsize *  12)
  image(t(zmat), axes=FALSE,xlab="Plate column",ylab="Plate row",col=colpal,cex.lab=1.2)  
  axis(1,at=(0:(ncol(zmat)-1))/(max(2,ncol(zmat))-1),labels=colnames(zmat),cex.axis=0.8)
  axis(2,at=(0:(nrow(zmat)-1))/(max(2,nrow(zmat))-1),labels=rownames(zmat),cex.axis=0.8)
  if(!missing(neginfo)) markneg(negcol=negcol)
  addlegend(plotvar=as.numeric(zmat), main=vardesc)
  dev.off()

 if("subplate" %in% infonames) {
  subplateids <- sort(unique(plateinfo$subplate))
  subplatecols <- rainbow(length(subplateids),alpha=0.2)
  zmatc <- xtabs( match(plateinfo$subplate,subplateids) ~ plateinfo$row + plateinfo$column) 
  zmatc[which(is.na(zmat))] <- NA
  png(paste0("Subplate",sfx,".png"), pointsize = cex.pointsize *  12)
   image(t(zmat), axes=FALSE,xlab="Plate column",ylab="Plate row",col=colpalgray,cex.lab=1.2)
   image(t(zmatc), axes=FALSE,col=subplatecols, add=TRUE) 
   axis(1,at=(0:(ncol(zmat)-1))/(max(2,ncol(zmat))-1),labels=colnames(zmat),cex.axis=0.8)
   axis(2,at=(0:(nrow(zmat)-1))/(max(2,nrow(zmat))-1),labels=rownames(zmat),cex.axis=0.8)
   if(!missing(neginfo)) markneg(negcol=negcol)
   addlegend(plotvar=plotvar, main=vardesc, colpalette=colpalgray, subid=subplateids, subcol=subplatecols)
   dev.off()
  }
 invisible(NULL)
} #plateplot

calcG <- function(snpsubset, sfx = "", puse, indsubset, depth.min = 0, depth.max = Inf, npc = 0, calclevel = 9, cocall.thresh = 0, mdsplot=FALSE,
                  mindepth.idr = 0.1, withPlotly=FALSE, withHeatmaply=withPlotly, plotly.group=NULL, plotly.group2=NULL, samp.info=NULL,samptype="diploid") {
  # example calls:
  #   Gfull <- calcG(npc=4) 
  #   GHWdgm.05 <- calcG(which(HWdis > -0.05),'HWdgm.05') # recalculate using Hardy-Weinberg disequilibrium cut-off at -0.05
  # sfx is text to add to GGBS5 as graph name, puse is allele freqs (all snps) to use
  # calclevel: 1: G5 only, 2: G5 + reports using G5, 3: all G, 9: all
  cat("Calculating G matrix, analysis code:", sfx, "\n")
  if (missing(snpsubset))   snpsubset <- 1:nsnps
  if (missing(indsubset))   indsubset <- 1:nind
  if (missing(puse))        puse <- p
  d4i <- 1.001  # min depth to use in inbreeding calcs. Anything in (1,2] OK for normal data, but a low value used in case effective depth (non-integer) being used

  Gpool <- NULL
  samptype <- tolower(samptype)
  if (!samptype %in% c("diploid","pooled")) samptype <- "diploid"
  if (samptype=="pooled") if (!exists("alleles")) {
    cat("Allele matrix does not exist. Using diploid samptype\n")
    samptype <- "diploid"
    }
  ## Some checks if using plotly
  if(withPlotly){
    if(!require(plotly))
      withPlotly = FALSE
    else{
      # Some checks on the plotly group vectors
      if(!is.vector(plotly.group)){
        warning("The first group object for use in plotly is not a vector.\nNo groups will be ploted in plotly.")
        plotly.group = NULL
      }
      else if(length(plotly.group) != length(indsubset)){
        warning("The first group vector for use in plotly does not equal the number of individuals.\nNo groups will be ploted in plotly.")
        plotly.group = NULL
      }
      if(!is.vector(plotly.group2)){
        warning("The second group object for use in plotly is not a vector.\nNo groups will be ploted in plotly.")
        plotly.group2 = NULL
      }
      else if(length(plotly.group2) != length(indsubset)){
        warning("The second group vector for use in plotly does not equal the number of individuals.\nNo groups will be ploted in plotly.")
        plotly.group2 = NULL
      }
      addpixel = ifelse(any(!is.null(c(plotly.group,plotly.group2))),200,0) 
      # Some checks on the hove information for plotly
      if(missing(samp.info)){
        if(!exists(seqID))
          samp.info <- list("Sample ID"=1:length(indsubset))
        else
          samp.info <- list("Sample ID"=seqID[indsubset])
      } else if(!is.list(samp.info)){
        warning("Sample information object is not a list")
        samp.info <- list("Sample ID"=seqID[indsubset])
      } else if(any(unlist(lapply(samp.info,function(x) length(x)!=length(indsubset))))){
        stop("Missing or extra sample information. Check that the 'samp.info' object is correct.")
      }
      hover.info <- apply(sapply(1:length(samp.info),function(x) paste0(names(samp.info)[x],": ",samp.info[[x]]),simplify = TRUE),1,paste0,collapse="<br>")
    }
  }
  if(withHeatmaply){
    if(!require(heatmaply) || compareVersion(as.character(packageVersion("heatmaply")), "0.15.0")==-1)
      withHeatmaply = FALSE
    if(!exists("hover.info"))
      hover.info <- apply(sapply(1:length(samp.info),function(x) paste0(names(samp.info)[x],": ",samp.info[[x]]),simplify = TRUE),1,paste0,collapse="<br>")
  }
  nsnpsub <- length(snpsubset)
  nindsub <- length(indsubset)
  # puse - determine whether global or indiv specifc
  if(is.matrix(puse)) {  # allows matrix puse to be specified relative to all or subsetted data
   cat("Using individual allele frequencies\n")
   if(ncol(puse) == nsnps & nsnps > nsnpsub) puse <- puse[,snpsubset]
   if(nrow(puse) == nind & nind > nindsub) puse <- puse[indsubset,]
   } else {
   cat("Using global allele frequencies\n")
   if(length(puse) == nsnps & nsnps > nsnpsub) puse <- puse[snpsubset]
   puse <- matrix(puse,byrow=TRUE,nrow=nindsub,ncol=nsnpsub)
   }
  if(!nrow(puse) == nindsub | !ncol(puse) == nsnpsub) stop("Dimensions of puse are incorrect.")
  #check for NA allele freqs
  pusecheck <- colSums(puse)
  pusena <- which(is.na(pusecheck))
  npusena <- length(pusena)
  if(npusena > 0 ) {
   cat("Warning:", npusena,"SNPs with NA allele frequencies were removed\n")
   snpsubset <- snpsubset[-pusena]
   pusecheck <- pusecheck[-pusena]
   puse <- puse[,-pusena]
   }
  if(any(pusecheck == 0 | pusecheck == 1 ) ) cat("Warning: Some MAFs are 0\n")
  depthsub <- depth[indsubset, snpsubset]
  
  cat("# SNPs: ", nsnpsub, "\n")
  cat("# individuals: ", nindsub, "\n")
  genon0 <- genon[indsubset, snpsubset]
  usegeno <- !is.na(genon[indsubset, snpsubset])
  if(samptype=="pooled") raf <- 2 * alleles[indsubset, seq(1, 2 * nsnps - 1, 2)][,snpsubset] / depthsub  # ref allele freq x 2 (keep on 0-2 scale)
  if (depth.min > 1 | depth.max < Inf) {
    depthsub[depthsub < depth.min] <- 0
    depthsub[depthsub > depth.max] <- 0
    genon0[depthsub==0] <- NA
    usegeno[depthsub==0] <- FALSE
    if(samptype=="pooled") raf[depthsub==0] <- NA
    }
  cocall <- tcrossprod(usegeno)
  cat("Mean co-call rate (for sample pairs):", mean(upper.vec(cocall)/nsnpsub), "\n")
  cat("Min  co-call rate (for sample pairs):", min(upper.vec(cocall)/nsnpsub), "\n")
  png(paste0("Co-call-", sfx, ".png"), pointsize = cex.pointsize * 12)
   hist(upper.vec(cocall)/nsnpsub, breaks = 50, xlab = "Co-call rate (for sample pairs)", main="", col = "grey")
   dev.off()
  lowpairs <- which(cocall/nsnpsub <= cocall.thresh & upper.tri(cocall),arr.ind=TRUE)
  if (have_rcpp) {
    sampdepth.max <- rcpp_rowMaximums(depthsub, nThreads)
  }
  else {
    sampdepth.max <- apply(depthsub, 1, max)
  }
  samp.removed <- NULL
  if(cocall.thresh >= 0) {  # remove samples which wont get self-rel (unless -ve cocall.thresh)
   samp.removed <- which(sampdepth.max < 2)
   ulow <- which(lowpairs[,1] %in% samp.removed | lowpairs[,2] %in% samp.removed)
   if(length(ulow) > 0) lowpairs <- lowpairs[-ulow,,drop=FALSE]
   }   
  while(nrow(lowpairs) > 0) {
   lowsamptab <- table(as.vector(lowpairs))
   lowsamp <- as.numeric(names(which.max(lowsamptab)))
   samp.removed <- c(samp.removed,lowsamp)
   lowpairs <- lowpairs[-(which(lowpairs[,1] == lowsamp | lowpairs[,2] == lowsamp)),,drop=FALSE]
   }

  if ((nsnpsub < nsnps | depth.min > 1 | depth.max < Inf) & calclevel %in% c(2,9)) {
    naa <- colSums(genon0 == 2, na.rm = TRUE)
    nab <- colSums(genon0 == 1, na.rm = TRUE)
    nbb <- colSums(genon0 == 0, na.rm = TRUE)
    p1 = (naa + nab/2)/(naa + nab + nbb)
    maf <- ifelse(p1 > 0.5, 1 - p1, p1)
    png(paste0("MAF", sfx, ".png"), pointsize = cex.pointsize * 12)
    hist(maf, breaks = 50, xlab = "MAF", col = "grey")
    dev.off()
    }
  samplesOK <- exists("samples")
  if(samplesOK) if(nrow(samples) != nind | any(is.na(samples[depth>0]))) samplesOK <- FALSE
  if (!gform == "chip" & calclevel > 2 & outlevel > 7 & samplesOK) {
   samples0 <- samples[indsubset, snpsubset] - 2 * puse
   samples0[is.na(genon0)] <- 0
   }
  genon0 <- genon0 - 2 * puse
  genon0[is.na(genon0)] <- 0     # equivalent to using 2p for missing genos
  
  sampdepthsub <- rowMeans(depthsub)
  missrate <- sum(as.numeric(depthsub == 0))/nrow(depthsub)/ncol(depthsub)
  cat("Proportion of missing genotypes: ", missrate, "Callrate:", 1-missrate,"\n")
  cat("Mean sample depth:", mean(sampdepthsub), "\n")
  
  Q0 <- puse * (1-puse)
  Q01 <- rowSums(Q0)
  Q0[!usegeno] <- 0
  div0a <- 2 * tcrossprod(Q0, usegeno)
  div0 <- sqrt(div0a) * t(sqrt(div0a))
  
  if (!gform == "chip" & calclevel > 2 & outlevel > 7  & samplesOK) {
    GGBS3top <- tcrossprod(samples0)
    GGBS3bot <- (div0 + diag(diag(div0)))
    GGBS3 <- GGBS3top/GGBS3bot  # faster in 3 steps
    rm(GGBS3top, GGBS3bot)
  } else {
    GGBS3 <- NULL
  }
  
  GGBS4top <- tcrossprod(genon0)
  GGBS4 <- GGBS4top/div0
  GGBS1 <- GGBS4top/crossprod(matrix(sqrt(2*Q01),nrow=1))
  rm(GGBS4top)

  if (samptype == "pooled") {
   raf <- raf - 2 * puse  # centred ref allele freq
   raf[depthsub == 0] <- 0     # equivalent to using 2p for missing genos
   Gpooltop <- tcrossprod(raf)
   Gpool <- Gpooltop/div0
   rm(Gpooltop)
   }
  
  # G5 diag adjustment
  genon01 <- genon0
  P0 <- puse
  P1 <- 1 - P0
  if (have_rcpp) {
    rcpp_assignP0P1Genon01(P0, P1, genon01, usegeno, depthsub, d4i, nThreads)
  }
  else {
    genon01[depth[indsubset, snpsubset] < d4i] <- 0
    P0[!usegeno | depthsub < d4i] <- 0
    P1[!usegeno | depthsub < d4i] <- 0
  }

#  div0 <- 2 * diag(tcrossprod(P0, P1))  # rowSums version faster
  div0 <- 2 * rowSums(P0 * P1)
  depthsub[depthsub < d4i] <- d4i  # not using depthsub values <=1 after this so set to >1 to avoid 0 divisor 
  Kdepth <- depth2K(depthsub)
  GGBS5d <- 1 + rowSums((genon01^2 - 2 * P0 * P1 * (1 + 2*Kdepth))/(1 - 2*Kdepth))/div0

  if (samptype == "pooled") {
   raf1 <- raf
   raf1[depthsub < d4i] <- 0
   diag(Gpool) <- 1 + rowSums((raf1^2 - 2 * P0 * P1 * (1 + 2*Kdepth))/(1 - 2*Kdepth))/div0
   rm(raf1)
   cat("Mean self-relatedness pools (Gpool diagonal):", mean(diag(Gpool)), "\n")
   }

  rm(Kdepth, div0, P0, P1, genon01)
  GGBS5 <- GGBS4
  diag(GGBS5) <- GGBS5d
  cat("Mean self-relatedness (G5 diagonal):", mean(GGBS5d), "\n")
  if(calclevel >1) cat("Mean relatedness (G5 off-diagonal):", (nindsub*mean(GGBS5)-mean(GGBS5d))/(nindsub-1), "\n")
  uhirel <- which(GGBS5 > hirel.thresh & upper.tri(GGBS5), arr.ind = TRUE)
  if (nrow(uhirel) > 0 & nsnpsub >= 999) 
    write.csv(data.frame(Indiv1 = seqID[indsubset][uhirel[, 1]], Indiv2 = seqID[indsubset][uhirel[, 2]], G5rel = GGBS5[uhirel], SelfRel1 = GGBS5d[uhirel[,1]], SelfRel2 = GGBS5d[uhirel[,2]] ), 
              paste0("HighRelatedness", sfx, ".csv"), row.names = FALSE)
  if (!npc == 0 ) {
   # check for missing elements and subset to remove
   pcasamps <- 1:nindsub
   if(length(samp.removed) > 0) {
    pcasamps <- pcasamps[-samp.removed]
    cat("SeqIDs removed for PCA and/or heatmap\n"); print(seqID[indsubset][samp.removed])
    }
   pcacolo <- fcolo[indsubset[pcasamps]]
   npc <- sign(npc) * min(abs(npc),nsnpsub,length(pcasamps))
   if (npc > 0) {
     if (samptype=="pooled") {
      temp <- sqrt(Gpool[pcasamps,pcasamps] - min(Gpool[pcasamps,pcasamps], na.rm = TRUE))
      } else {
      temp <- sqrt(GGBS5[pcasamps,pcasamps] - min(GGBS5[pcasamps,pcasamps], na.rm = TRUE))
      }
     if(withHeatmaply){
       if(length(table(plotly.group)) > 1) {
         temp_p <- heatmaply(x=round(temp,3), symm=TRUE, colors=rev(heat.colors(50)), hide_colorbar=TRUE,
                             width=1000, height=1000, plot_method="plotly", margins=c(0,0,0,0), seriate="none",
                             labRow = paste0("<br>",hover.info), labCol = paste0("<br>",hover.info), showticklabels=FALSE,
                             label_names=c("<b>Row</b>","<b>Column</b>","<b>Relatedness value</b>"),
#                             ColSideColors=pcacolo, RowSideColors=pcacolo,
                             ColSideColors=plotly.group[pcasamps], RowSideColors=plotly.group[pcasamps],
                             file=paste0("Heatmap-G5", sfx, ".html")) %>% layout(width=1000,height=1000)
       } else{
         temp_p <- heatmaply(x=round(temp,3), symm=TRUE, colors=rev(heat.colors(50)), hide_colorbar=TRUE,
                             width=1000, height=1000, plot_method="plotly", margins=c(0,0,0,0), seriate="none",
                             labRow = paste0("<br>",hover.info), labCol = paste0("<br>",hover.info), showticklabels=FALSE,
                             label_names=c("<b>Row</b>","<b>Column</b>","<b>Relatedness value</b>"),
                             file=paste0("Heatmap-G5", sfx, ".html"))
       }
       
       #htmlwidgets::saveWidget(temp_p, paste0("Heatmap-G5", sfx, ".html"))
     } else {
       png(paste0("Heatmap-G5", sfx, ".png"), width = 2000, height = 2000, pointsize = cex.pointsize *  18)
       if (require(parallelDist)) {
         cat("Using parallelDist function in heatmap\n")
         distfun <- parDist
       }
       else {
         cat("Using normal dist function in heatmap\n")
         distfun <- dist
       }
       if(length(table(pcacolo)) > 1) {
         hmout <- heatmap(temp, col = rev(heat.colors(50)), ColSideColors=pcacolo, RowSideColors=pcacolo, symm=TRUE, revC=FALSE, distfun=distfun)
       } else {
         hmout <- heatmap(temp, col = rev(heat.colors(50)), symm=TRUE, revC=FALSE, distfun=distfun)
       }
       hmdat <- data.frame(rowInd=hmout$rowInd,seqIDInd=indsubset[pcasamps[hmout$rowInd]],seqID=seqID[indsubset[pcasamps[hmout$rowInd]]])
       write.csv(hmdat,paste0("HeatmapOrder", sfx, ".csv"),row.names=FALSE,quote=FALSE)
       dev.off()
       }
     }
   }
  if(samptype=="pooled") selfrel <- diag(Gpool) else selfrel <- diag(GGBS5)
  if (calclevel %in% c(2,9)) {
    if(withPlotly){
      temp_p <- plot_ly(y=diag(GGBS4), x=selfrel, hoverinfo="text", text=hover.info, mode="markers", type="scatter",
                        width=480 + addpixel, height=480, marker=list(size=cex.pointsize*6),
                        color=plotly.group, symbol=plotly.group2) %>%
        layout(title="Self-relatedness estimates",xaxis=list(title = ifelse(samptype=="pooled","Using Gpool","Using G5")), yaxis=list(title = "Using G4"))
      htmlwidgets::saveWidget(temp_p, paste0("G", sfx, "-diag.html"))
    }
    png(paste0("G", sfx, "-diag.png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
    if(samptype=="pooled") {
     pairs(cbind(diag(GGBS4), diag(GGBS5), selfrel), col = fcolo[indsubset], main = "Self-relatedness estimates", labels = c("G4","G5","Gpool"), 
           upper.panel=plotpanel, gap=0)
     } else {
     plot(diag(GGBS4) ~ selfrel, col = fcolo[indsubset], main = "Self-relatedness estimates", 
                        xlab = ifelse(samptype=="pooled","Using Gpool","Using G5"), ylab = "Using G4")
     abline(a=0,b=1,col="red",lwd=2)
     }
    dev.off()
   }
  if (max(sampdepthsub) < Inf) {
    ylabel <- "Self-relatedness estimate using G5"
    if(samptype=="pooled") ylabel <- "Self-relatedness estimate using Gpool"
    if(withPlotly){
       temp_p <- plot_ly(y=selfrel, x=sampdepthsub, hoverinfo="text", text=hover.info, mode="markers", type="scatter",
                         width=480 + addpixel, height=480, marker=list(size=cex.pointsize*6),
                         color=plotly.group, symbol=plotly.group2) %>%
        layout(title=ylabel,xaxis=list(title = "Sample depth (log scale)", zeroline=FALSE, type='log'),
               yaxis=list(title = ylabel, zeroline=FALSE))
      htmlwidgets::saveWidget(temp_p, paste0("G", sfx, "diagdepth.html"))
    }
    png(paste0("G", sfx, "diagdepth.png"), width = 480, height = 480, pointsize = cex.pointsize * 12)
      plot(selfrel ~ sampdepthsub, col = fcolo[indsubset], ylab = ylabel, xlab = "Sample depth (log scale)", log="x")
      #plot(diag(GGBS5) ~ I(sampdepthsub + 1), col = fcolo[indsubset], ylab = "Self-relatedness estimate using G5", xlab = "Sample depth +1", log="x")
      dev.off()
    #Slippery slope
    uss <- which(sampdepthsub >= mindepth.idr)
    if(length(uss) > 0) {
     rellm <- lm(selfrel[uss] ~  log(sampdepthsub[uss]))
     slipslope <- coef(rellm)[2]
     pval <- anova(rellm)[1,"Pr(>F)"]
     cat("Self-Rel vs log(depth) regression = ", signif(slipslope,3)," p = ",signif(pval,3)," (min depth = ",mindepth.idr,")\n",sep="")
     }
   }
  if (calclevel %in% c(2,9)) {
   png(paste0("Gcompare", sfx, ".png"), width = 960, height = 960, pointsize = cex.pointsize *  18)
   if(samplesOK & gform != "chip" & calclevel > 2) midpair <- upper.vec(GGBS3) else midpair <- NULL
   if(samptype=="pooled") lastpair <- cbind(upper.vec(GGBS5),upper.vec(Gpool)) else lastpair <- upper.vec(GGBS5)
   Glabels <- "1"
   if(samplesOK & gform != "chip") Glabels <- c(Glabels,"3")
   Glabels <- c(Glabels,"5")
   if(samptype=="pooled") Glabels <- c(Glabels,"pool")
   if (length(Glabels) == 2)  # only GGBS1 and GGBS5 to compare
     plot(upper.vec(GGBS1) ~ upper.vec(GGBS5), col = "#80808060", pch = 16, main = "Off-diagonal comparisons", xlab = "Using G5", ylab = "Using G1")
   if (length(Glabels) > 2 & calclevel > 2) { 
     pairs(cbind(upper.vec(GGBS1), midpair, lastpair), col = "#80808060", pch = 16, main = "Off-diagonal comparisons", 
           labels = paste0("Using G", Glabels), upper.panel=plotpanel, gap=0)
     }
   dev.off()
   }
  npc <- abs(npc)
  if (npc >= 1) {
    ### PCA analysis on GGBS5 or Gpooled
    pcasymbol <- 1; if(length(pcasamps) < 100) pcasymbol <- 16
    if(samptype=="pooled") {
     PC <- svd(Gpool[pcasamps,pcasamps] - matrix(colMeans(Gpool[pcasamps,pcasamps]), nrow = length(pcasamps), ncol = length(pcasamps), byrow = TRUE), nu = npc)
     } else {
     PC <- svd(GGBS5[pcasamps,pcasamps] - matrix(colMeans(GGBS5[pcasamps,pcasamps]), nrow = length(pcasamps), ncol = length(pcasamps), byrow = TRUE), nu = npc)
     } 
    eval <- sign(PC$d) * PC$d^2/sum(PC$d^2)
    PC$x <- PC$u %*% diag(PC$d[1:npc],nrow=npc)  # nrow to get correct behaviour when npc=1
    cat("minimum eigenvalue: ", min(eval), "\n") 
    neprint <- min(2*npc,length(eval))
    cat("first",neprint,"eigenvalues:",eval[1:neprint],"\n")
    #diagnostic plots using first PC
    png(paste0("PC1vInb", sfx, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
      plot(I(selfrel-1)[pcasamps] ~ PC$x[, 1], cex = 0.6, col = pcacolo, pch = pcasymbol, xlab = "Principal component 1", ylab = "Inbreeding")
      dev.off()
    if(max(sampdepthsub) < Inf) {
     png(paste0("PC1vDepth", sfx, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
      plot(sampdepthsub[pcasamps] ~ PC$x[, 1], cex = 0.6, col = pcacolo, pch = pcasymbol, xlab = "Principal component 1", ylab = "Mean Sample Depth")
      dev.off()  
     }
    if(npc > 1) {
     if(withPlotly){
        temp_p <- plot_ly(y=PC$x[, 2],x=PC$x[, 1], type="scatter", mode="markers",
                          hoverinfo="text", text=hover.info, width=640 + addpixel, height=640,
                          marker=list(size=cex.pointsize*6), color=plotly.group, symbol=plotly.group2) %>%
          layout(xaxis=list(title="Principal component 1",zeroline=FALSE), yaxis=list(title="Principal component 2",zeroline=FALSE))
        htmlwidgets::saveWidget(temp_p, paste0("PC1v2G5", sfx, ".html"))
      } 
      png(paste0("PC1v2G5", sfx, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
      if(npc > 1) {
        plot(PC$x[, 2] ~ PC$x[, 1], cex = 0.6, col = pcacolo, pch = pcasymbol, xlab = "Principal component 1", ylab = "Principal component 2")
      } else {
        hist(PC$x[, 1], 50)
      }
      dev.off()
    }
    if (npc > 2) {
      pdf(paste0("PCG5", sfx, ".pdf"), pointsize = cex.pointsize * 12)
      pairs(PC$x[,1:npc], cex=0.6, pch = pcasymbol, label=paste(dimnames(PC$x)[[2]],round(eval,3),sep="\n")[1:npc], col=pcacolo)
      dev.off()
      }
    if(mdsplot) {
      if(samptype=="pooled") {
       mdsout <- cmdscale(dist(Gpool))
       } else {
       mdsout <- cmdscale(dist(GGBS5))
       } 
      if(withPlotly){
        temp_p <- plot_ly(x=mdsout[,1],y=mdsout[,2], type="scatter", mode="markers", marker=list(size=cex.pointsize*6),
                          width=640 + addpixel, height=640 + addpixel,
                          hoverinfo="text", text=hover.info, color=plotly.group, symbol=plotly.group2) %>%
          layout(xaxis=list(title="MDS coordinate 1",zeroline=FALSE), yaxis=list(title="MDS coordinate 2",zeroline=FALSE))
        htmlwidgets::saveWidget(temp_p, paste0("MDS1v2G5", sfx, ".html"))
      }
      else{
        png(paste0("MDS1v2G5", sfx, ".png"), width = 640, height = 640, pointsize = cex.pointsize *  15)
        plot(mdsout, cex = 0.6, col = pcacolo, xlab = "MDS coordinate 1", ylab = "MDS coordinate 2")
        dev.off()
      }
     }
    list(G1 = GGBS1, G4d = diag(GGBS4), G5 = GGBS5, Gpool = Gpool, indsubset = indsubset, samp.removed = samp.removed, PC = PC)  # add G3=GGBS3, if needed
  } else {
    list(G1 = GGBS1, G4d = diag(GGBS4), G5 = GGBS5, Gpool = Gpool, indsubset = indsubset, samp.removed = samp.removed)  # add G3=GGBS3, if needed
  }
}

calcGdiag <- function(snpsubset, puse, indsubset, depth.min = 0, depth.max = Inf, quiet=FALSE ) {
  if (missing(snpsubset))   snpsubset <- 1:nsnps
  if (missing(indsubset))   indsubset <- 1:nind
  if (missing(puse))        puse <- p
  d4i <- 1.001  # min depth to use in inbreeding calcs. Anything in (1,2] OK for normal data, but a low value used in case effective depth (non-integer) being used
  nsnpsub <- length(snpsubset)
  nindsub <- length(indsubset)
  depthsub <- depth[indsubset, snpsubset]
  if(is.matrix(puse)) {  # allows matrix puse to be specified relative to all or subsetted data
   if(!quiet) cat("Using individual allele frequencies\n")
   if(ncol(puse) == nsnps & nsnps > nsnpsub) puse <- puse[,snpsubset]
   if(nrow(puse) == nind & nind > nindsub) puse <- puse[indsubset,]
   } else {
   if(!quiet) cat("Using global allele frequencies\n")
   if(length(puse) == nsnps & nsnps > nsnpsub) puse <- puse[snpsubset]
   puse <- matrix(puse,byrow=TRUE,nrow=nindsub,ncol=nsnpsub)
   }
  if(!nrow(puse) == nindsub | !ncol(puse) == nsnpsub) stop("Dimensions of puse are incorrect.")
  if(min(depth) < d4i) depth[depth < d4i] <- d4i     # these obs become irrelevant
  if(!quiet) cat("Calculating G diags\n")
  if(!quiet) cat("# SNPs: ", nsnpsub, "\n")
  if(!quiet) cat("# individuals: ", nindsub, "\n")
  genon0 <- genon[indsubset, snpsubset]
  usegeno <- !is.na(genon[indsubset, snpsubset])
  if (depth.min > 1 | depth.max < Inf) {
    depthsub[depthsub < depth.min] <- 0
    depthsub[depthsub > depth.max] <- 0
    genon0[depthsub==0] <- NA
    usegeno[depthsub==0] <- FALSE
    }
  genon0 <- genon0 -  2 * puse
  genon0[is.na(genon0)] <- 0     # equivalent to using 2p for missing genos
  genon01 <- genon0
  genon01[depthsub < d4i] <- 0
#  P0 <- matrix(puse[snpsubset], nrow = nindsub, ncol = nsnpsub, byrow = T)
  P0 <- puse
  P1 <- 1 - P0
  P0[!usegeno | depthsub < d4i] <- 0
  P1[!usegeno | depthsub < d4i] <- 0
  div0 <- 2 * rowSums(P0 * P1)
  Kdepth <- depth2K(depth[indsubset, snpsubset])
  GGBS5d <- 1 + rowSums((genon01^2 - 2 * P0 * P1 * (1 + 2*Kdepth))/(1 - 2*Kdepth))/div0
}

ssdInb <- function(dpar=Inf, dmodel="bb", Inbtarget,snpsubset,puse,indsubset,quiet=FALSE, quieti=TRUE) {  
 if (missing(snpsubset))   snpsubset <- 1:nsnps
 if (missing(indsubset))   indsubset <- 1:nind
 if (missing(puse))        puse <- p
 if (length(Inbtarget) != length(indsubset)) stop("Different lengths for Inbtarget and indsubset")
 depth2K <<- depth2Kchoose (dmodel=dmodel, param=dpar)
 NInb <- calcGdiag(snpsubset=snpsubset,puse=puse, indsubset=indsubset, quiet=quieti)-1
 Inbss <- sum((NInb-Inbtarget)^2, na.rm=TRUE)
 if(!quiet) cat(dmodel,"par = ",dpar,"ss = ",Inbss,"\n")
 Inbss
 }

GRMPCA <- function(Guse, PCobj=NULL, npc=2, npcxtra=0, pcasymbol=0, plotname="PC", plotsfx="", pcacolo, plotord=NULL, legendf = NULL, 
                   legendpos = "", move.factor=0.05, cex.plotsize=1, softI=FALSE, ...) {
  # need to add withPlotly part
  npc <- abs(npc)
  ellipsislist <- list(...)
  callinfo <- sys.call()
  if("cex" %in% names(callinfo)) {
   cex.plotsize <- 1 # revert
   cexpos <- which(names(callinfo)=="cex")
   ellipsislist <- c(ellipsislist,list(cex=as.numeric(as.character(callinfo)[cexpos])))
   }
    if(is.null(PCobj)) {
     nsamps <- nrow(Guse)
     if(softI) {  # scale function to centre on column means (NA removed)
      PC <- softImpute(scale(Guse, scale=FALSE), rank.max = npc + npcxtra)
      } else {
      PC <- svd(scale(Guse, scale=FALSE), nu = npc)
      }
     PC$x <- PC$u[,1:npc] %*% diag(PC$d[1:npc],nrow=npc)  # nrow to get correct behaviour when npc=1
     } else {
     PC <- PCobj
     nsamps <- nrow(PC$x)
     }
  if (missing(pcacolo)) pcacolo <- rep("black",nsamps)
  if (missing(plotord)) plotord <- 1:nsamps
  if (!setequal(plotord,1:nsamps)) {cat("Warning: plotord is not an ordering (ignoring)\n"); plotord <- 1:nsamps }
  if(is.list(pcacolo)) { # interpret as a colourby object
   colobj <- pcacolo
   if(!is.numeric(colobj$collabels)) colobj$collabels[is.na(colobj$collabels)] <- "NA"
   pcacolo <- colobj$sampcol
   if (length(unique(colobj$symblist))==1 & all(pcasymbol==0)) pcasymbol <- colobj$symblist[1] 
   if(! is.null(colobj$symblabels) & all(pcasymbol==0)) pcasymbol <- 1   #default to 1 if symbols are for another factor
   if( length(colobj$symblist) >0 & all(pcasymbol==0)) pcasymbol <- colobj$sampsym
   legendsym <- pcasymbol   # used for colour legend
   if (length(unique(colobj$symblist))==1) legendsym <- colobj$symblist[1] else
   if(length(colobj$symblist) == length(colobj$collist)) if(
       length(colobj$collist) == qr(model.matrix(~ colobj$sampcol +  as.factor(colobj$sampsymb)))$rank) legendsym <- colobj$symblist
       # symbols and colours match
   if(! is.null(colobj$symblabels)) pcasymbol <- colobj$sampsym
   legname1 <- NULL; if(!is.null(colobj$col.name)) legname1 <- colobj$col.name
   legname2 <- NULL; if(!is.null(colobj$symb.name)) legname2 <- colobj$symb.name
   ncolumns <- 1; if(length(colobj$collist) > 6 & !is.numeric(colobj$collabels)) ncolumns <- 2
   ncolumns2 <- 1; if(length(colobj$symblist) > 6) ncolumns2 <- 2
   if(legendpos != "") legendf <- function() {  ## To do: use collegend instead 
     if(!is.numeric(colobj$collabels)) leginfo <- legend(legendpos, legend=colobj$collabels, col=colobj$collist, pch=legendsym, ncol=ncolumns, cex=0.9, title=legname1)$rect
     if(is.numeric(colobj$collabels)) { #raster instead of legend
      leginfo <- legend(legendpos,legend=rep("test",3),plot=FALSE,title=legname1)$rect # find where a 3 element legend would go
      topadj <- 1.5 * strheight(legname1)
      leginfo$top <- leginfo$top - topadj; leginfo$h <- leginfo$h - topadj  #find legend placement adj for title
      legend_image <- as.raster(matrix(rev(colobj$collist), ncol = 1))
      rasterImage(legend_image, leginfo$left+0.05*leginfo$w, leginfo$top-0.95*leginfo$h, 
                                 leginfo$left+0.55*leginfo$w, leginfo$top-0.05*leginfo$h)
      rticks <- pretty(colobj$collabels)
      rticks <- rticks[2:(length(rticks)-1)]
      nbars <- length(colobj$collist)
      rasth <- 0.9*leginfo$h
      rast0 <- leginfo$top-0.95*leginfo$h+0.5*rasth/nbars
      rasth <- rasth*(nbars-1)/nbars  # height between end mid-pints
      lines(x=rep(leginfo$left+0.6*leginfo$w,2),y=c(rast0,rast0+rasth))
      rtickss <- (rticks - min(colobj$collabels))/diff(range(colobj$collabels))
      for(itick in seq_along(rticks)){
       lines(x=leginfo$left+c(0.6,0.63)*leginfo$w,y=rep(rast0+rtickss[itick]*rasth,2))
       text(x=leginfo$left+0.65*leginfo$w, y= rast0+rtickss[itick]*rasth, labels=format(rticks)[itick],adj=c(0,0.5), cex=0.7)
       }
      text(x=leginfo$left+0.5*leginfo$w,y=leginfo$top, labels=legname1, xpd=NA, adj=c(0.5,0), cex=0.9)
      } #raster
     ncovered <- length(which(PC$x[,1] > leginfo$left & PC$x[,1] < leginfo$left + leginfo$w & PC$x[,2] < leginfo$top & PC$x[,2] > leginfo$top - leginfo$h))
     if(ncovered > 0) cat("Warning: ",ncovered,"data points covered by colour legend\n")
     if(!is.null(colobj$symblabels)) {
      leginfo2 <- legend(legendpos, legend=colobj$symblabels, ncol=ncolumns2, cex=0.9, title=legname2, plot=FALSE)$rect
      if(leginfo$left < mean(PC$x[,1])) x2 <- leginfo$left + leginfo$w *1.05 else x2 <- leginfo$left - leginfo2$w *1.05 
      if(leginfo$top < mean(PC$x[,2])) y2 <- min(PC$x[,2]) + leginfo2$h else y2 <- min(leginfo$top,max(PC$x[,2]))
      leginfo2 <- legend(x=x2, y=y2, legend=colobj$symblabels, col="black", pch=colobj$symblist, ncol=ncolumns2, cex=0.9, title=legname2)$rect
      ncovered <- length(which(PC$x[,1] > leginfo2$left & PC$x[,1] < leginfo2$left + leginfo2$w & PC$x[,2] < leginfo2$top & PC$x[,2] > leginfo2$top - leginfo2$h))
      if(ncovered > 0) cat("Warning: ",ncovered,"data points covered by symbol legend\n")
      }
     } #legendf
   }
   
   if (npc >= 1) {
    ### PCA analysis on GGBS5
    if(any(pcasymbol == 0)) {
     pcasymbol <- 1
     if(nsamps < 100) pcasymbol <- 16
     }
   if(length(pcasymbol)==1) pcasymbol <- rep(pcasymbol,nsamps)
   xlab0 <- "Principal component 1"
   ylab0 <- "Principal component 2"
   if(class(PC)=="lda") {   xlab0 <- "LDA component 1"; ylab0 <- "LDA component 2"; PC$d <- PC$svd}
   xtratxt <- ""; if(softI & is.null(PCobj)) xtratxt <- paste(" of",npc + npcxtra,"components")
 
   if("d" %in% names(PC)) {
    eval <- sign(PC$d) * PC$d^2/sum(PC$d^2)
    cat("minimum eigenvalue: ", min(eval), "\n") 
    neprint <- min(2*npc,length(eval))
    cat("first",neprint,"eigenvalues:",eval[1:neprint],"\n")
    xlab0 <- paste0(xlab0," (",round(100*eval[1],1),"%",xtratxt,")")
    ylab0 <- paste0(ylab0," (",round(100*eval[2],1),"%",xtratxt,")")
    }
   plotargs <- list(cex = 0.6, xlab = xlab0, ylab = ylab0)
   plotargs[names(ellipsislist)] <- ellipsislist
   png(paste0(plotname,"1v2", plotsfx, ".png"), width = 640*cex.plotsize, height = 640*cex.plotsize, pointsize = cex.pointsize *cex.plotsize *  15)
   if(npc > 1) {
#     plot(PC$x[, 2] ~ PC$x[, 1], cex = 0.6, col = pcacolo, pch = pcasymbol, ... ,
#          xlab = xlab0, ylab = ylab0)
     do.call(plot, c(list(x=PC$x[plotord, 1:2], col = pcacolo[plotord], pch = pcasymbol[plotord]),plotargs) )
     if(!is.null(legendf)) legendf()
   } else {
     hist(PC$x[, 1], 50, ...)
   }
   dev.off()
   if (npc > 2) {
     pdf(paste0(plotname,plotsfx, ".pdf"), pointsize = cex.pointsize * 12)
     pairs(PC$x[,1:npc], cex=0.6, pch = pcasymbol, label=paste(dimnames(PC$x)[[2]],round(eval,3),sep="\n")[1:npc], col=pcacolo, gap=1/10)
     dev.off()
     }
    }
   if(exists("colobj") & npc > 1) {
    if(!is.numeric(colobj$collabels)) {
     if( sum(duplicated(colobj$collist)) > 0) cat("Warning: ambiguous colours - redefine colour object\n")
     sampgroups <- colobj$collabels[match(colobj$sampcol,colobj$collist)]
     groupmean <- aggregate(PC$x[,1:2],list(group=sampgroups),mean)  # group, V1, V2
     groupmean$colour <- colobj$collist[match(groupmean$group,colobj$collabels)]
     #groupmean <- groupmean[match(colobj$collabels,groupmean$group),]
     png(paste0(plotname,"1v2", plotsfx, "-groups.png"), width = 640*cex.plotsize, height = 640*cex.plotsize, pointsize = cex.pointsize *cex.plotsize *  15)
      plot(PC$x[, 2] ~ PC$x[, 1], cex = 0.6, col = "white", ..., xlab=xlab0,ylab=ylab0)
#          xlab = paste0("Principal component 1 (",round(100*eval[1],1),"%",xtratxt,")"),
#          ylab = paste0("Principal component 2 (",round(100*eval[2],1),"%",xtratxt,")"))
      for(igroup in 1:nrow(groupmean)) {
       thislab <- groupmean$group[igroup]
       ugroup <- which(sampgroups == thislab)
       thisn <- length(ugroup)
       linecoords <- rbind(PC$x[ugroup,1:2],matrix(rep(unlist(groupmean[which(groupmean$group==thislab),2:3]),thisn),ncol=2,byrow=TRUE))
       linecoords <- linecoords[order(c(1:thisn,1:thisn)),]
       lines(linecoords[,1],linecoords[,2], col=groupmean$colour[igroup])
       }
      points(groupmean[,2:3],pch=21, col="grey", bg=groupmean$colour)
      #find closeness of labels and separate close ones
      grouppos <- groupmean[,2:3] 
      if(nrow(groupmean) > 1) {
       meandists <- as.matrix(dist(grouppos))
       mindist <- apply(meandists+max(meandists)*diag(nrow(meandists)),2,min)
       mindistpos <- apply(meandists+max(meandists)*diag(nrow(meandists)),2,which.min)
       overalldist <- sqrt((max(PC$x[,1])-min(PC$x[,1]))^2 + (max(PC$x[,2])-min(PC$x[,2]))^2)
       umove <- which(mindist/overalldist < move.factor)
       groupadj <- move.factor*overalldist
       grouppos[umove,] <- grouppos[umove,] + groupadj*((groupmean[,2:3]-groupmean[mindistpos,2:3])/mindist)[umove,]
       }
      for(igroup in 1:nrow(groupmean)) legend(grouppos[igroup,],legend=groupmean[igroup,1], bg="#FFFFFF80", box.col=groupmean$colour[igroup], cex=0.4, xjust=0.5, yjust=0.5, adj=0.2, xpd=TRUE)
      #if(!is.null(legendf)) legendf()
      dev.off()
    }
   }
  invisible(PC)
  }

labelpos <- function(x,y,xrange=range(x),yrange=range(y),gapp=0.02,neighbourp=0.1,xpd=FALSE) {
 #gapp is gap proportion from x,y (to allow a symbol at x,y)
 #neighbourp is the prop of x and y directions to define the neighbourhood
 x0<- xrange[1]; x1 <- xrange[2]
 y0<- xrange[1]; y1 <- xrange[2]
 xstd <- (x-x0)/(x1-x0)
 ystd <- (y-y0)/(y1-y0)
 labelxy <- data.frame(x,y) # start at the points but will overwrite
 for (i in 1:length(x)) {
  uneighbours <- which( sqrt((xstd-xstd[i])^2 + (ystd-ystd[i])^2) < neighbourp)
  if(length(uneighbours)<2) uneighbours <- 1:length(x) # sparse region - plot away from overall mean
  nmeanx <- mean(xstd[uneighbours],na.rm=TRUE)
  nmeany <- mean(ystd[uneighbours],na.rm=TRUE)
  gapx <- gapp/abs(xstd[i]-nmeanx)
  gapy <- gapp/abs(ystd[i]-nmeany)
  labelxy$x[i] <- x0+ (x1-x0) * ((1+gapx)*xstd[i] - gapx * nmeanx)
  labelxy$y[i] <- y0+ (y1-y0) * ((1+gapy)*ystd[i] - gapy * nmeany)
  }
 if(xpd){labelxy$x <- pmin(labelxy$x,x1); labelxy$y <- pmin(labelxy$y,y1); labelxy$x <- pmax(labelxy$x,x0); labelxy$y <- pmax(labelxy$y,y0)}
 labelxy
 }


writeG <- function(Guse, outname, outtype=0, indsubset,IDuse, metadf=NULL ) { # IDuse, metadf is for only those samples in Guse
 Gname <- deparse(substitute(Guse))
 samp.removed <- integer(0)
 if(isTRUE(outtype==0)) cat("Warning, no output will be produced \nSpecify outtype(s):\n 1 R datasets\n 2 G matrix\n 3 long G\n 4 inbreeding\n 5 t-SNE\n 6 PCA\n")
 if (is.list(Guse)) {
  if("PC" %in% names(Guse)) {PCtemp <- Guse$PC; nindout <- nrow(PCtemp$x)}
  if("samp.removed" %in% names(Guse)) {samp.removed <- Guse$samp.removed; nindout <- nindout+length(samp.removed)}
  if(any(outtype != 6)) { # ignore following if only PCs
   if(!"G5" %in% names(Guse)) stop("Guse is a list without a G5")
   Guse <- Guse$G5
   Gname <- "G5"
   nindout <- nrow(Guse)
   } 
  } else {  
   nindout <- nrow(Guse)
   }   
 if (missing(indsubset))   indsubset <- 1:nindout
 if (missing(outname))   outname <- "GBS-Gmatrix"
 IDname <- as.character(deparse(substitute(IDuse)))
 charpos <- regexpr("[",IDname,fixed=TRUE) ; if (charpos>0) IDname <- substr(IDname,1,charpos-1)
 charpos <- regexpr("$",IDname,fixed=TRUE) ; if (charpos>0) IDname <- substr(IDname,charpos+1,nchar(IDname))
 if (missing(IDuse))  { IDuse <- seqID; IDname <- "seqID" }
 IDuse <- IDuse[indsubset]
 if(any(outtype != 6)) Guse <- Guse[indsubset,indsubset]
 if(1 %in% outtype) {
  savelist <- list(Guse=Guse,IDuse=IDuse)
  charpos <- regexpr("[",Gname,fixed=TRUE) ; if (charpos>0) Gname <- substr(Gname,1,charpos-1)
  charpos <- regexpr("$",Gname,fixed=TRUE) ; if (charpos>0) Gname <- substr(Gname,charpos+1,nchar(Gname))
  names(savelist) <- c(Gname,IDname)
  save(list=names(savelist),file= paste0(outname,".RData"), envir=list2env(savelist))
  }
 if(2 %in% outtype) {
  colnames(Guse) <- IDuse
  Gout <- cbind(IDuse, Guse)
  colnames(Gout)[1] <- IDname
  write.csv(Gout, paste0(outname,".csv"),row.names=FALSE,quote=FALSE)
  }
 if(3 %in% outtype) {
  upperRel <- upper.vec(Guse)
  ID1 <- upper.vec(matrix(IDuse,nrow=length(IDuse),ncol=length(IDuse),byrow=FALSE))
  ID2 <- upper.vec(matrix(IDuse,nrow=length(IDuse),ncol=length(IDuse),byrow=TRUE))
  df.out <- data.frame(id1=ID1,id2=ID2,rel=upperRel)
  write.csv(df.out,paste0(outname,"-long.csv"),row.names=FALSE,quote=FALSE)
  }
 if(4 %in% outtype) {
  Inbreeding <- diag(Guse) - 1
  Gout <- cbind(IDuse,Inbreeding)
  colnames(Gout)[1] <- IDname
  write.csv(Gout, paste0(outname,"-Inbreeding.csv"),row.names=FALSE,quote=FALSE)
  }
 if(5 %in% outtype) { # input for projector.tensorflow.org
  colnames(Guse) <- IDuse
  write.table(Guse, paste0(outname,"-pca_vectors.tsv"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  if(is.null(metadf)) write.table(IDuse, paste0(outname,"-pca_metadata.tsv"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  if(!is.null(metadf)) write.table(metadf, paste0(outname,"-pca_metadata.tsv"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  }
 if(6 %in% outtype & !exists("PCtemp")) cat("Warning: PC output specified but no PC object found\n")
 if(6 %in% outtype & exists("PCtemp")) {
  colnames(PCtemp$x) <- paste0("PC",1:ncol(PCtemp$x))
  PCout <- data.frame(IDuse); colnames(PCout) <- IDname
  if(!is.null(metadf)) PCout <- metadf   #cbind(PCout,metadf)
  if(length(samp.removed) > 0) PCout <- PCout[-samp.removed,,drop=FALSE]
  PCout <- cbind(PCout,PCtemp$x)
  write.csv(PCout,paste0(outname,"-PC.csv"),row.names=FALSE,quote=FALSE)
  }
}

G5toDAWG <- function(Guse) {
 #Guse should be a G5 matrix calculated using allele frequencies of 0.5 for all SNPs
 #See Bilton thesis Eqn 6.20
 Z <- 0.5 + Guse/4
 Ztilde <- mean(upper.vec(Z),na.rm=TRUE)
 GWG <- 2* (Z - Ztilde) / (1-Ztilde)
 }

##### functions for G matrix comparison 
regpanel <- function(x,y,nvars=3, ...) {
 usr <- par("usr"); on.exit(par(usr))
 par(usr = c(0, 1, 0, 1))
 if(nvars==2)  par(usr = c(-1.5, 1, 0, 2.5))  # use lower RHS
 if(length(na.omit(x*y)) > 1 & var(x,na.rm=TRUE) > 0 & var(y,na.rm = TRUE)) {
  regn <- summary(lm(x ~ y)) # reverse order to match plot
  if(nrow(regn$coefficients)==1) regn$coefficients <- rbind(regn$coefficients,NA)  ## no slope
  regneqn <- paste0("y = ",signif(regn$coefficients[1,1],3),"\n (SE ",signif(regn$coefficients[1,2],2),")\n + ",
                           signif(regn$coefficients[2,1],3),"x\n (SE ",signif(regn$coefficients[2,2],2),")")
  regnsign <- sign(regn$coefficients[2,1])
  reqnn <- paste0("n = ",sum(regn$df[1:2]))
  reqn <- paste0("r = ",signif(regnsign * sqrt(regn$r.squared),3))
  textcex <- 0.3 + 2/(nvars^.75)   # need to make smaller when more panels
  text(0.5,0.7,labels=regneqn,cex=textcex)
  text(0.5,0.22,labels=reqnn,cex=textcex)
  text(0.5,0.1,labels=reqn,cex=textcex)
  }
 }

plotpanel <- function(x,y, ...) {
  points(x,y,...)
  abline(a=0,b=1,col="red",lwd=2)
}

GCompare <- function(Glist,IDlist,Gnames = paste0("G.",1:length(Glist)), plotname = "", whichplot="both", doBA=FALSE, ...) {
 #whichplot can be "both", "diag" or "off"
 regoutput <- function() {  # function to return regression results in a data frame
  regresults <- data.frame(x=character(0),y=character(0),Intercept=numeric(0),InterceptSE=numeric(0),
                            slope=numeric(0),slopeSE=numeric(0),n=integer(0),correlation=numeric(0))
  for(iG in 1:(nG-1)) for(jG in (iG+1):nG) {
    regn <- summary(lm(gcompare[,jG] ~ gcompare[,iG])) # reverse order to match plot
    if(nrow(regn$coefficients)==1) regn$coefficients <- rbind(regn$coefficients,NA)  ## no slope
    regresults <- rbind(regresults,data.frame(x=Gnames[iG],y=Gnames[jG],
                        Intercept=regn$coefficients[1,1], InterceptSE=regn$coefficients[1,2],
                        slope=regn$coefficients[2,1],slopeSE=regn$coefficients[2,2],n=sum(regn$df[1:2]),
                        correlation= sign(regn$coefficients[2,1])*sqrt(regn$r.squared)) )
   }
   regresults
  }
 if(doBA) doBA <- require(MethComp)
 nG <- length(Glist)
 if (length(IDlist) != nG) stop("ID list different length to G list")
 if (nG < 2) stop("Nothing to compare")
 allID <- unique(unlist(IDlist))
 for(iG in 1:nG) if(sum(duplicated(IDlist[[iG]])) > 0) cat("Warning: Duplicated IDs for ",Gnames[iG],"\n")
 if(whichplot %in% c("both","diag")) {
  gcompare <- diag(Glist[[1]])[match(allID,IDlist[[1]])]
  for(iG in 2:nG)  gcompare <- cbind(gcompare,diag(Glist[[iG]])[match(allID,IDlist[[iG]])])
  cat("Diagonals\n")
  print(regoutput())
  png(paste0("Gcompare-", plotname, "-diag.png"), width = 960, height = 960, pointsize = 18)
   if (nG>2)    pairs(gcompare,labels=Gnames,upper.panel=regpanel,lower.panel=plotpanel,main="Diagonal comparisons", nvars=nG, ...) else {
    plot(gcompare,xlab=Gnames[1],ylab=Gnames[2],main="Diagonal comparisons", ...)
    abline(a=0,b=1,col="red",lwd=2)
    if(!is.null(regpanel)) regpanel(gcompare[,2],gcompare[,1],nvars=2)
    }
   dev.off()
  if(doBA) {
   diffmax <- max(apply(gcompare,1,max,na.rm=TRUE)-apply(gcompare,1,min,na.rm=TRUE))
   gcompareBA <- data.frame(meth=Gnames[1],item=1:length(allID),y=diag(Glist[[1]])[match(allID,IDlist[[1]])])
   for(iG in 2:nG)  gcompareBA <- rbind(gcompareBA,data.frame(meth=Gnames[iG],item=1:length(allID),y=diag(Glist[[iG]])[match(allID,IDlist[[iG]])]))
   gcompareBA$meth <- factor(gcompareBA$meth,levels=Gnames)   # levels specified to keep order
   png(paste0("GcompareBA-", plotname, "-diag.png"), width = 960, height = 960, pointsize = 18)
    plot(Meth(gcompareBA),diff.range=diffmax,main="Diagonal comparisons")
    dev.off()
   }
  }
 if(whichplot %in% c("both","off")) {
  gcompare <- upper.vec(Glist[[1]][match(allID,IDlist[[1]]),match(allID,IDlist[[1]])])
  ncompare <- length(gcompare)
  for(iG in 2:nG)  gcompare <- cbind(gcompare,upper.vec(Glist[[iG]][match(allID,IDlist[[iG]]),match(allID,IDlist[[iG]])]))
  cat("Off-diagonals\n")
  print(regoutput())
  png(paste0("Gcompare-", plotname, "-offdiag.png"), width = 960, height = 960, pointsize = 18)
   if (nG>2) pairs(gcompare,labels=Gnames,upper.panel=regpanel,lower.panel=plotpanel,main="Off-diagonal comparisons", nvars=nG, ...) else {
    plot(gcompare,xlab=Gnames[1],ylab=Gnames[2],main="Off-diagonal comparisons", ...)
    abline(a=0,b=1,col="red",lwd=2)
    if(!is.null(regpanel)) regpanel(gcompare[,2],gcompare[,1],nvars=2)
    }
   dev.off()
  if(doBA) {
   diffmax <- max(apply(gcompare,1,max,na.rm=TRUE)-apply(gcompare,1,min,na.rm=TRUE))
   gcompareBA <- data.frame(meth=Gnames[1],item=1:ncompare,y= upper.vec(Glist[[1]][match(allID,IDlist[[1]]),match(allID,IDlist[[1]])]))
   for(iG in 2:nG)  gcompareBA <- rbind(gcompareBA,data.frame(meth=Gnames[iG],item=1:ncompare,y=upper.vec(Glist[[iG]][match(allID,IDlist[[iG]]),match(allID,IDlist[[iG]])])))
   gcompareBA$meth <- factor(gcompareBA$meth,levels=Gnames)
   png(paste0("GcompareBA-", plotname, "-offdiag.png"), width = 960, height = 960, pointsize = 18)
    plot(Meth(gcompareBA),diff.range=diffmax,main="Off-diagonal comparisons")
    dev.off()
   }
  }
 invisible(NULL)
 }


Gbend <- function(GRM,mineval=0.001, doplot=TRUE, sfx="", evalsum="free") {
  #function to bend GRM to make it positive definite. Use at your own risk!
  #evalsum="const" will rescale to original sum
  GEig <- eigen(GRM) 
  nevalset <- length(which(GEig$values < mineval))
  cat(nevalset,"eigenvalues modified\n")
  evalbent <- pmax(GEig$values,rep(mineval,nrow(GRM)))
  if(evalsum=="const") evalbent <- evalbent*sum(GEig$values)/sum(evalbent)
  GRMbent <- GEig$vectors %*% diag(evalbent) %*% t(GEig$vectors )
  if(doplot) {
  eigcol <- c("black","blue")
   png(paste0("Eigenvalues",sfx,".png"),  width = 960, height = 960,  pointsize = 18)
    plot(GEig$values,pch=16, cex=0.75, col=eigcol[(GEig$values<0)+1.0],xlab="Component", ylab="Eigenvalue",main=paste0("Eigenvalues",sfx))
    abline(h=mineval,col="grey",lwd=2)
    legend("topright",legend=c(">=0","<0"),col=eigcol,pch=16)
    dev.off()
   png(paste0("Self-Bending",sfx,".png"),  width = 960, height = 960,  pointsize = 18)
    plot(diag(GRMbent) ~ diag(GRM), xlab="Using original GRM",ylab="Using bent GRM", main=paste0("Self-relatedness",sfx))
    abline(a=0,b=1,col="red",lwd=3)
    dev.off()
   png(paste0("Rel-Bending",sfx,".png"),  width = 960, height = 960,  pointsize = 18)
    plot(upper.vec (GRMbent) ~ upper.vec (GRM), xlab="Using original GRM",ylab="Using bent GRM", main=paste0("Relatedness",sfx))
    abline(a=0,b=1,col="red",lwd=3)
    dev.off()
  }
  GRMbent
}

## function to use in writeVCF
genostring <- function(vec) {  #vec has gt, paa, pab ,pbb, llaa, llab, llbb, ref, alt
  nobj <- length(vec)/9
  outstr <-  gsub("NA",".",gsub("NA,","",paste(vec[1:nobj],paste(vec[nobj+(1:nobj)],vec[2*nobj+(1:nobj)],vec[3*nobj+(1:nobj)],sep=","), 
                   paste(vec[4*nobj+(1:nobj)],vec[5*nobj+(1:nobj)],vec[6*nobj+(1:nobj)],sep=","), 
                   paste(vec[7*nobj+(1:nobj)],vec[8*nobj+(1:nobj)],sep=","), sep=":")))
  }

## Write KGD back to VCF file
writeVCF <- function(indsubset, snpsubset, outname=NULL, ep=0.001, puse = p, IDuse, keepgt=TRUE, mindepth=0, allele.ref="C", allele.alt="G", 
                     verlabel="4.3",usePL=FALSE, contig.meta=FALSE, CHROM=NULL, POS=NULL){
  if (is.null(outname)) outname <- "GBSdata"
  filename <- paste0(outname,".vcf")
  if(gform=="chip" & max(depth) == Inf & !exists("alleles")) {  # create alleles for chip data
   ref <- alt <- matrix(0,nrow=nind, ncol=nsnps)
   ref[genon >= 1] <- Inf
   alt[genon <= 1] <- Inf
   }
  else if(!exists("alleles"))
    stop("Allele matrix does not exist. Change the 'alleles.keep' argument to TRUE and rerun KGD")
  else if(is.null(alleles))
    stop("Allele matrix object `alleles` is set to NULL.")
  else{
    ref <- alleles[, seq(1, 2 * nsnps - 1, 2)]
    alt <- alleles[, seq(2, 2 * nsnps, 2)]
  }
  ## subset data if required
  if(missing(indsubset))
    indsubset <- 1:nind
  if(missing(snpsubset))
    snpsubset <- 1:nsnps
  if (missing(IDuse)) IDuse <- seqID[indsubset]
  if(length(IDuse) == nind & length(indsubset) < nind) IDuse <- IDuse[indsubset] # assume given for all samples instead of subset
  if(length(IDuse) != length(indsubset) )
    stop("Lengths of IDuse and indsubset arguments does not match")
  is.big <- (as.numeric(length(indsubset)) * length(snpsubset) > 2^31-1 )
  ref <- ref[indsubset, snpsubset]
  alt <- alt[indsubset, snpsubset]
  genon0 <- genon[indsubset, snpsubset]
  pmat <- matrix(puse[snpsubset], nrow=length(indsubset), ncol=length(snpsubset), byrow=TRUE)
  if(length(allele.ref) == nsnps) allele.ref <- allele.ref[snpsubset]
  if(length(allele.alt) == nsnps) allele.alt <- allele.alt[snpsubset]
  
  # Meta information
  metalik <-  '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood">'
  if(usePL)   metalik <-  '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled Genotype Likelihood (rounded)">'
  metaInfo <- paste(paste0('##fileformat=VCFv',verlabel),paste0("##fileDate=",Sys.Date()),"##source=KGDpackage",
                    '##INFO=<ID=INDEL,Number=0,Type=Flag,Description="dummy">',
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probability">', 
                    metalik, 
                    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Read Counts">\n',sep="\n")
  cat(metaInfo, file=filename)
  if(contig.meta) write.table(cbind("##contig=<ID=",SNP_Names[snpsubset],">"),file=filename,sep="",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
  ## colnames:
  cat(c("#CHROM", "POS", "ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT", IDuse), file=filename, append=TRUE, sep="\t")
  cat("\n", file=filename, append=TRUE, sep="\t")
  ## set up the matrix to be written
  out <- matrix(nrow=length(snpsubset),ncol=9+length(indsubset))
  
  temp <- options()$scipen
  options(scipen=10)  #needed for formating
  ## Compute the Data line fields
  if(gform == "Tassel"){
    out[,1] <- chrom[snpsubset]
    out[,2] <- pos[snpsubset]
  }
  else{
    if(is.null(CHROM)) CHROM <- SNP_Names
    if(is.null(POS)) POS <- 1:nind
    out[,1] <- CHROM[snpsubset]
    out[,2] <- POS[snpsubset]
  }
  out[,3] <- SNP_Names[snpsubset]
  out[,4] <- allele.ref
  out[,5] <- allele.alt
  out[,6] <- rep(".", length(snpsubset))
  out[,7] <- rep(".", length(snpsubset))
  out[,8] <- rep(".", length(snpsubset))
  out[,9] <- rep("GT:GP:GL:AD", length(snpsubset))
  if(usePL)  out[,9] <- rep("GT:GP:PL:AD", length(snpsubset))

  
  ## compute probs
  paa <- (1-ep)^ref*ep^alt*pmat^2
  pab <- (1/2)^(ref+alt)*2*pmat*(1-pmat)
  pbb <- ep^ref*(1-ep)^alt*(1-pmat)^2
  psum <- paa + pab + pbb
  paa <- round(paa/psum,4)
  pab <- round(pab/psum,4)
  pbb <- round(pbb/psum,4)
  ## compute likelihood values
  compLike <- function(x) 1/2^(ref+alt)*((2-x)*ep + x*(1-ep))^ref * ((2-x)*(1-ep) + x*ep)^alt
  llaa <- log10(compLike(2))
  llab <- log10(compLike(1))
  llbb <- log10(compLike(0))
  llaa[is.infinite(llaa)] <- -1000
  llab[is.infinite(llab)] <- -1000
  llbb[is.infinite(llbb)] <- -1000
  #phred-scaled ...
  minll <- pmax(llaa,llab,llbb)
  plaa <- -10 * round(llaa - minll,0)
  plab <- -10 * round(llab - minll,0)
  plbb <- -10 * round(llbb - minll,0)
  llaa <- round(llaa,4)
  llab <- round(llab,4)
  llbb <- round(llbb,4)
  if(usePL) { # replace ll with pl versions
   llaa <- plaa
   llab <- plab
   llbb <- plbb
   }
  ## Create the data part
  ## Compute the genotype fields
  depthsub <- ref+alt
  is.na(genon) <- is.na(paa)  <- is.na(pab)  <- is.na(pbb)  <- is.na(llaa)  <- is.na(llab)  <- is.na(llbb) <- (depthsub < mindepth)
  rm(depthsub)
  genon0[is.na(genon0)] <- -1
  if(!keepgt) genon0[] <- -1  # set all elements to missing
  if (is.big) {
   gt <- apply(genon0+2, 2, function(x) c("./.","1/1","0/1","0/0")[x])
   out2[,-c(1:9)] <- apply(cbind(gt,paa,pab,pbb,llaa,llab,llbb,ref,alt),1,genostring)
   } else {
   gt <- sapply(as.vector(genon0), function(x) switch(x+2,"./.","1/1","0/1","0/0"))
   out[,-c(1:9)] <- matrix(gsub("NA",".",gsub("NA,","",paste(gt,paste(paa,pab,pbb,sep=","),paste(llaa,llab,llbb,sep=","), paste(ref,alt,sep=","), sep=":"))), 
                           nrow=length(snpsubset), ncol=length(indsubset), byrow=TRUE)
      # for missings, first change NA, to empty so that any set of NA,NA,...,NA changes to NA, then can set that to . which is vcf missing for the whole field
   }
  options(scipen=temp)

  outord <- 1:nrow(out)
  if(gform == "Tassel" | (!is.null(CHROM) & !is.null(POS))) outord <- order(out[,1],as.numeric(out[,2]))
  ## fwrite is much faster
  if(require(data.table))
    fwrite(split(t(out[outord,]), 1:(length(indsubset)+9)), file=filename, sep="\t", append=TRUE, nThread = 1)
  else
    write.table(out[outord,], file = filename, append = T, sep="\t", col.names = F, row.names=FALSE, quote=FALSE)
  return(invisible(NULL))
}

writeGBS <- function(indsubset,snpsubset,outname="HapMap.hmc.txt",outformat=gform,seqIDuse=seqID) {
 written <- FALSE
  if(!exists("alleles"))
    stop("Allele matrix does not exist. Change the 'alleles.keep' argument to TRUE and rerun KGD")
  else if(is.null(alleles))
    stop("Allele matrix object `alleles` is set to NULL.")
  if(nrow(alleles) != nind | ncol(alleles) != 2* nsnps) stop("Allele matrix does not correspond to genotype matrix")
  if(missing(indsubset)) indsubset <- 1:nind
  if(missing(snpsubset)) snpsubset <- 1:nsnps
  ref <- alleles[indsubset, seq(1, 2 * nsnps - 1, 2)[snpsubset]]
  alt <- alleles[indsubset, seq(2, 2 * nsnps, 2)[snpsubset]]
 if(outformat == "uneak") {
  depthsub <- ref + alt
  genonsub <- ref / depthsub
  genonsub <-  trunc(2*genonsub-1)+1
  HetCount_allele1 <- colSums(ref * (genonsub == 1), na.rm=TRUE)
  HetCount_allele2 <- colSums(alt * (genonsub == 1), na.rm=TRUE)
  allelespos <- seq(2, 2 * nsnps, 2)[snpsubset]
  allelespos <- sort(c(allelespos,allelespos-1))
  acounts <- colSums(alleles[indsubset,allelespos])
  Count_allele1 <- colSums(ref)
  Count_allele2 <- colSums(alt)
  psub <- Count_allele1/(Count_allele1 + Count_allele2)
  outmtx <- matrix(paste(t(ref),t(alt),sep="|"),nrow=length(snpsubset),ncol=length(indsubset))
  colnames(outmtx) <- seqIDuse[indsubset]
#  colnames(outmtx) <- c("rs",seqIDuse[indsubset],"hetc1","hetc2","acount1","acount2","p")
  write.table(cbind(rs=SNP_Names[snpsubset],outmtx,HetCount_allele1,HetCount_allele2,Count_allele1,Count_allele2,Frequency=round(psub,3)),
       outname,row.names=FALSE,quote=FALSE,sep="\t")
  written <- TRUE
  } 
 if(!written) cat("No file written output format not yet available\n")
 invisible(NULL)
 }

### functions for gender testing ####
upperboundary <- function(x){ 20*pmax(rep(0,length(x)),x)^2+0.2} 
lowerboundary <- function(x){ 0.1 + x} 

genderassign <- function(ped.df, index_Y_SNPs, index_X_SNPs, sfx="", hetgamsex = "M", homgamsex = "F", hetchrom = "Y", homchrom = "X", dojitter=FALSE) {
 #ped.df is a dataframe of individuals for gender prediction, as if read from a pedigree file
 #   optionally contains variables Sex (with values M, F or U for male, female, unknown)
 #              and Relationship (character, e.g. "progeny", "sire" or "dam")
 # uses upperboundary and lowerboundary vector functions. upperboundary can be nonlinear. 
 cat("Gender Prediction\n")
 cat(length(index_Y_SNPs), hetchrom,"chromosome SNPs\n")
 cat(length(index_X_SNPs), homchrom,"chromosome SNPs\n")
 genopos <- match(ped.df$seqID,seqID)

 #proporton of SNPs an individual has on the Y chromosome
 proportion_SNPs_Y <- apply(genon[genopos,index_Y_SNPs],1, function(x) sum(!is.na(x))/length(index_Y_SNPs)) #sums number non-missing SNPs #over total SNPs on Y
 #proportion of heterozygosity on the X chromosome
 proportion_Hetero_X <- apply(genon[genopos,index_X_SNPs],1, function(x) sum(x=="1")/sum(!is.na(x))) #######the genon matrix may have to be changed to hetero_X_genon  
 K_matrix <- 1-2*depth2K(depth[genopos,index_X_SNPs])
 K_matrix[K_matrix == -1]<- NaN    #turns values that are -1 (so depth of zero) to NaN
 row_sum_K_matrix<-rowSums(K_matrix,na.rm=TRUE) #sums up the rows 
 total_points_X<-apply(genon[genopos,index_X_SNPs],1,function(x) sum(!is.na(x))) #works out the number of non-missing points on the X chromosome
 E<-row_sum_K_matrix/total_points_X  #divides total K of a row by the number of SNPs on that chromosome
 new_prop_X<-proportion_Hetero_X/E #new depth adjusted proportion of heterozygosity used in subsequent analysis 

 gender_prediction <- rep("Uncertain",nrow(ped.df))
 gender_prediction[which(proportion_SNPs_Y > upperboundary(new_prop_X))] <- hetgamsex
 gender_prediction[which(proportion_SNPs_Y < lowerboundary(new_prop_X))] <- homgamsex
 females <- sum(gender_prediction==homgamsex) #counts number of females #found
 males <- sum(gender_prediction == hetgamsex) #counts numbers of males #found 
 uncertains <-sum(gender_prediction == "Uncertain") #counts number #of uncertains found 
 cat(females, homgamsex,"predicted\n") #prints number of females #predicted 
 cat(males, hetgamsex, "predicted\n") #prints number of males predicted
 cat(uncertains, "not able to be predicted\n") #prints number of #uncertains 

 if(!"Sex" %in% colnames(ped.df)) ped.df$Sex <- "U" 
 print(addmargins(table(Recorded=ped.df$Sex,gender_prediction,useNA="ifany")))
 gender_output <- data.frame(ped.df, gender_prediction,new_prop_X,proportion_SNPs_Y,sampdepth=sampdepth[genopos],stringsAsFactors=FALSE)
 write.csv(gender_output,file=paste0("gender_prediction",sfx,".csv"),row.names=FALSE)

 maxx <- max(0.37,new_prop_X,na.rm=TRUE) 
 maxy <- length(index_Y_SNPs)
 plotch <- rep(19,nrow(ped.df))
 if("Relationship" %in% colnames(ped.df)) {
  rellevels <- unique(ped.df$Relationship)
  plotcharset <- c(20,17,18,15,1:14,65:90)[1:length(rellevels)]
  plotch <- plotcharset[match(ped.df$Relationship,rellevels)]
  }
 gendercol<- c("blue","red","grey","grey")[match(ped.df$Sex ,c("M","F","U", NA))]
 maxxupper <- optimise(function(x) abs(upperboundary(x)-1),lower=0,upper=maxx)$minimum
 yval2plot <- maxy*proportion_SNPs_Y
 if(dojitter) yval2plot <- jitter(yval2plot)
 jtext <- ""; if(dojitter) jtext <- "(jittered)"
 png(file=paste0("GenderPlot",sfx,".png"), height=640, width=640, pointsize = cex.pointsize *  15)
  plot(new_prop_X, yval2plot, xlab=paste0("Heterozygosity (",homchrom," chr)"), ylab=paste(hetchrom,"chr SNPs",jtext), col=gendercol, 
        pch=plotch, cex.axis=1.2, cex.lab=1.2, ylim=c(0,maxy), xlim=c(0,maxx))
  ucheck <-  which(ped.df$Sex != gender_prediction)
  points(new_prop_X[ucheck], yval2plot[ucheck],col=gendercol[ucheck],pch=plotch[ucheck]) #highlight discrepancies
  edges <- par("usr"); #minx maxx miny maxy plotting area
#  poly1 <- data.frame(x1 = c(edges[1],edges[2],edges[2],edges[1],edges[1]), 
#                      y1 = c(lowerboundary(edges[1])*,maxy,lowerboundary(edges[2])*maxy,edges[3],edges[3],lowerboundary(edges[1])*maxy))
  poly1 <- data.frame(x1 = c(seq(edges[1],edges[2]+.01,0.01),edges[2],edges[1],edges[1]), 
                      y1 = c(lowerboundary(seq(edges[1],edges[2]+.01,0.01))*maxy,edges[3],edges[3],lowerboundary(edges[1])*maxy))
  #poly2 <- data.frame(x2 = c(seq(0,0.2,0.01),0,0), y2=c(upperboundary(seq(0,0.2,0.01)),1,0)*maxy) 
  poly2 <- data.frame(x2 = c(seq(edges[1],maxxupper,0.01),maxxupper,edges[1],edges[1]), 
                      y2 = c(maxy*upperboundary(seq(edges[1],maxxupper,0.01)),edges[4],edges[4],maxy*upperboundary(edges[1])))
  polygon(poly1,col = rgb(homgamsex=="F",0,homgamsex=="M",alpha=0.1),border=NA) 
  polygon(poly2,col = rgb(hetgamsex=="F",0,hetgamsex=="M",alpha=0.1),border=NA)
  legend(x=0.8*maxx,y=0.99*edges[4],legend=c("Male","Female","Unknown"),pch=19,col=c("blue","red","grey"),title="Recorded", cex=0.9)
  legend(x=0.55*maxx,y=0.99*edges[4],legend=c("Male","Female"),fill=c(rgb(0,0,1,alpha=0.2),rgb(1,0,0,alpha=0.2)),border=NA,title="Assigned", bg="white",cex=0.9)
  if("Relationship" %in% colnames(ped.df)) legend(x=0.8*maxx,y=0.75*maxy,legend=rellevels,pch=plotcharset,title="Type",cex=0.9)   
  dev.off()
 gender_output 
 }

