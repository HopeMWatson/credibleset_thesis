#Intro data with 100 simulations

#Intro data with 1 simulation - then need to loop over

## install packages - only need to do this first time
## install.packages("devtools")
## install_github("chr1swallace/coloc") # has the finemap.abf function
## install_github("chr1swallace/simGWAS") # simulated GWAS results

## each time you start, need to load these packages
library(devtools)
library(coloc)
library(simGWAS)

##CREATING P VALUES FOR 1 SIMULATION 
set.seed(42)

## functions
simref <- function(nsnps=100,nhaps=1000,lag=5) {
    maf <- runif(nsnps+lag,0.05,0.5) # common SNPs
    laghaps <- do.call("cbind", lapply(maf, function(f) rbinom(nhaps,1,f)))
    haps <- laghaps[,1:nsnps]
    for(j in 1:lag) 
        haps <- haps + laghaps[,(1:nsnps)+j]
    haps <- round(haps/(lag+1))
    
    snps <- colnames(haps) <- paste0("s",1:nsnps)
    freq <- as.data.frame(haps+1)
    freq$Probability <- 1/nrow(freq)
    freq
}

simdata <- function(freq,OR=1.5,nrep=10,N0=1000,N1=1000) {
    snps <- colnames(freq)[-ncol(freq)] # snp names
    CV=sample(snps,1) # random - different each time

    FP <- make_GenoProbList(snps=snps,W=CV,freq=freq)
    zsim <- simulated_z_score(N0=N0, # number of controls
                              N1=N1, # number of cases
                              snps=snps, # column names in freq of SNPs for which Z scores should be generated
                              W=CV, # causal variants, subset of snps
                              gamma.W=log(OR), # log odds ratios
                              freq=freq, # reference haplotypes
                              nrep=nrep)
    p <- 2*pnorm(-abs(zsim))

    ##GENERATING POSTERIOR PROBABILITIES USING FINEMAP
    MAF <- colMeans(freq[,snps]-1)
    PP <- matrix(0,nrow(p),ncol(p)) # this will hold the posterior probs
    results <- lapply(1:nrep, function(i) {
        tmp <- subset(finemap.abf(dataset=list(pvalues=p[i,], N=1000, MAF=MAF, s=0.5, type="cc"), p1=1e-04),snp!="null")
        tmp$SNP.PP <- tmp$SNP.PP/sum(tmp$SNP.PP)
        colnames(tmp) <- sub("\\.$","",colnames(tmp))
        colnames(tmp) <- sub("SNP.PP","PP",colnames(tmp))
        tmp$CV <- sub("SNP.","",tmp$snp)==sub("s","",CV)
        tmp[,c("snp","pvalues","MAF","PP","CV")]
    })
    results
}

credset <- function(pp, causal, thr=0.9,do.order=TRUE) {
  o <- if(do.order) {
    order(pp,decreasing=TRUE)
  } else {
    sample(seq_along(pp))
  }
  cumpp <- cumsum(pp[o])
  
  ## which is the first snp sums to cumsum > thr
  wh <- which(cumpp > thr)[1]
  wh
  
  ## these are the variants in the cred set
  credset <- o[1:wh]
  
  ## check their pp should sum to
  ## message("credset size: ",sum(pp[ credset ]))
  ## message("causal variant in credset? ", causal %in% credset)
  
  return(list(credset=credset,
              thr=thr,
              size=sum( pp[credset] ),
              contained=causal %in% credset))
}


## use the functions for one scenario
wrapper <- function(thr=0.9,...) {
    data <- simdata(freq,...) ## data is a list of data.frames
    class(data)
    length(data)
    class(data[[1]])
    ## look at the first one
    head(data[[1]])
    summary(data[[1]])
    
    ## work out the credsets
    cs <- lapply(data, function(d) {
        tmp.ord <- credset(d$PP, which(d$CV), thr=thr)
        tmp.noord <- credset(d$PP, which(d$CV), do.order=FALSE,thr=thr)
        data.frame(order=c(TRUE,FALSE),
                   thr=tmp.ord$thr,
                   size=c(tmp.ord$size,tmp.noord$size),
                   nvar=c(length(tmp.ord$credset),
                          length(tmp.noord$credset)),
                   covered=c(tmp.ord$contained,
                             tmp.noord$contained))
    })
    cs <- do.call("rbind",cs)
    cs.ord <- subset(cs,cs$order==TRUE)
    cs.noord <- subset(cs,cs$order==FALSE)
    
    ## what was the coverage?
    c(size.ord=mean(cs.ord$size), cov.ord=mean(cs.ord$covered),
      size.noord=mean(cs.noord$size), cov.noord=mean(cs.noord$covered))
}

#Ordering/No-Ordering
## only need one reference
freq <- simref()
## run the simulation 20 times
ord_noord <- replicate(200,wrapper())
## look at the average result - how well do size.ord and cov.ord match?  how about size.noord and cov.noord?
rowMeans(ord_noord)
se <- apply(ord_noord,1,function(x) sd(x)/sqrt(length(x)))
mn <- rowMeans(ord_noord)
Ord_noord <- cbind(average=mn, lci=mn - 1.96*se, uci=mn + 1.96*se)

#Different Thresholds 
##Target credible set sizes 
####thr parameter
#Holding arguements at OR=1.1 N=1000 s=0.5 
#Varying thr values: 0.5, 0.9, 0.99

#THRESHOLD 0.5
freq <- simref()
## run the simulation 20 times
## different target credible set size
thr_0.5 <- replicate(200,wrapper(OR=1.1,thr=0.5))
rowMeans(thr_0.5)

## can put confidence intervals on like this (see https://en.wikipedia.org/wiki/Standard_error)
se <- apply(thr_0.5,1,function(x) sd(x)/sqrt(length(x)))
mn <- rowMeans(thr_0.5)
Thr_0.5 <- cbind(average=mn, lci=mn - 1.96*se, uci=mn + 1.96*se)

#THRESHOLD 0.9
freq <- simref()
## run the simulation 20 times
## different target credible set size
thr_0.9 <- replicate(200,wrapper(OR=1.1,thr=0.9))
rowMeans(thr_0.9)

## can put confidence intervals on like this (see https://en.wikipedia.org/wiki/Standard_error)
se <- apply(thr_0.9,1,function(x) sd(x)/sqrt(length(x)))
mn <- rowMeans(thr_0.9)
Thr_0.9 <- cbind(average=mn, lci=mn - 1.96*se, uci=mn + 1.96*se)

#THRESHOLD 0.99
freq <- simref()
## run the simulation 20 times
## different target credible set size
thr_0.99 <- replicate(200,wrapper(OR=1.1,thr=0.99))
rowMeans(thr_0.99)

## can put confidence intervals on like this (see https://en.wikipedia.org/wiki/Standard_error)
se <- apply(thr_0.99,1,function(x) sd(x)/sqrt(length(x)))
mn <- rowMeans(thr_0.99)
Thr_0.99 <- cbind(average=mn, lci=mn - 1.96*se, uci=mn + 1.96*se)


##or=1.1
freq <- simref()
## run the simulation 20 times
## different target credible set size
or_1.1 <- replicate(200,wrapper(OR=1.1))
rowMeans(or_1.1)

## can put confidence intervals on like this (see https://en.wikipedia.org/wiki/Standard_error)
se <- apply(or_1.1,1,function(x) sd(x)/sqrt(length(x)))
mn <- rowMeans(or_1.1)
OR_1.1 <- cbind(average=mn, lci=mn - 1.96*se, uci=mn + 1.96*se)

##or=1.2
freq <- simref()
## run the simulation 20 times
## different target credible set size
or_1.2 <- replicate(200,wrapper(OR=1.2))
rowMeans(or_1.2)

## can put confidence intervals on like this (see https://en.wikipedia.org/wiki/Standard_error)
se <- apply(or_1.2,1,function(x) sd(x)/sqrt(length(x)))
mn <- rowMeans(or_1.2)
OR_1.2 <- cbind(average=mn, lci=mn - 1.96*se, uci=mn + 1.96*se)


##or=1.5
freq <- simref()
## run the simulation 20 times
## different target credible set size
or_1.5 <- replicate(200,wrapper(OR=1.5))
rowMeans(or_1.5)

## can put confidence intervals on like this (see https://en.wikipedia.org/wiki/Standard_error)
se <- apply(or_1.5,1,function(x) sd(x)/sqrt(length(x)))
mn <- rowMeans(or_1.5)
OR_1.5 <- cbind(average=mn, lci=mn - 1.96*se, uci=mn + 1.96*se)


rm(dummy)









