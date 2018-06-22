## I don't know how much R you know already.  Start at the beginning,
## use online R tutorials, the native help system, and books, and work
## through the steps below.  I've tried to write enough that if you
## already know R you won't be bored, so if you don't already know R,
## please don't assume you *have* to complete these tasks!

## install packages - only need to do this first time
install.packages("devtools")
library(devtools)
install_github("chr1swallace/coloc") # has the finemap.abf function
install_github("chr1swallace/simGWAS") # simulated GWAS results

## go read the vignettes on simGWAS and test the code there - see https://github.com/chr1swallace/simGWAS/blob/master/vignettes/intro.Rmd

## go read sections 1-3 of the coloc vignette, which tells you about fine mapping (a little) and test the code
library(coloc)
vignette("vignette",package="coloc")

## each time you start, need to load these packages
library(coloc)
library(simGWAS)

## THIS IS WHERE YOU write your own code!

## 1. first, write your own code to create a credible set.
## e.g. for a random PP
pp <- runif(10) 
pp <- pp/sum(pp) # must sum to 1

## That is - A sort the snps by decreasing posterior probability
## hint: look at the sort function, and the option "decreasing"
?sort
dcrpp <- pp[order(-pp)]

## B find how many SNPs are needed to get reach specific summed PP
## hint: look at the cumsum function
?cumsum

target <- 0.9

snps_needed <- length(which(cumsum(pp2) >= 0.9))
##answer gives 8L which is 8 SNPS needed. Visibly checked as 9th value is 0.937 cumprob

## 2. use some code from the simGWAS vignette to simulate some p values
##Only code added by me was the transformation of z values to p values 
library(simGWAS)


set.seed(42)
nsnps <- 100
nhaps <- 1000
lag <- 5 # genotypes are correlated between neighbouring variants
maf <- runif(nsnps+lag,0.05,0.5) # common SNPs
laghaps <- do.call("cbind", lapply(maf, function(f) rbinom(nhaps,1,f)))
haps <- laghaps[,1:nsnps]
for(j in 1:lag) 
  haps <- haps + laghaps[,(1:nsnps)+j]
haps <- round(haps/(lag+1))

snps <- colnames(haps) <- paste0("s",1:nsnps)
freq <- as.data.frame(haps+1)
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)


CV=sample(snps,1)
g1 <- c(1.4)

FP <- make_GenoProbList(snps=snps,W=CV,freq=freq)
z <- expected_z_score(N0=1000, # number of controls
                      N1=1000, # number of cases
                      snps=snps, # column names in freq of SNPs for which Z scores should be generated
                      W=CV, # causal variants, subset of snps
                      gamma.W=g1, # odds ratios
                      freq=freq, # reference haplotypes
                      GenoProbList=FP) # FP above
zsim <- simulated_z_score(N0=3000, # number of controls
              N1=2000, # number of cases
              snps=snps, # column names in freq of SNPs for which Z scores should be generated
              W=CV, # causal variants, subset of snps
              gamma.W=g1, # log odds ratios
              freq=freq, # reference haplotypes
          nrep=3)
                                   
                                   
p <- 2*pnorm(-abs(z))
                                   typeof(p)
listp <- list(p)
pp <- finemap.abf(listp)
## 3. run the result through finemap.abf() to get posterior probabilities
                                   library(coloc)
setClass("simdata",
         representation(df1="data.frame",df2="data.frame"))
setValidity("simdata", function(object) {
  n <- nrow(object@df1)
  if(nrow(object@df2)!=n)
    return("nrow of '@df1' should equal nrow of '@df2'")
})
setMethod("show", signature="simdata", function(object) {
  cat("pair of simulated datasets, with",ncol(object@df1)-1,"SNPs and",nrow(object@df1),"samples.\n")
})

sim.data <- function(nsnps=50,nsamples=200,causals=1:2,nsim=1) {
  cat("Generate",nsim,"small sets of data\n")
  ntotal <- nsnps * nsamples * nsim
  X1 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y1 <- rnorm(nsamples,rowSums(X1[,causals]),2)
  X2 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y2 <- rnorm(nsamples,rowSums(X2[,causals]),2)
  colnames(X1) <- colnames(X2) <- paste("s",1:nsnps,sep="")
  df1 <- cbind(Y=Y1,X1)
  df2 <- cbind(Y=Y2,X2)
  if(nsim==1) {
    return(new("simdata",
               df1=as.data.frame(df1),
               df2=as.data.frame(df2)))
  } else {
    index <- split(1:(nsamples * nsim), rep(1:nsim, nsamples))
    objects <- lapply(index, function(i) new("simdata", df1=as.data.frame(df1[i,]),
                                             df2=as.data.frame(df2[i,])))
    return(objects)
  }
}

## simulate some data and load the coloc library
set.seed(46411)
data <- sim.data(nsamples=1000,nsim=1)
library(coloc)

## 4. use your code from 1 to find the snps in the credible set for these data.
## Is the true (simulated) causal variant in the set?
## Can you wrap this code in a function, that will take a PP vector, target, and which snp is true causal variant, and return TRUE/FALSE according to whether the causal variant is in the set?

## For bonus points, can you extent the above function to take a *vector* of possible targets and return a vector of TRUE/FALSE?
target <- c(0.5,0.9,0.99) # vector of targets
