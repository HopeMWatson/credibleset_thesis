#Set of First taks - intro document 

## install packages - only need to do this first time
install.packages("devtools")
install_github("chr1swallace/coloc") # has the finemap.abf function
install_github("chr1swallace/simGWAS") # simulated GWAS results

## each time you start, need to load these packages
library(devtools)
library(coloc)
library(simGWAS)


## 1. first, write your own code to create a credible set.
## e.g. for a random PP
pp <- runif(10) 
pp <- pp/sum(pp) # must sum to 1

## That is - A sort the snps by decreasing posterior probability
## hint: look at the sort function, and the option "decreasing"

o <- order(pp,decreasing=TRUE)


##HOW MANY SNPS NEEDED TO REACH THRESHOLD

cumpp <- cumsum(pp[o])

thr <- 0.9
wh <- which(cumpp > thr)[1]
wh


## these are the variants in the cred set
o[1:wh]

## check their pp should sum to
sum(pp[ o[1:wh] ])

## what is the causal?
causal %in% o[1:wh]


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
zsim <- simulated_z_score(N0=1000, # number of controls
                          N1=1000, # number of cases
                          snps=snps, # column names in freq of SNPs for which Z scores should be generated
                          W=CV, # causal variants, subset of snps
                          gamma.W=g1, # log odds ratios
                          freq=freq, # reference haplotypes
                          nrep=3)


p <- 2*pnorm(-abs(zsim))

## 3. run the result through finemap.abf() to get posterior probabilities
#Need to assign necessary values to create a dataset that finemap.abf will recognise and run 
pvalues <- p 
p1 <- pvalues[1,]
p1.t <- t(p1)
p2 <- pvalues[2,]
p3 <- pvalues[3,]

MAF <- colMeans(haps)


help("finemap.abf")

##example code syntax from vignette 
##my.res <- finemap.abf(dataset=list(beta=b1$beta, varbeta=b1$varbeta, N=nrow(X1),sdY=sd(Y1),type="quant"))


my.res <- finemap.abf(dataset=list(pvalues=p, N=1000, MAF=MAF, s=0.5, type="cc"), p1=1e-04)
head(my.res)
tail(my.res)



## 4. use your code from 1 to find the snps in the credible set for these data.
## Is the true (simulated) causal variant in the set?
## Can you wrap this code in a function, that will take a PP vector, target, and which snp is true causal variant, and return TRUE/FALSE according to whether the causal variant is in the set?
#CREATING CREDSET FUNCTION
## wrap in a function
credset <- function(pp, causal, thr=0.9,do.order=TRUE) {
  o <- if(do.order) {
    order(pp,decreasing=TRUE)
  } else {
    seq_along(pp)
  }
  
  cumpp <- cumsum(pp[o])
  
  ## which is the first snp s t cumsum > 0.9
  wh <- which(cumpp > thr)[1]
  wh
  
  ## these are the variants in the cred set
  credset <- o[1:wh]
  
  ## check their pp should sum to
  message("credset size: ",sum(pp[ credset ]))
  message("causal variant in credset? ", causal %in% credset)
  
  return(list(credset=credset,
              thr=thr,
              size=sum( pp[credset] ),
              contained=causal %in% credset))
}


value <- credset(pp, causal)
credset(pp, causal,do.order=FALSE)
credset(pp, causal)

value$credset
value$thr
value$size
value$contained

## For bonus points, can you extent the above function to take a *vector* of possible targets and return a vector of TRUE/FALSE?
target <- c(0.5,0.9,0.99) # vector of targets
thrs <- c(0.5,0.9,0.99)
                            