#Intro data with 100 simulations

#Intro data with 1 simulation - then need to loop over

## install packages - only need to do this first time
install.packages("devtools")
install_github("chr1swallace/coloc") # has the finemap.abf function
install_github("chr1swallace/simGWAS") # simulated GWAS results

## each time you start, need to load these packages
library(devtools)
library(coloc)
library(simGWAS)

##CREATING P VALUES FOR 1 SIMULATION 
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
#Taking substring so function will recognise format of snps
causal=substr(CV,2,4)

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
                          nrep=100)


p <- 2*pnorm(-abs(zsim))

##GENERATING POSTERIOR PROBABILITIES USING FINEMAP
pvalues <- p 
MAF <- colMeans(haps)

#create matrix split (masplit) to create a list of 100 pvalues per each SNP that finemap.abf will recognise in arguements for pvalues
masplit <- split(p, col(p)) 

#Looping each SNP through finemap.abf to obtain posterior probability 
result = list(length(masplit))
for(i in seq_along(masplit)) {
  result[[i]] <- finemap.abf(dataset=list(pvalues=masplit[[i]], N=1000, MAF=MAF, s=0.5, type="cc"), p1=1e-04)
}

##Remove snp110 which is not relevant in simulation with only one causal variant specified to be in the set 

result <- lapply(result,subset,snp!="null")

##Adjust pp now that pp of snp110 of no causual snp 
lapply(result, function(x) {x$SNP.PP/sum(x$SNP.PP)})
#my.res$SNP.PP <- my.res$SNP.PP/sum(my.res$SNP.PP)
#result$SNP.PP <- lapply(result, result$SNP.PP/sum(result$SNP.PP))
#result$SNP.PP <- lapply(result, function(x) x$SNP.PP/sum($SNP.PP))         

#check to see snp110 was removed and pp were adjusted properly
head(result)
tail(result)

## no loop over and run credset; 
#result <- do.call("rbind", result)
##RUNNING POSTERIOR PROBABILITES THROUGH CREDSET FUNCTION

#select only SNP.PP from list 
pp <- lapply(result, function(x) x$SNP.PP)


credset <- function(pp, causal, thr=0.9,do.order=TRUE) {
  o <- if(do.order) {
    order(pp,decreasing=TRUE)
  } else {
    seq_along(pp)
  }
  
  

  cumpp <- cumsum(pp[o])
  
  ## which is the first snp sums to cumsum > 0.9
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

##should not do this since credset violates sensible parameters.
pp <- as.numeric(unlist(pp))


credset_data = list(pp)
for(i in seq_along(pp)) {
  result[[i]] <- credset(pp = pp, causal = causal, thr = 0.9, do.order = TRUE)
}

#Need to attempt to write to file in directory   write.csv(credset_data, file = '/home/hope/Documents/Thesis documents' )  

#credset_data <- lapply(pp,credset)
#pp <- as.data.frame(lapply(pp, unlist))

value <- credset(pp, causal)
credset(pp, causal,do.order=FALSE)
credset(pp, causal)

value$credset
value$thr
value$size
value$contained






