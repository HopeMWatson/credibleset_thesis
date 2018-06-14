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
dcrpp <- pp2[order(-pp2)]

## B find how many SNPs are needed to get reach specific summed PP
## hint: look at the cumsum function
?cumsum

target <- 0.9

snps_needed <- length(which(cumsum(pp2) <= 0.9))
##answer gives 8L which is 8 SNPS needed. Visibly checked as 9th value is 0.937 cumprob

## 2. use some code from the simGWAS vignette to simulate some p values

## 3. run the result through finemap.abf() to get posterior probabilities

## 4. use your code from 1 to find the snps in the credible set for these data.
## Is the true (simulated) causal variant in the set?
## Can you wrap this code in a function, that will take a PP vector, target, and which snp is true causal variant, and return TRUE/FALSE according to whether the causal variant is in the set?

## For bonus points, can you extent the above function to take a *vector* of possible targets and return a vector of TRUE/FALSE?
target <- c(0.5,0.9,0.99) # vector of targets
