#!/home/cew54/localc/bin/Rscript
source("dirs.rb")
library(randomFunctions)
library(magrittr)
library(abind)
## args=getArgs(defaults=list(chr=2, begin=43082452, end=44086664), numeric=c("chr","begin","end"))
args=getArgs(defaults=list(block="chr2_block0"))

################################################################################
## get reference haplotype data to swap alleles - use 1000GP hg38
message("reading haplotypes for block ",args$block)
tmp.imp=tempfile()
paste("./make_ref_block.rb ",args$block,tmp.imp) %>% system()

library(data.table)
leg=fread(paste0(tmp.imp, ".impute.legend"))
haps=scan(paste0(tmp.imp, ".impute.hap"),what=0) %>% matrix(., ncol=nrow(leg))
colnames(haps)=leg$ID

if(ncol(haps)>1000)
  haps=haps[,1:1000]
if(interactive() & ncol(haps)>100)
  haps=haps[,1:100] # for testing

################################################################################

library(simGWAS)
freq=(haps + 1) %>% as.data.frame()
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)

### store MAF, LD
snps <- colnames(haps)
LD=cor(haps)
MAF=colMeans(haps) %>% pmin(., 1-.)

## functions

## simulate log odds ratio
simg=function(n,sd=0.05) {
  g=rnorm(100*n,sd=sd) # sample larger values on purpose
  g=g[order(abs(g),decreasing=TRUE)]
  ## g[sample(1:n)] # sample to make sure A effect not systematically larger
  g[1:n] # don't sample to make sure A effect IS systematically larger
}

simall=function(cv) {
  g <- simg(length(cv)) # log OR
  z <- f.z(cv,g)
  ## abs(z) %>% apply(., 1, max) %>% summary()
  ## (apply(abs(z), 1, max) > 4.89) %>% mean
  vbeta <- f.vb(cv,g)
  beta <- z * sqrt(vbeta)
  keep=apply(abs(z[,match(cv,snps),drop=FALSE]),1,min) > 4.89
  if(!any(keep))
    return(NULL)
  list(beta=beta[keep,],vbeta=vbeta[keep,],g=rep(g,sum(keep)))
}

repsim=function(cv) {
  rep=replicate(100, simall(cv), simplify=FALSE)
  beta=lapply(rep, "[[", "beta") %>% do.call("rbind",.)
  vbeta=lapply(rep, "[[", "vbeta") %>% do.call("rbind",.)
  g=lapply(rep, "[[", "g") %>% unlist()
  list(beta=beta,vbeta=vbeta,g=g)
}

f.z=function(cv,g) {
            simulated_z_score(N0=10000, # number of controls
                            N1=10000, # number of cases
                            snps=snps, # column names in freq of SNPs for which Z scores should be generated
                            W=cv, # causal variants, subset of snps
                            gamma.W=g, # log odds ratios
                            freq=freq, # reference haplotypes
                            nrep=20)
}
f.vb=function(cv,g) {
            simulated_vbeta(N0=10000, # number of controls
                            N1=10000, # number of cases
                            snps=snps, # column names in freq of SNPs for which Z scores should be generated
                            W=cv, # causal variants, subset of snps
                            gamma.W=g, # log odds ratios
                            freq=freq, # reference haplotypes
                            nrep=20)
}

## scenarios
CV=sample(snps,2)

## A only
dA=repsim(CV[1])
z=dA$beta/sqrt(dA$vbeta)
dim(z)
abs(z) %>% apply(., 1, max) %>% summary()

## B only
dB=repsim(CV[2])

## A+B
dAB=repsim(CV)

save(dA,dB,dAB,MAF,LD,CV,
     file=tempfile(pattern=paste0(args$block,"_"),tmpdir=SIMDATA,fileext=".RData"))
