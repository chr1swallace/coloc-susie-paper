#!/home/cew54/localc/bin/Rscript
source("dirs.rb")
library(randomFunctions)
library(magrittr)
afile=list.files(SIMDATA,full=TRUE) %>% sample(.,1)
print(afile)
args=getArgs(defaults=list(file=afile,Ashared=1,Bshared=-1,nsim=10,ncopies=1),
             numeric=c("Ashared","Bshared","nsim","ncopies"))
(load(args$file)) # this contains: dA=datasets with A causal, dB=datasets with B causal, dAB=datasets with AB causal

library(mvtnorm)
if(args$ncopies>1) {
  library(Matrix)
  n0=1000*(args$ncopies-1)
  LD0=LD
  if(args$ncopies==3)
    LD0=bdiag(LD,LD) %>% as.matrix
  Z0=rmvnorm(args$nsim, mean=rep(0,n0), sigma=LD0)
  colnames(Z0)=paste0("null",1:n0)
  LD=bdiag(LD,LD0) %>% as.matrix()
}

devtools::load_all("~/RP/coloc")

## sample coloc data
sample_from_D=function(D) {
     rows=sample(1:nrow(D$beta),1)
     beta=D$beta[rows,]
     vbeta=D$vbeta[rows,]
     names(beta) <- names(vbeta) <- names(MAF)
     list(beta=beta,varbeta=vbeta,z=beta/sqrt(vbeta),snp=names(MAF),MAF=MAF,LD=LD,type="cc",
          N=20000,
          s=0.5)
}

zthr=c(0,0,0.5,1,1.5)
run_susie=function(D,i,zthr=c(0,0,0.5,1,1.5)) {
  times <- ret <- result <- vector("list",length(zthr))
  for(i in seq_along(zthr))
    times[[i]]=system.time(ret[[i]] <- runsusie(D,trimz=zthr[i]))[["elapsed"]]
  a=head(ret[[1]]$pip,-1)
  anull=tail(ret[[1]]$pip,1)
  for(i in 2:length(ret)) {
    b=ret[[i]]$pip[names(a)]
    bnull=tail(ret[[i]]$pip,1)
    result[[i]]= data.table(base_pnull=anull,
             pnull_diff=unname(bnull-anull),
             base_time=times[[1]],
             time_diff=times[[i]]-times[[1]],
             maxz=max(abs(D$z)),
             cv.z=D$z[ CV ] %>% signif(., 3) %>% paste(., collapse="/"),
             base_cv.pip=ret[[1]]$pip[ CV ] %>% signif(., 3) %>% paste(., collapse="/"),
             cv.pip=ret[[i]]$pip[ CV ] %>% signif(., 3) %>% paste(., collapse="/"),
             same_top=which.max(a) == which.max(b),
             max_abs_diff=max(abs(a-b)),
             mse=mean((a-b)^2),
             pipcr=cor(a,b),
             zthr=zthr[i])

  }
  rbindlist(result[2:(length(zthr))])
  ## c(time.trim=time.trim["elapsed"],
  ##   time.notrim=time.notrim["elapsed"],
  ##   nsnps=length(D$snp),
  ##   ncv=if(args$Bshare==1) { 2 } else { 1 },
  ##   pnull.trim=unname(a[length(a)]),
  ##   pnull.notrim=unname(b[length(b)]),
  ##   cor=cor(a[D$snp],b[D$snp]),
  ##   MSE=mean( (a[D$snp]-b[D$snp])^2 ))
}

library(Matrix)
run_sim=function(i) {
  ## D=lapply(1:args$ncopies, function(i) {
  D= if(args$Bshared==1) { message("dAB"); dAB } else { message("dA"); dA }
  D %<>% sample_from_D()
  ## })
  if(args$ncopies>1) {
    D=list(snp=c(D$snp, colnames(Z0)),
           beta=c(D$beta, Z0[i,] * rep(sqrt(D$varbeta),args$ncopies-1)),
           varbeta=rep(D$varbeta, args$ncopies),
           z=c(D$z,Z0[i,]),
           MAF=rep(D$MAF, args$ncopies),
           LD=D$LD, # already set to correct dimension at beginning
           type=D$type,
           N=D$N,
           s=D$s)
    D$snp %<>% make.unique()
    names(D$MAF) <- names(D$z) <- colnames(D$LD) <- rownames(D$LD) <- D$snp
  }
  tmp=run_susie(D,i)
}

RESULT=lapply(1:args$nsim, function(i) run_sim(i))
result=rbindlist(RESULT) #lapply(RESULT,function(x) as.list(x) %>% as.data.table()) %>%
  ## rbindlist()
result[,ncv:=if(args$Bshared==1) { 2 } else { 1 }]
result[,ncopies:=args$ncopies]

save(result,file=tempfile(pattern=paste0("approx_",args$Bshared),
                          tmpdir=file.path(SIMCOLOC,paste0("approx_",args$ncopies)),
                          fileext=".RData"))
