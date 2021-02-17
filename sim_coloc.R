#!/home/cew54/localc/bin/Rscript
source("dirs.rb")
library(randomFunctions)
library(magrittr)
afile=list.files(SIMDATA,full=TRUE) %>% sample(.,1)
print(afile)
args=getArgs(defaults=list(file=afile,Ashared=1,Bshared=-1,nsim=20,ncopies=2),
             numeric=c("Ashared","Bshared","nsim","ncopies"))
(load(args$file)) # this contains: dA=datasets with A causal, dB=datasets with B causal, dAB=datasets with AB causal

devtools::load_all("~/RP/coloc")


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
  dimnames(LD)=list(make.unique(rep(names(MAF),args$ncopies)),
                    make.unique(rep(names(MAF),args$ncopies)))
}

## hold results
SINGLE <- SUSIE_0 <- SUSIE_.5 <- SUSIE_1 <- COND <- MASK <- vector("list",args$nsim)

## sample coloc data
sample_from_D=function(D,i) {
  rows=sample(1:nrow(D$beta),1)
  beta=D$beta[rows,]
  vbeta=D$vbeta[rows,]
  names(beta) <- names(vbeta) <- names(MAF)
  D=list(beta=beta,varbeta=vbeta,z=beta/sqrt(vbeta),snp=names(MAF),MAF=MAF,LD=LD,type="cc",
         N=20000,
         s=0.5)
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
  D
}

run_single=function(D1,D2,i) {
  ret <- coloc:::coloc.abf(D1, D2, p12=5e-6)
  if(is.null(ret$summary))
    return(NULL)
  ret=ret$summary %<>% as.list %>% as.data.frame(., stringsAsFactors=FALSE)
  ret %>% cbind(., data.table(hit1=D1$snp[which.max(abs(D1$z))],
                              hit2=D2$snp[which.max(abs(D2$z))],
                              idx1=1,
                              idx2=1,
                              method="single",
                              run=i))
}

run_susie=function(D1,D2,zthr,i) {
  ret <- coloc:::coloc.susie(D1, D2, p12=5e-6,susie.args=list(trimz=zthr))
  ret %>% cbind(., data.table(method=paste("susie",zthr,sep="_"),
                              run=i))
}

run_condmask=function(D1,D2,i,method=c("mask","cond")) {
  method=match.arg(method)
  ret <- coloc:::coloc.signals(D1, D2, p12=5e-6,method=method)$summary
  ret %>% cbind(., data.table(method=method,
                              run=i))
}


run_sim=function(i) {
  D1={ if(args$Bshared==1) {
         message("dAB"); dAB
       } else {
         message("dA"); dA }
  } %>% sample_from_D(.,i)
  D2={ if(args$Ashared) {
         if(args$Bshared==-1) { message("dA"); dA } else { message("dAB"); dAB }
       } else {
         message("dB"); dB }
  } %>% sample_from_D(.,i)

  ## D1={ if(args$Bshared) { dAB } else { dA } } %>% sample_from_D()
  ## D2={ if(args$Ashared) {
  ##        if(args$Bshared==-1) { dA } else { dAB }
  ##      } else { dB } } %>% sample_from_D()

  RESULT=list(run_single(D1,D2,i),
              run_susie(D1,D2,zthr=0,i),
              run_susie(D1,D2,zthr=0.5,i),
              run_susie(D1,D2,zthr=1,i),
              run_susie(D1,D2,zthr=1.5,i),
              run_condmask(D1,D2,i,"mask"),
              run_condmask(D1,D2,i,"cond"))

  result=rbindlist(RESULT,fill=TRUE)
  result[,hit1.margz:=D1$z[ as.character(hit1) ]]
  result[,hit2.margz:=D2$z[ as.character(hit2) ]]
  result[!is.na(hit1),rA.1:=LD[ as.character(hit1), CV[1] ]]
  result[!is.na(hit1),rB.1:=LD[ as.character(hit1), CV[2] ]]
  result[!is.na(hit1),rA.2:=LD[ as.character(hit2), CV[1] ]]
  result[!is.na(hit1),rB.2:=LD[ as.character(hit2), CV[2] ]]
  result[,rho:=LD[ CV[1], CV[2] ]]
  result
}

RESULT=lapply(1:args$nsim, function(i) run_sim(i))
result=rbindlist(RESULT)

tmpdir=file.path(SIMCOLOC,paste("sim_ncopies",args$ncopies,"A",args$Ashared,"B",args$Bshared,sep="_"))
if(!file.exists(tmpdir))
  dir.create(tmpdir)

save(result,file=tempfile(pattern=paste0("sim_"), #,args$Ashared,"_",args$Bshared),
                          tmpdir=tmpdir,
                          fileext=".RData"))
