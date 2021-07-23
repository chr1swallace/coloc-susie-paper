#!/home/cew54/localc/bin/Rscript
source("dirs.rb")
library(randomFunctions)
library(magrittr)
library(data.table)
ok=FALSE
while(!ok) {
  afile=list.files(SIMDATA,full=TRUE) %>% sample(.,1)
  (load(afile)) # this contains: dA=datasets with A causal, dB=datasets with B causal, dAB=datasets with AB causal
  ok=length(dB$beta) && length(dAB$beta) && length(dA$beta)
}
print(afile)
args=getArgs(defaults=list(file=afile,Ashared=1,Bshared=-1,nsim=10,ncopies=3),
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

## devtools::load_all("~/RP/coloc")
library(coloc)

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

getpip=function(x) {
  colSums(x$alpha[ x$sets$cs_index,,drop=FALSE])
}

zthr=c(0,0.5,1,1.5)
  D= if(args$Bshared==1) { message("dAB"); dAB } else { message("dA"); dA }
  D %<>% sample_from_D()

run_susie=function(D,i,zthr=c(0,0.5,1,1.5)) {
  times <- ret <- result <- vector("list",length(zthr))
  for(j in seq_along(zthr))
    times[[j]]=system.time(ret[[j]] <- runsusie(D,trimz=zthr[j]))[["elapsed"]]
  sapply(ret, function(x) sum(getpip(x)))
  a=getpip(ret[[1]]) #$alpha[ret[[1]]$sets$cs_index,,drop=FALSE]) #/sum(ret[[1]]$pip)
  ## a=head(pip,-1)
  ## anull=tail(pip,1)
  NSETS=lapply(ret, function(x) {
    if(is.null(x$sets$cs)) { 0 } else { length(x$sets$cs) }
  })
  NSETS.in=lapply(ret, function(x) {
    if(is.null(x$sets$cs)) { 0 } else {
                                 A=sum(sapply(x$sets$cs,function(set) any(CV %in% names(set))))
                                 ## A=outer(CV, lapply(x$sets$cs,names), "%in%") # row=cv, col=set
                                 ## apply(A,2,any) %>% sum() # number of sets with >=1 CV
                               }
  })
  NOCV=setdiff(names(a),CV)
  for(j in 1:length(ret)) {
    b=getpip(ret[[j]])[names(a)]
    ## b=pip[names(a)]
    ## bnull=tail(pip,1)
    result[[j]]= data.table(#base_pnull=anull,
                                        #pnull_diff=unname(bnull-anull),
      base_tot=sum(unname(a)),
                            pcv_diff=sum(unname(b[CV]-a[CV])),
                            pnocv_diff=sum(unname(b[NOCV]-a[NOCV])),
                            tot_diff=sum(unname(b - a)),
      base_time=times[[1]],
             time_diff=times[[j]]-times[[1]],
             maxz=max(abs(D$z)),
             cv.z=D$z[ CV ] %>% signif(., 3) %>% paste(., collapse="/"),
             base_cv.pip=ret[[1]]$pip[ CV ] %>% signif(., 3) %>% paste(., collapse="/"),
             cv.pip=ret[[j]]$pip[ CV ] %>% signif(., 3) %>% paste(., collapse="/"),
             same_top=which.max(a) == which.max(b),
             max_pos_diff.noncv=max(b[NOCV]-a[NOCV]),
             min_neg_diff.noncv=min((b[NOCV]-a[NOCV])),
             mse.noncv=mean((a[NOCV]-b[NOCV])^2),
                          pipcr=cor(a,b),
      bias.noncv=mean(b[NOCV]-a[NOCV]))
  }
  result=rbindlist(result[2:(length(zthr))])
  cbind(result,
               data.table(base_nsets=NSETS[[1]],
                          base_nsets.in=NSETS.in[[1]],
                          nsets=unlist(NSETS[-1]),
                          nsets.in=unlist(NSETS.in[-1]),
                          zthr=zthr[-1])
               )

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

  D= if(args$Bshared==1) { message("dAB"); dAB } else { message("dA"); dA }
  D %<>% sample_from_D()
  ## tmp=run_susie(D,i)
  ## tmp

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
