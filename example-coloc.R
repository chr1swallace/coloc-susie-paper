setwd("~/Projects/coloc-susie")
source("dirs.rb")
library(randomFunctions)
library(magrittr)
## args=getArgs(defaults=list(chr=2, begin=43082452, end=44086664), numeric=c("chr","begin","end"))
source("~/E/sim_colours.R")
library(ggplot2)
library(ggrastr)
library(data.table)
## (load(file="~/coloc-susie.RData"))

## redefine colours to match plot_dataset
colA="dodgerblue2"
colB="green4"


## ## find a group 3 AB example with high H4
## RSQ_THR=0.5
## NPER_THR=10000 # because every dataset has at least 9660
## results=rbindlist(DATA)
## with(results, table(scenario,tested.cv))
## results[scenario=="sim_1_-1" & tested.cv!="AA",tested.cv:="?"]
## results[scenario=="sim_0_0" & tested.cv %in% c("BA","BB"),tested.cv:="?"]
## results[scenario=="sim_1_0" & tested.cv %in% c("BA"),tested.cv:="?"]
## results[grepl("NA",tested.cv),tested.cv:="?"]
## results[,tested.cv:=factor(tested.cv,levels=c("AA","BB","AB","BA","?"))]
## with(results, table(scenario,tested.cv))
## setnames(results, gsub("PP.|.abf","",names(results)))
## results[,rsq:=cut(rho^2,c(0,RSQ_THR,1),include.lowest=TRUE)]
## hist(results$rho^2); abline(v=RSQ_THR,col="red")
## results[,nper:=length(unique(index)),by=c("scenario","rsq")]
## results[,indexper:=cumsum(!duplicated(index)),by=c("scenario","rsq")]
## results[as.numeric(rsq)==1,.(number_run=max(indexper)),by=c("scenario","rsq")]

## results=results[indexper < NPER_THR]
## results[,number_run:=max(indexper),by=c("scenario","rsq")]
## unique(results[as.numeric(rsq)==1,.(scenario,rsq,number_run,nper)])

## results=results[scenario=="sim_ncopies_2_A_1_B_0"]
## idx=results[tested.cv=="AB" & H4 > 0.9]$index[1]
## results[index==idx & method=="cond"]
## results[index==idx & method=="susie_0"]
## idx

## (load(sub("_1$","",idx)))
## head(result)

source("dirs.rb")
library(randomFunctions)
library(magrittr)

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
  ret <- coloc:::coloc.susie(D1, D2, p12=5e-6,susie.args=list(trimz=zthr),nref=503)
  cbind(ret$summary, data.table(method=paste("susie",zthr,sep="_"),
                              run=i))
}

run_condmask=function(D1,D2,i,method=c("mask","cond"),mode=c("iterative","allbutone")) {
  method=match.arg(method)
  mode=match.arg(mode)
  ret <- coloc:::coloc.signals(D1, D2, p12=5e-6,method=method,mode=mode)$summary
  ret %>% cbind(., data.table(method=paste(method,mode,sep="_"),
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
  RESULT=list(## run_single(D1,D2,i),
              run_susie(D1,D2,zthr=0,i),
              ## run_susie(D1,D2,zthr=0.5,i),
              ## run_susie(D1,D2,zthr=1,i),
              ## run_susie(D1,D2,zthr=1.5,i),
              ## run_condmask(D1,D2,i,"mask"),
              run_condmask(D1,D2,i,"cond",mode="iterative"),
              run_condmask(D1,D2,i,"cond",mode="allbutone"))
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

devtools::load_all("~/RP/coloc")
library(mvtnorm)
  library(Matrix)
  args=getArgs(defaults=list(file="",Ashared=1,Bshared=0,nsim=1,ncopies=1),
               numeric=c("Ashared","Bshared","nsim","ncopies"))

## find examples
if(!interactive()) {
  all_files=list.files(SIMDATA,full=TRUE) %>% sample()
  found=list()
  for(i in seq_along(all_files)) {
    afile=all_files[i] # %>% sample(.,1)
    (load(afile)) # this contains: dA=datasets with A causal, dB=datasets with B causal, dAB=datasets with AB causal
    if(!length(dB$beta) || !length(dAB$beta) || !length(dA$beta))
      next
    for(j in 1:100) {
      ## set.seed(42)
      seed=sample(1:10000000,1) #.Random.seed
      set.seed(seed)
      result=invisible(run_sim(1))
      ## print(result[,.(method,PP.H4.abf,rB.2)])
      bad_cond_it=which(result$method=="cond_iterative" & result$PP.H4.abf>0.8 & result$rB.2 > 0.8)
      bad_cond_abo=which(result$method=="cond_allbutone" & result$PP.H4.abf>0.8 & result$rB.2 > 0.8)
      if(length(bad_cond_it) | length(bad_cond_abo)) {
        message("\t !!!\tfound\t",afile,"\t",seed)
        fwrite(result[,.(afile=afile,seed=seed,method,PP.H3.abf,PP.H4.abf,rA.2,rB.2)],
               file="found4.txt",
               append=TRUE,
               col.names=FALSE)
        ## cat(c(afile,seed,nrow(result),length(bad_cond_it),length(bad_cond_abo),"\n"),file="found3.txt",sep="\t",append=TRUE)
      }
    }
  }
  stop()
}

## load examples
found=fread("found4.txt",header=FALSE)
head(found)
setnames(found, c("afile","seed","method","H3","H4","rA.2","rB.2"))
summ=found[,.(max4=max(H4),min4=min(H4),nrow=.N),by=c("afile","seed","method")]
head(summ)
summ=summ[nrow==2]
head(isumm <- summ[min4>0.8 & method=="cond_iterative"])
head(asumm <- summ[min4>0.8 & method=="cond_allbutone"])


################################################################################

### ITERATIVE EXAMPLE

################################################################################

k=3
afile=isumm$afile[k]
seed=isumm$seed[k]
  (load(afile)) # this contains: dA=datasets with A causal, dB=datasets with B causal, dAB=datasets with AB causal
  set.seed(seed)
result=run_sim(1)
print(result[,.(method,PP.H3.abf,PP.H4.abf,rB.2)])

  set.seed(seed)
# dA
  D1=dA %>% sample_from_D(.,i)
# dAB
  D2=dAB %>% sample_from_D(.,i)

## run_condmask(D1,D2,1,"cond",mode="allbutone")
pos=sub("chr[0-9]+-","",names(D1$beta)) %>% as.numeric()
S1=runsusie(D1,null_weight=0)
S2=runsusie(D2,null_weight=0)
D1$position=pos
D2$position=pos
poscv=sub("chr[0-9]+-","",CV) %>% as.numeric()
par(mfrow=c(2,1))
plot_dataset(D1,S1); abline(v=poscv[1],col=colA)
plot_dataset(D2,S2); abline(v=poscv,col=c(colA,colB))

LD[CV,CV]^2

## step through coloc_signals, extract BF and plot
dataset1=D1
dataset2=D2
method="cond"
  fm1 <- finemap.signals(dataset1,
                           method=method,
                           maxhits=2,r2thr=0.01,pthr=1e-6)
    fm2 <- finemap.signals(dataset2,
                           method=method,
                           maxhits=2,r2thr=0.01,pthr=1e-4)

LD[names(fm2),CV]
# first is B, second is A

highA=rownames(LD)[ LD[,CV[1]] > 0.8 ] %>% sub("chr.*-","",.) %>% as.numeric()
highB=rownames(LD)[ LD[,CV[2]] > 0.8 ] %>% sub("chr.*-","",.) %>% as.numeric()
wA=which( LD[,CV[1]] > 0.8 )
wB=which( LD[,CV[2]] > 0.8 )

## method cond, mode iterative
cond1 <- est_all_cond(dataset1,fm1,mode="iterative")
X1 <- dataset1[ intersect(names(dataset1), c("N","sdY","type","s","position")) ]
cond2 <- est_all_cond(dataset2,fm2,mode="iterative")
cond3 <- est_all_cond(dataset2,fm2,mode="allbutone")
X2 <- dataset2[ intersect(names(dataset2), c("N","sdY","type","s","position")) ]
## cond1 <- est_all_cond(dataset1,fm1,mode="allbutone")
## X1 <- dataset1[ intersect(names(dataset1), c("N","sdY","type","s","position")) ]
## cond2 <- est_all_cond(dataset2,fm2,mode="allbutone")
cond2[[2]][,z:=beta/sqrt(varbeta)]
cond2[[1]][,z:=beta/sqrt(varbeta)]
cond3[[2]][,z:=beta/sqrt(varbeta)]
cond3[[1]][,z:=beta/sqrt(varbeta)]
X2 <- dataset2[ intersect(names(dataset2), c("N","sdY","type","s","position")) ]

annot=function(what=c("AB","A")) {
  what=match.arg(what)
  if(what=="A") {
    abline(v=poscv[1],col=c(colA)); text(x=poscv[1],y=5,col=c(colA),labels=c("A"),adj=0)
  } else {
    abline(v=poscv,col=c(colA,colB)); text(x=poscv,y=5,col=c(colA,colB),labels=c("A","B"),adj=0)
  }
}
colpoints=function(y,w,col) {
  points(D1$position[w],y[w], col=col, pch=16)
}



png("figure-example-iterative.png",width=8,height=6,units="in",pointsize=10,res=300)

iplot=1
par(mfcol=c(3,2),mar=c(4,3,2,1),mgp=c(2,1,0),bty="L",oma=c(0,0.2,0,0))
plot_dataset(D1); title(main="trait 1, observed data");
colpoints(y=-log10(pnorm(-abs(D1$z))*2), w=wA, col=colA)
annot("A")
label_plot()

plot_dataset(D2); title(main="trait 2, observed data")
colpoints(y=-log10(pnorm(-abs(D2$z))*2), w=wA, col=colA)
colpoints(y=-log10(pnorm(-abs(D2$z))*2), w=wB, col=colB)
annot()
legend("topright",legend=paste0("coloc with trait 1: PP H4=",signif(result[method=="cond_iterative"]$PP.H4.abf[1],2)))
label_plot()

X=X2; C=cond2[[2]]
plot_dataset(c(C,X)); title(main="trait 2, conditioned on proxy for B")
colpoints(y=-log10(pnorm(-abs(C$z))*2), w=wA, col=colA)
colpoints(y=-log10(pnorm(-abs(C$z))*2), w=wB, col=colB)
annot()
legend("topright",legend=paste0("coloc with trait 1: PP H4=",signif(result[method=="cond_iterative"]$PP.H4.abf[2],2)))
label_plot()

## dev.off()

## png("example-iterative+susie.png",width=8,height=4,units="in",pointsize=12,res=300)

## par(mfrow=c(3,1))
y=S1$lbf_variable[S1$sets$cs_index[1],]
plot_dataset(D1,alty=y,ylab="log BF",show_legend=FALSE,main="trait 1, SuSiE signal 1"); annot()
colpoints(y,wA,colA)
label_plot()
for(k in 1:2) {
  y=S2$lbf_variable[S2$sets$cs_index[k],]
  plot_dataset(D2,alty=y,ylab="log BF",show_legend=FALSE,main=paste("trait 2, SuSiE signal",k)); annot()
  colpoints(y,wA,colA)
  colpoints(y,wB,colB)
  legend("topright",legend=paste0("coloc with trait 1: PP H4=",signif(result[method=="susie_0"]$PP.H4.abf[k],2)))
label_plot()
}

dev.off()

## ## method cond, mode iterative
## cond1 <- est_all_cond(dataset1,fm1,mode="iterative")
## X1 <- dataset1[ intersect(names(dataset1), c("N","sdY","type","s","position")) ]
## cond2 <- est_all_cond(dataset2,fm2,mode="iterative")
## X2 <- dataset2[ intersect(names(dataset2), c("N","sdY","type","s","position")) ]
## par(mfrow=c(3,1))
## plot_dataset(D1); abline(v=poscv[1],col=colA); title(main="trait 1"); text(x=poscv[1],y=6,col=colA,labels="A",adj=0)
## plot_dataset(D2); abline(v=poscv,col=c(colA,colB)); title(main="trait 2")
## ## abline(v=sub("chr.*-","",names(fm1)[1]) %>% as.numeric(),col=c(colA),lty=2); #text(x=poscv,y=6,col=c(colA,colB),labels=c("A","B"),adj=0)
## ## abline(v=sub("chr.*-","",names(fm2)[1]) %>% as.numeric(),col=c(colB),lty=2); #text(x=poscv,y=6,col=c(colA,colB),labels=c("A","B"),adj=0)
## ## plot_dataset(D2); abline(v=poscv,col=c(colA,colB)); title(main="trait 2")
## text(x=poscv,y=6,col=c(colA,colB),labels=c("A","B"),adj=0)
## plot_dataset(c(cond2[[2]],X2)); abline(v=poscv,col=c(colA,colB)); title(main="trait 2, conditioned on B")
## text(x=poscv,y=6,col=c(colA,colB),labels=c("A","B"),adj=0)

## ## step through coloc_susie, extract BF and plot
## plot_dataset(D1,S1,what="pip")
## plot_dataset(D1,S1,what="alpha")
## plot_dataset(D2,S2,what="alpha")
## plot_dataset(D2,S2,what="alpha",alpha.row=2)

## par(mfrow=c(3,1))
## plot_dataset(D1,S1,what="alpha"); abline(v=poscv[1],col=colA); title(main="trait 1, only signal"); text(x=poscv[1],y=6,col=colA,labels="A",adj=0)
## plot_dataset(D2,S2,what="alpha"); abline(v=poscv,col=c(colA,colB)); title(main="trait 2, first signal")
## text(x=poscv,y=6,col=c(colA,colB),labels=c("A","B"),adj=0)
## plot_dataset(D2,S2,what="alpha",alpha.row=2); abline(v=poscv,col=c(colA,colB)); title(main="trait 2, second signal")
## text(x=poscv,y=6,col=c(colA,colB),labels=c("A","B"),adj=0)

################################################################################

## ABO example

################################################################################

label_plot=function() {
  mtext(letters[iplot],side=3,outer=FALSE,adj=-0.1,font=2)
  iplot <<- iplot + 1
}

k=1
afile=asumm$afile[k]
seed=asumm$seed[k]
  (load(afile)) # this contains: dA=datasets with A causal, dB=datasets with B causal, dAB=datasets with AB causal
  set.seed(seed)
result=run_sim(1)
print(result[,.(method,PP.H3.abf,PP.H4.abf,rB.2)])


  set.seed(seed)
# dA
  D1=dA %>% sample_from_D(.,i)
# dAB
  D2=dAB %>% sample_from_D(.,i)

## run_condmask(D1,D2,1,"cond",mode="allbutone")
pos=sub("chr[0-9]+-","",names(D1$beta)) %>% as.numeric()
S1=runsusie(D1,null_weight=0)
S2=runsusie(D2,null_weight=0)
D1$position=pos
D2$position=pos
poscv=sub("chr[0-9]+-","",CV) %>% as.numeric()
par(mfrow=c(2,1))
plot_dataset(D1,S1); abline(v=poscv[1],col=colA)
plot_dataset(D2,S2); abline(v=poscv,col=c(colA,colB))

LD[CV,CV]^2

## step through coloc_signals, extract BF and plot
dataset1=D1
dataset2=D2
method="cond"
  fm1 <- finemap.signals(dataset1,
                           method=method,
                           maxhits=2,r2thr=0.01,pthr=1e-6)
    fm2 <- finemap.signals(dataset2,
                           method=method,
                           maxhits=2,r2thr=0.01,pthr=1e-4)

LD[names(fm2),CV] # first is B, second is A

highA=rownames(LD)[ LD[,CV[1]] > 0.8 ] %>% sub("chr.*-","",.) %>% as.numeric()
highB=rownames(LD)[ LD[,CV[2]] > 0.8 ] %>% sub("chr.*-","",.) %>% as.numeric()
wA=which( LD[,CV[1]] > 0.8 )
wB=which( LD[,CV[2]] > 0.8 )

## method cond, mode ABO
cond1 <- est_all_cond(dataset1,fm1,mode="iterative")
X1 <- dataset1[ intersect(names(dataset1), c("N","sdY","type","s","position")) ]
cond2 <- est_all_cond(dataset2,fm2,mode="iterative")
cond3 <- est_all_cond(dataset2,fm2,mode="allbutone")
X2 <- dataset2[ intersect(names(dataset2), c("N","sdY","type","s","position")) ]
## cond1 <- est_all_cond(dataset1,fm1,mode="allbutone")
## X1 <- dataset1[ intersect(names(dataset1), c("N","sdY","type","s","position")) ]
## cond2 <- est_all_cond(dataset2,fm2,mode="allbutone")
cond2[[2]][,z:=beta/sqrt(varbeta)]
cond2[[1]][,z:=beta/sqrt(varbeta)]
cond3[[2]][,z:=beta/sqrt(varbeta)]
cond3[[1]][,z:=beta/sqrt(varbeta)]
X2 <- dataset2[ intersect(names(dataset2), c("N","sdY","type","s","position")) ]

annot=function(what=c("AB","A")) {
  what=match.arg(what)
  if(what=="A") {
    abline(v=poscv[1],col=c(colA)); text(x=poscv[1],y=5,col=c(colA),labels=c("A"),adj=0)
  } else {
    abline(v=poscv,col=c(colA,colB)); text(x=poscv,y=5,col=c(colA,colB),labels=c("A","B"),adj=0)
  }
}

png("figure-example-abo.png",width=8,height=8,units="in",pointsize=10,res=300)

iplot=1
par(mfcol=c(4,2),mar=c(4,3,2,1),mgp=c(2,1,0),bty="L",oma=c(0,0.2,0,0))
plot_dataset(D1); title(main="trait 1, observed data");
points(D1$position[wA],-log10(pnorm(-abs(D1$z[wA]))*2), col=colA, pch=16)
annot("A")
label_plot()

plot_dataset(D2); title(main="trait 2, observed data")
points(D2$position[wA],-log10(pnorm(-abs(D2$z[wA]))*2), col=colA, pch=16)
points(D2$position[wB],-log10(pnorm(-abs(D2$z[wB]))*2), col=colB, pch=16)
annot()
label_plot()
## abline(v=names(fm2) %>% sub("chr.*-","",.) %>% as.numeric(),
##        col=c(colB,colA),
##        lty=2)

X=X2; C=cond3[[1]]
plot_dataset(c(C,X));  title(main="trait 2, conditioned on proxy for A")
points(D2$position[wA],-log10(pnorm(-abs(C$z[wA]))*2), col=colA, pch=16)
points(D2$position[wB],-log10(pnorm(-abs(C$z[wB]))*2), col=colB, pch=16)
annot()
legend("topright",legend=paste0("coloc with trait 1: PP H4=",signif(result[method=="cond_allbutone"]$PP.H4.abf[1],2)))
label_plot()

X=X2; C=cond3[[2]]
plot_dataset(c(C,X)); title(main="trait 2, conditioned on proxy for B")
points(X1$position[wA],-log10(pnorm(-abs(C$z[wA]))*2), col=colA, pch=16)
points(X1$position[wB],-log10(pnorm(-abs(C$z[wB]))*2), col=colB, pch=16)
## abline(v=sub("chr.*-","",names(fm2)[1]) %>% as.numeric(),col=c(colB),lty=2); #text(x=poscv,y=6,col=c(colA,colB),labels=c("A","B"),adj=0)
annot()
legend("topright",legend=paste0("coloc with trait 1: PP H4=",signif(result[method=="cond_allbutone"]$PP.H4.abf[2],2)))
label_plot()

## par(mfrow=c(3,1))
y=S1$lbf_variable[S1$sets$cs_index[1],]
plot_dataset(D1,alty=y,ylab="log BF",show_legend=FALSE,main="trait 1, SuSiE signal 1"); annot()
colpoints(y,wA,colA)
label_plot()

par(mfg=c(3,2))
for(k in 1:2) {
  y=S2$lbf_variable[S2$sets$cs_index[k],]
  plot_dataset(D2,alty=y,ylab="log BF",show_legend=FALSE,main=paste("trait 2, SuSiE signal",k)); annot()
  colpoints(y,wA,colA)
  colpoints(y,wB,colB)
  legend("topright",legend=paste0("coloc with trait 1: PP H4=",signif(result[method=="susie_0"]$PP.H4.abf[k],2)))
  label_plot()
}

dev.off()
