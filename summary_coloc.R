#!/home/cew54/localc/bin/Rscript
setwd("~/Projects/coloc-susie")
source("dirs.rb")
library(randomFunctions)
library(magrittr)
## args=getArgs(defaults=list(chr=2, begin=43082452, end=44086664), numeric=c("chr","begin","end"))
source("~/E/sim_colours.R")
NCOPIES_TO_SHOW=1000

## example with high LD between CVs (0.0.99)
## f="data/sim/coloc/_223c973ddcaff_255f719b075c4.RData"
library(pbapply)
dirs=list.files(SIMCOLOC,full=TRUE,pattern="sim_")
d=dirs[1]
for(d in dirs) {
  files=list.files(d,full=TRUE) %>% grep("summary",.,value=TRUE,invert=TRUE)
  message(basename(d),"\tfiles found: ",length(files))
}

library(ggplot2)
library(ggrastr)
library(data.table)

proc_file=function(f,d,rsq=0.5) {
  (load(f))
  result[,scenario:=basename(d)][,index:=paste(f,run,sep="_")]
  result[,tested.1:="NA"]
  result[,tested.2:="NA"]
  result[rA.1^2 > rsq & rA.1^2 > rB.1^2, tested.1:="A"]
  result[rB.1^2 > rsq & rB.1^2 > rA.1^2, tested.1:="B"]
  result[rA.2^2 > rsq & rA.2^2 > rB.2^2, tested.2:="A"]
  result[rB.2^2 > rsq & rB.2^2 > rA.2^2, tested.2:="B"]
  result[,tested.cv:=paste0(tested.1,tested.2)]
  ## result[,tested.cv:=paste0(ifelse(rA.1^2 > rB.1^2, "A", "B"), ifelse(rA.2^2 > rB.2^2, "A", "B")) ]
  ## ss=split(result, result[,.(index,method)])
  ## result[,.(PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,method,hit1.margz,hit2.margz,scenario,index,conclusion)]
  result
}

d=dirs[1]
DATA=vector("list",length(dirs))
names(DATA)=basename(dirs)
for(d in dirs) {
  files=list.files(d,full=TRUE) %>% grep("summary",.,value=TRUE,invert=TRUE)
  f=files[1]
  message(basename(d),"\tfiles found: ",length(files))
  if(!length(files))
    next
  if(length(files)>1000)
    files=sample(files,1000)
  DATA[[ basename(d) ]] = lapply(files,proc_file,d=d) %>% rbindlist()
  ## scores[is.nan(spec.B),spec.B:=0]
  ## scores[is.nan(spec.A),spec.A:=0]
}
save(DATA,file="~/coloc-susie-11-42.RData")

################################################################################

## preparation

(load(file="~/coloc-susie-11-42.RData"))

RSQ_THR=0.5
NPER_THR=10000 # because every dataset has at least 9660
results=rbindlist(DATA)
results=results[!grepl("mask",method)]
results[,method:=factor(as.character(method))]
with(results, table(scenario,tested.cv))
results[scenario=="sim_1_-1" & tested.cv!="AA",tested.cv:="?"]
results[scenario=="sim_0_0" & tested.cv %in% c("BA","BB"),tested.cv:="?"]
results[scenario=="sim_1_0" & tested.cv %in% c("BA"),tested.cv:="?"]
results[grepl("NA",tested.cv),tested.cv:="?"]
results[,tested.cv:=factor(tested.cv,levels=c("AA","BB","AB","BA","?"))]
with(results, table(scenario,tested.cv))
setnames(results, gsub("PP.|.abf","",names(results)))
results[,rsq:=cut(rho^2,c(0,RSQ_THR,1),include.lowest=TRUE)]
hist(results$rho^2); abline(v=RSQ_THR,col="red")
results[,nper:=length(unique(index)),by=c("scenario","rsq")]
results[,indexper:=cumsum(!duplicated(index)),by=c("scenario","rsq")]
results[as.numeric(rsq)==1,.(number_run=max(indexper)),by=c("scenario","rsq")]
results[,method:=sub("iterative","it",method)]
results[,method:=sub("allbutone","abo",method)]

results=results[indexper < NPER_THR]
results[,number_run:=max(indexper),by=c("scenario","rsq")]
unique(results[as.numeric(rsq)==1,.(scenario,rsq,number_run,nper)])

results[,method:=sub("susie_0","susie",method)]
head(results)
## ggplot(results[method %in% c("cond_it","cond_abo")],
##        aes(x=H4)) +
##   geom_histogram() + facet_grid(scenario~method)

################################################################################

library(cowplot); theme_set(theme_minimal())
summ.h3=results[scenario=="sim_ncopies_1_A_0_B_0",.(single=max(H3*(method=="single")),
                                                    susie=max(H3*(method=="susie")),
                                                    z1=max(abs(hit1.margz)),
                                                    z2=max(abs(hit2.margz))),
             by=c("index")]
summ.h4=results[scenario=="sim_ncopies_1_A_1_B_-1",.(single=max(H4*(method=="single")),
                                                     susie=max(H4*(method=="susie")),
                                                     z1=max(abs(hit1.margz)),
                                                     z2=max(abs(hit2.margz))),
             by=c("index")]
summ=rbind(summ.h3,summ.h4)
summ[,p:=pnorm(pmin(z1,z2),lower=FALSE)*2]
summ[,class:=paste(#ifelse(single>0,"single call","single nocall"),
                   ifelse(susie>0,"one or more","none"))]
ggplot(summ,aes(x=class,y=-log10(p))) + geom_boxplot(varwidth=TRUE) +
  labs(x="number of coloc-SuSiE comparisons",y="min (over traits) { max (over snps) -log10(p) }") +
  geom_hline(yintercept=-log10(5e-8),linetype="dashed") +
  ## scale_y_continuous(breaks=-log10(c(1e-8,1e-10,1e-15,1e-20)),labels=c(8,10,15,20)) +
  background_grid(major="y")
ggsave("figure-min-max-z.png",height=4,width=4)

## ggplot(summ,aes(x=pmin(z1,z2),y=pmax(z1,z2),col=susie==0)) + geom_point() + geom_density2d() +
##   facet_wrap(~ (susie==0) )

## fill in single when susie doesn't run

results.x=results[method=="susie"]
results.x=rbind(results.x,
                results[method=="single" & !(index %in% results.x$index)])
results.x$method="hybrid"
results=rbind(results,results.x)

################################################################################


library(viridis)
library(ggplot2)
library(ggpubr)
library(cowplot)

## ## plot average posterior by n
## kk <- results[!is.na(H0) & n<=4, .(H0=sum(H0),H1=sum(H1),H2=sum(H2),H3=sum(H3),H4=sum(H4)),
##               by=c("scenario","method","n","rsq")]
## m <- melt(kk,c("scenario","method","n","rsq"))
## m[,method:=factor(method,levels=c("single","cond","mask","susie"))]
## p.r <- ggplot(m,aes(x=n,y=value/(NPER*n),fill=variable)) + geom_col() +
##   facet_grid(method ~ scenario + rsq,space="free_x",scales="free_x") + theme_pubr() +
##   scale_fill_viridis_d("hypoth.") +
##   ylab("Avg. posterior") +
##   xlab("Tested variants") +
##   theme(legend.position="right",
##         panel.border=element_rect(linetype="solid",fill=NA)
##         ##, strip.text.y = element_blank()
##         )
## p.rr
kk <- results[!is.na(H0) & as.numeric(rsq)==1,
              .(H0=sum(H0,na.rm=TRUE),H1=sum(H1,na.rm=TRUE),H2=sum(H2,na.rm=TRUE),H3=sum(H3,na.rm=TRUE),H4=sum(H4,na.rm=TRUE)),
              by=c("scenario","method","tested.cv","rsq","number_run")]
m <- melt(kk,c("scenario","method","tested.cv","rsq","number_run"))
## m[,method:=factor(method,levels=c("single","cond","mask","susie_0","susie_0.5","susie_1","susie_1.5"))]
m[,method:=factor(method,levels=c("single","cond_it","cond_abo",#"mask_it","mask_abo",
                                  "susie","hybrid"))]
p.r <- ggplot(m,aes(x=tested.cv,y=value/number_run,fill=variable)) + geom_col() +
  facet_grid(method ~ scenario,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hypoth.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid(major="y") +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA)
        ##, strip.text.y = element_blank()
        )
## p.r
 ## facet_grid(vore ~ conservation, labeller = labeller(vore = capitalize))
m[,nsnps:=paste0(gsub("sim_ncopies_|_A.*","",scenario),"000")]
m[,scen:=sub(".*_A","A",scenario)]
m[,meth.nsnps:=paste(method,nsnps)]
m[,scen.nsnps:=paste(scen,nsnps)]
bottomline <- ggplot(m[rsq=="[0,0.5]" & nsnps==NCOPIES_TO_SHOW],
                     aes(x=tested.cv,y=value/number_run,fill=variable)) +
  geom_col() +
  facet_grid(scen ~ method,space="free_x",scales="free_x") +
  scale_fill_viridis_d("hypoth.") +
  scale_y_continuous("Avg. posterior",breaks=c(0,0.25,0.5,0.75,1),labels=c("0","","","","1"),
                     limits=c(0,1)) +
  xlab("Tested variant pairs") +
  theme_pubr(base_size=10) +
  background_grid(major="y") +
  theme(legend.position="bottom",
        panel.border=element_rect(linetype="solid",fill=NA,size=0.1),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        axis.text.y=element_text(),
        strip.text.y = element_blank()
        )
label_custom=function(labels,multi_line=TRUE) {
  substr(labels,nchar(labels)-4,nchar(labels))
}
bottomline2 <- ggplot(m[rsq=="[0,0.5]"&  nsnps==3000],aes(x=tested.cv,y=value/number_run,fill=variable)) +
  geom_col() +
  facet_grid(scen ~ method,space="free_x",scales="free_x") + #,
             ## labeller=labeller(scen.nsnps=label_custom)) +
  scale_fill_viridis_d("hypoth.") +
  scale_y_continuous("Avg. posterior",breaks=c(0,0.25,0.5,0.75,1),labels=c("0","","","","1"),
                     limits=c(0,1)) +
  xlab("Tested variant pairs") +
  theme_pubr(base_size=10) +
  background_grid(major="y") +
  theme(legend.position="bottom",
        panel.border=element_rect(linetype="solid",fill=NA,size=0.1),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        axis.text.y=element_text()
        ## strip.text.y = element_blank()
        )
  dt=dcast(
    m[ as.numeric(rsq)==1,.(scenario=scen,nsnps_in_region=nsnps,method,inferred_cv_pair=tested.cv,variable,avg_posterior=value/number_run)]
   ,scenario+nsnps_in_region+method+inferred_cv_pair ~ variable, value.var="avg_posterior")
dt[scenario=="A_0_B_0",scenario:="A-B"]
dt[scenario=="A_1_B_-1",scenario:="A-A"]
dt[scenario=="A_1_B_0",scenario:="A-AB"]
dt[scenario=="A_1_B_1",scenario:="AB-AB"]
fwrite(dt,file="table-coloc-sims.csv")

## bottomline

## plot number of tests per result
kk=results[,.(n=.N),by=c("index","scenario","method")]
table(kk$n)
kk[n>4,n:=5]
bottom_n=ggplot(kk, aes(x=n)) +
  geom_histogram(binwidth=1,fill="grey40",colour="grey80") +
  scale_x_continuous(breaks=c(1,2,3,4,5),labels=c("1","2","3","4","5+")) +
  facet_grid(method ~ scenario) + theme_pubr() +
  background_grid(major="y") +
  theme(legend.position="right",
        ## panel.border=element_rect(linetype="solid",fill=na)
        , strip.text.x = element_blank()
        )



################################################################################

## four plots
RELATIONS=list(sim_0_0=data.table(to=c("trait 1","trait 2"),from=c("A","B")),
  "sim_1_-1"=data.table(to=c("trait 1","trait 2"),from=c("A","A")),
               sim_1_0=data.table(to=c("trait 1","trait 2","trait 2"),from=c("A","A","B")),
               sim_1_1=data.table(to=c("trait 1","trait 2","trait 1","trait 2"),
                                  from=c("A","A","B","B")))
relations=RELATIONS[[3]]

library(cowplot)
theme_set(theme_cowplot(font_size=28))
## theme_set(theme_minimal())
library(seaborn)
simCols <- seaborn:::SEABORN_PALETTES$pastel[c(6,8,1)] # 8 is grey
library(igraph)
nodes=data.table(name=c("A","B","trait 1","trait 2"),
                 x=c(-1,-1,1,1),
                 y=c(1,-1,1,-1),
                 col=simCols[c(1,1,3,3)], #cols[c(2,10,10,10,10,2,4)],
                 class=c("variant","variant","trait","trait"),
                 StringsAsFactors=FALSE)

plotter=function(tag) {
  relations=RELATIONS[[tag]]
  rdf=merge(relations,nodes[,.(name,x,y)],by.x="from",by.y="name")
  rdf <- merge(rdf,nodes[,.(name,x,y,col)],by.x="to",by.y="name",suffixes=c(".from",".to"))
  rdf[,col:=simCols[2]][,row:="row"][,column:="column"]
  cscale <- structure(unique(c(nodes$col,rdf$col)),names=unique(c(nodes$col,rdf$col)))
  ggplot(nodes, aes(x=x,y=y)) +
    geom_point(aes(colour=col,size=class),pch=20) +
    geom_segment(aes(x=x.from+0.3,y=y.from,xend=x.to-0.3,yend=y.to,col=col),data=rdf,
                 size=0.2,col="black",
                 arrow=arrow(length = unit(0.05, "npc"), type="closed")) +
    geom_text(aes(label=name),size=2) +
    ## xlim(0,3) + ylim(1.9,3.5) +
    scale_colour_manual(values=cscale) +
    scale_size_manual(values=c(trait=12,variant=8)) +
    xlim(-1.6,1.6) + ylim(-1.8,1.8) +
    facet_grid(row ~ column) +
    ## ggtitle(patterns$match[i]) +
    ## theme_pubr() +
    theme(legend.position="none",axis.line=element_blank(),
          axis.title=element_blank(),
        panel.border=element_rect(linetype="solid",fill=NA),
         strip.text.x = element_blank(),
         strip.text.y = element_blank(),
        axis.ticks=element_blank(), axis.text=element_blank())
    ## coord_flip()
}

topplots=lapply(names(RELATIONS),plotter)
w=3
topline=plot_grid(plotlist=list(NULL,topplots[[1]],NULL,topplots[[2]],NULL,topplots[[3]],NULL,topplots[[4]],NULL),
                  ncol=1,
                  ## labs=c("","a","","b","","c","","d",""),
                  rel_heights=c(1,w,1,w,1,w,1,w,3))
## plot_grid(topline,bottomline,bottom_n,ncol=1,rel_heights=c(.2,.8,.5))
both=plot_grid(topline,bottomline + theme(strip.background=element_blank(),strip.text=element_text(face="bold")),rel_widths=c(.2,.8),nrow=1)
both2=plot_grid(topline,bottomline2 + theme(strip.background=element_blank(),strip.text=element_text(face="bold")),rel_widths=c(.2,.8),nrow=1)
## topline

library(grid)
lines <- lapply(c(0.356,0.561,0.757), function(Y) {
  linesGrob(
  x = unit(c(0,1),"npc"),
  y = unit(c(Y,Y),"npc"),
  gp = gpar(lty=2,col="grey")
  )
})

## bad1=results[method=="cond_allbutone" & H4>0.8 & rB.2>0.8 & scenario=="sim_ncopies_1_A_1_B_0"]$index
## bad2=results[method=="cond_iterative" & H4>0.8 & rB.2>0.8 & scenario=="sim_ncopies_1_A_1_B_0"]$index
## intersect(bad1,bad2)

final=ggdraw(both) #+
  ## draw_grob(lines[[1]]) +
  ## draw_grob(lines[[2]]) +
  ## draw_grob(lines[[3]])
## both
## final
ggsave("figure-coloc.png",plot=final,height=6,width=8)

final2=ggdraw(both2) # +
  ## draw_grob(lines[[1]]) +
  ## draw_grob(lines[[2]]) +
  ## draw_grob(lines[[3]])
## both
ggsave("figure-coloc-supp.png",plot=final2,height=8,width=8)

################################################################################

## PP
f=function(x) {
  l=sapply(x, length)
  if(any(l>1))
    stop("l>1! Panic!")
  ret=numeric(length(x))
  ret[l==0]=NA
  ret[l==1]=unlist(x[l==1])
  ret
}
results[,trait1.pp.A:=f(trait1.pp.A)]
results[,trait2.pp.A:=f(trait2.pp.A)]
results[,trait1.pp.B:=f(trait1.pp.B)]
results[,trait2.pp.B:=f(trait2.pp.B)]
results=results[order(method,tested.cv,coloc.pp.A)]
results[H4>0.9,xA:=1:.N,by=c("method","tested.cv")]
results=results[order(method,tested.cv,coloc.pp.B)]
results[H4>0.9,xB:=1:.N,by=c("method","tested.cv")]
results[,lower.pp.A:=pmin(trait1.pp.A,trait2.pp.A)]
results[,higher.pp.A:=pmax(trait1.pp.A,trait2.pp.A)]
results[,lower.pp.B:=pmin(trait1.pp.B,trait2.pp.B)]
results[,higher.pp.B:=pmax(trait1.pp.B,trait2.pp.B)]
results[,max.trait.pp:=pmax(higher.pp.A,higher.pp.B)]
results[,max.coloc.pp:=pmax(coloc.pp.A,coloc.pp.B)]

library(ggplot2)
library(cowplot)
library(ggnewscale)
theme_set(theme_minimal())
## dt[scenario=="A_0_B_0",scenario:="A-B"]
## dt[scenario=="A_1_B_-1",scenario:="A-A"]
## dt[scenario=="A_1_B_0",scenario:="A-AB"]
## dt[scenario=="A_1_B_1",scenario:="AB-AB"]

thr=0.9
rA=results[H4>thr & nsnps==1000 ,
           .(coloc.pp=max(coloc.pp.A),max.trait.pp=max(higher.pp.A)),
           by=c("index","nsnps","method","scenario")]
rA[,shared:=!grepl("A_0_B_0",scenario)]
rB=results[H4>thr & nsnps==1000 ,
           .(coloc.pp=max(coloc.pp.B),max.trait.pp=max(higher.pp.B)),
           by=c("index","nsnps","method","scenario")]
rB[,shared:=grepl("A_1_B_1",scenario)]
rs=rbind(rA,rB)
rs[,num_shared:="shared CV: 1"]
rs[grepl("A_0_B_0",scenario),num_shared:="shared CV: 0"]
rs[grepl("A_1_B_1",scenario),num_shared:="shared CV: 2"]
rs[,method:=factor(method,levels=c("single","cond_abo","cond_it","susie","hybrid"))]
rsumm=rs[shared==TRUE,
         .(pc=paste0(signif(100*mean(coloc.pp >= max.trait.pp,na.rm=TRUE),3),"%")),
                     ## "avg. diff: ",signif(mean((max.coloc.pp - max.trait.pp)/max.trait.pp,na.rm=TRUE),3))),
         by=c("shared","method","num_shared")]

ggplot(rs[shared==TRUE], aes(x=max.trait.pp,y=coloc.pp)) +
  geom_point(size=0.5) +
  geom_abline() +
  ## geom_smooth(se=FALSE) +
  geom_label(aes(label=pc),x=0.15,y=0.85,data=rsumm,col="red",size=3.5) +
  ## geom_boxplot(varwidth=TRUE) +
  ## scale_colour_manual("source",values=c(lower="grey",higher="lightblue",coloc="red")) +
  ## labs(x="Tested pair",y="Posterior probability") +
  ## new_scale("color") +
  ## geom_smooth(se=FALSE) +
  ## scale_colour_manual("lines",values=c(lower="lightblue",higher="lightblue",coloc="#ff000000")) +
  scale_x_continuous(breaks=c(0,1),labels=c("0","1")) +
  scale_y_continuous(breaks=c(0,.5,1)) +
  facet_grid( num_shared ~ method) +
  theme(legend.position="bottom",strip.text=element_text(face="bold")) +
  background_grid(major="y")
ggsave("fig-pp.png",height=2*8/5,width=8)

## rs=results[H4>0.5 & nsnps==1000,.(max.coloc.pp=max(max.coloc.pp),max.trait.pp=max(max.trait.pp)),by=c("index","nsnps","method")]
## rsumm=rs[,.(pc=paste0(signif(100*mean(max.coloc.pp >= max.trait.pp,na.rm=TRUE),3),"%")),by="method"]
## ggplot(rs, aes(x=max.trait.pp,y=max.coloc.pp)) +
##   geom_point(size=1) +
##   geom_abline() +
##   ## geom_smooth(se=FALSE) +
##   geom_label(aes(label=pc),x=0.25,y=0.75,data=rsumm,col="red") +
##   ## geom_boxplot(varwidth=TRUE) +
##   ## scale_colour_manual("source",values=c(lower="grey",higher="lightblue",coloc="red")) +
##   ## labs(x="Tested pair",y="Posterior probability") +
##   ## new_scale("color") +
##   ## geom_smooth(se=FALSE) +
##   ## scale_colour_manual("lines",values=c(lower="lightblue",higher="lightblue",coloc="#ff000000")) +
##   facet_wrap(~ method) +
##   theme(legend.position="bottom") +
##   background_grid(major="y")

## mA=melt(results[nsnps %in% c(1000,3000) & H4>0.9 & tested.cv %in% c("AA","AB","BA","BB")],
##         id.vars=c("nsnps","H4","tested.cv","method"),
##         measure.vars=c("lower.pp.A","higher.pp.A","coloc.pp.A"))
## mA[,variable:=c("lower","higher","coloc")[variable]]
## mA[,"H4.accurate":=tested.cv=="AA"]
## mB=melt(results[nsnps %in% c(1000,3000) & H4>0.9 & tested.cv %in% c("BB","AB","BA","AA")],
##         id.vars=c("nsnps","H4","tested.cv","method"),
##         measure.vars=c("lower.pp.B","higher.pp.B","coloc.pp.B"))
## mB[,variable:=c("lower","higher","coloc")[variable]]
## mB[,"H4.accurate":=tested.cv=="BB"]
## mA[,cv:="A"]
## mB[,cv:="B"]
## m=rbind(mA,mB)
## m[,tested.pair:=ifelse(H4.accurate, "AA or BB", "AB or BA")]
## ## m[method=="susie_0",method:="susie"]
## head(m)

## ## library(ggplot2)
## ## library(cowplot)
## ## library(ggnewscale)
## ## theme_set(theme_cowplot())
## ## ggplot(m, aes(x=xA,y=pp.A,col=factor(variable))) +
## ##   geom_point(size=1,alpha=0.1) +
## ##   scale_colour_manual("points",values=c(lower="grey",higher="grey",coloc="red")) +
## ##   ## new_scale("color") +
## ##   geom_smooth(se=FALSE) +
## ##   ## scale_colour_manual("lines",values=c(lower="lightblue",higher="lightblue",coloc="#ff000000")) +
## ##   facet_grid(method ~ tested.cv,scales="free_x",space="free_x")

## library(ggplot2)
## library(cowplot)
## library(ggnewscale)
## theme_set(theme_cowplot())
## ggplot(m[nsnps==1000], aes(x=tested.cv,y=value,col=factor(variable))) +
##   geom_boxplot(varwidth=TRUE) +
##   scale_colour_manual("source",values=c(lower="grey",higher="lightblue",coloc="red")) +
##   labs(x="Tested pair",y="Posterior probability") +
##   ## new_scale("color") +
##   ## geom_smooth(se=FALSE) +
##   ## scale_colour_manual("lines",values=c(lower="lightblue",higher="lightblue",coloc="#ff000000")) +
##   facet_grid(cv~ method) +
##   theme(legend.position="bottom") +
##   background_grid(major="y")

## ggsave("figure-pp.png",height=8,width=8)

## library(ggplot2)
## library(cowplot)
## library(ggnewscale)
## theme_set(theme_cowplot())
## ggplot(m, aes(x=tested.pair,y=value,col=factor(variable))) +
##   geom_boxplot(varwidth=TRUE) +
##   scale_colour_manual("source",values=c(lower="grey",higher="lightblue",coloc="red")) +
##   labs(x="Tested pair",y="Posterior probability") +
##   ## new_scale("color") +
##   ## geom_smooth(se=FALSE) +
##   ## scale_colour_manual("lines",values=c(lower="lightblue",higher="lightblue",coloc="#ff000000")) +
##   facet_grid(nsnps~ method) +
##   theme(legend.position="bottom") +
##   background_grid(major="y")
## ggsave("figure-pp.png",height=8,width=8)
