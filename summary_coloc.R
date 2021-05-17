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
save(DATA,file="~/coloc-susie.RData")

load(file="~/coloc-susie.RData")

RSQ_THR=0.5
NPER_THR=10000 # because every dataset has at least 9660
results=rbindlist(DATA)
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

head(results)
ggplot(results[method %in% c("cond_it","cond_abo")],
       aes(x=H4)) +
  geom_histogram() + facet_grid(scenario~method)

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
              .(H0=sum(H0),H1=sum(H1),H2=sum(H2),H3=sum(H3),H4=sum(H4)),
              by=c("scenario","method","tested.cv","rsq","number_run")]
m <- melt(kk,c("scenario","method","tested.cv","rsq","number_run"))
## m[,method:=factor(method,levels=c("single","cond","mask","susie_0","susie_0.5","susie_1","susie_1.5"))]
m[,method:=factor(method,levels=c("single","cond_it","cond_abo","mask_it","mask_abo","susie_0","susie_0.5","susie_1","susie_1.5"))]
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
bottomline <- ggplot(m[rsq=="[0,0.5]"& nsnps==NCOPIES_TO_SHOW],aes(x=tested.cv,y=value/number_run,fill=variable)) +
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
bottomline2 <- ggplot(m[rsq=="[0,0.5]"& nsnps!=NCOPIES_TO_SHOW],aes(x=tested.cv,y=value/number_run,fill=variable)) +
  geom_col() +
  facet_grid(scen.nsnps ~ method,space="free_x",scales="free_x",
             labeller=labeller(scen.nsnps=label_custom)) +
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
    m[as.numeric(rsq)==1,.(scenario=scen,nsnps_in_region=nsnps,method,inferred_cv_pair=tested.cv,variable,avg_posterior=value/number_run)]
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
both=plot_grid(topline,bottomline,rel_widths=c(.2,.8),nrow=1)
both2=plot_grid(topline,bottomline2,rel_widths=c(.2,.8),nrow=1)
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

final2=ggdraw(both2) +
  draw_grob(lines[[1]]) +
  draw_grob(lines[[2]]) +
  draw_grob(lines[[3]])
## both
ggsave("figure-coloc-supp.png",plot=final2,height=8,width=8)
