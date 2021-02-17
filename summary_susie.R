#!/home/cew54/localc/bin/Rscript
setwd("~/Projects/coloc-susie")
source("dirs.rb")
library(randomFunctions)
library(magrittr)
## args=getArgs(defaults=list(chr=2, begin=43082452, end=44086664), numeric=c("chr","begin","end"))
source("~/E/sim_colours.R")

## example with high LD between CVs (0.0.99)
## f="data/sim/coloc/_223c973ddcaff_255f719b075c4.RData"
library(pbapply)
dirs=list.files(SIMCOLOC,pattern="approx_",full=TRUE)
d=dirs[1]
for(d in dirs) {
  files=list.files(d,full=TRUE) %>% grep("summary",.,value=TRUE,invert=TRUE)
  message(basename(d),"\tfiles found: ",length(files))
}

library(ggplot2)
library(ggrastr)
library(data.table)

proc_file=function(f,d) {
  (load(f))
  result[,scenario:=basename(d)][,index:=paste(f,1:.N,sep="_")]
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


NPER_THR=1000
results=rbindlist(DATA)
with(results, table(scenario))
results[,nsnp:=as.numeric(sub("approx_","",scenario))]
head(results)


results[,pip.cv1:=sub("/.*","",cv.pip) %>% as.numeric()]
results[ncv==2,pip.cv2:=sub(".*/","",cv.pip) %>% as.numeric()]
results[,base_pip.cv1:=sub("/.*","",base_cv.pip) %>% as.numeric()]
results[ncv==2,base_pip.cv2:=sub(".*/","",base_cv.pip) %>% as.numeric()]

results[,z.cv1:=sub("/.*","",cv.z) %>% as.numeric()]
results[nsnp==1,z.cv2:=sub(".*/","",cv.z) %>% as.numeric()]
results[,.(n=.N,prop.same=mean(same_top)),by=c("ncv","nsnp","zthr")]

results[nsnp < 1000, nsnp:=1000 * nsnp]

m=melt(results,c("nsnp","ncv","zthr"),list(c("pip.cv1","pip.cv2"),c("z.cv1","z.cv2"),c("base_pip.cv1","base_pip.cv2")))
setnames(m,c("variable","value1","value2","value3"),c("cv","pip","z","base_pip"))
with(m,table(ncv,cv))
m=m[!is.na(pip) & as.numeric(cv) <= as.numeric(ncv)]
m[,pip_diff:=pip-base_pip]
head(m)

m[,cvscen:=ifelse(ncv==1, "single CV", paste("two CV, variant",cv))]
m[,ncvscen:=ifelse(ncv==1, "single CV", "two CV")]
results[,ncvscen:=ifelse(ncv==1, "single CV", "two CV")]
m[,snpscen:=paste(nsnp," SNPs")]
results[,snpscen:=paste(nsnp," SNPs")]
## ggplot(m,aes(x=base_pip,y=pip,col=factor(zthr),size=abs(z))) +
##   geom_point() +
##   geom_smooth() +
##   facet_grid(nsnp ~ ncv + cv, labeller=label_both)


## ggplot(results,aes(x=base_pnull,y=base_pnull+pnull_diff,col=factor(zthr))) +
##   geom_point() +
##   geom_smooth() +
##   facet_grid(nsnp ~ ncv, labeller=label_both)
library(cowplot); theme_set(theme_cowplot(font_size=10))
library(seaborn)
library(ggridges)
plots=list(
  pnull=ggplot(results, aes(y=factor(zthr),x=pnull_diff,fill=factor(zthr),col=factor(zthr))) +
    geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
  labs(x="Pr(null) in trimmed - Pr(null) in full data")
,time= ggplot(results, aes(y=factor(zthr),x=time_diff+base_time,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
   labs(x="Time to run SuSiE")
,max_abs_diff=ggplot(results, aes(y=factor(zthr),x=max_abs_diff,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
   labs(x="Maximum absolute pip difference")
 ,pipcr=ggplot(results, aes(y=factor(zthr),x=pipcr,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=1,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
    labs(x="Correlation between PIP")
 ,pip=ggplot(m, aes(y=factor(zthr),x=pip_diff,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
    labs(x="Difference in PIP at causal variants"))


plotsplus = lapply(plots,function(p)
  p + geom_density_ridges(stat = "binline", scale=0.95,bins=20, draw_baseline = FALSE) +
  scale_fill_seaborn("|Z| threshold for trimming") +
  scale_colour_seaborn("|Z| threshold for trimming") +
  background_grid(major="y") +
  ylab("|Z| threshold") +
  theme(legend.position="bottom") +
  facet_grid(ncvscen ~ snpscen,space="free"))

plotsplus$pip = plotsplus$pip +
  facet_grid(cvscen ~ snpscen, space="free") +
  scale_x_continuous(breaks=c(-1,0,1))

## stand alone
## plotsplus$time
plotsplus$pip
ggsave("~/Projects/coloc-susie/figure-pipdiff.png",plot=plotsplus$pip,height=6,width=6)

## other diagnostics
plotlist=plotsplus[c("pnull","max_abs_diff","pipcr")]
plotlist[1:2] %<>% lapply(., function(p) p + theme(legend.position="none"))
plot_grid(plotlist=plotlist[c(1,3)],ncol=1)
ggsave("~/Projects/coloc-susie/figure-otherdiagnostics.png",height=8,width=6)

## alternative time
## ggplot(results, aes(x=factor(nsnp),y=base_time + time_diff,col=factor(zthr))) +
##   geom_boxplot() +
##   facet_wrap(~ncv,labeller=label_both) +
##   background_grid(major="y")

## alternative time
mtime=results[,.(y=median(base_time+time_diff),ymax=quantile(base_time+time_diff,.75),ymin=quantile(base_time+time_diff,.25)),
              by=c("zthr","nsnp","ncvscen")]
ggplot(mtime, aes(x=nsnp,y=y,ymin=ymin,ymax=ymax,col=factor(zthr))) +
  geom_pointrange() +
  geom_path(aes(group=factor(zthr))) +
  facet_wrap(~ncvscen) +
  labs(x="Number of SNPs in region",y="Time (seconds)") +
  scale_colour_seaborn("|Z| threshold for trimming") +
  theme(legend.position="bottom") +
  background_grid(major="y")
ggsave("~/Projects/coloc-susie/figure-time.png",height=4,width=6)

################################################################################

## all abandonned below here


library(ggridges)
ggplot(results, aes(y=factor(zthr),x=max_abs_diff,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
  geom_density_ridges(stat = "binline", scale=0.95,bins=20, draw_baseline = FALSE,
                      ## point_size = 0.4, point_alpha = 1,
                      ## position = position_raincloud(adjust_vlines = TRUE)) +
                      ) +
  scale_fill_seaborn("|Z| threshold for trimming") +
  scale_colour_seaborn("|Z| threshold for trimming") +
  background_grid(major="y") +
  labs(x="Maximum absolute difference in ppi",y="|Z|") +
  background_grid(major="y") +
  facet_grid(nsnp ~ ncv, labeller=label_both,space="free",scales="free_x")

library(ggridges)

ggplot(results, aes(x=factor(nsnp*1000),y=max_abs_diff,group=factor(nsnp))) +
  geom_boxplot() +
  facet_grid(ncv ~ zthr, labeller=label_both)

ggplot(results, aes(x=base_pnull,y=pnull_diff+base_pnull)) +
  geom_abline(),
  geom_point(alpha=0.1) + geom_density2d() +
  facet_grid(ncv + nsnp ~ zthr, labeller=label_both)

ggplot(results, aes(x=base_time,y=time_diff+base_time)) +
  geom_abline() +
  geom_point(alpha=0.1) + geom_density2d() +
  facet_grid(ncv + nsnp ~ zthr, labeller=label_both)


m=melt(results, measure.vars=list(pnull=c("pnull.trim","pnull.notrim"),
                                  time=c("time.trim.elapsed","time.notrim.elapsed")))
m[,trim:=c("trim","notrim")[variable]]
head(m)
tail(m)

library(viridis)
library(ggplot2)
library(ggpubr)
library(cowplot); theme_set(theme_cowplot(font_size=10,line_size=0.3))
library(seaborn)

ggplot(m, aes(x=factor(nsnps),y=time,col=trim)) +
  geom_boxplot() +
  labs(x="Number of SNPs in region",y="Elapsed time (seconds)")

ggplot(results, aes(x=pnull.notrim,y=pnull.trim,col=factor(ncv))) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~factor(nsnps)) +
  labs(x="P(null), full data",y="P(null), trimmed data")

ggplot(results, aes(x=as.factor(ncv),y=cor,col=factor(nsnps))) +
  geom_boxplot() +
  scale_colour_seaborn("Number of SNPs in region") +
  labs(x="Number of causal SNPs in region",y="Correlation between PPI") +
  theme(legend.position="bottom")

ggplot(results, aes(x=as.factor(ncv),y=MSE,col=factor(nsnps))) +
  geom_boxplot(outlier.size=0.5) +
  scale_colour_seaborn("Number of SNPs in region") +
  labs(x="Number of causal SNPs in region",y="Mean sq error between PPI") +
  theme(legend.position="bottom")

ggplot(results, aes(x=as.factor(ncv),y=nsnps*MSE,col=factor(nsnps))) +
  geom_boxplot(outlier.size=0.5) +
  scale_colour_seaborn("Number of SNPs in region") +
  labs(x="Number of causal SNPs in region",y="Sum of sq error between PPI") +
  theme(legend.position="bottom")

ggsave("~/susie_approx.png",width=6,height=6)

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
kk <- results[!is.na(H0), .(H0=sum(H0),H1=sum(H1),H2=sum(H2),H3=sum(H3),H4=sum(H4)),
              by=c("scenario","method","tested.cv","rsq")]
m <- melt(kk,c("scenario","method","tested.cv","rsq"))
m[,method:=factor(method,levels=c("single","cond","mask","susie"))]
p.r <- ggplot(m,aes(x=tested.cv,y=value/NPER,fill=variable)) + geom_col() +
  facet_grid(method ~ scenario + rsq,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hypoth.") +
  ylab("Avg. posterior") +
  xlab("Tested variants") +
  background_grid(major="y") +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA)
        ##, strip.text.y = element_blank()
        )
## p.r

bottomline <- ggplot(m[rsq=="[0,0.5]"],aes(x=tested.cv,y=value/NPER,fill=variable)) + geom_col() +
  facet_grid(method ~ scenario,space="free_x",scales="free_x") + theme_pubr() +
  scale_fill_viridis_d("hypoth.") +
  ylab("Avg. posterior") + ylim(0,1) +
  xlab("Tested variant pairs") +
  background_grid(major="y") +
  theme(legend.position="right",
        panel.border=element_rect(linetype="solid",fill=NA)
        , strip.text.x = element_blank()
        )
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
        panel.border=element_rect(linetype="solid",fill=NA)
        , strip.text.x = element_blank()
        )


################################################################################

## four plots
RELATIONS=list(sim_0_0=data.table(from=c("trait 1","trait 2"),to=c("A","B")),
  "sim_1_-1"=data.table(from=c("trait 1","trait 2"),to=c("A","A")),
               sim_1_0=data.table(from=c("trait 1","trait 2","trait 2"),to=c("A","B","A")),
               sim_1_1=data.table(from=c("trait 1","trait 2","trait 1","trait 2"),
                                  to=c("A","A","B","B")))

library(cowplot)
theme_set(theme_cowplot(font_size=28))
library(seaborn)
simCols <- seaborn:::SEABORN_PALETTES$pastel[c(6,8,1)] # 8 is grey
library(igraph)
nodes=data.table(name=c("A","B","trait 1","trait 2"),
                 x=c(-1,1,-1,1),
                 y=c(-1,-1,1,1),
                 col=simCols[c(1,1,3,3)], #cols[c(2,10,10,10,10,2,4)],
                 class=c("variant","variant","trait","trait"),
                 stringsAsFactors=FALSE)

plotter=function(tag) {
  relations=RELATIONS[[tag]]
  rdf=merge(relations,nodes[,.(name,x,y)],by.x="from",by.y="name")
  rdf <- merge(rdf,nodes[,.(name,x,y,col)],by.x="to",by.y="name",suffixes=c(".from",".to"))
  rdf[,col:=simCols[2]][,row:="row"][,column:="column"]
  cscale <- structure(unique(c(nodes$col,rdf$col)),names=unique(c(nodes$col,rdf$col)))
  ggplot(nodes, aes(x=x,y=y)) +
    geom_segment(aes(x=x.from,y=y.from,xend=x.to,yend=y.to,col=col),data=rdf,size=1) +
    geom_point(aes(colour=col,size=class),pch=20) +
    geom_text(aes(label=name)) +
    ## xlim(0,3) + ylim(1.9,3.5) +
    scale_colour_manual(values=cscale) +
    scale_size_manual(values=c(trait=22,variant=12)) +
    xlim(-1.6,1.6) + ylim(-1.2,1.5) +
    facet_grid(row ~ column) +
    ## ggtitle(patterns$match[i]) +
    theme(legend.position="none",axis.line=element_blank(),
          axis.title=element_blank(),
        panel.border=element_rect(linetype="solid",fill=NA),
         strip.text.x = element_blank(),
         strip.text.y = element_blank(),
          axis.ticks=element_blank(), axis.text=element_blank())
}

topplots=lapply(names(RELATIONS),plotter)
w=1.5
topline=plot_grid(plotlist=list(NULL,topplots[[1]],NULL,topplots[[2]],NULL,topplots[[3]],NULL,topplots[[4]],NULL),
                  nrow=1,rel_widths=c(1,w,1,w,1,w,2,w,2.5))
## plot_grid(topline,bottomline,bottom_n,ncol=1,rel_heights=c(.2,.8,.5))
both=plot_grid(topline,bottomline,rel_heights=c(.2,.8),ncol=1)
## both
ggsave("~/coloc_susie.png",both,height=7,width=14)
