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
}


d=dirs[1]
DATA=vector("list",length(dirs))
names(DATA)=basename(dirs)
for(d in dirs) {
  files=list.files(d,full=TRUE) %>% grep("summary",.,value=TRUE,invert=TRUE)
  ## if(length(files)>200)
  ##   files=files[1:200]
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

results[,base_nsets.out:=base_nsets-base_nsets.in]
results[,nsets.out:=nsets-nsets.in]
results[,diff_nsets.in:=nsets.in-base_nsets.in]
results[,diff_nsets.out:=nsets.out-base_nsets.out]
results[,pip.cv1:=sub("/.*","",cv.pip) %>% as.numeric()]
results[ncv==2,pip.cv2:=sub(".*/","",cv.pip) %>% as.numeric()]
results[,base_pip.cv1:=sub("/.*","",base_cv.pip) %>% as.numeric()]
results[ncv==2,base_pip.cv2:=sub(".*/","",base_cv.pip) %>% as.numeric()]

results[,z.cv1:=sub("/.*","",cv.z) %>% as.numeric()]
results[nsnp==1,z.cv2:=sub(".*/","",cv.z) %>% as.numeric()]
results[,.(n=.N,prop.same=mean(same_top)),by=c("ncv","nsnp","zthr")]
results[nsnp < 1000, nsnp:=1000 * nsnp]
results[,ncvscen:=ifelse(ncv==1, "single CV", "two CV")]
results[,snpscen:=paste(nsnp," SNPs")]

library(cowplot)
theme_set(theme_cowplot(font_size=10))
results[,sets_same:=nsets.in==base_nsets.in & nsets.out==nsets.out]

with(results, table(snpscen,ncvscen,zthr))

ggplot(results, aes(x=pipcr,fill=factor(zthr),col=factor(zthr))) +
  geom_histogram(position="dodge",binwidth=0.05) +
  scale_fill_seaborn("|Z| threshold for trimming") +
  scale_colour_seaborn("|Z| threshold for trimming") +
  background_grid(major="y") +
  labs(x="Correlation between PIP before and after trimming",y="Count") +
  theme(legend.position="bottom") +
  ## scale_y_log10() +
  ## scale_colour_seaborn("|Z| threshold for trimming") +
  facet_grid(snpscen ~ ncvscen)
ggsave("~/Projects/coloc-susie/figure-pipcr.png",height=4,width=6)

## alternative time
mtime=results[,.(y=median(base_time+time_diff),ymax=quantile(base_time+time_diff,.75),ymin=quantile(base_time+time_diff,.25)),
              by=c("zthr","nsnp","ncvscen")]
mtime=rbind(mtime, results[zthr==0.5,.(y=median(base_time), ymax=quantile(base_time,.75),ymin=quantile(base_time,.25),zthr=0),
                           by=c("nsnp","ncvscen")])
ggplot(mtime, aes(x=nsnp,y=y,ymin=ymin,ymax=ymax,col=factor(zthr))) +
  geom_pointrange() +
  geom_path(aes(group=factor(zthr))) +
  facet_wrap(~ncvscen) +
  scale_x_continuous(breaks=c(1000,2000,3000)) +
  labs(x="Number of SNPs in region",y="Time (seconds)") +
  scale_colour_seaborn("|Z| threshold for trimming") +
  theme(legend.position="bottom") +
  background_grid(major="y")
ggsave("~/Projects/coloc-susie/figure-time.png",height=4,width=6)

################################################################################

## all abandonned below here

m=rbind(results[,.(base_pip=base_pip.cv1,pip=pip.cv1,zthr,pipcr,snpscen,
                   ncvscen=ifelse(ncvscen=="single CV", ncvscen, paste(ncvscen, "(CV 1)")))],
                results[ncvscen!="single CV",.(base_pip=base_pip.cv2,pip=pip.cv2,zthr,pipcr,snpscen,
                   ncvscen=ifelse(ncvscen=="single CV", ncvscen, paste(ncvscen, "(CV 2)")))])

with(m, table(snpscen,ncvscen,zthr))


ggplot(m, aes(x=base_pip,y=pip)) +
  geom_abline() +
  geom_point(alpha=0.1) +
  facet_grid(factor(zthr) ~ ncvscen + snpscen)

## ggplot(m, aes(x=base_pip,y=pip)) +
##   geom_abline() +
##   geom_hex() +
##   facet_grid(factor(zthr) ~ ncvscen + snpscen)

## ggplot(results, aes(x=base_pip.cv1,y=pip.cv1)) + #,col=factor(nsets.in==base_nsets.in & nsets.out==nsets.out)
##   geom_abline(intercept=0,slope=1) +
##   geom_point(alpha=0.1) +
##   facet_grid(factor(zthr) ~ ncvscen + snpscen)

## ggplot(results, aes(x=base_pip.cv2,y=pip.cv2)) + #,col=factor(nsets.in==base_nsets.in & nsets.out==nsets.out)
##   geom_abline(intercept=0,slope=1) +
##   geom_point(alpha=0.1) +
##   facet_grid(factor(zthr) ~ ncvscen + snpscen)

## ggplot(results, aes(x=max_pos_diff.noncv,y=min_neg_diff.noncv)) +#,col=factor(nsets.in==base_nsets.in & nsets.out==nsets.out))) +
##   geom_abline() +
##   geom_point(alpha=0.1) +
##   facet_grid(factor(zthr) + factor(sets_same) ~ ncvscen + snpscen)

## results[,.(mse1=mean((pip.cv1-base_pip.cv1)^2),
##            mse2=mean((pip.cv2-base_pip.cv2)^2)), by=c("zthr","ncvscen","snpscen")]

## ggplot(results,aes(x=base_tot+tot_diff)) + geom_histogram() + facet_wrap(~zthr)

## ggplot(results,aes(x=base_tot,y=tot_diff)) + geom_hex() + geom_smooth() +
##   facet_wrap(~zthr)

## fill barchart
## library(cowplot)
## library(seaborn)
## m=melt(results, measure.vars=c("diff_nsets.in","diff_nsets.out"))
## head(m)
## ggplot(m, aes(x=value,fill=as.factor(zthr))) +
##   geom_bar(position="dodge") +
##   scale_fill_seaborn("|Z| threshold for trimming") +
##   scale_colour_seaborn("|Z| threshold for trimming") +
##   background_grid(major="y") +
##   facet_grid(ncvscen + variable ~ snpscen)

## ggplot(m, aes(x=value,fill=variable)) +
##   geom_bar(position="dodge") +
##   scale_fill_seaborn("|Z| threshold for trimming") +
##   scale_colour_seaborn("|Z| threshold for trimming") +
##   background_grid(major="y") +
##   facet_grid(ncvscen + factor(zthr) ~ snpscen)

## ggplot(results, aes(x=pip.cv1-base_pip.cv1,fill=factor(zthr))) +
##   geom_histogram(position="dodge",aes(y=..density..)) +
##   scale_fill_seaborn("|Z| threshold for trimming") +
##   scale_colour_seaborn("|Z| threshold for trimming") +
##   background_grid(major="y") +
##   facet_grid(ncvscen + factor(diff_nsets.in==0) ~ snpscen)


## bivariate plot
library(cowplot)
library(ggExtra)
library(viridis)
p=ggplot(results[zthr>0], aes(x=pcv_diff,y=pnocv_diff)) + #,col=pnocv_diff)) +
  ## geom_hex(binwidth=c(0.01/21,2/11),colour="grey") +
  ## geom_bin2d(binwidth=c(0.018/17,4/19),colour="grey") +
  ## geom_hline(yintercept=0) + geom_vline(xintercept=0,width=0.1) +
  geom_jitter(alpha=0.01) +
  ## geom_point(colour="transparent") +
  ## geom_density2d(alpha=0.2) +
  background_grid() +
  geom_rug(aes(x=pcv_diff,y=pcv_diff),alpha=0.01,data=results) +
  scale_colour_viridis() +
  scale_fill_viridis() +
  facet_grid(ncvscen + zthr ~ snpscen)
p


## https://www.r-bloggers.com/2014/05/finding-the-midpoint-when-creating-intervals/
midpoints <- function(x){
lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
lower+(upper-lower)/2
}

min(results$pnocv_diff)
max(results$pnocv_diff)
min(results$pcv_diff)
max(results$pcv_diff)
## results[,xbin:=cut(pnocv_diff,breaks=seq(-0.009,0.009,length.out=40),include.lowest=TRUE)]
results[,xbin:=cut(pnocv_diff,breaks=seq(-20,20,length.out=80),include.lowest=TRUE)]
results[,ybin:=cut(pcv_diff,breaks=seq(-2,2,length.out=20),include.lowest=TRUE)]
m=results[zthr > 0,.(n=.N),by=c("xbin","ybin","ncvscen","zthr","snpscen")]
m[,pc:=100*n/sum(n),by=c("ncvscen","zthr","snpscen")]
m[,x:=midpoints(xbin)]
m[,y:=midpoints(ybin)]
msub=m[ncvscen=="two CV" & snpscen=="1000  SNPs"]
ggplot(msub[zthr>0]) + #,col=pnocv_diff)) +
  ## geom_hex(binwidth=c(0.01/21,2/11),colour="grey") +
  geom_tile(aes(x=x,y=y,fill=pc),colour="grey") +
  geom_text(aes(x=x,y=y,label=round(pc),colour=100-pc),data=msub[zthr>0]) +
  ## geom_hline(yintercept=0) + geom_vline(xintercept=0,width=0.1) +
  ## geom_jitter(alpha=0.01) +
  ## geom_point(colour="transparent") +
  ## geom_density2d(alpha=0.2) +
  background_grid() +
  ## geom_rug(aes(x=pnull_diff,y=pcv_diff),alpha=0.02,data=results[zthr>0]) +
  scale_fill_viridis("percent") +
  scale_colour_viridis("percent",guide=FALSE) +
  facet_grid(ncvscen + zthr ~ snpscen) +
  theme(legend.position="bottom")

for(isnp in c(1000,2000,3000)) {
plots=lapply(c("single CV","two CV"), function(incv) {
  lapply(c(0.5,1,1.5), function(izthr) {
    msub=m[ncvscen==incv & snpscen==paste0(isnp,"  SNPs")]
    show_legend = izthr==0.5
    show_xlab = izthr==1.5
    show_ylab = incv=="single CV"
    p=ggplot(msub[zthr==izthr]) + #,col=pnocv_diff)) +
      geom_tile(aes(x=x,y=y,fill=pc),colour="grey") +
      geom_text(aes(x=x,y=y,label=round(pc,1),colour=100-pc)) +
      geom_point(aes(x=pnocv_diff,y=pcv_diff),colour="transparent",
                 data=results[ncvscen=="two CV" & snpscen=="1000  SNPs" & zthr==izthr]) +
      scale_fill_viridis("percent") +
      scale_colour_viridis("percent",guide=FALSE) +
      labs(x=if(show_xlab) { "Change in total PP at non-causal variants" } else { "" },
           y=if(show_ylab) { "Change in total PP at causal variants" } else { "" }) +
      ggtitle(paste0(incv,", |Z| threshold = ",izthr)) +
      theme_cowplot(font_size=10) +
      background_grid() +
      theme(legend.position="none",#if(show_legend) {c(1,1)} else { "none" },
            legend.justification=c("left","bottom"),
            legend.direction="horizontal",
            plot.title=element_text(face="bold"))
    ggMarginal(p,size=10)
  })
}) %>% do.call("c",.)
pgrid=plot_grid(plotlist=plots,ncol=2,align="hv",axis="xy",byrow=FALSE)
## save_plot(paste0("pipgrid-",isnp,".png"),plot=pgrid,base_height=10,base_width=8)
ggsave(paste0("figure-pipgrid-",isnp,".png"),plot=pgrid,height=10,width=8,units="in")
}



## univariate plots
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

results[,ptot_diff:=pnocv_diff+pcv_diff+pnull_diff]

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
,max_pos_diff=ggplot(results, aes(y=factor(zthr),x=max_pos_diff.noncv,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
   labs(x="Max PIP difference excluding CV")
,min_neg_diff=ggplot(results, aes(y=factor(zthr),x=min_neg_diff.noncv,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
   labs(x="Min PIP difference excluding CV")
,totdiff=ggplot(results, aes(y=factor(zthr),x=ptot_diff,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
   labs(x="Bias in PIP excluding CV")
,nocvdiff=ggplot(results, aes(y=factor(zthr),x=pnocv_diff,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
   labs(x="Bias in PIP excluding CV")
,cvdiff=ggplot(results, aes(y=factor(zthr),x=pcv_diff,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
   labs(x="Bias in PIP excluding CV")
,nulldiff=ggplot(results, aes(y=factor(zthr),x=pnull_diff,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
   labs(x="Bias in PIP excluding CV")
 ,pipcr=ggplot(results, aes(y=factor(zthr),x=pipcr,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=1,linetype="dashed") +
  ## geom_density_ridges(stat = "density", aes(height = stat(density)), col="grey") +
    labs(x="Correlation between PIP")
 ,pip=ggplot(m, aes(y=factor(zthr),x=pip_diff,fill=factor(zthr),col=factor(zthr))) +
  geom_vline(xintercept=0,linetype="dashed") +
    labs(x="Difference in PIP at causal variants"))


plotsplus = lapply(plots,function(p)
  p + geom_density_ridges(stat = "binline", scale=0.95,bins=21, draw_baseline = FALSE) +
  scale_fill_seaborn("|Z| threshold for trimming") +
  scale_colour_seaborn("|Z| threshold for trimming") +
  background_grid(major="y") +
  ylab("|Z| threshold") +
  theme(legend.position="bottom") +
  facet_grid(ncvscen ~ snpscen,space="free"))

msum=m[,.(mid=paste0(100*round(mean(abs(pip_diff)<0.1),2),"%")),
          by=c("cvscen","snpscen","zthr")]
plotsplus$pip =
  plotsplus$pip +
  geom_text(aes(label=mid,y=zthr*2+1.5),x=0.5,data=msum) +
  facet_grid(cvscen ~ snpscen, space="free") +
  scale_x_continuous(breaks=c(-1,0,1))

## stand alone
## plotsplus$time
plotsplus$pip

ggsave("~/Projects/coloc-susie/figure-pipdiff.png",plot=plotsplus$pip,height=6,width=6)

## other diagnostics
plotlist=plotsplus[c("pnull","pipcr")]
plotlist[1] %<>% lapply(., function(p) p + theme(legend.position="none"))
plot_grid(plotlist=plotlist[],ncol=1)

ggsave("~/Projects/coloc-susie/figure-otherdiagnostics.png",height=8,width=6)

plotlist=plotsplus[c("max_pos_diff","min_neg_diff","bias")]
plotlist=plotsplus[c("totdiff","cvdiff","nocvdiff","nulldiff")]
plot_grid(plotlist=plotlist[],ncol=1)

library(MASS)
library(ggplot2)
library(viridis)
## theme_set(theme_bw(base_size = 16))

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

results[zthr>0,dens:=get_density(x=pnull_diff, y=pcv_diff, n=100), by=c("ncvscen","zthr","snpscen")]
## results[zthr>0,dens:=dens/max(dens), by=c("ncvscen","zthr","snpscen")]
ggplot(results[zthr>0][order(dens)], aes(x=pnull_diff,y=pcv_diff,col=dens)) +
  ## geom_hex() +
  ## geom_hline(yintercept=0) + geom_vline(xintercept=0,width=0.1) +
  geom_point(alpha=1,size=5) +
  ## geom_density2d(alpha=0.2) +
  background_grid() +
  scale_colour_viridis() +
  facet_grid(ncvscen + zthr ~ snpscen, space="free")

## alternative time
## ggplot(results, aes(x=factor(nsnp),y=base_time + time_diff,col=factor(zthr))) +
##   geom_boxplot() +
##   facet_wrap(~ncv,labeller=label_both) +
##   background_grid(major="y")


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
