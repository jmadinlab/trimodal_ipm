
#source("R/data_prep.R")
#source("R/params.R")


#######################################
# TRADEOFFS 
#######################################
head(params)

# survival
params$p_mort<-inv.logit((params$p.slp*log10(0.01))+params$p.int)
params$av.surv<-aggregate(pred~spp, s.pred,mean)$pred
params$survcm<-aggregate(pred~spp, s.pred[s.pred$area<log10(pi*(11/100/2)^2),], FUN=mean)$pred
max.surv<-do.call(rbind, lapply(split(s.pred,as.factor(s.pred$spp)), function(x) {return(x[which.max(x$pred),])}))
params$safe.size<-max.surv$area

# growth
gdat$growth<-gdat$area_cm2_next-gdat$area_cm2	
params$av.growth<-aggregate(growth~spp, gdat[gdat$growth>0,], mean)$growth

# reproduction
params$f.cm2<-aggregate(f.cm2~spp, fec,mean)$f.cm2
params$en.cm2<-params$f.cm2*params$eggC
params$min.r<-aggregate(area_cm2~spp, fec[fec$reproductive==1,], min)$area_cm2

#params2<-params
#params2<-rbind(params2, rep(NA, ncol(params)))
#params2[12,"spp"]<-"CR"
#params2[12,"r.int"]<-0.2/1000
#params2[12,"min.r"]<-pi*((2/10/2)^2) # 1.2-2 mm diam


colnames(params)
pca<-prcomp(params[,c("r.int", "f.cm2","survcm","av.surv","p_mort")], scale=T, center=T)
#pca<-prcomp(params[,c("r.int", "f.cm2","av.surv", "p_mort")], scale=T, center=T)
pcdat<-data.frame(spp=params$spp, x=pca$x[,1], y=pca$x[,2], abun=params$abundance_05)
rot<-data.frame(spp=rownames(pca$rotation), x=pca$rotation[,1], y=pca$rotation[,2])
exp<-round(c(summary(pca)[[1]][1]^2/sum(summary(pca)[[1]]^2),summary(pca)[[1]][2]^2/sum(summary(pca)[[1]]^2)),3)*100


#pdf("figs/p7_tradeoffs.pdf")
ggplot()+
geom_segment(data=rot, aes(x=0,y=0, xend=x, yend=y), arrow=arrow(length=unit(1, "mm")))+
geom_text(data=rot, aes(x*1.2, y*1.2, label=spp))+
geom_point(data=pcdat, aes(x,y, col=spp, size=abun))+
geom_text(data=pcdat, aes(x,y, label=spp), size=3)+
scale_colour_manual(values=c(cols))+
guides(colour="none", size="none")+scale_radius(range=c(5,12))+
labs(x=paste("PC1 (",exp[1],"%)", sep=""),y=paste("PC2 (",exp[2],"%)", sep=""))+
theme_bw()
#dev.off()




# linear models..

lm.plot<-function(x, y, n){
	dat<-params
	dat$x<-params[,x]
	dat$y<-params[,y]
    ggplot(dat, aes(x,y))+	
    geom_smooth(method="lm", formula=y~poly(x,n))+
	geom_point(aes(fill=spp), shape=21, size=3.5)+
	geom_text(aes(label=spp), size=1.5)+
	labs(x=x, y=y)+
	scale_fill_manual(values=paste(cols))+guides(fill="none")+
	theme(axis.text=element_text(size=6),axis.title=element_text(size=10))
}

lm.plots<-plot_grid(
lm.plot("r.int","min.r",1)+scale_x_log10()+scale_y_log10(),
lm.plot("r.int","p_mort",1)+scale_x_log10(), 
lm.plot("r.int","safe.size",1)+scale_x_log10(),
lm.plot("r.int","survcm",2)+scale_x_log10(),
lm.plot("r.int","av.surv",2)+scale_x_log10(),
lm.plot("r.int","f.cm2",2)+scale_x_log10(),
lm.plot("r.int","en.cm2",1)+scale_x_log10(),
align="hv")
  

