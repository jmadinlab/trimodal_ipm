






# size at maximum survival 
max.surv<-do.call(rbind, lapply(split(s.pred, as.factor(s.pred$spp)), function(x) {return(x[which.max(x$pred),])}))
params$safe.size2<-max.surv$area
rec.cm <- 11
params$survcm2<-aggregate(pred~spp, s.pred[s.pred$area<log10(pi*(rec.cm/100/2)^2),], FUN=mean)$pred

# recruit survival
min.S <- aggregate(area~morph, aggregate(area~spp+morph, s.pred, min), max)
s.pred$min.s <- min.S$area[match(s.pred$morph, min.S$morph)]
max.S <- aggregate(area~morph, aggregate(area~spp+morph, s.pred, max), min)
s.pred$max.s <- max.S$area[match(s.pred$morph, min.S$morph)]
quarts <- aggregate(area~morph, s.pred, FUN=function(x){quantile(x, 0.2 )})
s.pred$quarts <- quarts$area[match(s.pred$morph, quarts$morph)]
s.pred$less <-ifelse(s.pred$area < s.pred$quarts & s.pred$area>s.pred$min.s & s.pred$area < s.pred$max.s, "y", "n")
head(s.pred)
ggplot()+
geom_line(data=s.pred, aes(area, pred, group=spp), size=0.1, col="grey")+
geom_line(data=s.pred[s.pred$less=="y",], aes(area, pred, col=spp))+scale_colour_manual(values=cols)+guides(col="none")
rec.surv<- do.call(rbind, lapply(split(s.pred, as.factor(s.pred$spp)), function(x) {
	x[x$less=="y",]
	}))
head(rec.surv)
avrec<-aggregate(pred~spp, rec.surv, mean)
params$rec.surv <- avrec$pred[match(params$spp, avrec$spp)]

# minimum at reproductive maturity
params$min.r2<-aggregate(area_cm2~spp, fec[fec$reproductive==1,], min)$area_cm2

# size 50% maturity (how delayed)
mat50 <- aggregate(area~spp, m.pred[m.pred$pred>0.5,], min)
params$mat50<-mat50$area[match(params$spp, mat50$spp)]

# barplots
bars <- function(t){
	params$x <- params[,t]
ggplot(params, aes(x=reorder(spp, -x), y=x, fill=spp))+guides(fill="none")+
geom_bar(stat="identity")+
ylab(t)+ggtitle(t)+
theme(axis.title.x=element_blank())+
facet_wrap(~morphology, scales="free_x", nrow=1)+
scale_fill_manual(values=cols)}

plot_grid(
bars("eggC"),
bars("f.int")+coord_cartesian(ylim=c(15,16.2)),
bars("f.slp")+coord_cartesian(ylim=c(2.3,2.8)),
ncol=1)

# lms
lm.plot<-function(x, y, n){
	dat<-params
	dat$x<-params[,x]
	dat$y<-params[,y]
    ggplot(dat, aes(x,y))+	
    geom_smooth(method="lm", formula=y~poly(x,n), col="black", size=0.3)+
    geom_path(aes(group=morphology), linetype="dotted")+
	geom_point(aes(fill=spp), shape=21, size=3.5, stroke=0.1)+
	geom_text(aes(label=spp), size=1.5)+
	labs(x=x, y=y)+
	scale_fill_manual(values=cols)+guides(fill="none")+	theme(axis.text=element_text(size=6),axis.title=element_text(size=10), plot.title=element_text(size=7))+theme_classic()
}

plot_grid(
lm.plot("mat50","f.slp", 1)+
ggtitle("delayed maturity = delayed eggs")+
ylab("higher slope = more delayed"),
lm.plot("f.slp", "rec.surv", 2)+scale_y_log10()+
ggtitle("delayed eggs = juvenile survival"),
lm.plot("eggC", "f.slp", 1),
lm.plot("f.slp", "f.cm2", 1),
lm.plot("r.int", "f.slp", 1)
)


# pca
head(params)
params$safe.sizex<-10^params$safe.size
pca2<-prcomp(params[,c("rsize.gr", "f.slp", "eggC", "f.int", "f.colony")], scale=T, center=T)
#biplot(pca2)
params$Rpca1 <- pca2$x[,1][match(params$spp, rownames(pca2$x))]
params$Rpca2 <- pca2$x[,2][match(params$spp, rownames(pca2$x))]
Rvecs <- data.frame(pca2$rotation, lab=rownames(pca2$rotation))



ggplot()+
geom_path(data=params, aes(Rpca1, Rpca2, group=morphology))+
geom_point(data=params, aes(Rpca1, Rpca2, fill=spp, size=abundance_05), shape=21)+
geom_segment(data=Rvecs, aes(0,0, xend=PC1, yend=PC2))+
geom_text(data=Rvecs, aes(PC1*1.2, PC2*1.2, label=lab))+
scale_fill_manual(values=cols)+guides(col="none", size="none", fill="none")




gdat$growth <- 10^gdat$area_next - 10^gdat$area
params$av.growth <- aggregate(growth~spp, gdat, mean)$growth
params$av.pmort <- aggregate(p_mort~spp, gdat, mean)$p_mort
gr.plots<-plot_grid(
lm.plot("r.int","min.r",1)+scale_x_log10()+scale_y_log10()+labs(x="Max colony growth rate", y="Min size at maturity"),
lm.plot("r.int","p_mort",1)+scale_x_log10()+scale_y_log10()+labs(x="Max colony growth rate", y="Partial mortality rate"),
lm.plot("r.int","f.colony",1)+scale_x_log10()+scale_y_log10()+scale_y_log10()+labs(x="Max colony growth rate", y="Colony fecundity"),
lm.plot("r.int","survcm",2)+scale_x_log10()+scale_y_log10()+labs(x="Max colony growth rate", y="Recruit survival"),
lm.plot("r.int","av.surv",2)+scale_x_log10()+scale_y_log10()+labs(x="Max colony growth rate", y="Adult survival"),
lm.plot("r.int","f.cm2",2)+scale_x_log10()+scale_y_log10()+labs(x="Max colony growth rate", y="Fecundity per unit area"),
align="hv")
gr.plots






