
rectheme <- theme(axis.text=element_text(size=7), axis.title=element_text(size=8), axis.line=element_line(size=0.25), axis.ticks=element_line(size=0.25), plot.title=element_text(size=8, face="bold", hjust=0.5), 
plot.background=element_blank(), panel.background=element_blank())


#######################################
# RECRUITMENT PLOT
#######################################
tab<-readPNG("data/coral_silhouettes/tabularC.png")
tab<-rasterGrob(tab, interpolate=TRUE)
mas<-readPNG("data/coral_silhouettes/massiveC.png")
mas<-rasterGrob(mas, interpolate=TRUE)
cor<-readPNG("data/coral_silhouettes/corymboseC.png")
cor<-rasterGrob(cor, interpolate=TRUE)
cor2<-readPNG("data/coral_silhouettes/corymboseC2.png")
cor2<-rasterGrob(cor2, interpolate=TRUE)
dig<-readPNG("data/coral_silhouettes/digitateC.png")
dig<-rasterGrob(dig, interpolate=TRUE)
brn<-readPNG("data/coral_silhouettes/branchingC.png")
brn<-rasterGrob(brn, interpolate=TRUE)

lsize <- 0.5

p.plot <- ggplot(tiles, aes(Family, p.rec))+
geom_segment(data=NULL, inherit.aes=F, aes(x="Acroporidae",xend=-Inf, y=rec.fam[1,2], yend=rec.fam[1,2]), col="grey", linetype="dotted")+
geom_segment(data=NULL, inherit.aes=F, aes(x="Merulinidae",xend=-Inf, y=rec.fam[2,2], yend=rec.fam[2,2]), col="grey", linetype="dotted")+
geom_boxplot(size=0.2, width=0.7, fill="grey90", outlier.size=0)+
geom_text(data=aggregate(p.rec~Family, tiles, max), aes(Family, p.rec*6.5,label=Family), size=1.75, fontface="italic")+
scale_y_log10(limits=c(min(pairs$rec), 0.1),breaks=c(10^-6, 10^-4, 10^-2))+
coord_flip()+
theme_classic()+ 
ggtitle("Recruitment success")+
rectheme+theme(axis.line=element_blank(), axis.ticks=element_blank(),axis.title=element_blank(), axis.text=element_blank())
p.plot

 xmin <- log10(0.01)
 xmax <- log10(0.05)
 xmin2 <- log10(0.013)
 xmax2 <- log10(0.037)
 ymin <- pairs[pairs$rec==0.01,c("morph", "logdiff")]
 
 
 ylims <- c(min(c(pairs$logdiff, pairs2$logdiff))*1.1, max(c(pairs$logdiff, pairs2$logdiff))*1.1)
 
pairplot <- ggplot(pairs,aes(rec, logdiff, col=morph))+
geom_segment(data=NULL, inherit.aes=F, aes(y=Inf, yend=-Inf, x=rec.fam[1,2], xend=rec.fam[1,2]), col="grey", linetype="dotted")+
geom_segment(data=NULL, inherit.aes=F, aes(y=Inf,yend=-Inf, x=rec.fam[2,2], xend=rec.fam[2,2]), col="grey", linetype="dotted")+
		geom_line(size=lsize )+
		ylim(ylims)+
		geom_point(data=comp,aes(rec, logdiff))+
		scale_x_log10(breaks=c(10^-6, 10^-4,10^-2), labels=c(
		expression(10^-6),
		expression(10^-4),
		expression(10^-2)))+
		coord_cartesian(xlim=c(min(pairs$rec), 0.1))+
		geom_hline(yintercept=0,  size=0.25)+
 annotation_custom(brn, xmin, xmax, ymin[5,2]-0.05, ymin[5,2]+0.15)+
 annotation_custom(tab, xmin, xmax, ymin[6,2]+0.05, ymin[6,2]+0.25)+
 annotation_custom(mas, xmin2, xmax2, ymin[4,2]-0.1, ymin[4,2]+0.1)+
 annotation_custom(cor, xmin2, xmax2, ymin[2,2], ymin[2,2]+0.2)+
 annotation_custom(dig, xmin2, xmax2, ymin[3,2]-0.2, ymin[3,2])+
 annotation_custom(cor2, xmin2, xmax2, ymin[1,2]-0.1, ymin[1,2]+0.1)+
		labs(x="Probability of larval settlement", y=expression(log*(lambda)[common]~-~log*(lambda)[rare]))+
scale_colour_manual(values=colsC)+guides(col="none")+
theme_classic()+rectheme+theme(plot.title=element_blank())
pairplot


#######################################
# REC SIZE PLOT
#######################################

store2$rsize_cm <- sqrt(((10^store2$rec.size)*10000)/pi)*2
pairs2$rsize_cm <- sqrt(((10^pairs2$rec.size)*10000)/pi)*2		
minrec <- min(store2$rsize_cm)
maxrec <- max(store2$rsize_cm)
rszizlim <- c(minrec, 12.5)

s.rec <- rsize
s.rec$s.rec <- s.rec$recsize_1yr

sizebar <- ggplot(rsize, aes(Family, recsize_1yr))+
geom_segment(data=NULL, inherit.aes=F, aes(x="Acroporidae",xend=-Inf, y=rsize.fam[1,2], yend=rsize.fam[1,2]), col="grey", linetype="dotted")+
geom_segment(data=NULL, inherit.aes=F, aes(x="Merulinidae",xend=-Inf, y=rsize.fam[2,2], yend=rsize.fam[2,2]), col="grey", linetype="dotted")+geom_boxplot(size=0.2, width=0.7, outlier.size=0.1, fill="grey90")+
coord_flip(ylim=rszizlim)+
geom_text(data=aggregate(recsize_1yr~Family, rsize, max), aes(Family, recsize_1yr+1.5,label=Family), size=1.75, fontface="italic")+
		guides(fill="none")+
		ggtitle("Recruit size")+
		theme_classic()+rectheme+theme(axis.line= element_blank(),axis.ticks=element_blank(),axis.title=element_blank(), axis.text=element_blank())
		sizebar
						
  xmin <- 11.9
  xmax <-12.6
  xmin2 <- 11.7
  xmax2 <- 12.9
  
 ymin <- pairs2[pairs2$rec.size==max(pairs2$rec.size),c("morph", "logdiff")]

sizediff <-	ggplot(pairs2,aes(rsize_cm, logdiff, col=morph))+
geom_segment(data=NULL, inherit.aes=F, aes(y=Inf,yend=-Inf, x=rsize.fam[1,2], xend=rsize.fam[1,2]), col="grey", linetype="dotted")+
geom_segment(data=NULL, inherit.aes=F, aes(y=Inf,yend=-Inf, x=rsize.fam[2,2], xend=rsize.fam[2,2]), col="grey", linetype="dotted")+
#geom_boxplot(size=0.2, width=0.7, outlier.size=0.1, fill="grey90")+
		geom_line(size=lsize )+
		geom_point(data=comp,aes(rec.size.cm, logdiff))+
		#coord_cartesian(xlim=rszizlim)+
		scale_colour_manual(values=colsC)+
		scale_x_continuous(limits=rszizlim)+
		ylim(ylims)+
		guides(col="none")+
		geom_hline(yintercept=0,  size=0.25)+
		ggtitle("common vs. rare")+
		labs(x="Diameter at 1 year (cm)")+
 annotation_custom(brn, xmin, xmax, ymin[5,2]-0.05, ymin[5,2]+0.01)+
 annotation_custom(mas, xmin, xmax, ymin[4,2]-0.06, ymin[4,2]+0.04)+
 annotation_custom(cor, xmin, xmax, ymin[2,2]-0.045, ymin[2,2]+0.06)+
 annotation_custom(dig, xmin, xmax, ymin[3,2]-0.01, ymin[3,2]+0.09)+
 annotation_custom(cor2, xmin, xmax, ymin[1,2]-0.04, ymin[1,2]+0.06)+
  annotation_custom(tab, xmin2, xmax2, ymin[6,2]-0.05, ymin[6,2]+0.1)+
theme_classic()+rectheme+theme(plot.title=element_blank(), axis.title.y=element_blank())
sizediff 
		
#######################################
# CONTOURS
#######################################

scapes <- pairs3[!pairs3$morph=="corymbose",]
scapes$posi <- ifelse(scapes$X1 > 0, "common > rare", "rare > common")
scapes$morph <- ifelse(scapes$morph=="corymbose_2", "corymbose", scapes$morph)
scapes$morph <- ifelse(scapes$morph=="staghorn", "arborescent", scapes$morph)
scapes$morph <- ifelse(scapes$morph=="massive", "boulder", scapes$morph)
scapes$morph <- factor(scapes$morph, levels=rev(c("tabular", "arborescent", "corymbose","digitate", "boulder")))

colsC2 <- colsC
names(colsC2)[6]<-"boulder"
names(colsC2)[3]<-"arborescent"

sp.points <- params2[,c("morph","rec.size", "spp")]
sp.points <- sp.points[!sp.points$morph=="corymbose",]
sp.points$rec <- ifelse(sp.points$morph=="massive", 10^-4, 10^-3)
sp.points$morph <- ifelse(sp.points$morph=="corymbose_2", "corymbose", sp.points$morph)
sp.points$morph <- ifelse(sp.points$morph=="staghorn", "arborescent", sp.points$morph)
sp.points$morph <- ifelse(sp.points$morph=="massive", "boulder", sp.points$morph)
sp.points$morph <- factor(sp.points$morph, levels=rev(c("tabular", "arborescent", "corymbose","digitate", "boulder")))

blueish <- colsC[6]

legend <- ggplot()+
geom_point(data=scapes, aes(x=id, y=sqrt(((10^time)*10000)), fill=posi), shape=22)+
theme_classic()+
scale_fill_manual(values=c(blueish, "white"))+
guides(col="none",  fill=guide_legend(direction="vertical", title=expression(Delta~log*(lambda)), title.hjust=0.5, title.vjust=-1, label.hjust=100))+
theme(legend.text=element_text(size=6),
	#legend.title=element_text(size=6, face="bold"),
	legend.title=element_blank(),
	plot.margin=unit(c(1, 4, 1, 1), "cm"),
	legend.key.height=unit(1,"mm"),
	legend.position=c(0.9,0 ),
	legend.background=element_blank(),
	legend.key.width=unit(5,"mm"))
	
	sp.points$eggC <- params$eggC[match(sp.points$spp, params$spp)]
	sp.points$recEn <- sp.points$rec * (10^(scale(sp.points$eggC)/max(abs(scale(sp.points$eggC))))) 



contours <- ggplot()+
geom_contour_filled(data=scapes, aes(x=id, y=sqrt(((10^time)*10000)), z=X1, fill=morph, alpha=..level..),breaks=c(0,0.1, 0.2, 0.4, 0.8, 1.6), col="grey", size=lsize)+
geom_contour(data=scapes, aes(x=id, y=sqrt(((10^time)*10000)), z=X1, col=morph), breaks=c(0,0.1, 0.2, 0.4, 0.8,1.6),  size=lsize)+
geom_segment(data=sp.points, aes(x=rec, xend=recEn, y=4.5, yend=sqrt(((10^rec.size)*10000))), col="black", arrow=arrow(length=unit(1, "mm")))+
geom_point(data=sp.points, inherit.aes=F, aes(x=rec,y=4.5), shape=21, col="black", size=0.5)+
#geom_point(data=sp.points, inherit.aes=F, aes(x=rec,y=sqrt(((10^rec.size)*10000))), shape=3, stroke=1, size=0.25, col="black")+
facet_wrap(~morph, nrow=1)+
labs(x="Probability of larval settlement", y="Recruit diameter")+
scale_fill_manual(values=colsC2)+
scale_colour_manual(values=colsC2)+
guides(fill="none", alpha="none", col="none")+
ggtitle(expression(bold(Delta*fitness~landscapes)))+
	#scale_y_log10(expand=c(0,0))+
	scale_y_continuous(expand=c(0,0))+
	scale_x_log10(breaks=c(10^-5, 10^-4,10^-3), labels=c(
		expression(10^-5),
		expression(10^-4),
		expression(10^-3)), expand=c(0,0))+
	theme_bw()+
	theme(strip.background=element_blank(),  
	strip.text=element_text(size=7, vjust=-2),
	plot.background=element_blank(),
	panel.background=element_rect(fill="white"),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	plot.title=element_text(size=8, hjust=0.5, vjust=-7, face="bold"),
	axis.text.x=element_text(size=5), 
	axis.text.y=element_text(size=5, angle=90, hjust=0.5), 
	axis.title=element_text(size=7))#+coord_flip()
contours

#######################################
# FIG 2
#######################################

fig2 <- plot_grid(
plot_grid(
plot_grid(p.plot, ggplot()+theme_void(),pairplot, ncol=1, align="v", rel_heights=c(0.25,-0.08,1)), 
plot_grid(sizebar,ggplot()+theme_void(), sizediff, ncol=1, align="v", rel_heights=c(0.25,-0.08,1)), labels=c("A","B"), label_size=9),
contours, 
NULL,
get_legend(legend),
NULL,
ncol=1, rel_heights=c(1,0.75, -0.08,0.03,0.05), labels=c("","C"), label_size=9)
fig2



