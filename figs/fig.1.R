
#######################################
# IMAGES
#######################################
tab<-readPNG("data/coral_silhouettes/tabularG.png")
tab<-rasterGrob(tab, interpolate=TRUE)
mas<-readPNG("data/coral_silhouettes/massiveG.png")
mas<-rasterGrob(mas, interpolate=TRUE)
cor<-readPNG("data/coral_silhouettes/corymboseG.png")
cor<-rasterGrob(cor, interpolate=TRUE)
cor<-readPNG("data/coral_silhouettes/corymboseG.png")
cor<-rasterGrob(cor, interpolate=TRUE)
dig<-readPNG("data/coral_silhouettes/digitateG.png")
dig<-rasterGrob(dig, interpolate=TRUE)
brn<-readPNG("data/coral_silhouettes/branchingG.png")
brn<-rasterGrob(brn, interpolate=TRUE)
#######################################

pcdat<-data.frame(spp=params$spp, x=pca$x[,1], y=pca$x[,2], abun=params$abundance_05)
rot<-data.frame(spp=rownames(pca$rotation), x=pca$rotation[,1], y=pca$rotation[,2])

rot$num <- c(2,4,6,7,3,5,1)
rot

group_name <- c(
expression(1*.~Fecundity[colony]),
expression(2*.~Growth[max]), 
expression(3*.~Partial~mortality[100]), 
expression(4*.~Fecundity[area]), 
#expression(5*.~Fecundity[intercept]),
expression(5*.~Mature~size[min]), 
expression(6*.~Survival[juvenile]), 
expression(7*.~Survival[adult]))

rot$x2 <- rot$x
rot$x2[rot$spp=="r.int"] <- rot$x[rot$spp=="r.int"]+0.03

pcdat$spp <- factor(pcdat$spp, levels=order)

pcplot <- ggplot()+
geom_segment(data=rot, aes(x=0,y=0, xend=x*1.3, yend=y*1.3), arrow=arrow(length=unit(1, "mm")))+
geom_text(data=NULL, aes(x=-3.47, y=seq(1.1,2.3,length.out=7)), label=rev(group_name),size=2, fontface="bold", hjust=0)+
geom_text_repel(data=rot, aes(x2*1.6, y*1.5,label=num),  size=2, fontface="bold", force=0.001)+
#geom_point(data=pcdat, aes(x,y, col=spp, size=abun))+
geom_point(data=pcdat, aes(x,y, fill=spp), size=4, shape=21, stroke=0.15)+
geom_text_repel(data=pcdat, aes(x+0.38,y-0, label=spp), size=2, force=0.0005)+
#scale_colour_manual(values=c(cols), labels=labs)+
scale_fill_manual(values=c(cols), labels=labs)+
scale_radius(range=c(5,10))+
xlim(c(min(pcdat$x)*1.1,max(pcdat$x)*1.2))+
ylim(c(min(pcdat$y)*1.2,max(pcdat$y)*1.1))+
labs(x=paste("PC1 (",exp[1],"%)", sep=""),y=paste("PC2 (",exp[2],"%)", sep=""))+
annotation_custom(tab, xmin=-3.5, xmax=-2, ymin=-1.7, ymax=-0.7)+
annotation_custom(brn, xmin=-3, xmax=-0.5, ymin=-0.5, ymax=0.3)+
annotation_custom(cor, xmin=-0.6, xmax=0.5, ymin=0.9, ymax=1.4)+
#annotation_custom(dig, xmin=2, xmax=3, ymin=-0.9, ymax=0.5)+
annotation_custom(dig, xmin=0.5, xmax=1.5, ymin=0, ymax=0.5)+
annotation_custom(mas, xmin=2.5, xmax=3.85, ymin=-1.1, ymax=-0.7)+
ggtitle("Demographic trade-offs")+
theme_bw()+theme(legend.title=element_blank(), legend.text=element_text(size=7, face="italic"), legend.key.size=unit(1, "mm"), 
#legend.position=c(0.18, 0.83), 
panel.grid.minor=element_line(color="grey98"),
panel.grid.major=element_line(color="grey98"),
legend.background=element_blank(), 
axis.text=element_text(size=6), axis.title=element_text(size=7),
plot.title=element_text(face="bold", size=8, hjust=0.5))
#dev.off()
pcplot


abunplot<- ggplot(params, aes(x=reorder(spp, -abundance_05), y=abundance_05, fill=spp))+
guides(fill="none")+
scale_y_continuous(expand=c(0,0), limits=c(0, 4500))+
geom_bar(stat="identity", col="black", size=0.1, width=0.75)+
scale_fill_manual(values=cols)+
theme_classic()+
geom_text(aes(y=abundance_05+450, label=spp), size=2)+
ggtitle("Abundance")+labs(y="N colonies")+
theme(plot.title=element_text(size=8, hjust=0.5, face="bold"), 
axis.title.x=element_text(size=8), 
#axis.ticks.y=element_blank(),
#axis.line.y=element_blank(),
axis.title.y=element_blank(),
axis.text.x=element_text(size=5), 
strip.text=element_blank(),
strip.background=element_blank(),
axis.text.y=element_blank()
)+coord_flip()+facet_wrap(~morphology, ncol=1, scales="free_y")
#abunplot


size.av <- size.av[order(size.av$area, decreasing=TRUE),]
size.av$spp <- factor(size.av$spp, levels=size.av$spp)
ss$spp <- factor(ss$spp, levels=size.av$spp)
dat$spp <- factor(dat$spp, levels=size.av$spp)


size <- ggplot()+
#geom_density(data=dat[!is.na(dat$spp),], aes(x=area_cm2/10000), fill=NA, col="grey")+
geom_density(data=ss[!is.na(ss$spp),], aes(x=10^area, fill=spp), col="black", size=0.15)+
geom_text(data=size.av, aes(x=0.001, y=0.75, label=spp),size=2)+
geom_segment(data=size.av, aes(x=10^area,xend=10^area, y=Inf, yend=-Inf), col="black", size=0.15)+
geom_segment(data=aggregate(area_cm2~spp, dat[!is.na(dat$spp),], mean), aes(x=area_cm2/10000,xend=area_cm2/10000, y=Inf, yend=-Inf), col="slategrey", size=0.15)+
facet_wrap(~spp, ncol=1, strip.position="left")+
scale_y_continuous(breaks=c(0.5))+
scale_fill_manual(values=cols)+
scale_colour_manual(values=cols)+
scale_x_log10()+
guides(fill="none", col="none")+
labs(x=expression(area~(m^2)))+
ggtitle("Size")+
theme_classic()+
theme(strip.background=element_blank(), 
strip.text=element_blank(), 
axis.title.y=element_blank(),
plot.background=element_blank(),
axis.line.y=element_blank(), 
axis.ticks.y=element_blank(),
#axis.text.y=element_text(size=5),
axis.title.x=element_text(size=8),
axis.text.x=element_text(size=5, angle=30, hjust=1),
axis.text.y=element_blank(),
plot.title=element_text(size=8, face="bold", hjust=0.5))
size


#limits=c(-100, max(sad$abundance)+100), expand=c(0,0)
sadplot <- ggplot(sad, aes(x=abundance))+geom_histogram(bins=20, col="white", fill="grey", size=0.3)+
scale_colour_manual(values=cols)+
scale_y_sqrt(expand=c(0,0), breaks=c(100))+
scale_x_continuous(breaks=c(0,2500,5000),)+
guides(col="none")+
labs(x="N colonies", y="N species")+
theme_classic()+
theme(panel.background=element_blank(), plot.background=element_blank(), 
#axis.line=element_blank(),
#axis.text.y=element_blank(),
axis.text.y=element_text(size=5, margin=margin(0,-5,0,0)),
axis.line.y=element_blank(),
axis.ticks.y=element_blank(),
axis.text.x=element_text(size=5), axis.title=element_text(size=8))
sadplot 


fig.1 <- plot_grid(size,
pcplot+guides(fill="none"), 
plot_grid(abunplot+theme(axis.title.x=element_blank())+scale_y_continuous(limits=c(0,max(sad$abundance+100)), breaks=c(0,2500,5000)), NULL, 
plot_grid(NULL, sadplot+theme(axis.title.y=element_blank()), rel_widths=c(-0.035, 1))
, 
NULL, get_legend(pcplot+guides(fill = guide_legend(override.aes = list(size=1.2)))), ncol=1, rel_heights=c(1, -0.05, 0.3,-0.05, 0.6)),
rel_widths=c(0.25,1, 0.3), labels=c("A", "B", "C"), label_size=8, hjust=c(0,-5, 0), nrow=1)
fig.1



fig.1 <- plot_grid(size,
pcplot+guides(fill="none"), 
plot_grid(abunplot, NULL, get_legend(pcplot+guides(fill = guide_legend(override.aes = list(size=1.2)))), ncol=1, rel_heights=c(1, -0.1,0.6)),
rel_widths=c(0.25,1, 0.3), labels=c("A", "B", "C"), label_size=8, hjust=c(0,-5, 0), nrow=1)
fig.1


