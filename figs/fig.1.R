
#######################################
# IMAGES
#######################################
tab<-readPNG("figs/coral_silhouettes/tabularG.png")
tab<-rasterGrob(tab, interpolate=TRUE)
mas<-readPNG("figs/coral_silhouettes/massiveG.png")
mas<-rasterGrob(mas, interpolate=TRUE)
cor<-readPNG("figs/coral_silhouettes/corymboseG.png")
cor<-rasterGrob(cor, interpolate=TRUE)
cor<-readPNG("figs/coral_silhouettes/corymboseG.png")
cor<-rasterGrob(cor, interpolate=TRUE)
dig<-readPNG("figs/coral_silhouettes/digitateG.png")
dig<-rasterGrob(dig, interpolate=TRUE)
brn<-readPNG("figs/coral_silhouettes/branchingG.png")
brn<-rasterGrob(brn, interpolate=TRUE)

#######################################
# SIZE
#######################################

ss$spp <- params$spp[match(ss$species, params$species)]
ord.size <- params[order(params$size.ss, decreasing=TRUE),"spp"]
params$spp <- factor(params$spp, levels=ord.size)
ss$spp <- factor(ss$spp, levels=ord.size)

size <- ggplot()+
#geom_density(data=dat[!is.na(dat$spp),], aes(x=area_cm2/10000), fill=NA, col="grey")+
geom_density(data=ss[!is.na(ss$spp),], aes(x=10^area, fill=spp), col="black", size=0.15)+
geom_text(data=params, aes(x=0.001, y=0.75, label=spp),size=2)+
geom_segment(data=params, aes(x=10^size.ss,xend=10^size.ss, y=Inf, yend=-Inf), col="black", size=0.15)+
geom_segment(data=params, aes(x=10^size.dat,xend=10^size.dat, y=Inf, yend=-Inf), col="slategrey", size=0.15)+
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

params$spp <- factor(params$spp, levels=order)
ss$spp <- factor(ss$spp, levels=order)

#######################################
# PCA
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
expression(5*.~Mature~size[min]), 
expression(6*.~Survival[juvenile]), 
expression(7*.~Survival[adult]))

rot$x2 <- rot$x
rot$x2[rot$spp=="r.int"] <- rot$x[rot$spp=="r.int"]+0.03

pcdat$spp <- factor(pcdat$spp, levels=order)

labs2 <- paste(labs, " (", order, ")", sep="")

pcplot <- ggplot()+
geom_segment(data=rot, aes(x=0,y=0, xend=x*1.3, yend=y*1.3), arrow=arrow(length=unit(1, "mm")))+
geom_text(data=NULL, aes(x=-3.47, y=seq(1.1,2.3,length.out=7)), label=rev(group_name),size=2, fontface="bold", hjust=0)+
geom_text_repel(data=rot, aes(x2*1.6, y*1.5,label=num),  size=2, fontface="bold", force=0.001)+
geom_point(data=pcdat, aes(x,y, fill=spp), size=4, shape=21, stroke=0.15)+
geom_text_repel(data=pcdat, aes(x+0.38,y-0, label=spp), size=2, force=0.0005)+
#scale_colour_manual(values=c(cols), labels=labs)+
scale_fill_manual(values=cols, labels=labs2)+
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
guides(fill = guide_legend(override.aes = list(size=1.2)))+
ggtitle("Demographic trade-offs")+
theme_bw()+theme(legend.title=element_blank(), legend.text=element_text(size=7, face="italic"), legend.key.size=unit(0.1, "mm"), 
panel.grid.minor=element_line(color="grey98"),
panel.grid.major=element_line(color="grey98"),
legend.background=element_blank(), 
axis.text=element_text(size=6), axis.title=element_text(size=7),
plot.title=element_text(face="bold", size=8, hjust=0.5))
pcplot

#######################################
# ABUNDANCE
#######################################

abunplot<-plot_grid(NULL,
plot_grid(
ggplot(data=abun.BT, aes(reorder(spp, -abundance_05), N/10, fill=spp))+
#geom_bar(data=tri2.av, aes(reorder(species, -abun05), N/10, fill=species),stat="identity")+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
facet_grid(.~morphology, scales="free_x", space="free_x")+
guides(fill="none")+
labs(subtitle="2005",y="Colonies per m")+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
ggtitle("Abundance")+
theme(axis.title=element_blank(), 
axis.text.x=element_blank(),
axis.text.y=element_text(size=7),
plot.subtitle=element_text(size=7, hjust=0.5),
panel.spacing.x=unit(1,"mm"),
plot.title=element_text(size=8, hjust=0.5, face="bold"), 
strip.text=element_blank(), strip.background=element_blank())
,
ggplot()+theme_void(),
ggplot(data=abun.LIT[abun.LIT$year==2011,], aes(reorder(spp, -abundance_05), N/10, fill=spp))+
#geom_bar(data=tri2.av, aes(reorder(species, -abun05), N/10, fill=species),stat="identity")+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
facet_grid(.~morphology, scales="free_x", space="free_x")+
guides(fill="none")+
labs(subtitle="2011",y="Colonies per m")+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
theme(axis.title=element_blank(), 
axis.text.x=element_blank(),
axis.text.y=element_text(size=7),
panel.spacing.x=unit(1,"mm"),
plot.subtitle=element_text(size=7, hjust=0.5),
plot.title=element_text(size=8, hjust=0.5, face="bold"), 
strip.text=element_blank(), strip.background=element_blank())
,
ggplot()+theme_void(),
ggplot(data=abun.LIT[abun.LIT$year==2014,], aes(reorder(spp, -abundance_05), N/10, fill=spp))+
#geom_bar(data=tri2.av, aes(reorder(species, -abun05), N/10, fill=species),stat="identity")+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
facet_grid(.~morphology, scales="free_x", space="free_x")+
guides(fill="none")+
labs(subtitle="2014", y="Colonies per m")+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
theme(axis.title=element_blank(), 
axis.text.x=element_text(size=6, angle=90, vjust=0.5),
axis.text.y=element_text(size=7),
panel.spacing.x=unit(1,"mm"),
plot.subtitle=element_text(size=7, hjust=0.5),
plot.title=element_text(size=8, hjust=0.5, face="bold"), 
strip.text=element_blank(), strip.background=element_blank()),
ncol=1, rel_heights=c(1.1, -0.05,1,-0.05,1.1),
align="hv", axis="lr"),
rel_widths=c(0.03,1))
abunplot



fig.1 <- plot_grid(size,
plot_grid(pcplot+guides(fill="none"),NULL,
get_legend(pcplot+guides(fill = guide_legend(ncol=2, override.aes = list(size=1.5)))), ncol=1, rel_heights=c(1,-0.02, 0.25)), abunplot,
rel_widths=c(0.28,1, 0.47), labels=c("a", "b", "c"), label_size=8, hjust=c(0,-5, 0), nrow=1)+
draw_label(expression(Colonies~per~m^2), 0.737, 0.78, size=6, angle=90)+
draw_label("Colonies per m", 0.737, 0.5, size=6, angle=90)+
draw_label("Colonies per m", 0.737, 0.18, size=6, angle=90)
fig.1
