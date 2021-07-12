
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

tab2<-readPNG("data/coral_silhouettes/tabularG.png")
tab2<-rasterGrob(tab2, interpolate=TRUE)
mas2<-readPNG("data/coral_silhouettes/massiveG.png")
mas2<-rasterGrob(mas2, interpolate=TRUE)
cor2<-readPNG("data/coral_silhouettes/corymboseG.png")
cor2<-rasterGrob(cor2, interpolate=TRUE)
cor2<-readPNG("data/coral_silhouettes/corymboseG.png")
cor2<-rasterGrob(cor2, interpolate=TRUE)
dig2<-readPNG("data/coral_silhouettes/digitateG.png")
dig2<-rasterGrob(dig2, interpolate=TRUE)
brn2<-readPNG("data/coral_silhouettes/branchingG.png")
brn2<-rasterGrob(brn2, interpolate=TRUE)
#######################################


library("viridis")

plots <- list()

for (sp in spp) {
	#sp <- "Acy"
	y <- y.list[[sp]]
	
	ipm.p <- reshape(data.frame(t(ipm.p.list[[sp]])),
	direction="long", varying=list(c(1:100)), ids=y, times=y,
	new.row.names=c(1:10000))
	ipm.r <- reshape(data.frame(t(ipm.r.list[[sp]])),
	direction="long", varying=list(c(1:100)), ids=y, times=y,
	new.row.names=c(1:10000))
	
	sub<-gdat[gdat$spp==sp,]
	#sub<-subset(sub, area>=min(b) & area_next>=min(b))

plot<-ggplot()+
	geom_raster(data=ipm.p, aes(id, time,fill=X1))+
	geom_point(data=ipm.r[ipm.r$X1>0,], aes(id, time, col=X1, alpha=log(X1)), size=0.1)+
	ggtitle(paste(params$species[params$spp==sp]))+
	geom_point(data=sub, aes(x=area, y=area_next), col="grey", size=0.01)+
	scale_fill_viridis(n.breaks=3, trans="sqrt")+
	scale_colour_viridis(option="C", trans="log", breaks=c(10^3, 10^5, 10^7), labels=c(expression(~10^3),expression(~10^5), expression(~10^7)))+
	#scale_colour_viridis(option="C", trans="log", breaks=c(100, 1000,10000,round(max(ipm.r$X1),digits=-4)))+
	theme_bw()+
	guides(alpha="none", fill=guide_colourbar(order=1), colour=guide_colourbar(order=2))+
	geom_abline(col="black", slope=1, size=0.1, linetype="dashed")+
	scale_x_continuous(expand=c(0,0), limits=c(min(ipm.p$id), max(ipm.p$id)), labels=function(x){10^x}, breaks=c(-3,-2,-1,0,1))+
	scale_y_continuous(expand=c(0,0), limits=c(min(ipm.p$id), max(ipm.p$id)), labels=function(x){10^x}, breaks=c(-3,-2,-1,0,1) )+
	xlab(expression(area[~t]~(m^2)))+
	ylab(expression(area[~t~+1]~(m^2)))+
	theme(legend.title=element_blank(), legend.key.size=unit(2,"mm"), 
	legend.text=element_text(size=5, colour="white"), axis.title.x=element_text(size=8,color="white"), 
	legend.spacing.x = unit(0.05, "cm"),
	legend.margin = margin(0,0.1,0,0, unit="cm"),
	axis.title.y=element_blank(), 
	legend.key.height=unit(1.4,"mm"), legend.key.width=unit(1,"mm"),
	plot.title=element_text(size=7, face="bold.italic", hjust=0.5), legend.position=c(0.22,0.91), 
	legend.background=element_blank(), 
	legend.box = "horizontal", 
	plot.background=element_blank(), 
	axis.text.y=element_text(size=6, angle=90, hjust=0.5),
	axis.text.x=element_text(size=6))
	plot
	
	plots[[sp]] <- plot
	}

#fig.s2 <- plot_grid(plotlist=plots)
#fig.s2

fig.s2 <- plot_grid(

plot_grid(plots[["Ahy"]]+annotation_custom(tab,0,0.9,-1.5,-1), NULL, plots[["Acy"]]+annotation_custom(tab2,0,0.9,-1.5,-1)+
theme(axis.title.y=element_text(size=8, angle=90, hjust=1.2)), rel_heights=c(1,-0.12,1), ncol=1, align="v"),

plot_grid(plots[["Ain"]]+annotation_custom(brn,0,0.9,-1.7,-1), NULL, plots[["Aro"]]+annotation_custom(brn2,0,0.9,-1.7,-1), rel_heights=c(1,-0.12,1), ncol=1),

plot_grid(plots[["Ana"]]+annotation_custom(cor,0,0.9,-2,-1.5), NULL, plots[["Asp"]]+annotation_custom(cor2,0,0.9,-2,-1.5)+
theme(axis.title.x=element_text(colour="black",size=8)), rel_heights=c(1,-0.12,1), ncol=1),

#plot_grid(plots[["AN"]], NULL, plots[["AL"]], NULL, plots[["AM"]],rel_heights=c(1,-0.1,1, -0.1,1), ncol=1),

plot_grid(plots[["Adi"]]+annotation_custom(dig,0,0.9,-2.5,-2), NULL, plots[["Ahu"]]+annotation_custom(dig2,0,0.9,-2.5,-2), rel_heights=c(1,-0.12,1), ncol=1),

plot_grid(plots[["Gre"]]+annotation_custom(mas,0,0.9,-3,-2.2), NULL, plots[["Gpe"]]+annotation_custom(mas2,0,0.9,-3,-2.2), rel_heights=c(1,-0.12,1), ncol=1),

align="h",
nrow=1, rel_widths=c(1.1,1,1,1,1))

fig.s2

ggsave("figs/fig.s2.png", fig.s2, width=23, height=10.5, units="cm", dpi = 300)

