
library("reshape2")

plotsE <- list()
for (sp in spp) {
	#sp <- "Ahy"
	y <- y.list[[sp]]
	sens <- reshape(data.frame(t(s.list[[sp]])),
	direction="long", varying=list(c(1:100)), ids=y, times=y,
	new.row.names=c(1:10000))
	e.k <- reshape(data.frame(t(eK.list[[sp]])),
	direction="long", varying=list(c(1:100)), ids=y, times=y,
	new.row.names=c(1:10000))
	e.p <- reshape(data.frame(t(eP.list[[sp]])),
	direction="long", varying=list(c(1:100)), ids=y, times=y,
	new.row.names=c(1:10000))
	e.r <- reshape(data.frame(t(eR.list[[sp]])),
	direction="long", varying=list(c(1:100)), ids=y, times=y,
	new.row.names=c(1:10000))
	
	sub<-gdat[gdat$spp==sp,]
	#sub<-subset(sub, area>=min(b) & area_next>=min(b))
	
	#limits=c(min(ipm.p$id), max(ipm.p$id)), 
	young <- log10(pi*((10/2)/100)^2)
	
plotE<-ggplot()+
	geom_raster(data=e.k, aes(id, time,fill=X1))+
		geom_point(data=sub, aes(x=area, y=area_next), col="grey", size=0.01)+
	scale_fill_viridis(n.breaks=3, trans="sqrt")+
	#scale_fill_viridis(n.breaks=3)+
	theme_bw()+
	guides(alpha="none", fill=guide_colourbar(order=1), colour=guide_colourbar(order=2))+
	geom_abline(col="black", slope=1, size=0.1, linetype="dashed")+
	geom_point(data=NULL, aes(x=young, y=young), shape=4, col="white")+
	scale_x_continuous(expand=c(0,0), labels=function(x){10^x}, breaks=c(-3,-2,-1,0,1))+
	scale_y_continuous(expand=c(0,0), labels=function(x){10^x}, breaks=c(-3,-2,-1,0,1) )+
	xlab(expression(area[~t]~(m^2)))+
	ylab(expression(area[~t~+1]~(m^2)))+
	ggtitle(paste(sp))+
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
	
	plotsE[[sp]] <- plotE
	}


plot_grid(plotlist=plotsE)

spp

eplot <- plot_grid(

plot_grid(plotsE[["Ahy"]]+annotation_custom(tab,0,0.9,-1.5,-1), NULL, plotsE[["Acy"]]+annotation_custom(tab2,0,0.9,-1.5,-1)+
theme(axis.title.y=element_text(size=8, angle=90, hjust=1.2)), rel_heights=c(1,-0.12,1), ncol=1, align="v"),

plot_grid(plotsE[["Ain"]]+annotation_custom(brn,0,0.9,-1.7,-1), NULL, plotsE[["Aro"]]+annotation_custom(brn2,0,0.9,-1.7,-1), rel_heights=c(1,-0.12,1), ncol=1),

plot_grid(plotsE[["Ami"]]+annotation_custom(cor,0,0.9,-2,-1.5), NULL, plotsE[["Asp"]]+annotation_custom(cor2,0,0.9,-2,-1.5)+
theme(axis.title.x=element_text(colour="black",size=8)), rel_heights=c(1,-0.12,1), ncol=1),

#plot_grid(plots[["AN"]], NULL, plots[["AL"]], NULL, plots[["AM"]],rel_heights=c(1,-0.1,1, -0.1,1), ncol=1),

plot_grid(plotsE[["Adi"]]+annotation_custom(dig,0,0.9,-2.5,-2), NULL, plotsE[["Ahu"]]+annotation_custom(dig2,0,0.9,-2.5,-2), rel_heights=c(1,-0.12,1), ncol=1),

plot_grid(plotsE[["Gre"]]+annotation_custom(mas,0,0.9,-3,-2.2), NULL, plotsE[["Gpe"]]+annotation_custom(mas2,0,0.9,-3,-2.2), rel_heights=c(1,-0.12,1), ncol=1),

align="h",
nrow=1, rel_widths=c(1.1,1,1,1,1))

eplot



