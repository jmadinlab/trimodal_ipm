
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

ipm.k.list <- list()
ipm.p.list <- list()
ipm.r.list <- list()
y.list <- list()

for (sp in spp) {
	rec <- 1
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	#sub<-gdat[gdat$spp==sp,]
    mod <- bigmatrix()
	y.list[[sp]] <- y 
	ipm.k.list[[sp]] <- mod$P
	ipm.p.list[[sp]] <- mod$P
	ipm.r.list[[sp]] <- mod$R
	}

library("viridis")

plots <- list()
for (sp in spp) {
	#sp <- "Acy"
	ydat <- data.frame(n=c(1:n), y=y.list[[sp]])
	ipm.p <- melt(ipm.p.list[[sp]])
	ipm.p$y1 <- ydat$y[match(ipm.p$Var1, ydat$n)]
	ipm.p$y2 <- ydat$y[match(ipm.p$Var2, ydat$n)]
	ipm.r <- melt(ipm.r.list[[sp]])
	ipm.r$y1 <- ydat$y[match(ipm.r$Var1, ydat$n)]
	ipm.r$y2 <- ydat$y[match(ipm.r$Var2, ydat$n)]
	#sub <- gdat[gdat$spp==sp,]
plot<-ggplot()+
	geom_raster(data=ipm.p, aes(y2, y1,fill=value))+
	geom_point(data=ipm.r[ipm.r$value>0,], aes(y2, y1, col=value, alpha=log(value)), size=0.1)+
	ggtitle(paste(params$species[params$spp==sp]))+
	#geom_point(data=sub, aes(x=area, y=area_next), col="grey", size=0.01)+
	scale_fill_gradient(low="white", high="grey",n.breaks=3)+
	scale_colour_viridis(trans="log", breaks=c(10^3, 10^5, 10^7), labels=c(expression(~10^3),expression(~10^5), expression(~10^7)))+
	theme_bw()+
	guides(alpha="none", fill=guide_colourbar(order=1), colour=guide_colourbar(order=2))+
	geom_abline(col="black", slope=1, size=0.1, linetype="dashed")+
	scale_x_continuous(expand=c(0,0), limits=c(min(ipm.p$y2), max(ipm.p$y2)), labels=function(x){10^x}, breaks=c(-3,-2,-1,0,1))+
	scale_y_continuous(expand=c(0,0), limits=c(min(ipm.p$y1), max(ipm.p$y1)), labels=function(x){10^x}, breaks=c(-3,-2,-1,0,1) )+
	xlab(expression(area[~t]~(m^2)))+
	ylab(expression(area[~t~+1]~(m^2)))+
	theme(legend.title=element_blank(), legend.key.size=unit(2,"mm"), 
	legend.text=element_text(size=5), axis.title.x=element_text(size=8,color="white"), 
	legend.spacing.x = unit(0.05, "cm"),
	legend.margin = margin(0,0.1,0,0, unit="cm"),
	axis.title.y=element_blank(), 
	legend.key.height=unit(1.4,"mm"), legend.key.width=unit(1,"mm"),
	plot.title=element_text(size=7, face="bold.italic", hjust=0.5), legend.position=c(0.22,0.91), 
	legend.background=element_blank(), 
	legend.box = "horizontal", 
	plot.background=element_blank(), 
	axis.text.y=element_text(size=4.5, angle=90, hjust=0.5),
	axis.text.x=element_text(size=4.5))

	plots[[sp]] <- plot
	}


fig.s5 <- plot_grid(
plot_grid(plots[["Ahy"]]+annotation_custom(tab,-0.5,0.2,-2,-1.5), NULL, plots[["Acy"]]+annotation_custom(tab,-0.5,0.2,-2,-1.3)+
theme(axis.title.y=element_text(size=8, angle=90, hjust=1.2)), rel_heights=c(1,-0.12,1), ncol=1, align="v"),
plot_grid(plots[["Ain"]]+annotation_custom(brn,-0.5,0.5,-2,-1.3), NULL, plots[["Aro"]]+annotation_custom(brn,-0.5,0.5,-2,-1.3), rel_heights=c(1,-0.12,1), ncol=1),
plot_grid(plots[["Asp"]]+annotation_custom(cor,-1,0.5,-2.5,-2), NULL, plots[["Ami"]]+annotation_custom(cor,-1,0.5,-2.3,-1.8)+
theme(axis.title.x=element_text(colour="black",size=8)), rel_heights=c(1,-0.12,1), ncol=1),
plot_grid(plots[["Adi"]]+annotation_custom(dig,-1,0.3,-2.8,-2.2), NULL, plots[["Ahu"]]+annotation_custom(dig,-1,0.3,-2.8,-2.2), rel_heights=c(1,-0.12,1), ncol=1),
plot_grid(plots[["Gre"]]+annotation_custom(mas,-0.8,0,-3,-2.2), NULL, plots[["Gpe"]]+annotation_custom(mas,-0.8,0,-3,-2.2), rel_heights=c(1,-0.12,1), ncol=1),
align="h",
nrow=1, rel_widths=c(1.1,1,1,1,1))
fig.s5



