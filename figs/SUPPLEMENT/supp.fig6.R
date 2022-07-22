
lsize <- 0.5


#######################################
# CONTOURS
#######################################

scapes <- pairs3
head(pairs3)
scapes$morph <- factor(scapes$morph, levels=rev(c("tabular", "staghorn", "corymbose","corymbose_2","digitate", "massive")))

sp.points <- params2[,c("morph","rec.size", "spp")]
sp.points$rec <- params$rec[match(sp.points$morph, params2$morph)]
sp.points$morph <- factor(sp.points$morph, levels=rev(c("tabular", "staghorn", "corymbose","corymbose_2","digitate", "massive")))

blueish <- colsC[2]

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
sp.points$rec.const <- params$rec.const[match(sp.points$spp, params$spp)]

library("metR")

sp.points<-aggregate(.~morph, subset(sp.points, select=-spp), mean)

fig.s6 <- ggplot(data=scapes, aes(x=rec, y=sqrt(((10^recsize)*10000)), z=value))+
geom_contour_filled(aes(fill=morph, alpha=..level..),breaks=c(0,0.05,0.1, 0.2, 0.4, 0.8), col="grey", size=lsize)+
geom_contour(aes(col=morph), breaks=c(0,0.05,0.1, 0.2, 0.4, 0.8,1.6),size=lsize)+
geom_text_contour(breaks = c(0,0.05, 0.1, 0.2, 0.4), size=1.5, skip=0, stroke=0.15)+
geom_point(data=sp.points, inherit.aes=F, aes(x=rec,y=sqrt(((10^rec.size)*10000))), shape=3, stroke=1, size=0.25, col="black")+
facet_wrap(~morph, nrow=2)+
labs(x="Probability of larval settlement", y="Recruit diameter")+
scale_fill_manual(values=colsC)+
scale_colour_manual(values=colsC)+
guides(fill="none", alpha="none", col="none")+
ggtitle(expression(bold(Delta*fitness~landscapes)))+
	scale_y_log10(expand=c(0,0))+
	#scale_y_continuous(expand=c(0,0))+
	scale_x_log10(breaks=c(10^-5, 10^-4,10^-3), labels=c(
		expression(10^-5),
		expression(10^-4),
		expression(10^-3)), expand=c(0,0))+
	theme_bw()+
	theme(strip.background=element_blank(),  
	strip.text=element_text(size=7, vjust=-2),
	panel.background=element_rect(fill="white"),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	plot.title=element_text(size=8, hjust=0.5, vjust=-7, face="bold"),
	axis.text.x=element_text(size=5), 
	axis.text.y=element_text(size=5, angle=90, hjust=0.5), 
	axis.title=element_text(size=7))#+coord_flip()
fig.s6

