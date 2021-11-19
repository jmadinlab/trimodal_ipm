
rectheme <- theme(axis.text=element_text(size=7), axis.title=element_text(size=8), axis.line=element_line(size=0.25), axis.ticks=element_line(size=0.25), plot.title=element_text(size=8, face="bold", hjust=0.5),  panel.background=element_blank())

comp$maxgen <- maxdat$gen[match(comp$morph, maxdat$morph)]

projplot <- ggplot()+
geom_line(data=proj2, aes(gen, diff, col=morph), linetype="dotted", size=0.2)+
geom_line(data=proj2[proj2$gen <= proj2$maxgen,], aes(gen, diff, col=morph), size=0.75)+
scale_y_log10()+
#scale_x_log10()+
#geom_segment(data=proj4[proj4$gen==proj4$max,], aes(xend=gen, x=gen,y=Inf, yend=-Inf), linetype="dotted", col="grey")+
#geom_text(data=maxdat, aes(x, y, label=lab),  col="grey", angle=90, size=2, hjust=0)+
#geom_segment(aes(x=-Inf, xend=Inf, y=real, yend=real, col=morph), linetype="dotted")+
#geom_point(data=comp, inherit.aes=F, aes(difftime, diff, col=morph))+
#geom_point(aes(max,real, col=morph))+
#geom_segment(data=dmaxdat, aes(y=diff, yend=diff, x=norm, xend=gen), arrow=arrow(length=unit(1,"mm")))+
#geom_segment(data=NULL, aes(y=5.769612, yend=5.769612, x=2, xend=1), arrow=arrow(length=unit(1,"mm")))+
geom_point(data=proj2[proj2$gen==proj2$maxgen,], aes(gen, diff, fill=morph), shape=21, stroke=0.1, size=2)+
scale_colour_manual(values=colsC)+
scale_fill_manual(values=colsC)+
theme_classic()+
theme(axis.line=element_line(size=0.1), axis.title=element_text(size=8), 
plot.title=element_text(size=8, hjust=0.5, face="bold"))+
coord_cartesian(ylim=c(1,55), xlim=c(1,85))+
guides(colour="none", fill="none")+
#scale_x_continuous(expand=c(0,0))+
labs(x="Time (years)", y=expression(N[" common "]/N[" rare"]))+
ggtitle("Time taken to project \n observed gaps in abundance")
projplot


#######################################
# RECRUITMENT PLOT
#######################################
tab<-readPNG("data/coral_silhouettes/tabularG.png")
tab<-rasterGrob(tab, interpolate=TRUE)
mas<-readPNG("data/coral_silhouettes/massiveG.png")
mas<-rasterGrob(mas, interpolate=TRUE)
cor<-readPNG("data/coral_silhouettes/corymboseG.png")
cor<-rasterGrob(cor, interpolate=TRUE)
cor2<-readPNG("data/coral_silhouettes/corymboseG.png")
cor2<-rasterGrob(cor2, interpolate=TRUE)
dig<-readPNG("data/coral_silhouettes/digitateG.png")
dig<-rasterGrob(dig, interpolate=TRUE)
brn<-readPNG("data/coral_silhouettes/branchingGt.png")
brn<-rasterGrob(brn, interpolate=TRUE)


comp$labels <- paste(comp$Common, "-", comp$Rare, sep="")

mas_x <- comp$doubdiff[comp$morph=="massive"]
tab_x <- comp$doubdiff[comp$morph=="tabular"]
brn_x <- comp$doubdiff[comp$morph=="staghorn"]
dig_x <- comp$doubdiff[comp$morph=="digitate"]
cor_x <- comp$doubdiff[comp$morph=="corymbose"]
cor2_x <- comp$doubdiff[comp$morph=="corymbose_2"]

comp$rank<-NA
comp$rank[order(comp$lamdiff)] <- 1:nrow(comp)
mas_y <- comp$rank[comp$morph=="massive"]
tab_y <- comp$rank[comp$morph=="tabular"]
brn_y <- comp$rank[comp$morph=="staghorn"]
dig_y <- comp$rank[comp$morph=="digitate"]
cor_y <- comp$rank[comp$morph=="corymbose"]
cor2_y <- comp$rank[comp$morph=="corymbose_2"]

lamcomp2 <- ggplot()+
geom_bar(data=comp, aes(x=lamdiff, y=reorder(labels, lamdiff), fill=morph), stat="Identity",  col="black", size=0.1, width=0.75)+
labs(x=expression(Difference~"in"~lambda))+
geom_segment(data=comp, aes(y=reorder(labels, -lamdiff), yend=reorder(labels, -lamdiff), x=lamdiff, xend=doubdiff), arrow=arrow(length=unit(1,"mm")))+
guides(fill="none")+scale_fill_manual(values=colsC)+
scale_x_continuous(expand=c(0,0), limits=c(0,1.1))+
ggtitle("Difference in fitness \n(common-rare)")+
geom_segment(data=NULL, aes(y=2, yend=2, x=0.8, xend=0.9), arrow=arrow(length=unit(1,"mm")))+
geom_text(data=NULL,  aes(x=0.85, y=1.2, label="recruitment x2 in \nall taxa"), size=2)+
annotation_custom(mas, mas_x+0.01, mas_x+0.13, mas_y-0.5, mas_y+0.5)+
annotation_custom(brn, brn_x+0.01, brn_x+0.14, brn_y-0.5, brn_y+0.5)+
annotation_custom(tab, tab_x+0.01, tab_x+0.18, tab_y-0.7, tab_y+0.7)+
annotation_custom(cor, cor_x+0.01, cor_x+0.15, cor_y-0.5, cor_y+0.5)+
annotation_custom(cor2, cor2_x+0.01, cor2_x+0.15, cor2_y-0.5, cor2_y+0.5)+
annotation_custom(dig, dig_x+0.01, dig_x+0.15, dig_y-0.6, dig_y+0.6)+
theme_classic()+theme(axis.title.y=element_blank(), axis.title.x=element_text(size=7), 
plot.title=element_text(size=8, hjust=0.5, face="bold"),
axis.text=element_text(size=7))
lamcomp2








lsize <- 0.5


#######################################
# CONTOURS
#######################################

#scapes <- pairs3[!pairs3$morph=="corymbose",]
scapes <- pairs3
head(pairs3)
scapes$posi <- ifelse(scapes$X1 > 0, "common > rare", "rare > common")
#scapes$morph <- ifelse(scapes$morph=="corymbose_2", "corymbose", scapes$morph)
scapes$morph <- ifelse(scapes$morph=="staghorn", "arborescent", scapes$morph)
#scapes$morph <- ifelse(scapes$morph=="massive", "boulder", scapes$morph)
scapes$morph <- factor(scapes$morph, levels=rev(c("tabular", "arborescent", "corymbose","corymbose_2","digitate", "massive")))

colsC2 <- colsC
names(colsC2)[6]<-"massive"
names(colsC2)[3]<-"arborescent"

sp.points <- params2[,c("morph","rec.size", "spp")]
sp.points$rec <- params$rec[match(sp.points$morph, params2$morph)]

#sp.points <- sp.points[!sp.points$morph=="corymbose",]
#sp.points$rec <- ifelse(sp.points$morph=="massive", 10^-4, 10^-3)
#sp.points$morph <- ifelse(sp.points$morph=="corymbose_2", "corymbose", sp.points$morph)
sp.points$morph <- ifelse(sp.points$morph=="staghorn", "arborescent", sp.points$morph)
#sp.points$morph <- ifelse(sp.points$morph=="massive", "boulder", sp.points$morph)
sp.points$morph <- factor(sp.points$morph, levels=rev(c("tabular", "arborescent", "corymbose","corymbose_2","digitate", "massive")))


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

contours <- ggplot(data=scapes, aes(x=rec, y=sqrt(((10^recsize)*10000)), z=value))+
geom_contour_filled(aes(fill=morph, alpha=..level..),breaks=c(0,0.05,0.1, 0.2, 0.4, 0.8), col="grey", size=lsize)+
geom_contour(aes(col=morph), breaks=c(0,0.05,0.1, 0.2, 0.4, 0.8,1.6),size=lsize)+
geom_text_contour(breaks = c(0,0.05, 0.1, 0.2, 0.4), size=1.5, skip=0, stroke=0.15)+
#geom_segment(data=sp.points, inherit.aes=F, aes(x=5*10^-4, xend=rec, y=5, yend=sqrt(((10^rec.size)*10000))), col="black", arrow=arrow(length=unit(1, "mm")))+
#geom_point(data=sp.points, inherit.aes=F, aes(x=5*10^-4,y=5), shape=21, col="black", size=0.5)+
geom_point(data=sp.points, inherit.aes=F, aes(x=rec,y=sqrt(((10^rec.size)*10000))), shape=3, stroke=1, size=0.25, col="black")+
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
	panel.background=element_rect(fill="white"),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	plot.title=element_text(size=8, hjust=0.5, vjust=-7, face="bold"),
	axis.text.x=element_text(size=5), 
	axis.text.y=element_text(size=5, angle=90, hjust=0.5), 
	axis.title=element_text(size=7))#+coord_flip()
contours




fig3<-plot_grid(
plot_grid(lamcomp2, projplot, nrow=1, rel_widths=c(0.7,1), labels=c("A","B"), label_size=9), 
NULL,
contours,
ncol=1, rel_heights=c(1,-0.0,0.9),labels=c("","","C"), label_size=9)
fig3



