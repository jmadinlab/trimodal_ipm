#######################################
# MORE REC FUNCTIONS FOR PLOTTING
#######################################

rec.xDET <- 10^(seq(-5,-2, speed))

storeDET<-data.frame()
for (sp in spp) {
	for (rec in rec.xDET){
	  rec.size <- params$rec.size[params$spp==sp]
	  h <- mesh()$h
	  y <- mesh()$y
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  storeDET<-rbind(storeDET, data.frame(sp=sp, rec=rec, lam=bigmatrix()$lam))
   } 
  }

# model lam based on tiles estimates
storeTIL<-data.frame()
for (sp in spp) {
		#sp <- "Ahy"
		fam <- params$family[params$spp==sp]
	  rec.size <- params$rec.size[params$spp==sp]
	  rec.lit <- tiles$p.rec[tiles$Family==fam]
	  q25 <- quantile(rec.lit, 0.25)
	  q75 <- quantile(rec.lit, 0.75)
	  rec.lit<-seq(from=q25, to=q75, length.out=10)
	  for (rec in rec.lit){
	  h <- mesh()$h
	  y <- mesh()$y
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  storeTIL<-rbind(storeTIL, data.frame(sp=sp, fam=fam, rec=rec, lam=bigmatrix()$lam))
   } 
  }

    storeTIL$genus<-ifelse(storeTIL$fam=="Acroporidae","Acropora", "Goniastrea")
      storeDET$genus<-ifelse(storeDET$sp %in% c("Gre","Gpe"),"Goniastrea", "Acropora")
  storeLIM <- subset(storeDET, rec > 1*10^-5 & lam <2.2)
  

  storeCI<- aggregate(lam~sp+fam, storeTIL, mean)
  storeCI$lam<-params$lam.est[match(storeCI$sp, params$spp)]
  storeCI$sd<-aggregate(lam~sp, storeTIL, sd)$lam
  storeCI$n<-aggregate(lam~sp, storeTIL, length)$lam
  storeCI$se<- storeCI$sd/sqrt(storeCI$n)
  storeCI$upper <- storeCI$lam+(1.96*storeCI$se)
  storeCI$lower <- storeCI$lam-(1.96*storeCI$se)

rectheme <- theme(axis.text=element_text(size=7), axis.title=element_text(size=8), axis.line=element_line(size=0.25), axis.ticks=element_line(size=0.25), plot.title=element_text(size=8, face="bold", hjust=0.5), panel.background=element_blank())

lsize <- 0.5

tiles$genus<-ifelse(tiles$Family=="Acroporidae","Acropora", "Goniastrea")

p.plot <- ggplot(tiles, aes(p.rec, genus))+
geom_boxplot(size=0.2, width=0.7, fill="grey90", outlier.size=0)+
geom_text(data=aggregate(p.rec~genus, tiles[tiles$Family=="Acroporidae",], min), aes(y=genus, x=p.rec/2.1,label=genus), size=1.9, fontface="italic", col="grey")+
geom_text(data=aggregate(p.rec~genus, tiles[tiles$Family=="Merulinidae",], min), aes(y=genus, x=p.rec*12,label=genus), size=1.9, fontface="italic", col="grey")+
#scale_y_log10(limits=c(min(pairs$rec), 0.1),breaks=c(10^-6, 10^-4, 10^-2))+
scale_x_log10(limits=c(10^-7, 10^-2))+
#xlim(c(10^-5, 10^-2))+
coord_cartesian( xlim=c(10^-5, 10^-2))+
theme_classic()+ 
ggtitle("Recruitment vs. fitness")+
rectheme+theme(axis.line=element_blank(), axis.ticks=element_blank(),axis.title=element_blank(), axis.text=element_blank())
p.plot



  segs <- data.frame(meds=c(median(tiles$p.rec[tiles$Family==
 "Acroporidae"]),median(tiles$p.rec[tiles$Family==
 "Merulinidae"])), genus=c("Acropora", "Goniastrea"), max=c(2.2,1.4))
 
   recplot <- ggplot(storeTIL, aes(x=rec, y=lam, col=sp))+
geom_segment(data=segs , inherit.aes=F, aes(y=max, yend=-Inf, x=meds, xend=meds), col="grey", linetype="dotted")+
geom_hline(yintercept=1, size=0.2)+
geom_line(data=storeLIM, size=0.5)+
geom_line(size=1.25)+
scale_x_log10(breaks=c(10^-5, 10^-4,10^-3,10^-2), labels=c(
		expression(10^-5),
		expression(10^-4),
		expression(10^-3),
		expression(10^-2)))+
scale_colour_manual(values=cols)+
facet_grid(genus~., space="free_y", scales="free_y")+
 labs(x=expression(Probability~of~recruitment~(P[rec])), y=expression(Fitness~(lambda)))+
coord_cartesian(xlim=c(10^-5, 10^-2))+
 theme_classic()+rectheme+theme(strip.background=element_blank(), strip.text=element_blank())+guides(colour="none")
 recplot
  
  
  estplot2<-ggplot(storeCI, aes(x=reorder(sp, -lam), y=lam, fill=sp))+
geom_hline(yintercept=1)+
geom_segment(data=storeCI, aes(x=reorder(sp, -lam), xend=reorder(sp, -lam), y=lower, yend=upper), size=0.25)+
#geom_boxplot(size=0.1,outlier.size=0)+
geom_point(shape=21, size=2, stroke=0.1)+
#stat_summary(geom="point", fun="mean", shape=21, size=2.5, stroke=0.2)+
ggtitle("Fitness estimates")+
scale_fill_manual(values=cols)+
 labs(y=expression(lambda))+
 theme_classic()+guides(colour="none", fill="none")+rectheme+
 theme(axis.text.x=element_blank(), axis.title.x=element_blank())
 estplot2

change<-ggplot(diffs, aes(reorder(spp, -lam), diff, fill=spp))+
stat_summary(geom="bar", fun="mean", col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
#geom_bar(stat="identity", position="dodge")+
#facet_wrap(~type, ncol=1)+
#geom_segment(data=diff[diff$type=="cover",], aes(x=spp, xend=spp, y=0, yend=diff, col=spp), arrow=arrow(length=unit(1,"mm")))+
#labs(y="Change in abundance", x="species", title="bars=individuals", subtitle="arrows=cover")+
#stat_summary(data=diffs, aes(reorder(spp, -lam), cover2, fill=spp), geom="point", fun="mean", size=1, shape=21)+
scale_fill_manual(values=cols)+guides(fill="none", col="none")+
#ggtitle("bars=individuals, arrows=cover")+
geom_hline(yintercept=0)+rectheme+
scale_colour_manual(values=cols)+
theme_classic()+
ggtitle("Observed change in \nabundance (2011-14)")+
ylab("Relative change in N")+
rectheme+theme(axis.text.x=element_text(angle=90, vjust=0.5, size=7), axis.title.x=element_blank(), plot.title=element_text(size=8, hjust=0.5, face="bold"))
change

examples2 <- 
plot_grid(NULL, plot_grid(NULL, 
plots[["Ahy"]]+annotation_custom(tab,-0.5, 0.2,-2,-1.4)+theme(plot.title=element_blank())+guides(col="none",fill="none"),
NULL,
plots[["Gre"]]+annotation_custom(mas,-0.7, 0.2,-3.5,-2.3)+theme(plot.title=element_blank())+guides(col="none",fill="none"),
ncol=1, rel_heights=c(0.1,1,-0.15,1))+
draw_label(expression(size[" "*t]~(m^2)), 0.7, 0.05, size=8, hjust=1)+
draw_label(expression(size[" "*t+1]~(m^2)), -0.05, 0.55, size=8, angle=90)+
draw_label("IPM matrices (examples)", 0.55, 0.97, fontface="bold", size=8),
nrow=1, rel_widths=c(0.12,1))
examples2


get_legend(plots[["Gre"]])


 recplotX<-plot_grid(p.plot, ggplot()+theme_void(),recplot, ncol=1, align="v", rel_heights=c(0.25,-0.05,1), axis="lr")

LEGENDS <- plot_grid(get_legend(plots[["Ahy"]]+theme(legend.text=element_text(colour="black"),legend.box = "vertical", legend.direction = "vertical",legend.position="left")+guides(fill=guide_colourbar(frame.colour = "black", frame.linewidth = 0.1, ticks.colour="black",ticks.linewidth = 1), colour=guide_colourbar(frame.colour = "black", frame.linewidth = 0.1, ticks.colour="black",ticks.linewidth = 1)))
,NULL,
get_legend(plots[["Gre"]]+theme(legend.text=element_text(colour="black"), legend.direction = "vertical",legend.box = "vertical",legend.position="left")+guides(fill=guide_colourbar(frame.colour = "black", frame.linewidth = 0.1, ticks.colour="black",ticks.linewidth = 1), colour=guide_colourbar(frame.colour = "black", frame.linewidth = 0.1, ticks.colour="black",ticks.linewidth = 1))),
NULL,
rel_heights=c(1,-0.1,1,0.2),ncol=1)
#LEGENDS

fig2<-plot_grid(examples2,NULL,LEGENDS, recplotX, plot_grid(estplot2, change, ncol=1, align="v", rel_heights=c(0.5,1),labels=c("C","D"), label_size=9), nrow=1, rel_widths=c(1,-0.04,0.2, 1.5,1.1), labels=c("A","","","B",""), label_size=9)+
draw_label("Acropora spp.", 0.4, 0.78, size=8, fontface="italic", hjust=0)+
draw_label("Goniastrea spp.", 0.4, 0.3, size=8, fontface="italic", hjust=0)
fig2


