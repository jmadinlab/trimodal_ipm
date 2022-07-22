
source("figs/SUPPLEMENT/supp.fig5.R")

examples <- 
plot_grid(NULL,plot_grid(NULL,
plots[["Ahy"]]+ggtitle("Ahy")+theme(legend.position="right", plot.title=element_blank(),legend.box="vertical",legend.spacing.y=unit(0.5,"mm"), legend.margin=margin(0,0,0,-7)),
plots[["Ahu"]]+ggtitle("Adi")+theme(legend.position="right", plot.title=element_blank(),legend.box="vertical", legend.spacing.y=unit(0.5,"mm"), legend.margin=margin(0,0,0,-7)),
nrow=1, rel_widths=c(0.2,1,1))+
draw_label(expression(size[" "*t]~(m^2)), 0.6, 0.14, size=8, hjust=1)+
draw_label(expression(size[" "*t+1]~(m^2)), 0.06, 0.65, size=8, angle=90)+
draw_label(expression(bold(Example~IPMs~(before~bolditalic(P[rec])))), 0.55, 1.05, fontface="bold", size=8), nrow=2, rel_heights=c(0.15,1))
examples

params$lam.est3 <- aggregate(lam~spp, boot, median)$lam
params$log.lam3 <- log(params$lam.est)
arrows2 <- dcast(params[!params$spp=="Ana", ], morphology~abundance_pair,  value.var="log.lam3")
arrows2$Common2 <- arrows2$Common-0.03
arrows2$Rare2 <- arrows2$Rare+0.02
params$AB <- ifelse(params$abundance_pair=="Rare","R","C")
params$AB[params$spp=="Ana"]<-"D"


boot$morphology <- params$morphology[match(boot$spp, params$spp)]
morph.ord <- c("massive","digitate","corymbose","staghorn","tabular")
boot$morphology <- factor(boot$morphology, levels=morph.ord)
arrows2$morphology <- factor(arrows2$morphology, levels=morph.ord)
params$morphology <- factor(params$morphology, levels=morph.ord)

arrows2$yaxis <- c(38, 30,15,13,13)
params$yaxis <- arrows2$yaxis[match(params$morphology, arrows2$morphology)]
arrows2

arrows2$diff <- arrows2$Common - arrows2$Rare

lab.point <- data.frame(morphology="massive", x=c(0.2, 0.2, 0.2), y=c(30, 20, 10), lab = c('C = Common','R = Rare','D = Declined'))
lab.point$morphology <- factor(lab.point$morphology, levels=morph.ord)

params$log.lam3[params$spp=="Ana"]<-
params$log.lam3[params$spp=="Ana"]+0.02
params$log.lam3[params$spp=="Ami"]<-
params$log.lam3[params$spp=="Ami"]-0.02


lam.dens <- plot_grid(NULL, ggplot(boot)+
geom_segment(data=arrows2, aes(x=-Inf, xend=Inf, y=yaxis*1.2, yend=yaxis*1.2), col="white")+
geom_segment(data=arrows2, aes(y=yaxis, yend=yaxis, x=Rare2, xend=Common2), arrow=arrow(length=unit(1, "mm")))+
geom_text(data=lab.point, aes(x,y, label=lab), size=1.8, hjust=0, nudge_y=0.3)+
geom_density(aes(x=log(lam), col=spp, fill=spp), alpha=0.5)+
geom_vline(xintercept=0, linetype="dotted")+
geom_point(data=params, aes(log.lam3, yaxis,col=spp), shape=21, size=3, stroke=0.7, fill="white")+
geom_point(data=params, aes(log.lam3, yaxis, col=spp), size=3, alpha=0.5)+
geom_text(data=params[!params$spp=="Ana",], aes(log.lam3, yaxis, label=AB), size=1.8, fontface="bold")+
geom_text(data=params[params$spp=="Ana",], aes(log.lam3, yaxis, label="D"), size=1.8, fontface="bold")+
#geom_hline(yintercept=0, size=1, col="white")+
geom_hline(yintercept=-0.1, size=1, col="black")+
#scale_y_sqrt(expand=c(0,0))+
scale_y_continuous(expand=c(0,0))+
#lims(x=c(min(boot$lam), max(boot$lam)))+
ggtitle(expression(bold(Population~growth)))+
xlim(c(min(log(boot$lam)), max(log(boot$lam))-0.02))+
coord_cartesian(xlim=c(-0.5, max(log(boot$lam))))+
#scale_x_log10(limits=c(min(boot$lam), max(boot$lam)))+
facet_wrap(~morphology, ncol=1, scales="free_y", strip.position="left")+
labs(y="Density", x=expression(log(lambda)))+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
theme_classic()+theme(strip.text=element_blank(), strip.background=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=8), axis.title.x=element_text(size=9),plot.title=element_text(size=8, hjust=0.5, face="bold")), nrow=1, rel_widths=c(0.1,1))
lam.dens

#r.long$label <-ifelse(r.long$variable=="rec.tran", "Transects","Constant")
#r.long$label <- factor(r.long$label, levels=c("Fitted","Constant"))
params$genus <- ifelse(params$family=="Acroporidae", "Acr","Gon")


arrows <- dcast(params[!params$spp=="Ana", ], morphology~abundance_pair,  value.var="rec.one")
arrows$Common2 <- 10^(log10(arrows$Common)+0.2)
arrows$Rare2 <- arrows$Rare
arrows[arrows$morphology %in% c("staghorn", "tabular"),c("Rare2","Common2")] <- arrows[arrows$morphology %in% c("staghorn", "tabular"),c("Rare","Common")]

arrows$diff <-   arrows$Rare/arrows$Common


storeDET2 <- storeDET[storeDET$rec>1*10^-5 & storeDET$rec<2*10^-1,]

rechange <- ggplot()+geom_line(data=storeDET2, aes(rec, log(lam), col=spp))+
scale_x_log10(limits=c(min(storeDET2$rec),max(storeDET2$rec)), breaks=c(10^-5, 10^-4,10^-3,10^-2, 10^-1), labels=c(
		expression(10^-5), expression(10^-4), expression(10^-3),expression(10^-2), expression(10^-1)))+
geom_hline(yintercept=0, size=0.1)+
geom_point(data=params, aes(x=rec, y=log(lam.est)), size=1, col="white")+
geom_point(data=params, aes(x=rec, y=log(lam.est), col=spp, fill=spp), size=1, alpha=0.5)+
labs(x=expression(Recruitment~probablity~(P[rec])))+
geom_segment(data=params[params$family=="Acroporidae",], aes(x=rec, xend=rec,y=-0.52, yend=-0.51), arrow=arrow(length=unit(1, "mm") ))+
geom_segment(data=params[params$family=="Merulinidae",], aes(x=rec, xend=rec,y=-0.52, yend=-0.51), arrow=arrow(length=unit(1, "mm") ))+
geom_text(data=NULL, aes(x=params$rec[1], y=-0.46), label="Acr", size=2.5, col="grey")+
geom_text(data=NULL, aes(x=params$rec[11], y=-0.46), label="Gon", size=2.5,col="grey")+
guides(col="none", fill="none")+
theme_classic()+
labs(x=expression(Recruitment~probablity~(P[rec])), y=expression(log(lambda)))+
coord_cartesian(ylim=c(-0.55,0.65))+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
theme( plot.background=element_blank(), axis.text=element_text(size=8), axis.title=element_text(size=9), plot.title=element_text(size=8, hjust=0.5, face="bold"), panel.background=element_blank())
rechange 

arrows$Common2[arrows$morphology=="tabular"] <- 9*10^-05
arrows$Common2[arrows$morphology=="staghorn"] <- 5.5*10^-04 

growplot <- ggplot()+
geom_segment(data=arrows, aes(y=morphology, yend=morphology, x=1, xend=Common), col="grey95", size=1)+
geom_segment(data=arrows, aes(y=morphology, yend=morphology, x=Rare2, xend=Common2), arrow=arrow(length=unit(1, "mm")))+
#geom_point(aes(fill=spp), shape=21, size=3, stroke=0.1)+
geom_point(data=params, aes(rec.one, morphology,col=spp), shape=21, size=3, stroke=0.7, fill="white")+
geom_point(data=params, aes(rec.one, morphology, col=spp), size=3, alpha=0.5)+
scale_x_log10(limits=c(min(storeDET2$rec),1), breaks=c(10^-5, 10^-4,10^-3,10^-2, 10^-1), labels=c(
		expression(10^-5), expression(10^-4), expression(10^-3),expression(10^-2), expression(10^-1)))+
		coord_cartesian(xlim=c(min(storeDET2$rec),max(storeDET2$rec)))+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
theme_classic()+
ggtitle("Recruitment vs population growth")+
geom_text(data=params[!params$spp=="Ana",], aes(rec.one, morphology, label=AB), size=1.8, fontface="bold")+
geom_text(data=params[params$spp=="Ana",], aes(rec.one, morphology, label="D"), size=1.8, fontface="bold")+
scale_y_discrete(expand=c(0.3,0.3))+
labs(x=expression(Recruitment~probablity~(P[rec])))+ 
theme(plot.subtitle=element_text( size=7, hjust=1, vjust=-5), axis.title.y=element_blank(), axis.title=element_blank(),axis.text=element_blank(), plot.background=element_blank(), plot.title=element_text(size=8, hjust=0.5, face="bold"), panel.background=element_blank(), axis.line=element_blank(),axis.ticks=element_blank())
growplot



#######################################
#GENERATION TIME
#######################################	

sds<-aggregate(lam~spp, boot, function(x){ quantile(x, 0.95) })
params$lam.sd1 <- (sds$lam[match(params$spp, sds$spp)])
sds<-aggregate(lam~spp, boot, function(x){ quantile(x, 0.05) })
params$lam.sd2 <- (sds$lam[match(params$spp, sds$spp)])
sds<-aggregate(lam~spp, boot, median)
params$lam.mn <- (sds$lam[match(params$spp, sds$spp)])

sds<-aggregate(GT~spp, boot, function(x){ quantile(x, 0.95) })
params$GT.sd1 <- (sds$GT[match(params$spp, sds$spp)])
sds<-aggregate(GT~spp, boot, function(x){ quantile(x, 0.05) })
params$GT.sd2 <- (sds$GT[match(params$spp, sds$spp)])
sds<-aggregate(GT~spp, boot, median)
params$GT.mn <- (sds$GT[match(params$spp, sds$spp)])

ggplot()+geom_hline(yintercept=1, size=0.1)+
geom_segment(data=params, aes(x=GT, xend=GT, y=lam.sd1, yend=lam.sd2), size=0.2)+
geom_point(data=params,aes(GT, lam.mn, fill=spp), shape=21, stroke=0.2, size=3)+
scale_fill_manual(values=cols)+guides(fill="none")+
theme_classic()


GTplot<- ggplot(boot, aes(GT, log(lam)))+
geom_hline(yintercept=0, size=0.1)+
geom_segment(data=params, aes(x=GT, xend=GT, y=log(lam.sd1), yend=log(lam.sd2), col=spp), size=0.5)+
geom_segment(data=params, aes(x=GT.sd1, xend=GT.sd2, y=log(lam.mn), yend=log(lam.mn), col=spp), size=0.5)+
geom_point(data=params, aes(GT, log(lam.mn),col=spp), shape=21, size=3, stroke=0.7, fill="white")+
geom_text(data=params[params$spp=="Adi",], aes(GT, log(lam.mn), label=AB), size=1.8, fontface="bold")+
geom_point(data=params, aes(GT, log(lam.mn), col=spp), size=3, alpha=0.5)+
#geom_point(shape=21, aes(col=spp), alpha=0.2)+
#coord_cartesian(xlim=c(min(boot$GT),100))+
scale_colour_manual(values=cols)+
scale_fill_manual(values=cols)+
scale_x_log10()+
geom_text(data=params[!params$spp=="Adi",], aes(GT, log(lam.mn), label=AB), size=1.8, fontface="bold")+
guides(col="none", fill="none")+
theme_classic()+
ggtitle("Generation times")+
labs(x="Generation Time (years)", y=expression(log(lambda)))+theme(plot.title=element_text(size=8, hjust=0.5, face="bold"), axis.title=element_text(size=9), axis.text=element_text(size=8))
GTplot



figCD <- plot_grid(growplot,ggplot()+theme_void(),rechange, GTplot,
ncol=1, rel_heights=c(0.44,-0.07,1, 1), align="v", axis="lr", labels=c("b","","","d"), label_size=9)
#figCD

fig.2 <- plot_grid(plot_grid(examples, lam.dens, nrow=2, rel_heights=c(0.45, 1), labels=c("a","c"), label_size=9),figCD, nrow=1, rel_widths=c(1,1))+
annotation_custom(dig,0.29,0.38,0.905,0.93)+
annotation_custom(tab,0.09,0.13,0.9,0.94)+
annotation_custom(mas,0.05,0.13,0.56,0.59)+
annotation_custom(dig,0.05,0.13,0.45,0.49)+
annotation_custom(cor,0.05,0.13,0.34,0.38)+
annotation_custom(brn,0.04,0.14,0.24,0.29)+
annotation_custom(tab,0.06,0.12,0.125,0.155)+
#draw_text("Recruitment needed\nfor pop. growth", 0.55, 0.895, hjust=0, size=6.5)+
draw_label(x=0.63, y=0.89, label=expression(italic(P[rec])~at~lambda>=1), hjust=0, size=7)

fig.2
