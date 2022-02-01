
examples1 <- 
plot_grid(NULL, plot_grid(NULL, 
plots[["Ahy"]]+theme(legend.position="right", plot.title=element_blank(),legend.box="vertical",legend.spacing.y=unit(0.5,"mm"), legend.margin=margin(0,0,0,-7)),
NULL,
plots[["Gpe"]]+theme(legend.position="right", plot.title=element_blank(),legend.box="vertical",legend.spacing.y=unit(0.5,"mm"), legend.margin=margin(0,0,0,-7)),NULL,
#guides(col="none",fill="none")
ncol=1, rel_heights=c(0.35,1,-0.1,1,0.1))+
draw_label(expression(size[" "*t]~(m^2)), 0.65, 0.08, size=8, hjust=1)+
draw_label(expression(size[" "*t+1]~(m^2)), -0.05, 0.55, size=8, angle=90)+
draw_label("Example IPMs", 0.45, 0.95, fontface="bold", size=8)+
draw_label(expression("("*bold(before)~bolditalic(P[rec])*")"), 0.45, 0.9, fontface="bold", size=8),
nrow=1, rel_widths=c(0.12,1))
examples1

examples2 <- 
plot_grid(NULL,plot_grid(NULL,
plots[["Ahy"]]+ggtitle("Ahy")+theme(legend.position="right", plot.title=element_blank(),legend.box="vertical",legend.spacing.y=unit(0.5,"mm"), legend.margin=margin(0,0,0,-7)),
plots[["Gpe"]]+ggtitle("Adi")+theme(legend.position="right", plot.title=element_blank(),legend.box="vertical", legend.spacing.y=unit(0.5,"mm"), legend.margin=margin(0,0,0,-7)),
nrow=1, rel_widths=c(0.1,1,1))+
draw_label(expression(size[" "*t]~(m^2)), 0.6, 0.17, size=8, hjust=1)+
draw_label(expression(size[" "*t+1]~(m^2)), 0.03, 0.65, size=8, angle=90)+
draw_label(expression(bold(Example~IPMs~(before~bolditalic(P[rec])))), 0.55, 1.05, fontface="bold", size=8), nrow=2, rel_heights=c(0.15,1))
examples2

arrows2 <- dcast(params[!params$spp=="Ana", ], morphology~abundance_pair,  value.var="lam.est")
arrows2$Common2 <- arrows2$Common*0.95
arrows2$Rare2 <- arrows2$Rare*1.05
params$AB <- ifelse(params$abundance_pair=="Rare","R","C")

arrows2$yaxis <- c(30,38,35,14,11)
params$yaxis <- arrows2$yaxis[match(params$morphology, arrows2$morphology)]
arrows2

arrows2$diff <- arrows2$Common - arrows2$Rare

boot$morphology <- params$morphology[match(boot$spp, params$spp)]
morph.ord <- c("massive","digitate","corymbose","staghorn","tabular")
boot$morphology <- factor(boot$morphology, levels=morph.ord)
arrows2$morphology <- factor(arrows2$morphology, levels=morph.ord)
params$morphology <- factor(params$morphology, levels=morph.ord)

lab.point <- data.frame(morphology="massive", x=c(1.5, 1.5, 1.5), y=c(30, 20, 10), lab = c('C = Common','R = Rare','D = Declined'))
lab.point$morphology <- factor(lab.point$morphology, levels=morph.ord)

lam.dens <- plot_grid(NULL, ggplot(boot)+
geom_rect(data=lab.point, aes(xmin=1.4, xmax=1.9, ymin=4, ymax=37), fill=NA, col="grey", size=0.5)+
geom_segment(data=arrows2, aes(x=-Inf, xend=Inf, y=yaxis*1.2, yend=yaxis*1.2), col="white")+
geom_segment(data=arrows2, aes(y=yaxis, yend=yaxis, x=Rare2, xend=Common2), arrow=arrow(length=unit(1, "mm")))+
geom_text(data=lab.point, aes(x,y, label=lab), size=1.8, hjust=0, nudge_y=0.3)+
geom_density(aes(x=lam, col=spp, fill=spp), alpha=0.5)+
geom_vline(xintercept=1, linetype="dotted")+
geom_point(data=params, aes(lam.est, yaxis,col=spp), shape=21, size=3, stroke=0.7, fill="white")+
geom_point(data=params, aes(lam.est, yaxis, col=spp), size=3, alpha=0.5)+
geom_text(data=params[!params$spp=="Ana",], aes(lam.est, yaxis, label=AB), size=1.8, fontface="bold")+
geom_text(data=params[params$spp=="Ana",], aes(lam.est, yaxis, label="D"), size=1.8, fontface="bold")+
geom_hline(yintercept=0, size=1, col="white")+
geom_hline(yintercept=-0.1, size=1, col="black")+
#scale_y_sqrt(expand=c(0,0))+
scale_y_continuous(expand=c(0,0))+
lims(x=c(min(boot$lam), max(boot$lam)))+
ggtitle(expression(bold(Fitness~(lambda))))+
xlim(c(0.5, 2))+
#scale_x_log10()+
facet_wrap(~morphology, ncol=1, scales="free_y", strip.position="left")+
labs(y="Density", x=expression(lambda))+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
theme_classic()+theme(strip.text=element_blank(), strip.background=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.line.y=element_blank(), axis.ticks.y=element_blank(), plot.title=element_text(size=8, hjust=0.5, face="bold")), nrow=1, rel_widths=c(0.1,1))
lam.dens

plot_grid(examples1, lam.dens, nrow=1, rel_widths=c(0.5, 1))




r.long$label <-ifelse(r.long$variable=="rec.tran", "Transects","Constant")
#r.long$label <- factor(r.long$label, levels=c("Fitted","Constant"))
params$genus <- ifelse(params$family=="Acroporidae", "Acr","Gon")

params$AB <- ifelse(params$abundance_pair=="Rare","R","C")

arrows <- dcast(params[!params$spp=="Ana", ], morphology~abundance_pair,  value.var="rec.one")
arrows$Common2 <- arrows$Common*1.5
arrows$Rare2 <- arrows$Rare*0.6
arrows[arrows$morphology %in% c("staghorn", "tabular"),c("Rare2","Common2")] <- arrows[arrows$morphology %in% c("staghorn", "tabular"),c("Rare","Common")]

arrows$diff <-   arrows$Rare/arrows$Common

fitplot <- ggplot()+
geom_path(data=r.long, aes(label, value, group=spp), size=0.1, col="grey")+
#geom_point(data=r.long, aes(label, value, fill=spp), shape=21, stroke=0.1, size=2.7)+
geom_point(data=r.long, aes(label, value, col=spp), shape=4,  size=1, stroke=1)+
#geom_point(data=params, aes(x="Constant",rec.mean), fill="grey", shape=21, stroke=0.5, size=2)+
geom_point(data=params, aes(x="Constant",rec.mean), fill="grey", shape=4,  size=1, stroke=1)+
geom_text(data=params[params$spp %in% c("Gre","Ahy"),], aes(x="Constant",rec.mean, label=genus), nudge_y=0, nudge_x=-0.4, angle=0, hjust=0.5, size=2.2, fontface="italic")+
theme_classic()+
scale_y_log10(limits=c(min(storeDET$rec),max(storeDET$rec)), breaks=c(10^-5, 10^-4,10^-3,10^-2, 10^-1), labels=c(
		expression(10^-5), expression(10^-4), expression(10^-3),expression(10^-2), expression(10^-1)))+
scale_fill_manual(values=cols)+scale_colour_manual(values=cols)+
coord_flip()+
guides(col="none", fill="none")+ggtitle("Recruitment vs fitness")+theme_void()+theme(plot.title=element_text(size=8, face="bold", hjust=0.5))
fitplot



rechange <- ggplot(storeDET, aes(rec, log(lam)))+geom_line(aes(col=spp))+
scale_x_log10(limits=c(min(storeDET$rec),max(storeDET$rec)), breaks=c(10^-5, 10^-4,10^-3,10^-2, 10^-1), labels=c(
		expression(10^-5), expression(10^-4), expression(10^-3),expression(10^-2), expression(10^-1)))+
geom_hline(yintercept=0, size=0.1)+
guides(col="none")+
theme_classic()+
#ggtitle("Recruitment vs fitness")+
labs(x=expression(Recruitment~probablity~(P[rec])), y=expression(log(lambda)))+
coord_cartesian(ylim=c(-0.5,1))+
#coord_cartesian(ylim=c(0.6,2))+
scale_colour_manual(values=cols)+
theme(axis.title.x=element_blank(), plot.background=element_blank(), axis.text=element_text(size=8), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), plot.title=element_text(size=8, hjust=0.5, face="bold"))
rechange 


growplot <- ggplot()+
geom_segment(data=arrows, aes(y=morphology, yend=morphology, x=Rare2, xend=Common2), arrow=arrow(length=unit(1, "mm")))+
#geom_point(aes(fill=spp), shape=21, size=3, stroke=0.1)+
geom_point(data=params, aes(rec.one, morphology,col=spp), shape=21, size=3, stroke=0.7, fill="white")+
geom_point(data=params, aes(rec.one, morphology, col=spp), size=3, alpha=0.5)+
scale_x_log10(limits=c(min(storeDET$rec),max(storeDET$rec)), breaks=c(10^-5, 10^-4,10^-3,10^-2, 10^-1), labels=c(
		expression(10^-5), expression(10^-4), expression(10^-3),expression(10^-2), expression(10^-1)))+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
#geom_text(data=NULL, aes(y="tabular", x=10^-2, label='C = Common'), size=1.8, hjust=0, nudge_y=0.3)+
#geom_text(data=NULL, aes(y="tabular", x=10^-2, label='R = Rare'), size=1.8, hjust=0, nudge_y=-0.3)+
#geom_text(data=NULL, aes(y="tabular", x=10^-2, label='D = Declined'), size=1.8, hjust=0, nudge_y=-0.9)+
theme_classic()+
geom_text(data=params[!params$spp=="Ana",], aes(rec.one, morphology, label=AB), size=1.8, fontface="bold")+
geom_text(data=params[params$spp=="Ana",], aes(rec.one, morphology, label="D"), size=1.8, fontface="bold")+
labs(x=expression(Recruitment~probablity~(P[rec])))+ 
ggtitle("Recruitment needed\nfor pop. growth")+
theme(plot.title=element_text( size=7, hjust=1, vjust=-10), axis.title.y=element_blank(), axis.text=element_text(size=8), axis.title.x=element_text(size=9), plot.background=element_blank())
growplot


figCD <- plot_grid(fitplot, ggplot()+theme_void(), rechange, ggplot()+theme_void(), growplot,
ncol=1, rel_heights=c(0.4,-0.05, 1,-0.15,  1), align="v", axis="lr", labels=c("C","","",""), label_size=9)
figCD

fig.2 <- plot_grid(plot_grid(examples2, lam.dens, nrow=2, rel_heights=c(0.45, 1), labels=c("A","C"), label_size=9),figCD, nrow=1, rel_widths=c(1,1))+
annotation_custom(dig,0.28,0.37,0.91,0.93)+
annotation_custom(tab,0.07,0.11,0.91,0.94)+
annotation_custom(mas,0.05,0.13,0.54,0.57)+
annotation_custom(dig,0.05,0.13,0.43,0.47)+
annotation_custom(cor,0.05,0.13,0.32,0.36)+
annotation_custom(brn,0.04,0.14,0.22,0.27)+
annotation_custom(tab,0.06,0.12,0.105,0.135)
fig.2






fig.2b <- plot_grid(examples1, figCD, lam.dens, nrow=1, rel_widths=c(0.55, 1, 0.75), labels=c("A","B","C"), label_size=9)+
#annotation_custom(dig,0.28,0.37,0.91,0.93)+
#annotation_custom(tab,0.07,0.11,0.91,0.94)+
#annotation_custom(mas,0.05,0.13,0.54,0.57)+
#annotation_custom(dig,0.05,0.13,0.43,0.47)+
#annotation_custom(cor,0.05,0.13,0.32,0.36)+
#annotation_custom(brn,0.04,0.14,0.22,0.27)+
annotation_custom(tab,0.06,0.12,0.105,0.135)
fig.2b



fig.2c <- plot_grid(plot_grid(plot_grid(NULL, examples2, rel_widths=c(0.12, 1)), figCD,  nrow=2, rel_heights=c(0.45, 1), labels=c("A","B"), label_size=9),lam.dens, nrow=1, rel_widths=c(1,0.8))+
#annotation_custom(dig,0.28,0.37,0.91,0.93)+
#annotation_custom(tab,0.07,0.11,0.91,0.94)+
#annotation_custom(mas,0.05,0.13,0.54,0.57)+
#annotation_custom(dig,0.05,0.13,0.43,0.47)+
#annotation_custom(cor,0.05,0.13,0.32,0.36)+
#annotation_custom(brn,0.04,0.14,0.22,0.27)+
annotation_custom(tab,0.06,0.12,0.105,0.135)
fig.2c

