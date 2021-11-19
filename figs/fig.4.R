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
brn<-readPNG("data/coral_silhouettes/branchingG.png")
brn<-rasterGrob(brn, interpolate=TRUE)



params$AB <- ifelse(params$abundance_pair=="Rare","R","C")
params2$AB <- ifelse(params2$abundance_pair=="Rare","R","C")


ek.av$morphology<-factor(ek.av$morphology, levels=c("massive","digitate","corymbose","staghorn","tabular"))
sdat$morphology<-factor(sdat$morphology, levels=c("massive","digitate","corymbose","staghorn","tabular"))
s.avs <- aggregate(area~morphology, sdat, median)
segs2<-data.frame(morphology="tabular")
segs2$morphology<-factor(segs2$morphology, levels=c("massive","digitate","corymbose","staghorn","tabular"))

size.elas<-ggplot(ek.av[!ek.av$spp=="Asp",], aes(size,X2))+
geom_bar(stat="identity", position=position_dodge(preserve = "single"), aes(fill=spp), col="black", size=0.1, width=0.33)+
#geom_line()+
#geom_point()+
facet_wrap(~morphology, ncol=1)+
scale_fill_manual(values=cols)+
scale_colour_manual(values=cols)+
geom_boxplot(data=sdat[!sdat$spp=="Asp",], aes(x=area, y=0.9, fill=spp), width=0.2, outlier.size=0.01, size=0.1, outlier.colour="grey")+
guides(fill="none", col="none")+
geom_segment(data=s.avs, aes(x=area, xend=area, y=0, yend=0.7), linetype="dotted")+
#geom_point(data=ek.av.av, aes(x=X3, y=0.6, col=spp))+
labs(x=expression(size~(m^2)), y="elasticity")+
scale_y_continuous(expand=c(0,0), breaks=c(0,0.5,1))+
theme_classic()+theme(strip.background=element_blank(), strip.text=element_text(size=8),
axis.text.y=element_text(size=6), axis.text.x=element_text(size=6, angle=45, hjust=1), axis.title=element_text(size=8), axis.line=element_line(size=0.1), plot.title=element_text(hjust=0.5, size=8, face="bold"))+
geom_segment(data=segs2, aes(x=0.5, xend=0.5, y=0.8, yend=1))+
geom_segment(data=segs2, aes(x=0.5, xend=0.4, y=0.8, yend=0.8))+
geom_segment(data=segs2, aes(x=0.5, xend=0.4, y=1, yend=1))+
geom_segment(data=segs2, aes(x=0.5, xend=0.7, y=0.9, yend=0.9))+
geom_segment(data=segs2, aes(x=0.7, xend=0.7, y=0.5, yend=0.9))+
geom_text(data=segs2,aes(x=1, y=0.35, label="colony \nsizes"), size=2, hjust=1)+
ggtitle(expression(bold(Effects~of~colony~sizes~on~lambda)))+
#ggtitle("Demographic sensitivity to size")+
scale_x_continuous(breaks=c(-4,-3,-2,-1,0), labels=c(expression(0.0001),expression(0.001), expression(0.01),  expression(0.1),  expression(1)))
size.elas


params$morphology <- factor(params$morphology, levels=c("staghorn", "tabular", "corymbose","digitate","massive"))
elasplot<-ggplot(params, aes(reorder(spp, eR), eR, fill=spp))+
geom_bar(stat="identity", size=0.1, col="black")+
scale_fill_manual(values=cols)+
guides(fill="none")+
facet_grid(.~morphology, scales="free_x", space="free_x")+
labs(y=expression(R[elasticity]))+
geom_text(aes(y=eR+0.02,label=AB), size=1.8)+
scale_y_continuous(expand=c(0,0), limits=c(0,0.45), breaks=c(0,0.1, 0.2, 0.3))+
ggtitle("Demographic sensitivity \nto reproduction")+
theme_classic()+theme(strip.text=element_blank(),
 strip.background=element_blank(), axis.title.x=element_blank(),  axis.text.y=element_text(size=7),axis.text.x=element_text(size=7, angle=90, vjust=0.5), axis.title.y=element_text(size=8), plot.title=element_text(size=8, hjust=0.5, face="bold"), axis.line.y=element_line(size=0.1))#+coord_flip()
elasplot


params$ec.se<-aggregate(Carbon_ug_corrected~spp+morph, ec, sd)$Carbon_ug_corrected


reproplot<-ggplot()+
#geom_path(data=params2, aes(f.int, eggC, group=morph, col=morphology), linetype="dotted", size=0.2)+
#geom_path(data=params2, aes(f.int, eggC, col=morph), arrow=arrow(type="closed", length=unit(3,"mm")),linetype="dotted", size=0.2)+
#geom_point(col="white", stroke=0.1, size=6)+
geom_smooth(data=params[params$family=="Acroporidae",], aes(fec1cm, eggC), method="lm", formula=y~poly(x,1), se=F, size=0.2, col="black")+
geom_segment(data=params,aes(x=fec1cm, xend=fec1cm, y=eggC-ec.se, yend=eggC+ec.se), size=0.2)+
geom_point(data=params,aes(fec1cm, eggC, fill=spp), shape=21, stroke=0.2, size=3)+
#geom_point(data=params[params$spp %in% c("Aro","Acy","Ahu","Ami","Gpe"),], aes(f.int, eggC, fill=spp), shape=21, stroke=0.1, size=3.5)+
geom_text(data=params, aes(fec1cm, eggC, label=AB), size=1.8)+
#geom_text(data=params[params$AB=="R",], aes(fec1cm, eggC, label=AB), size=1.8)+
#geom_text_repel(data=params, aes(f.int+0, eggC+0,label=spp), size=2, force=0.1)+
geom_segment(aes(x=335, xend=350, y=47, yend=43),col="grey",#colsC[3], 
arrow=arrow(type="closed", length=unit(0.8,"mm")), size=0.2)+
geom_segment(aes(x=405, xend=435, y=41.5, yend=40),col="grey",#colsC[2], 
arrow=arrow(type="closed", length=unit(0.8,"mm")), size=0.2)+
geom_segment(aes(x=700, xend=800, y=38, yend=36),col="grey",#colsC[1], 
arrow=arrow(type="closed", length=unit(0.8,"mm")), size=0.2)+
geom_segment(aes(x=900, xend=1000, y=28, yend=33),col="grey",#colsC[4], 
arrow=arrow(type="closed", length=unit(0.8,"mm")), size=0.2)+
geom_segment(aes(x=920, xend=1020, y=27, yend=27.5),col="grey",#colsC[5], 
arrow=arrow(type="closed", length=unit(0.8,"mm")), size=0.2)+
geom_segment(aes(x=880, xend=930, y=10, yend=10.5),col="grey",#colsC[6], 
arrow=arrow(type="closed", length=unit(0.8,"mm")), size=0.2)+
ggtitle("Reproductive investments")+
scale_x_log10(breaks=c(400, 600, 900, 1200))+
geom_text(data=NULL, aes(330, 10, label='C = Common'), size=1.8, hjust=0)+
geom_text(data=NULL, aes(330, 6, label='R = Rare'), size=1.8, hjust=0)+
#ylim(3,55)+
labs(x=expression(Egg~number~(eggs~cm^-2)), y= "Egg mass (g of Carbon)")+
scale_fill_manual(values=cols)+guides(fill="none", col="none")+
scale_colour_manual(values=colsC)+
theme_classic()+rectheme+
theme(plot.title=element_text(size=8, hjust=0.5))
reproplot


fig4 <- plot_grid(size.elas, plot_grid(elasplot,reproplot, ncol=1, rel_heights=c(1, 1.2),labels=c("B","C"), label_size=9), nrow=1, rel_widths=c(1,1),labels=c("A",""), label_size=9)+
draw_plot(mas, 0.92, 0.62, 0.06, 0.1)+
draw_plot(dig, 0.85, 0.65, 0.06, 0.1)+
draw_plot(cor, 0.79, 0.79, 0.06, 0.1)+
draw_plot(tab, 0.67, 0.81, 0.094, 0.1)+
draw_plot(brn, 0.62, 0.83, 0.035, 0.1)
fig4 


#fig4 <- plot_grid(size.elas, plot_grid(elasplot,reproplot, nrow=1, rel_widths=c(1, 1.2),labels=c("B","C"), label_size=9), ncol=1, rel_heights=c(1,1.1),labels=c("A",""), label_size=9)+
#draw_plot(mas, 0.39, 0.1, 0.06, 0.1)+
#draw_plot(dig, 0.33, 0.12, 0.06, 0.1)+
#draw_plot(cor, 0.26, 0.29, 0.06, 0.1)+
#draw_plot(tab, 0.16, 0.31, 0.094, 0.1)+
#draw_plot(brn, 0.12, 0.36, 0.035, 0.1)
#fig4 

#fig4 <- plot_grid(size.elas, plot_grid(elasplot,reproplot, ncol=1, rel_heights=c(1, 1.2),labels=c("B","C"), label_size=9), ncol=2, rel_widths=c(1.8,1),labels=c("A",""), label_size=9)+
#draw_plot(mas, 0.39, 0.1, 0.06, 0.1)+
#draw_plot(dig, 0.33, 0.12, 0.06, 0.1)+
#draw_plot(cor, 0.26, 0.29, 0.06, 0.1)+
#draw_plot(tab, 0.16, 0.31, 0.094, 0.1)+
#draw_plot(brn, 0.12, 0.36, 0.035, 0.1)
#fig4 
