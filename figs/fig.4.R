



params$AB <- ifelse(params$abundance_pair=="Common","C","R")
params$AB[params$spp=="Ana"]<-"D"



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
geom_text(data=NULL, aes(330, 8, label='R = Rare'), size=1.8, hjust=0)+
geom_text(data=NULL, aes(330, 6, label='D = Declined'), size=1.8, hjust=0)+
#ylim(3,55)+
labs(x=expression(Egg~number~(eggs~cm^-2)), y= "Egg mass (g of Carbon)")+
scale_fill_manual(values=cols)+guides(fill="none", col="none")+
scale_colour_manual(values=colsC)+
theme_classic()+
theme(plot.title=element_text(size=8, hjust=0.5, face="bold"), axis.title=element_text(size=10))
reproplot

arrows3 <- dcast(params[!params$spp=="Ana", ], morphology~abundance_pair,  value.var="lam.mn")
arrows3$GT.C <- dcast(params[!params$spp=="Ana", ], morphology~abundance_pair,  value.var="GT")$Common
arrows3$GT.R <- dcast(params[!params$spp=="Ana", ], morphology~abundance_pair,  value.var="GT")$Rare
arrows3$Common2 <- ifelse(arrows3$GT.R > 20, arrows3$Common * 0.98, arrows3$Common * 0.95)
arrows3$Rare2 <- ifelse(arrows3$GT.R > 20, arrows3$Rare * 1.02, arrows3$Rare * 1.1)
arrows3$GT.C2 <- ifelse(arrows3$GT.C > 20, arrows3$GT.C *1.1, arrows3$GT.C * 1)
arrows3$GT.R2 <- ifelse(arrows3$GT.R > 20, arrows3$GT.R *0.9, arrows3$GT.R * 1)

GTplot<- ggplot()+geom_hline(yintercept=1, size=0.1)+
geom_segment(data=params, aes(x=GT, xend=GT, y=lam.sd1, yend=lam.sd2), size=0.2)+
#geom_segment(data=arrows3, aes(x=GT.R2, xend=GT.C2, y=Rare2, yend=Common2), col="grey", arrow=arrow(type="closed", length=unit(0.8,"mm")))+
geom_point(data=params,aes(GT, lam.mn, fill=spp), shape=21, stroke=0.2, size=3)+
geom_text(data=params, aes(GT, lam.mn, label=AB), size=1.8)+
#scale_x_log10()+
#scale_y_log10()+
scale_fill_manual(values=cols)+guides(fill="none")+
theme_classic()+
ggtitle("Generation times")+
labs(x="Generation Time (years)", y=expression(Fitness~(lambda)))+theme(plot.title=element_text(size=8, hjust=0.5, face="bold"), axis.title=element_text(size=10))


fig.4 <- reproplot #plot_grid(GTplot, reproplot, label_size=9, labels=c("A", "B"))
