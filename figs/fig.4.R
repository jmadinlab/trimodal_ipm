

params$AB <- ifelse(params$abundance_pair=="Common","C","R")
params$AB[params$spp=="Ana"]<-"D"

reproplot<-ggplot()+
geom_smooth(data=params[params$family=="Acroporidae",], aes(fec1cm, eggC), method="lm", formula=y~poly(x,1), se=F, size=0.2, col="black")+
geom_segment(data=params,aes(x=fec1cm, xend=fec1cm, y=eggC-eggC.se, yend=eggC+eggC.se), size=0.2)+
geom_point(data=params,aes(fec1cm, eggC, fill=spp), shape=21, stroke=0.2, size=3)+
geom_text(data=params, aes(fec1cm, eggC, label=AB), size=1.8)+
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

fig.4 <- reproplot 