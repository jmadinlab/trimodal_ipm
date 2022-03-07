

diffs$diff2 <- ((diffs$X2014+0.05)/(diffs$X2011+0.05))^(1/3)
fig.s4A <- ggplot(diffs, aes(reorder(spp, -lam), diff, fill=spp))+
stat_summary(geom="bar", fun="mean", col="black", size=0.1)+
#stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
#stat_summary(aes(reorder(spp, -diff), cover, fill=spp), geom="point", fun="mean", shape=4)+
scale_fill_manual(values=cols)+guides(fill="none")+theme_classic()+theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=10),plot.title=element_text(size=8, hjust=0.5, face="bold"))+labs(y=expression(Delta*N))+
ggtitle("Change in abundance on \ntransects (2011-2014)")

recfit<-melt(params[,c("spp","rec.tran","rec.mean")])
recfit$variable2 <- ifelse(recfit$variable=="rec.mean","constant","fitted")
recfitplot <- ggplot()+geom_line(data=recfit, aes(value, variable2, group=spp), size=0.1)+geom_point(data=recfit, aes(value, variable2, fill=spp), shape=21, stroke=0.1, size=3)+geom_point(data=NULL, aes(x=params$rec.mean[1], y="constant"), size=3)+geom_point(data=NULL, aes(x=params$rec.mean[11], y="constant"), size=3)+scale_x_log10()+scale_fill_manual(values=cols)+guides(fill="none")+
#scale_y_discrete(limits=rev)+
theme_classic()+labs(x=expression(Recruitment~probability~(italic(P[rec]))))+
geom_text(data=NULL, aes(x=params$rec.mean[11], y=0.7, label="Gon"), size=3)+geom_text(data=NULL, aes(x=params$rec.mean[1], y=0.7, label="Acr"), size=3)+theme(axis.title.y=element_blank(),plot.title=element_text(size=8, hjust=0.5, face="bold"))+
ggtitle("Fitted vs. constant recruitment")
recfitplot 

fig.s4B <-ggplot(params, aes(reorder(spp, -lam.est), log(lam.est), fill=spp))+
geom_hline(yintercept=0)+
geom_bar(stat="identity", size=0.1, col="black")+
#geom_point(size=2, shape=21, stroke=0.1)+
scale_fill_manual(values=cols)+
guides(fill="none", col="none")+theme_classic()+theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5), axis.title.y=element_text(size=10), plot.title=element_text(size=8, hjust=0.5, face="bold"))+labs(y=expression(log(lambda)[~IPMs]))+ggtitle("Fitness (constant recruitment)")

fig.s4C <-ggplot(params[params$abundance_pair %in% c("Common","Rare"),], aes(log(lam.est), log(lam.tran)))+
geom_abline(slope=1, size=0.1)+geom_smooth(method="lm", linetype="dotted")+
geom_point(aes(col=spp), size=3.5)+
geom_text(aes(label=spp), size=2)+
#scale_y_log10()+scale_x_log10()+
theme_classic()+ggtitle("Fitness vs change \nin abundance")+
labs(y="Change in \nabundance (2011-14)", x=expression(lambda~(constant~recruitment)))+theme(plot.title=element_text(size=8, hjust=0.5, face="bold"))+
scale_colour_manual(values=cols)+guides(col="none")

fig.s4 <- plot_grid(plot_grid(fig.s4A, fig.s4B, ncol=1, rel_heights=c(1,0.9), align="v", labels=c("A","B"), label_size=12), fig.s4C,
labels=c("","C"), label_size=12)
fig.s4


fig.s4 <- plot_grid(plot_grid(fig.s4A, fig.s4B, ncol=1, rel_heights=c(1,1.1), align="v", labels=c("A","C"), label_size=12), plot_grid(recfitplot,fig.s4C,
labels=c("B","D"), ncol=1, label_size=12, rel_heights=c(0.6, 1)), rel_widths=c(0.9,1))
fig.s4

