
##### FIGURE S2

head(abun.BT)

figS2A<-plot_grid(
ggplot(abun.BT, aes(reorder(spp, -abundance_05), N/10, fill=spp))+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0)+
guides(fill="none")+
labs(y=expression(Colonies~per~m^2))+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
labs(title="N colonies (27 x BTs)",subtitle="2005")+
theme(axis.title.x=element_blank(), 
axis.text.x=element_blank(),
plot.title=element_text(size=8, hjust=0.5, face="bold"),
plot.subtitle=element_text(size=9, hjust=0.5),strip.background=element_blank())
,
ggplot(data=abun.LIT, aes(reorder(species, -abundance_05), N/10, fill=spp))+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0)+
facet_wrap(~year, ncol=1)+
guides(fill="none")+
labs(title="N colonies (6 x LITs)",y="Colonies per m")+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
theme(axis.title.x=element_blank(), 
plot.title=element_text(size=8, hjust=0.5, face="bold"),
axis.text.x=element_text(angle=55, hjust=1, face="italic"), strip.background=element_blank()),
ncol=1, rel_heights=c(1,3))
figS2A


figS2B <- ggplot(data=abun.LIT, aes(reorder(species, -abundance_05), cover, fill=spp))+
#geom_bar(data=tri2.av, aes(reorder(species, -abun05), N/10, fill=species),stat="identity")+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0)+
facet_wrap(~year, ncol=1)+
guides(fill="none")+
labs(y="% cover")+
ggtitle("Coral Cover")+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
theme(axis.title.x=element_blank(), 
plot.title=element_text(size=8, hjust=0.5, face="bold"),
axis.text.x=element_text(angle=55, hjust=1, face="italic"), strip.background=element_blank())
figS2B

#######################################
# SADs
#######################################

sad <- read.csv("data/BT_counts_2005.csv")
sad <- aggregate(Abundance~Species, sad, sum)

sad$Species[sad$Species=="Acropora fat dig"] <- "Acropora cf. digitifera"
sad$spp <- params$spp[match(sad$Species, params$species)]
sad$spp <- as.character(sad$spp)
sad$spp[is.na(sad$spp)]<-"Other"
cols2 <- c(cols, "Other"="Grey")
sad <- sad[order(sad$Abundance),]
sad$rank <- c(nrow(sad):1)
sad$morphology<-params$morph[match(sad$spp, params$spp)]

sadplot <- plot_grid(
ggplot(sad[!is.na(sad$morphology),], aes(y=morphology, x=rank))+
geom_line(aes(col=morphology), arrow=arrow(length=unit(1,"mm")))+
ggtitle("Abundance distribution")+
xlim(c(0,70))+
geom_text(data=sad[sad$spp %in% c("Ahy","Adi","Ana","Ain","Gre"),], aes(label=morphology), size=1.8, hjust=1)+
scale_colour_manual(values=colsC)+
theme_void()+
theme(plot.title=element_text(size=8, hjust=0.5, face="bold"))+
guides(col="none"),
ggplot(subset(sad, Abundance>50), aes(rank, Abundance/2700))+geom_segment(aes(rank, xend=rank, y=0, yend=Abundance/2700, col=spp), size=1)+
scale_colour_manual(values=cols2)+
labs(x="Rank", y="Colony density (2005)")+
geom_text_repel(data=sad[!sad$spp=="Other",], aes(rank, y=(Abundance/2700)+0.2, label=spp), angle=45, size=2, force=0.001, hjust=0)+
xlim(c(0,70))+
theme_classic()+
guides(colour="none"),
ncol=1, align="v", axis="lr", rel_heights=c(0.3,1))
sadplot


diffs$diff2 <- ((diffs$X2014+0.05)/(diffs$X2011+0.05))^(1/3)
diffs$spp <- params$spp[match(diffs$species, params$species)]
diffs$lam.est <- params$lam.est[match(diffs$spp, params$spp)]


figS2D <- ggplot(diffs, aes(reorder(spp, -lam.est), diff, fill=spp))+
stat_summary(geom="bar", fun="mean", col="black", size=0.1)+
scale_fill_manual(values=cols)+guides(fill="none")+theme_classic()+theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_text(size=10),plot.title=element_text(size=8, hjust=0.5, face="bold"))+labs(y=expression(Delta*N))+
ggtitle("Change in abundance on \ntransects (2011-2014)")

figS2E <-ggplot(params, aes(reorder(species, -lam.est), log(lam.est), fill=spp))+
geom_hline(yintercept=0)+
geom_bar(stat="identity", size=0.1, col="black")+
#geom_point(size=2, shape=21, stroke=0.1)+
scale_fill_manual(values=cols)+
guides(fill="none", col="none")+theme_classic()+theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=55, hjust=1, face="italic"), axis.title.y=element_text(size=10), plot.title=element_text(size=8, hjust=0.5, face="bold"))+labs(y=expression(log(lambda)[~IPMs]))+ggtitle("Fitness from IPMs")
figS2E

FigS2new <- plot_grid(figS2D, ggplot()+theme_void(), figS2E, align="v", axis="lr", ncol=1, rel_heights=c(1,0.2,1.4), labels=c("D","","E"), label_size=12)
FigS2new

fig.s2 <- plot_grid(figS2A, plot_grid(sadplot, plot_grid(figS2B,NULL, rel_widths=c(1,0.2)), ncol=1, rel_heights=c(0.7,1), labels=c("B","C"), label_size=12),FigS2new, nrow=1, rel_widths=c(1,1.2,1), labels=c("A","", ""), label_size=12)
fig.s2

