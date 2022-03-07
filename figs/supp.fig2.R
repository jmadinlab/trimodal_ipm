
head(d08B)

triplot<-plot_grid(
ggplot(d08B, aes(reorder(spp, -abundance_05), N/10, fill=spp))+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0)+
guides(fill="none")+
labs(y=expression(Colonies~per~m^2))+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
labs(title="27 x BTs",subtitle="2005")+
theme(axis.title.x=element_blank(), 
axis.text.x=element_blank(),
plot.title=element_text(size=8, hjust=0.5, face="bold"),
plot.subtitle=element_text(size=9, hjust=0.5),strip.background=element_blank())
,
ggplot(data=abun, aes(reorder(species, -abundance_05), N/10, fill=spp))+
#geom_bar(data=tri2.av, aes(reorder(species, -abun05), N/10, fill=species),stat="identity")+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0)+
facet_wrap(~year, ncol=1)+
guides(fill="none")+
labs(title="6 x LITs",y="Colonies per m")+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
theme(axis.title.x=element_blank(), 
plot.title=element_text(size=8, hjust=0.5, face="bold"),
axis.text.x=element_text(angle=55, hjust=1, face="italic"), strip.background=element_blank()),
ncol=1, rel_heights=c(1,3))
triplot


covplot <- ggplot(data=abun, aes(reorder(species, -abundance_05), cover, fill=spp))+
#geom_bar(data=tri2.av, aes(reorder(species, -abun05), N/10, fill=species),stat="identity")+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0)+
facet_wrap(~year, ncol=1)+
guides(fill="none")+
labs(y="% cover")+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
theme(axis.title.x=element_blank(), 
plot.title=element_text(size=8, hjust=0.5, face="bold"),
axis.text.x=element_text(angle=55, hjust=1, face="italic"), strip.background=element_blank())
covplot

#######################################
# SADs
#######################################

sad <- read.csv("data/abundance/dornelas2008.csv")

head(sad)
sad$species[sad$species=="Acropora fat dig"] <- "Acropora cf. digitifera"
sad$spp <- params$spp[match(sad$species, params$species)]
sad$spp <- as.character(sad$spp)
sad$spp[is.na(sad$spp)]<-"Other"
cols2 <- c(cols, "Other"="Grey")
sad
sad <- sad[order(sad$abundance),]
sad$rank <- c(nrow(sad):1)
sad$morphology<-params$morph[match(sad$spp, params$spp)]

sadplot <- plot_grid(
ggplot(sad[!is.na(sad$morphology),], aes(y=morphology, x=rank))+
geom_line(aes(col=morphology), arrow=arrow(length=unit(1,"mm")))+
xlim(c(0,70))+
geom_text(data=sad[sad$spp %in% c("Ahy","Adi","Ana","Ain","Gre"),], aes(label=morphology), size=1.8, hjust=1)+
scale_colour_manual(values=colsC)+
theme_void()+
#scale_colour_manual(values=colsC)+
guides(col="none"),
ggplot(subset(sad, abundance>50), aes(rank, abundance/2700))+geom_segment(aes(rank, xend=rank, y=0, yend=abundance/2700, col=spp), size=1)+
scale_colour_manual(values=cols2)+
labs(x="Species abundance rank", y="Colony density (2005)")+
geom_text_repel(data=sad[!sad$spp=="Other",], aes(rank, y=(abundance/2700)+0.2, label=spp), angle=45, size=2, force=0.001, hjust=0)+
xlim(c(0,70))+
theme_classic()+
guides(colour="none"),
ncol=1, align="v", axis="lr", rel_heights=c(0.3,1))

sadplot


fig.s2 <- plot_grid(triplot, plot_grid(sadplot, plot_grid(covplot,NULL, rel_widths=c(1,0.2)), ncol=1, rel_heights=c(0.7,1), labels=c("B","C"), label_size=12), nrow=1, rel_widths=c(1,1.2), labels=c("A",""), label_size=12)
fig.s2

