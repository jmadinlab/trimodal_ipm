
vars$type <- ifelse(vars$t %in% dems, "Abundance", ifelse(vars$t=="abundance_05", "Abundance", "Demographic traits"))
#vars$type <- factor(vars$type, levels=c("Demographic traits", "Abundance", "Abundance"))

abunvar <- subset(vars, t=="abundance_05")
vars2 <- vars[!vars$t=="abundance_05",]


vars2$t[order(vars2$within)]
group_name <- c(
expression(italic(P[rec])),
expression(Survival[juvenile]), 
expression(Partial~mortality),
expression(Growth[max]),  
expression(Fecundity[colony]),
expression(Mature~size),
expression(Survival[adult]),
"Egg Mass",
expression(Fecundity[area]), 
expression(Fitness~(lambda)),
expression(Abundance["5yrs"]),
expression(Abundance["10yrs"])
)

#pal[1:3]
library("RColorBrewer")
pal <- brewer.pal(n = 3, name = 'Greys')
#https://www.schemecolor.com/leather.php

withinplot <- ggplot(vars2, aes(y=reorder(t, -within), x=within))+
#geom_point(aes(col=type))+
geom_bar(stat="Identity", col="black", size=0.1, width=0.75, aes(fill=type))+
geom_bar(stat="Identity", col="black", size=0.1, width=0.75,fill="grey", alpha=0.5)+
scale_x_continuous(expand=c(0,0), limits=c(0,80))+
scale_y_discrete(labels=rev(group_name))+
labs(x="% variation within \nmorphological groups", y="parameter")+
guides(fill="none")+
#coord_cartesian(xlim=c(0,55))+
theme_classic()+
geom_point(data=abunvar, aes(y="Gen10", x=within), shape=4, stroke=0.7)+
geom_segment(data=NULL, aes(y="f.cm2", yend="survcm", x=20, xend=20))+
geom_segment(data=NULL, aes(y="f.cm2", yend="f.cm2", x=20, xend=19))+
geom_segment(data=NULL, aes(y="survcm", yend="survcm", x=20, xend=19))+
geom_text(data=NULL, aes(y="min.r", x=22, label="Demographic \ntraits"), size=2.5, hjust=0)+
geom_text(data=NULL, aes(y="lam.est", x=80, label="Abundances\nin 2005"), size=2.5, hjust=1)+
geom_segment(data=NULL, aes(y="Gen5", yend=1.5, x=77.3, xend=77.3))+
scale_fill_manual(values=c("#9F7159","#DF9D6C"))+ #pal[c(2:3)]
theme(axis.title.y=element_blank(), legend.title=element_blank(), legend.key.width=unit(2, "mm"), legend.key.height=unit(2, "mm"), legend.position=c(0.6, 0.9), 
legend.background=element_blank(),
panel.grid.major.x=element_line(),
panel.grid.minor.x=element_line(),
axis.line=element_line(size=0.1),
legend.text=element_text(size=8),
axis.title.x=element_text(size=8), axis.text.y=element_text(size=8))
withinplot


Fig4AB <- plot_grid(projplot+guides(col="none")+labs(x="Time (years)", y=expression(N[common]/N[rare]))+ggtitle("Projected abundance \ndifferences")+theme(plot.title=element_text(size=8, hjust=0.5, face="bold")), withinplot+ggtitle("Disassociation from \nmorphology")+theme(plot.title=element_text(size=8, hjust=0.5, face="bold")), rel_widths=c(0.8, 1), labels=c("A","B"), label_size=9)
	Fig4AB
