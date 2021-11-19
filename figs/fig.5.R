
vars$type <- ifelse(vars$t %in% dems, "Abundance", ifelse(vars$t=="abundance_05", "Abundance", "Demographic traits"))
#vars$type <- factor(vars$type, levels=c("Demographic traits", "Abundance", "Abundance"))

vars$t[order(vars$within)]
group_name <- c(
expression(Survival[juvenile]), 
expression(Growth[max]), 
expression(Partial~mortality), 
expression(Fecundity[colony]),
expression(Survival[adult]),
expression(Mature~size),
"Egg Mass",
expression(Fecundity[area]), 
expression(R[elasticity]),
expression(Fitness~(lambda[con])),
expression(Fitness~(lambda[rec])),
"Abundance"
)

#pal[1:3]
library("RColorBrewer")
pal <- brewer.pal(n = 3, name = 'Greys')

withinplot <- ggplot(vars, aes(y=reorder(t, -within), x=within))+
#geom_point(aes(col=type))+
geom_bar(stat="Identity", col="black", size=0.1, width=0.75, aes(fill=type))+
scale_x_continuous(expand=c(0,0), limits=c(0,80))+
scale_y_discrete(labels=rev(group_name))+
labs(x="% variation within \nmorphological groups", y="parameter")+
guides(fill="none")+
#coord_cartesian(xlim=c(0,55))+
theme_classic()+
geom_segment(data=NULL, aes(y="f.cm2", yend="survcm", x=20, xend=20))+
geom_segment(data=NULL, aes(y="f.cm2", yend="f.cm2", x=20, xend=19))+
geom_segment(data=NULL, aes(y="survcm", yend="survcm", x=20, xend=19))+
geom_text(data=NULL, aes(y="av.surv", x=22, label="Demographic \ntraits"), size=2.5, hjust=0)+
scale_fill_manual(values=pal[c(2:3)])+
theme(axis.title.y=element_blank(), legend.title=element_blank(), legend.key.width=unit(2, "mm"), legend.key.height=unit(2, "mm"), legend.position=c(0.6, 0.9), 
legend.background=element_blank(),
panel.grid.major.x=element_line(),
panel.grid.minor.x=element_line(),
axis.line=element_line(size=0.1),
legend.text=element_text(size=8),
axis.title.x=element_text(size=8), axis.text.y=element_text(size=8))
withinplot

fig.5 <- withinplot
