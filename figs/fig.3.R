library("RColorBrewer")

plotdat <- proj

projplot <- ggplot(plotdat, aes(gen, diff, col=morph))+
geom_line(linetype="dotted")+
#scale_x_log10()+scale_y_log10()+
geom_line(data=plotdat[plotdat$gen<=plotdat$maxgen,])+
geom_point(data=plotdat[plotdat$gen==plotdat$maxgen,])+
coord_cartesian(ylim=c(0.5,65), xlim=c(0.7,45))+
scale_y_sqrt(expand=c(0,0))+
scale_x_sqrt(expand=c(0,0))+
scale_colour_manual(values=colsC)+
theme_classic()+theme(legend.position=c(0.75, 0.2), legend.title=element_blank(), legend.background=element_blank(), legend.key.height=unit(1,"mm"), legend.text=element_text(size=7), axis.title=element_text(size=9))
projplot

####### PANEL B

vars$type <- ifelse(vars$t %in% dems, "Abundance", ifelse(vars$t=="abundance_05", "Abundance", "Demographic traits"))

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
expression(Pop.~growth~(lambda)),
expression(Abundance["5yrs"]),
expression(Abundance["10yrs"])
)


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
scale_fill_manual(values=c("#9F7159","#DF9D6C"))+ 
theme(axis.title.y=element_blank(), legend.title=element_blank(), legend.key.width=unit(2, "mm"), legend.key.height=unit(2, "mm"), legend.position=c(0.6, 0.9), 
legend.background=element_blank(),
panel.grid.major.x=element_line(),
panel.grid.minor.x=element_line(),
axis.line=element_line(size=0.1),
legend.text=element_text(size=8),
axis.title.x=element_text(size=8), axis.text.y=element_text(size=8))
withinplot


Fig4AB <- plot_grid(projplot+labs(x="Time (years)", y=expression(N[common]/N[rare]))+ggtitle("Projected abundance \ndifferences")+theme(plot.title=element_text(size=8, hjust=0.5, face="bold")), withinplot+ggtitle("Disassociation from \nmorphology")+theme(plot.title=element_text(size=8, hjust=0.5, face="bold")), rel_widths=c(0.8, 1), labels=c("a","b"), label_size=9)
	Fig4AB

# Panel C-D

pal <- brewer.pal(n = 3, name = 'Greys')

comp$labels <- paste(comp$Common, " - ", comp$Rare, sep="")

mas_x <- comp$lamdiff[comp$morph=="massive"]
tab_x <- comp$lamdiff[comp$morph=="tabular"]
brn_x <- comp$lamdiff[comp$morph=="staghorn"]
dig_x <- comp$lamdiff[comp$morph=="digitate"]
cor_x <- comp$lamdiff[comp$morph=="corymbose"]

comp$rank<-NA
comp$rank[order(comp$lamdiff)] <- nrow(comp):1
mas_y <- comp$rank[comp$morph=="massive"]
tab_y <- comp$rank[comp$morph=="tabular"]
brn_y <- comp$rank[comp$morph=="staghorn"]
dig_y <- comp$rank[comp$morph=="digitate"]
cor_y <- comp$rank[comp$morph=="corymbose"]
cor2_y <- comp$rank[comp$morph=="corymbose_2"]

m.diff$label <- comp$label[match(m.diff$morph, comp$morph)]


fig4CD <- plot_grid(ggplot(m.diff, aes(d, reorder(label, -d), fill=param))+
geom_bar(stat="Identity", position="stack", col="black", size=0.2, width=0.6)+
scale_fill_manual(values=c("white",pal[c(2:3)]))+
labs(x=expression(lambda~difference~(Common~-~Rare)))+
xlim(c(-0.02, 0.6))+
annotation_custom(mas, mas_x+0.01, mas_x+0.13, mas_y-0.4, mas_y+0.4)+
annotation_custom(brn, brn_x+0.05, brn_x+0.14, brn_y-0.5, brn_y+0.5)+
annotation_custom(tab, tab_x+0.03, tab_x+0.12, tab_y-0.5, tab_y+0.5)+
annotation_custom(cor, cor_x+0.01, cor_x+0.15, cor_y-0.5, cor_y+0.5)+
annotation_custom(dig, dig_x+0.01, dig_x+0.15, dig_y-0.5, dig_y+0.5)+
geom_vline(xintercept=0)+
theme_classic()+theme(axis.title.y=element_blank(), legend.title=element_blank(),  legend.key.width=unit(2,"mm"), legend.key.height=unit(1,"mm"), axis.line.y=element_blank(), legend.position=c(0.8,0.8),axis.text.x=element_text(size=9),axis.text.y=element_text(size=7, angle=30), axis.title=element_text(size=9)),
ggplot(sums, aes(lamdiff, d))+geom_abline(slope=1, size=0.1)+geom_point(aes(col=morph), size=2)+
labs(x=expression(~Original~lambda~difference), y = "Sum of \ndifferences")+guides(col="none")+scale_colour_manual(values=colsC)+theme_classic()+theme(axis.text=element_text(size=9),axis.title=element_text(size=9)), 
rel_widths=c(1,0.5), labels=c("c","d"), label_size=9)

fig.3<- plot_grid(plot_grid(NULL, Fig4AB, rel_widths=c(0.02,1)), NULL, fig4CD, ncol=1, rel_heights=c(1,0.12, 0.7))+
draw_label(expression(bold(Source~of~lambda~differences)), 0.55, 0.42, size=8)+
draw_label(expression(bold(within~groups)), 0.55, 0.39, size=8)+
draw_label("?", 0.63, 0.9, size=8, fontface="bold")
fig.3
