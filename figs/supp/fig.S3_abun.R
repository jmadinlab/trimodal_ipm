



triplot<-plot_grid(
ggplot(params, aes(reorder(spp, -abundance_05), abundance_05/2700, fill=spp))+geom_bar(stat="identity", width=0.8, col="black", size=0.1)+
guides(fill="none")+
labs(y=expression(Colonies~per~m^2))+
scale_fill_manual(values=cols)+
scale_y_continuous(expand=c(0,0))+
geom_hline(yintercept=0)+
theme_classic()+
labs(title="27 x BTs",subtitle="2004")+
theme(axis.title.x=element_blank(), 
axis.text.x=element_blank(),
plot.title=element_text(size=8, hjust=0.5, face="bold"),
plot.subtitle=element_text(size=9, hjust=0.5),strip.background=element_blank())
,
ggplot(data=tri2, aes(reorder(species, -abun05), N/10, fill=spp))+
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


covplot <- ggplot(data=tri2, aes(reorder(species, -abun05), cover, fill=spp))+
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

sad <- read.csv("data/dornelas2008.csv")

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
labs(x="Species abundance rank", y="Colony density (2004)")+
geom_text_repel(data=sad[!sad$spp=="Other",], aes(rank, y=(abundance/2700)+0.2, label=spp), angle=45, size=2, force=0.001, hjust=0)+
xlim(c(0,70))+
theme_classic()+
guides(colour="none"),
ncol=1, align="v", axis="lr", rel_heights=c(0.3,1))

sadplot






fig.s3 <- plot_grid(triplot, plot_grid(sadplot, plot_grid(covplot,NULL, rel_widths=c(1,0.2)), ncol=1, rel_heights=c(0.7,1)), nrow=1, rel_widths=c(1,1.2))
fig.s3


# OTHER ANALYSIS OF ABUNDANCE



#coral cover
###########################

non <- c("Sinul","Sarco","Helio", "Mille", "other", "anemo", "soft","Dead ","clam", "Palyt", "Rubbl", "Cladi", "Other","Xenia","Isis ", "zooan", "coral")
ab <- ab[!ab$genus %in% non, ]
unique(ab$genus)


cov <- aggregate(intercept_cm~ transect +site +campaign, ab, sum)
cov$cover <- cov$intercept_cm/(1000)*100

av_cov <- aggregate(cover~ site +campaign, cov, mean)
av_cov$sd<-aggregate(cover~ site +campaign, cov, sd)$cover

ggplot(av_cov, aes(campaign, cover))+geom_point()+geom_line(aes(group=site))+
geom_vline(xintercept=2015, linetype="dotted")+
geom_segment(aes(x=campaign, xend=campaign, y=cover-sd, yend=cover+sd))+facet_wrap(~site)

# genus composition
###########################

gen <-  aggregate(intercept_cm~ transect +site +campaign +genus, ab, sum)
gen$cover <- gen$intercept_cm/(1000)*100
head(gen)

av_gen <- aggregate(cover~ site +campaign+genus, gen, mean)
av_gen$sd<-aggregate(cover~ site +campaign+genus, gen, sd)$cover

ggplot(aggregate(cover~genus, av_gen, sum), aes(reorder(genus, -cover), cover))+geom_bar(stat="Identity")+
theme(axis.text.x=element_text(angle=55, hjust=1))

av_gen$genus2 <- ifelse(av_gen$genus=="Acrop", "Acropora", ifelse(av_gen$genus=="Porit", "Porites", ifelse(av_gen$genus=="Pocil", "Pocillopora", ifelse(av_gen$genus=="Monti", "Montipora", ifelse(av_gen$genus=="Gonia", "Goniastrea", ifelse(av_gen$genus=="Dipsa", "Dipsastrea", ifelse(av_gen$genus=="Platy", "Platygyra", ifelse(av_gen$genus=="Lobop", "Lobophyllia",ifelse(av_gen$genus=="Isopo", "Isopora",ifelse(av_gen$genus=="Diplo", "Diploastrea", ifelse(av_gen$genus=="Favit", "Favites", ifelse(av_gen$genus=="Echin", "Echinophyllia", "Other"))))))))))))

ggplot(av_gen, aes(campaign, cover))+
#geom_segment(aes(x=campaign, xend=campaign, y=cover-sd, yend=cover+sd))+
geom_point(size=0.1, col="grey")+geom_line(aes( group=site), col="grey")+
geom_vline(xintercept=2015, linetype="dotted")+facet_wrap(~genus2)+scale_y_log10()+stat_summary(geom="line", fun="mean")


# trimodal cover/comp
###########################
av_gen$genus3 <- ifelse(av_gen$genus=="Acrop", "Acropora", ifelse(av_gen$genus=="Porit", "Porites", ifelse(av_gen$genus=="Pocil", "Pocillopora", ifelse(av_gen$genus=="Monti", "Montipora", ifelse(av_gen$genus=="Gonia", "Goniastrea", ifelse(av_gen$genus=="Dipsa", "Dipsastrea", ifelse(av_gen$genus=="Platy", "Platygyra", ifelse(av_gen$genus=="Turbi", "Turbinarea",ifelse(av_gen$genus=="Isopo", "Isopora", ifelse(av_gen$genus=="Coelo", "Coeloseris", "Other Scleractinia"))))))))))
#av_gen$genus3 <- av_gen$genus

tri_cov <- subset(av_cov, site=="Trimodal")
tri_gen <- subset(av_gen, site=="Trimodal")

tri_gen$genus3[tri_gen$genus=="Coelo"]<-"Goniastrea"

#order <- subset(tri_gen, campaign==2011)
order <- aggregate(cover~genus3,tri_gen,mean)
ord <- order$genus3[order(order$cover, decreasing=T)]
#tri_gen$genus3 <- factor(tri_gen$genus3, levels=ord)
ord

# need 0s
tri_gen$ID <- paste(tri_gen$campaign, tri_gen$genus3)

all_combs <- data.frame(site ="Trimodal", campaign = rep(unique(tri_gen$campaign), length(unique(tri_gen$genus3))), genus3= rep(unique(tri_gen$genus3), each=length(unique(tri_gen$campaign))))
all_combs$ID <- paste(all_combs$campaign, all_combs$genus3)
all_combs$cover <- tri_gen$cover[match(all_combs$ID, tri_gen$ID)]
all_combs$cover[is.na(all_combs$cover)]<-0
all_combs$genus3 <- factor(all_combs$genus3, levels=ord)
unique(all_combs$genus3)
head(all_combs)

ggplot(all_combs, aes(campaign, cover, fill=genus3, group=genus3))+geom_area(position = 'stack', col="black", size=0.1)+
geom_vline(xintercept=2014+(4/12), linetype="dotted")+
geom_vline(xintercept=2015+(3/12), linetype="dotted")+
scale_y_continuous(expand=c(0,0))





 

