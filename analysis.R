
rm(list = ls())
library("ggrepel")
library("grid")
library("png")
library("MASS") 
library("reshape2")
library("ggplot2")
library("cowplot")
library("stats4")

#######################################
# LOAD DATA
#######################################
source("R/functions.R")
source("R/data_prep.R")

# remove tiny colonies (<3cm diameter)
diam <- 3
gdat <- subset(gdat, area > log10((pi*(diam/2)^2)/10000))
sdat <- subset(sdat, area > log10((pi*(diam/2)^2)/10000))
fec <- subset(fec, area > log10((pi*(diam/2)^2)/10000))

sdat<-sdat[!sdat$colony_id==287,]
gdat<-gdat[!gdat$colony_id==287,]

# ordering
#################################### 
order <- c("Ahy","Acy","Ain","Aro","Ana","Asp","Ami","Adi","Ahu", "Gre","Gpe")
params$spp <- factor(params$spp, levels=order)
labs <- params$species[order(params$spp)]
cols<-as.character(params$cols)
names(cols) <- params$spp
cols <- cols[order]

#######################################
# MODEL PARAMETERS
#######################################
source("R/params.R")

# Give Gre the growth parameters of Gpe
params[params$spp=="Gre", c("p.int", "p.slp","p.sig")] <- params[params$spp=="Gpe", c("p.int", "p.slp","p.sig")]
params[params$spp=="Gre", c("g.int", "g.slp","g.var")] <- params[params$spp=="Gpe", c("g.int", "g.slp","g.var")]

# source("figs/supp.fig1.R") # Fig S1 - demographic models
# fig.s1
# ggsave("figs/supp.fig1.jpg", fig.s1, width=23, height=11.5, units="cm")

#######################################
# MORPHOLOGICAL PAIRS
#######################################

# Corymbose is a group of 3
# Remove Ana for comparisons, because it declined
params2 <- params[!params$spp=="Ana",]
params2$morph <- params2$morphology
comp<-dcast(params2[!params2$spp=="Ana",], morph~abundance_pair, value.var="spp")
comp

colsC <- cols[names(cols) %in% comp$Common]
names(colsC)<-comp$morph[match(names(colsC), comp$Common)]

#######################################
# DEMOGRAPHIC SPACE
#######################################

# 1 - total fecundity
params$f.colony<-aggregate(fecundity~spp, fec,mean)$fecundity 
# 2 - growth (r.int)
# 3 - proportion mortality at 100cm2
params$p_mort<-inv.logit((params$p.slp*log10(0.01))+params$p.int)
# 4 - fecundity per area
params$f.cm2<-aggregate(f.cm2~spp, fec,mean)$f.cm2
# 5 - minimum at reproductive maturity
params$min.r<-1/aggregate(area_cm2~spp, fec[fec$reproductive==1,], min)$area_cm2
# 6 - recruit survival, rec.cm = 15
params$survcm<-aggregate(pred~spp, s.pred[s.pred$area<log10(pi*(15/100/2)^2),], FUN=mean)$pred
# 7 - average modelled survival
params$av.surv<-aggregate(pred~spp, s.pred,mean)$pred

# PCA of demographic parameters
colnames(params)
rownames(params) <- params$spp
traits <- c("r.int","f.cm2","survcm", "av.surv","p_mort", "min.r", "f.colony")
pca<-prcomp(params[,traits], scale=T, center=T)
biplot(pca)

# explained variation
exp<-round(c(summary(pca)[[1]][1]^2/sum(summary(pca)[[1]]^2),summary(pca)[[1]][2]^2/sum(summary(pca)[[1]]^2)),3)*100
exp

#######################################
# ABUNDANCE
#######################################

abun[,c("spp","morphology","abundance_05")] <- params[match(abun$species, params$species), c("spp","morphology","abundance_05")]

d08B[,c("spp","morphology","abundance_05")] <- params[match(d08B$Species, params$species), c("spp","morphology","abundance_05")]

# source("figs/supp.fig2.R") # Figure S2 - abundances	
# fig.s2
# ggsave("figs/supp.fig2.jpg", fig.s2, width=15, height=15, units="cm")

#######################################
# FIGURE 1
#######################################

source("figs/fig.1.R")
fig.1
#ggsave("figs/fig.1.jpg", fig.1, width=15, height=9.5, units="cm", dpi = 300)

#######################################
# IPM MESH AND BOUNDARIES
#######################################

# repeat function to generate mesh for each species
max.size <- 0.3 # max(ss$area)
n <- 100
mesh <- function(){
	min.size <- rec.size # 3.5
	b <- seq(min.size, max.size, length=n)
	h <- b[2] - b[1]
	b <- c(min(b)-(2*h), min(b)-h, b)
	y <- 0.5 * (b[1:n]+b[2:(n+1)])
	I <- y >= rec.size
	return(list(b=b,h=h,y=y,I=I, rec.size = rec.size))}
	
#######################################
# RECRUITMENT PARAMETERS
#######################################

# recruit size - based on max growth, mean within groups
agg <- aggregate(r.int~morphology, params, mean)
params$r.int2 <- agg$r.int[match(params$morphology, agg$morphology)]
params$rec.size <- log10(pi*(params$r.int2*0.85)^2)
#r.limit <- 12
#params$rec.size <- ifelse(params$rec.size>log10(pi*((r.limit/2)/100)^2),log10(pi*((r.limit/2)/100)^2), params$rec.size)
# alternatively: *0.8 recsize without limit.

# recruitment probability - same within families
params$rec <- ifelse(params$family=="Merulinidae",0.0003162278, 0.0011364637)

#######################################
# PLOT IPMS
#######################################

# source("figs/supp/supp.fig3.R")
# fig.s3
# ggsave("figs/supp.fig3.jpg", fig.s3, width=20, height=10.2, units="cm")

#######################################
# ESTIMATE LAM
#######################################

lam.est <- NULL
for (sp in spp) {
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
lam.est<-rbind(lam.est, data.frame(sp, lam=bigmatrix()$lam))
}
params$lam.est <- lam.est$lam

ggplot(params, aes(reorder(spp, -lam.est), lam.est))+geom_bar(stat="identity", aes(fill=spp))+scale_fill_manual(values=cols)+geom_hline(yintercept=1)

params[,c("spp","lam.est", "g.slp")]

#######################################
# CHANGE IN ABUNDANCE (TRANSECTS)
#######################################

head(abun)
diffs <- data.frame(dcast(abun, species+tran~year, value.var="N"))
diffs$diff <- (diffs$X2014-diffs$X2011)
diffs$spp<-params$spp[match(diffs$species, params$species)]
diffs$lam <- params$lam.est[match(diffs$spp, params$spp)]

ggplot(diffs, aes(reorder(spp, -lam), diff, fill=spp))+
stat_summary(geom="bar", fun="mean", col="black", size=0.1)+
scale_fill_manual(values=cols)+guides(fill="none")

d.mean<-aggregate(X2011~species, diffs, mean)
d.mean$X2014 <- aggregate(X2014~species, diffs, mean)$X2014
d.mean$X2014 <- ifelse(d.mean$X2014==0, 0.05, d.mean$X2014)
d.mean$X2011 <- ifelse(d.mean$X2011==0, 0.05, d.mean$X2011)
# give absent species a small abundance (0.1)^^
d.mean$lam.tran <- (d.mean$X2014/d.mean$X2011)^(1/3)
d.mean

params$lam.tran <- d.mean$lam.tran[match(params$species, d.mean$species)]

ggplot(params, aes(lam.est, lam.tran))+geom_text(aes(label=spp))+geom_abline(slope=1)+geom_smooth(method="lm", se=F, linetype="dotted")+scale_y_log10()+scale_x_log10()

summary(lm(log(lam.tran)~log(lam.est), params))

# source("figs/supp.fig4.R")
# fig.s4
# ggsave("figs/supp.fig4.jpg", fig.s4, width=15, height=13, units="cm")


#######################################
# RECRUITMENT SENSITIVITY I
#######################################	
rec.xDET <- 10^(seq(-5,-1, 0.1))

storeDET<-data.frame()
for (sp in spp) {
	for (rec in rec.xDET){
	  rec.size <- params$rec.size[params$spp==sp]
	  h <- mesh()$h
	  y <- mesh()$y
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  storeDET<-rbind(storeDET, data.frame(spp=sp, rec=rec, lam=bigmatrix()$lam))
   } 
  }

ggplot(storeDET, aes(rec, lam))+geom_line(aes(col=spp))+scale_x_log10()+scale_colour_manual(values=cols)+coord_cartesian(ylim=c(0.6,2))+guides(col="none")

#######################################
# REC AT TRANSECT LAMS
#######################################	

rec.transect <- NULL
for(sp in spp){
	#sp <- "Gre"
sub <- storeDET[storeDET$spp==sp,]
r.sp <- params$lam.tran[params$spp==sp]
n.rec <- which(abs(sub$lam-r.sp)==min(abs(sub$lam-r.sp)))
if(n.rec==1){
newlam<-	min(sub$lam)+0.01
n.rec2 <- which(abs(sub$lam-newlam)==min(abs(sub$lam-newlam)))
new <- sub[n.rec2,]
} else { new <- sub[n.rec,] }
rec.transect <- rbind(rec.transect, cbind(new, lam.r=r.sp))
}

rec.transect

params$rec.tran <- rec.transect$rec[match(params$spp, rec.transect$spp)]

# rec at 1... 
rec.one <- NULL
for(sp in spp){
sub <- storeDET[storeDET$spp==sp,]
n.rec <- which(abs(sub$lam-1)==min(abs(sub$lam-1)))
rec.one <- rbind(rec.one, cbind(sub[n.rec,], lam.r=r.sp))
}
rec.one
params$rec.one <- rec.one$rec[match(params$spp, rec.one$spp)]

#######################################
# REC AT TRANSECT LAMS
#######################################	

agg.rec <- aggregate(log10(rec.tran)~family, params, mean)
params$rec.mean <- 10^agg.rec[match(params$family, agg.rec$family),"log10(rec.tran)"]
#params$rec.mean <- 10^mean(log10(params$rec.tran))
#params$rec.mean <- median(params$rec.tran)

r.long <- melt(params[c("spp","rec.mean", "rec.tran")])
r.long

#r.long<-r.long[!(r.long$spp=="Gpe" & r.long$variable=="rec.tran"),] 


plot_grid(
ggplot(r.long, aes(variable, value))+
geom_path(aes(group=spp), size=0.1)+
#geom_jitter(aes(col=spp), height=0.05, width=0.05)+
geom_point(aes(col=spp))+
geom_point(data=params, aes(x="rec.mean",rec.mean))+
geom_text(data=params[params$spp %in% c("Gre","Ahy"),], aes(x="rec.mean",rec.mean, label=family), nudge_y=0.1, nudge_x=0.2, angle=15, hjust=0, size=2)+
scale_y_log10(limits=c(min(storeDET$rec),max(storeDET$rec)))+
scale_colour_manual(values=cols)+
#geom_path(aes(col=spp, group=spp))+
coord_flip()+
guides(col="none"), 
ggplot(storeDET, aes(rec, lam))+geom_line(aes(col=spp))+scale_x_log10(limits=c(min(storeDET$rec),max(storeDET$rec)))+
geom_hline(yintercept=1)+
#guides(col="none")+
scale_colour_manual(values=cols)+coord_cartesian(ylim=c(0.6,2)),
ncol=1, rel_heights=c(0.5,1), align="v", axis="lr")

#######################################
# ESTIMATE LAM
#######################################

lam.est <- NULL
for (sp in spp) {
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec.mean[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
lam.est<-rbind(lam.est, data.frame(sp, lam=bigmatrix()$lam))
}
params$lam.est2 <- lam.est$lam

#ggplot()

ggplot(params, aes(reorder(spp, -lam.est2), lam.est2))+geom_bar(stat="identity", aes(fill=spp))+scale_fill_manual(values=cols)+geom_hline(yintercept=1)

params[,c("spp","lam.est","lam.est2", "g.slp")]

params$lam.est <- params$lam.est2
  ######### NEW REC!!!!
params$rec <- params$rec.mean

#######################################
# BOOTSTRAP
#######################################	
gdat_orig <- gdat
sdat_orig <-sdat
fec_orig <-fec

#source("R/bootstap_lam.R")
boot <- read.csv("data/lam.range.csv")

ggplot(boot)+
geom_density(aes(x=log(lam), col=spp, fill=spp), alpha=0.5)+
geom_vline(xintercept=0, linetype="dotted")+
#geom_point(data=params, aes(lam.est, y=0), col="grey")+
#geom_point(data=params, aes(x=lam.est2, y=0))+
facet_wrap(~morphology, scale="free_y", ncol=1)+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
#xlim(c(0.5, 1.8))+
#scale_x_log10(limits=c(0.5, 2))+
theme_classic()

head(boot)
boot$GT <- log(boot$r)/log(boot$lam)
ggplot(boot, aes(GT, log(lam)))+
geom_hline(yintercept=0)+
geom_point(shape=21, aes(col=spp), alpha=0.2)+
coord_cartesian(xlim=c(min(boot$GT),100))+
scale_colour_manual(values=cols)+scale_x_log10()+
geom_point(data=params, aes(GT, log(lam.est)))+guides(col="none")+
theme_classic()


#######################################
# ELASTICITY & IPM MEASURES 
#######################################
 # eigen-things combined to get sensitivity/elasticity matrices.

s.list <- list()
eK.list <- list()
eR.list <- list()
eP.list <- list()
demovals <- NULL
sizevals <- NULL

for (sp in spp) {
	#sp <- "Ahy"
	rec <- params$rec[params$spp==sp]
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	K <- bigmatrix()$K
	P <- bigmatrix()$P
	R <- bigmatrix()$R
	#image(t(K^0.2))
	
	lam <- Re(eigen(K)$values[1]) # population growth
	w.eigen<-Re(eigen(K)$vectors[,1]) # right eigenvec
	stable.size <- w.eigen/sum(w.eigen) # sable size dist 
	v.eigen<-Re(eigen(t(K))$vectors[,1]) # left eigenvec
	repro.val <- v.eigen/v.eigen[1] # rel. reprodutive values 
	v.dot.w<-sum(stable.size*repro.val*h) # reproductive val * stable size
	sens<-outer(repro.val,stable.size,"*")/v.dot.w   # sensitivity matrix
	K.elas <- sens*(K/h)/lam *h^2 # elasticity matrices (h varies)
	P.elas<-(P/h)*sens/lam # survival elasticity matrix
	eP=sum(P.elas)*h^2 # total survival elasticity
	R.elas<-(R/h)*sens/lam # reproduction elasticity matrix
	eR=sum(R.elas)*h^2 # total reproduction elasticity
	
	# Net reproductive rate/Generation time from IPMbook monocarp
	N <- solve(diag(n)-P)
	R0 <- abs(eigen(R %*% N)$values[1])
	GT <- log(R0)/log(lam)
	
	demovals <- rbind(demovals, data.frame(spp=sp, eR, eP, R0, GT))
	sizevals <- rbind(sizevals, data.frame(spp=sp, area=y, stable.size, repro.val, v.dot.w))
		s.list[[sp]] <- sens
		eR.list[[sp]] <- R.elas
		eP.list[[sp]] <- P.elas
		eK.list[[sp]] <- K.elas
		}

params[,colnames(demovals)]<- demovals[match(demovals$spp, params$spp),]
lapply(eK.list, function(x){sum(x)})
rowSums(demovals[,c("eP","eR")]) # summing to 1

plot_grid(
ggplot(sizevals, aes(x=area, y=repro.val, col=spp))+geom_line()+
scale_colour_manual(values=cols)+
scale_y_log10()+guides(col="none"),
ggplot(sizevals, aes(x=area, y=stable.size, col=spp))+geom_line()+
scale_colour_manual(values=cols)+guides(col="none"),
ggplot(params, aes(reorder(spp, -eR), eR, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none"), 
ggplot(params, aes(reorder(spp, -R0), R0, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none")+scale_y_log10(),
ggplot(params, aes(reorder(spp, -GT), GT, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none")+scale_y_log10(),
ggplot(params, aes(r.int, eR))+geom_text(aes(label=spp))+scale_y_sqrt()+scale_x_log10())

#######################################
#GENERATION TIME
#######################################	

sds<-aggregate(lam~spp, boot, function(x){ quantile(x, 0.95) })
params$lam.sd1 <- (sds$lam[match(params$spp, sds$spp)])
sds<-aggregate(lam~spp, boot, function(x){ quantile(x, 0.05) })
params$lam.sd2 <- (sds$lam[match(params$spp, sds$spp)])
sds<-aggregate(lam~spp, boot, median)
params$lam.mn <- (sds$lam[match(params$spp, sds$spp)])

sds<-aggregate(GT~spp, boot, function(x){ quantile(x, 0.95) })
params$GT.sd1 <- (sds$GT[match(params$spp, sds$spp)])
sds<-aggregate(GT~spp, boot, function(x){ quantile(x, 0.05) })
params$GT.sd2 <- (sds$GT[match(params$spp, sds$spp)])
sds<-aggregate(GT~spp, boot, median)
params$GT.mn <- (sds$GT[match(params$spp, sds$spp)])

ggplot()+geom_hline(yintercept=1, size=0.1)+
geom_segment(data=params, aes(x=GT, xend=GT, y=lam.sd1, yend=lam.sd2), size=0.2)+
geom_point(data=params,aes(GT, lam.mn, fill=spp), shape=21, stroke=0.2, size=3)+
geom_text(data=params, aes(GT, lam.mn, label=AB), size=1.8)+
#scale_x_log10()+#scale_y_log10()+
scale_fill_manual(values=cols)+guides(fill="none")+
theme_classic()

#######################################
# FIG 2
#######################################	
 
source("figs/fig.2.R") # check ordering
fig.2
#ggsave("figs/fig.2.jpg", fig.2, width=15, height=12.7, units="cm", dpi = 300)

#######################################
# COMPARE LAMS
#######################################

params2$lam.est <- params$lam.est[match(params2$spp, params$spp)]

morph.diff <- function(x){
	#x <-"m.int"
	params2$x.temp <-params[match(params2$spp, params$spp), x]
	comp$x.common <- params2$x.temp[match(comp$Common, params2$spp)]
	comp$x.rare <- params2$x.temp[match(comp$Rare, params2$spp)]
	comp$x.diff <- comp$x.common-comp$x.rare
	comp$x.div <- comp$x.common/comp$x.rare
	comp
	}

comp$AC<-morph.diff("abundance_05")$x.common
comp$AR<-morph.diff("abundance_05")$x.rare
comp$abundiff <- morph.diff("abundance_05")$x.div

comp$lamC <- morph.diff("lam.est")$x.common
comp$lamR <- morph.diff("lam.est")$x.rare
comp$lamdiff <- morph.diff("lam.est")$x.diff
comp$logdiff <- log(comp$lamC) - log(comp$lamR)

# Double recruitment
doub_lam <- NULL
for(sp in spp){
	rec <- params$rec[params$spp==sp] * 2
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	doub_lam <- rbind(doub_lam, data.frame(spp=sp, lam=bigmatrix()$lam))
	}
comp$doubC <- doub_lam$lam[match(comp$Common, doub_lam$spp)]
comp$doubR <- doub_lam$lam[match(comp$Rare, doub_lam$spp)]
comp$doubdiff <- comp$doubC - comp$doubR

ggplot(comp)+
geom_bar(aes(y=reorder(morph, -lamdiff), x=lamdiff), stat="identity")+
geom_segment(aes(y=reorder(morph, -lamdiff), yend=reorder(morph, -lamdiff), x=lamdiff, xend=doubdiff), arrow=arrow(length=unit(1,"mm")))

#######################################
# TIME NEEDED TO GET DIFFERENCES
#######################################

# method 1 

comp$doub.time <- log(2)/log(1+comp$lamdiff)
comp$difftime <- log(comp$abundiff)/log(1+comp$lamdiff)
		
ggplot(comp[comp$difftime>0,], aes(x=difftime, y=reorder(morph, -difftime)))+
geom_bar(stat="Identity", aes(fill=morph))+
labs(x="Years to project \nabundance differences")+
guides(fill="none")+scale_fill_manual(values=colsC)

# method 2 - projecting lamda

n1 <- 100
ngen <- 100
pars <- params2[, c("spp","morph", "abundance_pair","lam.est")]
proj <- merge(pars, data.frame(gen=1:ngen))
proj$n <- n1*proj$lam.est^proj$gen

ggplot(proj, aes(gen, n))+geom_line(aes(group=spp, col=spp))+scale_colour_manual(values=cols)+facet_wrap(~morph, scales="free_y")

ssq <- NULL
for(i in seq(1, 50, 1)){
Gen5 <- subset(proj, gen==i )
params2$Gen5<-scale(Gen5$n[match(params2$spp, Gen5$spp)])
mod <- aov(Gen5~morphology, data=params2)
summary(mod)
within <- summary(mod)[[1]]["Residuals","Sum Sq"]
between <- summary(mod)[[1]]["morphology","Sum Sq"]
ssq <- rbind(ssq, data.frame(i, within1=within, between=between, within2=(within/(within+between))))
}

plot_grid(ggplot()+geom_point(data=ssq, aes(i, within1))+geom_point(data=ssq, aes(i, between), col="red")+scale_y_log10(),
ggplot(ssq, aes(i, within2))+geom_point() )

proj2 <- dcast(morph+gen~abundance_pair, data=proj, value.var="n")
proj2$diff <- proj2$Common/proj2$Rare
proj2$real <- comp$abundiff[match(proj2$morph, comp$morph)]
gens <- aggregate(gen~morph, proj2[proj2$diff<proj2$real,], max)
proj2$maxgen <- gens$gen[match(proj2$morph, gens$morph)] +1 #to cross the line
head(proj2)

ggplot(proj2, aes(gen, diff, col=morph))+geom_line(linetype="dotted")+coord_cartesian(ylim=c(0,300))+geom_line(data=proj2[proj2$gen<=proj2$maxgen,])+geom_point(data=proj2[proj2$gen==proj2$maxgen,])
	
# method 3 - projecting lamda from size dist.. 

proj.size <- NULL
for(sp in spp){
	#sp <-"Ahu"
	rec <- params$rec[params$spp==sp] 
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	lam=bigmatrix()$lam
	v=bigmatrix()$v
	w=bigmatrix()$w
	ngen <- 200
	npop <- 100
	v <- matrix(c(1, rep(0,99)))
pop  <- matrix(nrow=npop, ncol=ngen) # size classes by gens.
pop[,1] <- round(v * npop) # inital population size
#plot(y, pop[,1])
#for(j in 2:ngen){ pop[,j] <- pop[,j-1]*lam }
for(j in 2:ngen){ pop[,j] <- bigmatrix()$K %*% as.matrix(pop[,j-1]) }
pop.long <- melt(pop, value.name = "NperSize")
colnames(pop.long)[1]<-"Size"
colnames(pop.long)[2]<-"gen"
head(pop.long)
ggplot(data=pop.long[pop.long$gen %in% c(1,5,10,20),], aes(x=Size, y=NperSize,  group=Size))+geom_point(shape=21)+facet_wrap(~gen, scales="free_y")
psize<-	aggregate(NperSize~gen, data=pop.long, sum)
psize$n <- npop*lam^psize$gen # compare with simple method...
proj.size <- rbind(proj.size, cbind(psize, spp=sp, abun=params$abundance_pair[params$spp==sp], morph=params$morph[params$spp==sp]))
}

head(proj.size)
	ggplot()+
	geom_point(data=proj.size, aes(gen, NperSize, col=spp))+
	scale_y_log10()+scale_colour_manual(values=cols)+
	facet_wrap(~morph)
	
ssq <- NULL
for(i in seq(2, 50, 1)){
Gen5 <- subset(proj.size, gen==i )
params$Gen5<-scale(Gen5$NperSize[match(params$spp, Gen5$spp)])
mod <- aov(Gen5~morphology, data=params)
summary(mod)
within <- summary(mod)[[1]]["Residuals","Sum Sq"]
between <- summary(mod)[[1]]["morphology","Sum Sq"]
ssq <- rbind(ssq, data.frame(i, within1=within, between=between, within2=(within/(within+between))))
}

plot_grid(ggplot()+geom_point(data=ssq, aes(i, within1))+geom_point(data=ssq, aes(i, between), col="red")+scale_y_log10(),
ggplot(ssq, aes(i, within2))+geom_point() )

	
	
	
	
	hist(proj3$diff)

proj.size2 <- proj.size[!proj.size$spp=="Ana",]	
#proj.size$morph <- ifelse(proj.size$spp=="Ana", "corymbose_2",as.character(proj.size$morph))
#cory2 <- proj.size[proj.size$spp=="Ami",]
#cory2$morph <- "corymbose_2"
#proj.size2 <- rbind(proj.size, cory2)

proj3 <- dcast(morph+gen~abun, data=proj.size2, value.var="NperSize")
proj3$diff <- proj3$Common/proj3$Rare
proj3$real <- comp$abundiff[match(proj3$morph, comp$morph)]
gens <- aggregate(gen~morph, proj3[proj3$diff<proj3$real,], max)
proj3$maxgen <- gens$gen[match(proj3$morph, gens$morph)] +1 
head(proj3)	
gens 
tail(proj3)

aggregate(diff~morph, proj3, min)

projdat <- proj3[proj3$gen==proj3$maxgen,]
projdat

	
ggplot(proj3, aes(gen, diff, col=morph))+
geom_line(linetype="dotted")+
geom_line(data=proj3[proj3$gen<=proj3$maxgen,])+
geom_point(data=proj3[proj3$gen==proj3$maxgen,])+ 
#scale_y_log10()+#scale_x_log10()+
scale_y_sqrt()+scale_x_sqrt()+
coord_cartesian(ylim=c(0,60), xlim=c(0,100))

#######################################
# WITHIN VS BETWEEN
#######################################

Gen5 <- subset(proj.size, gen==5 )
params$Gen5<-Gen5$NperSize[match(params$spp, Gen5$spp)]

Gen10 <- subset(proj.size, gen==10)
params$Gen10<-Gen10$NperSize[match(params$spp, Gen10$spp)]

dems <- c("lam.est", "Gen5",  "Gen10", "abundance_05", "rec")

vars <- NULL
for (t in c(traits, "eggC",dems)){
	#t <- "Gen15"
	params$t <- scale(params[,t])
mod <- aov(t~morphology, data=params)
summary(mod)
within <- summary(mod)[[1]]["Residuals","Sum Sq"]
between <- summary(mod)[[1]]["morphology","Sum Sq"]
pval <- summary(mod)[[1]]["morphology","Pr(>F)"]
summary(mod)
vars <- rbind(vars, data.frame(t, within=100*(within/(within+between)), between=100*(between/(within+between)), pval))
}
ggplot(vars, aes(within,reorder(t, -within)))+geom_bar(stat="identity")

#######################################
# SENSITIVITY ANALYSIS!!  Within morph groups
#######################################

params_orig <- params # CAREFUL set original
	
source("R/params_morph.R")
head(p.morph)

pars <- c("m.int", "m.slp","f.int","f.slp","g.int","g.slp","g.var", "s.int", "s.slp", "s.slp.2")

p.types <- data.frame(pars, type = c("reproduction","reproduction","reproduction","reproduction","growth","growth","growth", "survival","survival","survival"))

types <- unique(p.types$type)

test <- NULL
test2 <- NULL
for (i in types) {
#i <- "m.int"
#i = "maturity"
params <- params_orig 
params[,pars] <- p.morph[match(params$morphology, p.morph$morph), pars]
pars2 <- p.types$pars[p.types$type==i]
params[,pars2] <- params_orig[match(params$spp, params_orig$spp), pars2]
params
for (sp in spp) {
	#sp <- "Ahy"
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
test<-rbind(test, data.frame(spp=sp, param=i, lam=bigmatrix()$lam, morph=params$morphology[params$spp==sp], rec="normal"))
rec <- params$rec.mean[params$spp==sp] * 2
test2<-rbind(test2, data.frame(spp=sp, param=i, lam=bigmatrix()$lam, morph=params$morphology[params$spp==sp], rec="double"))
 } 
}

diffs <- NULL
diffs2 <- NULL
for (i in types) {
	#i <- "maturity"
	params <- params_orig 
	sub <- test[test$param==i, ]
	params$x <- sub$lam[match(params$spp, sub$spp)]
	diffs <- rbind(diffs, data.frame(morph=comp$morph, d=(morph.diff("x")$x.common - morph.diff("x")$x.rare), param=i, rec="normal"))
	sub2 <- test2[test2$param==i, ]
	params$x <- sub2$lam[match(params$spp, sub2$spp)]
	diffs2 <- rbind(diffs2, data.frame(morph=comp$morph, d=(morph.diff("x")$x.common - morph.diff("x")$x.rare), param=i, rec="double"))
	}
head(diffs)


sums <- aggregate(d~morph+rec, diffs, sum)
sums$lamdiff <- comp$lamdiff[match(sums$morph, comp$morp)]

sums2 <- aggregate(d~morph+rec, diffs2, sum)
sums2$lamdiff <- comp$doubdiff[match(sums2$morph, comp$morp)]


ggplot(rbind(sums,sums), aes(lamdiff, d))+geom_point(aes(col=morph, shape=rec))+
scale_y_log10()+scale_x_log10()+
#scale_y_sqrt()+scale_x_sqrt()+
geom_abline(slope=1)+
labs(x="Real lamda difference", y = "Sum of differences")

ggplot()+
geom_bar(data=diffs, aes(d, reorder(morph, -d), fill=param), stat="identity", position="stack", col="black")+
geom_point(data=comp, aes(lamdiff, morph, col=morph), shape=4, stroke=1, size=0.5)+
#scale_fill_manual(values=c("white",pal[c(2:3)]))+
labs(x="Lamda difference (Common - Rare)")+
geom_vline(xintercept=0)+
#scale_x_sqrt()+
guides(col="none")+scale_colour_manual(values=colsC)+
theme_classic()+theme(axis.title.y=element_blank(), legend.title=element_blank(), legend.position=c(0.8, 0.8), legend.key.width=unit(2,"mm"), legend.key.height=unit(1,"mm"))


#######################################
# FIGURE 3
#######################################

source("figs/fig.3.R")
fig.3
#ggsave("figs/fig.3.jpg", fig.3, width=14, height=12, units="cm", dpi = 250)

#######################################
# REPRODUCTIVE TRADEOFF
#######################################

params <- params_orig 

ec$morph<-params$morphology[match(ec$spp, params$spp)]
ecplot <- ggplot()+
geom_boxplot(data=ec, aes(y=reorder(morph, -Carbon_ug_corrected), x=Carbon_ug_corrected, fill=spp), outlier.size=0.2, size=0.1,position = position_dodge(width = 1))+scale_fill_manual(values=cols)+guides(fill="none")

params$fec1cm <- exp(params$f.int +  params$f.slp * log10(1)) /10000
params2$fec1cm <- params$fec1cm[match(params2$spp, params$spp)]

plot_grid(ecplot,ggplot(params, aes(fec1cm, eggC))+geom_text(aes(label=spp))+geom_smooth(data=params[params$family=="Acroporidae",], method="lm"), nrow=1, rel_widths=c(1.15,1))


#######################################
# FIGUREEEEEEEEE
#######################################

source("figs/fig.4.R")
fig.4
#ggsave("figs/fig.4.jpg", fig.4, width=5.5, height=8, units="cm", dpi = 250)

#######################################
# RECRUITMENT SENSITIVITY II
#######################################	


# lambda at different rec/recsize combinations
rec.x <- 10^(seq(-6,-2,0.5))
recsize.x <- seq(min(params$rec.size)-0.1, max(params$rec.size)+0.1,0.1)
3
store3 <- list()
for (sp in spp) {
	temp<-matrix(NA,length(rec.x), length(recsize.x))
	for (i in 1:length(rec.x)){
		rec <- rec.x[i]
		for(j in 1:length(recsize.x)){
			rec.size <- recsize.x[j]
		h <- mesh()$h
		y <- mesh()$y	
   mod <- bigmatrix()
   temp[i,j] <- mod$lam
  }
 } 
 store3[[sp]] <- temp  
  }

minlam <- min(log(do.call("rbind", store3)))
maxlam <- max(log(do.call("rbind", store3)))

par(mfrow=c(3,4))
sp.conts <- NULL
for(sp in spp){
	#sp<-"Ahy"
  image(log10(rec.x), recsize.x, log(store3[[sp]]), main=sp, zlim=c(minlam, maxlam))
  abline(h=params[params$spp==sp, "rec.size"], lty=2)
  contour(log10(rec.x), recsize.x, log(store3[[sp]]), add=TRUE, levels=0)
  }
   contour(log10(rec.x), recsize.x, log(store3[[spp[1]]]), levels=0, col=cols[1], lwd=3)
 for(i in 1:length(spp)){
 	sp <- spp[i]
  contour(log10(rec.x), recsize.x, log(store3[[sp]]), add=TRUE, levels=0, col=cols[i], lwd=3)
  conts <- melt(store3[[spp[i]]])
  sp.conts <- rbind(sp.conts, cbind(conts,rec = rep(rec.x, length(recsize.x)), recsize = rep(recsize.x, each=length(rec.x)), spp=sp))
  }
head(sp.conts)

##### PLOT

params2$rec <- params$rec[match(params2$spp, params$spp)]
params2$rec.size <- params$rec.size[match(params2$spp, params$spp)]

comp$rec <- params2$rec[match(comp$morph, params2$morph)] 
comp$rec.size <- params2$rec.size[match(comp$morph, params2$morph)] 
comp$rec.size.cm <- sqrt(((10^comp$rec.size)*10000)/pi)*2	

# difference in lambda across rec/recsize combinations
pairs3 <- NULL
for(m in comp$morph){
	#m <- "tabular"
	spC <- comp$Common[comp$morph==m]
	spR <- comp$Rare[comp$morph==m]
	logdiff <- log(store3[[spC]]) - log(store3[[spR]])
	temp <- melt(logdiff)
	temp$rec <- rep(rec.x, length(recsize.x))
	temp$recsize <- rep(recsize.x, each=length(rec.x))
    pairs3 <- rbind(pairs3, cbind(temp, morph=m))
}
head(pairs3)

pairs3$lam2 <- ifelse(pairs3$value < 0, NA, pairs3$value)

ggplot()+
geom_raster(data=pairs3, aes(x=rec, y=recsize, fill=lam2))+
#geom_point(data=params2, aes(rec.const, rec.size.const ), col="white", shape=3, stroke=1, size=0.3)+
geom_point(data=params2, aes(rec, rec.size ), col="white", shape=3)+
	scale_y_continuous(expand=c(0,0))+
	scale_x_log10(expand=c(0,0))+
	scale_fill_distiller(palette="Spectral")+
	facet_wrap(~morph, scales="free")+
	labs(x="rec rate", y="rec size")+
	theme(strip.background=element_blank())


# source("figs/supp.fig5.R")	# Figure S5
#fig.s5
#ggsave("figs/supp.fig5.jpg", fig.s5, width=9, height=8, units="cm")
