
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
sdat<-sdat[!sdat$colony_id==287,]
gdat<-gdat[!gdat$colony_id==287,]

#######################################
# MODEL PARAMETERS
#######################################
source("R/params.R")
fig.s1
#ggsave("figs/supp/fig.s1.png", fig.s1, width=25, height=13, units="cm", dpi = 300)

#######################################
# DEMOGRAPHIC SPACE
#######################################
# size
ggplot(ss, aes(x=area, col=spp))+
geom_density()+scale_colour_manual(values=cols)

# proportion mortality at 100cm2
params$p_mort<-inv.logit((params$p.slp*log10(0.01))+params$p.int)

# average modelled survival
params$av.surv<-aggregate(pred~spp, s.pred,mean)$pred

# recruit survival
rec.cm <- 11
params$survcm<-aggregate(pred~spp, s.pred[s.pred$area<log10(pi*(rec.cm/100/2)^2),], FUN=mean)$pred

# size at maximum survival 
max.surv<-do.call(rbind, lapply(split(s.pred, as.factor(s.pred$spp)), function(x) {return(x[which.max(x$pred),])}))
params$safe.size<-max.surv$area

# fecundity per area
params$f.cm2<-aggregate(f.cm2~spp, fec,mean)$f.cm2
#ggplot(params, aes(f.cm2, f.int))+geom_text(aes(label=spp))

# total fecundity
head(fec)
params$f.colony<-aggregate(fecundity~spp, fec,mean)$fecundity

# average size
size.av <- aggregate(area~spp, ss[!is.na(ss$spp),], mean)
params$size.av <- size.av$area[match(params$spp, size.av$spp)]

# minimum at reproductive maturity
params$min.r<-1/aggregate(area_cm2~spp, fec[fec$reproductive==1,], min)$area_cm2

# PCA of demographic parameters
colnames(params)
rownames(params) <- params$spp
traits <- c("r.int","f.cm2","survcm", "av.surv","p_mort", "min.r", "f.colony")
pca<-prcomp(params[,traits], scale=T, center=T)
biplot(pca)

# explained variation
exp<-round(c(summary(pca)[[1]][1]^2/sum(summary(pca)[[1]]^2),summary(pca)[[1]][2]^2/sum(summary(pca)[[1]]^2)),3)*100
exp

# Figure S2 - fast/slow continuum
#source("figs/supp/fig.S2_lms.R")	
#fig.s2

#######################################
# ABUNDANCE
#######################################

head(abun)
ab2004 <- ggplot(params, aes(x=abundance_05/2700, y=species))+
geom_bar(stat="identity")+
ggtitle("2004")+xlab("colonies per m2")
ab2011_14 <- ggplot(data=abun, aes(x=N/10, y=species))+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
facet_grid(.~year, scales="free_x", space="free_x")+xlab("Colonies per m")

plot_grid(ab2004, ab2011_14, rel_widths=c(1,1.5))

# Figure S3 - abundances
#source("figs/supp/fig.S3_abun.R")	
#fig.s3

#######################################
# FIGURE 1
#######################################

source("figs/fig.1.R")
fig.1
#ggsave("figs/fig.1.png", fig.1, width=15, height=9.5, units="cm", dpi = 300)

#######################################
# IPM MESH AND BOUNDARIES
#######################################

max.size <- 0.3 # max(ss$area)
n <- 100
mesh <- function(){
	min.size <- rec.size
	b <- seq(min.size, max.size, length=n)
	h <- b[2] - b[1]
	b <- c(min(b)-(2*h), min(b)-h, b)
	y <- 0.5 * (b[1:n]+b[2:(n+1)])
	I <- y >= rec.size
	return(list(b=b,h=h,y=y,I=I, rec.size = rec.size))}
	
	#######################################
# RECRUITMENT LIT REVIEW
#######################################

#  Assume: (A) closed system, (B) 11 species = all acroporids/favids
density <- params$abundance_05/2700# 275*1*10m 
area.m2 <- log10(density*(10^(params$size.av)))
params$fec.m2 <- exp(params$f.int+params$f.slp*area.m2) 
fam <- aggregate(fec.m2~family, params, sum)
tiles$eggs<-fam$fec.m2[match(tiles$Family, fam$family)]
tiles$N_m2_year[tiles$N_m2_year==0] <-1
tiles$p.rec <- tiles$N_m2_year/tiles$eggs
rec.fam <- aggregate(p.rec~Family, tiles, median)

plot_grid(ggplot(tiles, aes(x=(N_m2_year/eggs), y=Family))+
geom_boxplot()+scale_x_log10()+ggtitle("Recruitment"),
ggplot(rsize, aes(x=recsize_1yr, y=Family))+geom_boxplot()+ggtitle("Recruit size"), ncol=1)


# recsize
rec.size.const <- log10(pi*((5/2)/100)^2) # 5cm diameter
params$rec.size <- log10(pi*(params$r.int*(10/12))^2)
r.limit <- log10(pi*((10/100)/2)^2)
params$rec.size <- ifelse(params$rec.size > r.limit,r.limit,params$rec.size)

# recruitment 
rec.const <- 10^-3
params$rec <- rec.fam$p.rec[match(params$family, rec.fam$Family)]
params$rec <- ifelse(params$family=="Merulinidae", 5.900139e-05, 1.5*10^-3)


#######################################
# GENERATE IPMS
#######################################

par(mfcol=c(2, 6))

ipm.k.list <- list()
ipm.p.list <- list()
ipm.r.list <- list()
y.list <- list()

for (sp in spp) {
	#sp<-"Ahy"
	rec <- 1
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	sub<-gdat[gdat$spp==sp,]
    mod <- bigmatrix()
   # plot
	image(y, y, t(mod$P)^0.3)   
	points(sub$area, sub$area_next, cex=0.25)
	title(sp, line=-1)
	abline(0, 1, lty=2)
	y.list[[sp]] <- y 
	ipm.k.list[[sp]] <- mod$P
	ipm.p.list[[sp]] <- mod$P
	ipm.r.list[[sp]] <- mod$R
	}

# Figure S4
source("figs/supp/fig.s4_ipms.R")	
fig.s4


#######################################
# ESTIMATE LAM
#######################################

lam.est <- NULL
lam.const <- NULL
for (sp in spp) {
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
lam.est<-rbind(lam.est, data.frame(sp, lam=bigmatrix()$lam))
rec.size <- rec.size.const
rec <- rec.const
h <- mesh()$h
y <- mesh()$y
lam.const<-rbind(lam.const, data.frame(sp, lam=bigmatrix()$lam))
}
params$lam.est <- lam.est$lam
params$lam.const <- lam.const$lam

ggplot(params, aes(reorder(spp, -lam.est), lam.est))+geom_bar(stat="identity", aes(fill=spp))+scale_fill_manual(values=cols)+geom_hline(yintercept=1)

#######################################
# INCREASE/DECREASE
#######################################

head(abun)
diffs <- data.frame(dcast(abun, species+tran~year, value.var="N"))
sp.avs<-aggregate(N~species, abun, mean)
diffs$orig <- sp.avs$N[match(diffs$species, sp.avs$species)]
diffs$diff <- (diffs$X2014-diffs$X2011)/diffs$orig

diffsC <- data.frame(dcast(abun, species+tran~year, value.var="cover"))
sp.avs<-aggregate(cover~species, abun, mean)
diffsC$orig <- sp.avs$cover[match(diffsC$species, sp.avs$species)]
diffs$cover <- (diffsC$X2014 - diffsC$X2011)/diffsC$orig

diffs$spp<-params$spp[match(diffs$species, params$species)]
diffs$lam <- params$lam.est[match(diffs$spp, params$spp)]

ggplot(diffs, aes(reorder(spp, -lam), diff, fill=spp))+
stat_summary(geom="bar", fun="mean", col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
stat_summary(aes(reorder(spp, -lam), cover, fill=spp), geom="point", fun="mean")

#######################################
# RECRUITMENT SENSITIVITY I
#######################################	

# lambda at different rec/recsize combinations
rec.x <- 10^(seq(-6,-2,0.5))
recsize.x <- seq(min(params$rec.size)-0.05, max(params$rec.size)+0.1,0.1)
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

#######################################
# FIG2
#######################################

#speed <- 0.1 # slow but final
speed <- 0.5 # fast
source("figs/fig.2.R")
fig2
#ggsave("figs/fig.2.png", fig2, width=15, height=8, units="cm", dpi = 250)

#######################################
# COMPARE MORPHOLOGIES
#######################################

params2<-rbind(params, params[6,]) 
params2$morph<-as.character(params2$morph)
params2$morph[c(7, 12)]<-c("corymbose_2","corymbose_2")#AN/AM
comp<-dcast(params2, morph~abundance_pair, value.var="spp")

colsC <- cols[names(cols) %in% comp$Common]
names(colsC)<-comp$morph[order(comp$Common)]

comp$AC<-params2$abundance_05[match(comp$Common, params2$spp)]
comp$AR<-params2$abundance_05[match(comp$Rare, params2$spp)]
comp$abundiff <- comp$AC/comp$AR	
comp$lamC <- params$lam.est[match(comp$Common, params$spp)]	
comp$lamR <- params$lam.est[match(comp$Rare, params$spp)]
comp$lamdiff <- comp$lamC - comp$lamR
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

proj2 <- dcast(morph+gen~abundance_pair, data=proj, value.var="n")
proj2$diff <- proj2$Common/proj2$Rare
proj2$real <- comp$abundiff[match(proj2$morph, comp$morph)]
gens <- aggregate(gen~morph, proj2[proj2$diff<proj2$real,], max)
proj2$maxgen <- gens$gen[match(proj2$morph, gens$morph)] +1 #to cross the line
head(proj2)

ggplot(proj2, aes(gen, diff, col=morph))+geom_line(linetype="dotted")+coord_cartesian(ylim=c(0,60))+geom_line(data=proj2[proj2$gen<=proj2$maxgen,])+geom_point(data=proj2[proj2$gen==proj2$maxgen,])
	
#######################################
# RECRUITMENT SENSITIVITY II
#######################################	
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
geom_point(data=params2, aes(rec.const, rec.size.const ), col="white", shape=3, stroke=1, size=0.3)+
geom_point(data=params2, aes(rec, rec.size ), col="white", shape=3)+
	scale_y_continuous(expand=c(0,0))+
	scale_x_log10(expand=c(0,0))+
	scale_fill_distiller(palette="Spectral")+
	facet_wrap(~morph, scales="free")+
	labs(x="rec rate", y="rec size")+
	theme(strip.background=element_blank())

# Figure S5
#source("figs/supp/fig.s5_rec.R")	
#fig.s5
	
#######################################
# FIG3
#######################################

# Figure 3
source("figs/fig.3.R")
fig3
#ggsave("figs/fig.3.png", fig3, width=14, height=9, units="cm", dpi = 250)

#######################################
# REPRODUCTIVE TRADEOFF
#######################################
ec$morph<-params$morphology[match(ec$spp, params$spp)]
ecplot <- ggplot()+
geom_boxplot(data=ec, aes(y=reorder(morph, -Carbon_ug_corrected), x=Carbon_ug_corrected, fill=spp), outlier.size=0.2, size=0.1,position = position_dodge(width = 1))+scale_fill_manual(values=cols)+guides(fill="none")

#mat<-inv.logit(params$m.int[params$spp==sp]+params$m.slp[params$spp==sp]*x)
params$fec1cm <- exp(params$f.int +  params$f.slp * log10(1)) /10000
params2$fec1cm <- params$fec1cm[match(params2$spp, params$spp)]

plot_grid(ecplot,ggplot(params, aes(fec1cm, eggC))+geom_text(aes(label=spp))+geom_smooth(data=params[params$family=="Acroporidae",], method="lm"), nrow=1, rel_widths=c(1.15,1))

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
# ELASTICITY ACROSS SIZES
#######################################

spp2 <- spp
ek.hist <- NULL
for (sp in spp2) {
#	sp <- "Gpe"
	y <- y.list[[sp]]
	e.k <- melt(eK.list[[sp]])
	e.k$y <- rep(y, length(y))
	e.k.sum <- aggregate(value~Var1+y, e.k, sum)
	ek.hist <- rbind(ek.hist, cbind(e.k.sum, spp=sp))
	}
head(ek.hist)

ek.hist$size <- round(ek.hist$y/0.45)*0.45
ek.av <- aggregate(value~size+spp, ek.hist, sum)
ek.av$morphology<-params$morphology[match(ek.av$spp, params$spp)]
ek.av$X2 <- ek.av$value#*h^2  #/ek.av$max

ggplot(ek.av[!ek.av$spp=="Asp",], aes(size,X2))+
geom_bar(stat="identity", position=position_dodge(preserve = "single"), aes(fill=spp), col="black", size=0.1, width=0.33)+
facet_wrap(~morphology, scales="free_x", nrow=2)

#######################################
# Figure 4
#######################################
source("figs/fig.4.R")
fig4 
#ggsave("figs/fig.4.png", fig4, width=10, height=11, units="cm", dpi = 250)

#######################################
# WITHIN VS BETWEEN
#######################################

params
dems <- c("lam.est", "eR","lam.const")
vars <- NULL
for (t in c(traits, "eggC", dems,  "abundance_05")){
	#t <- "f.cm2"
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

source("figs/fig.5.R")
fig.5
#ggsave("figs/fig.5.png", withinplot, width=6, height=6, units="cm", dpi = 250)

