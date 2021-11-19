
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
# TO DO: 
# bootstrapping
# fix goni
# density dependence
# storms

#######################################
# LOAD DATA
#######################################
source("R/functions.R")
source("R/data_prep.R")
#sdat<-sdat[!sdat$colony_id==287,]
#gdat<-gdat[!gdat$colony_id==287,]

#sdat<-sdat[!sdat$colony_id==358,]
#gdat<-gdat[!gdat$colony_id==358,]

#sdat<-sdat[!sdat$colony_id==337,]
#gdat<-gdat[!gdat$colony_id==337,]

#gdat<-gdat[!gdat$colony_id==345,]
#sdat<-sdat[!sdat$colony_id==345,]

#gdat<-gdat[!gdat$colony_id==341,]
#sdat<-sdat[!sdat$colony_id==341,]


#gdat <- gdat[!gdat$colony_id==356,]
#sdat <- sdat[!sdat$colony_id==356,]

#gdat <- gdat[!gdat$colony_id==331,]
#sdat <- sdat[!sdat$colony_id==331,]

#gdat <- gdat[!gdat$colony_id==351,]
#sdat <- sdat[!sdat$colony_id==351,]

#gdat <- gdat[!gdat$colony_id==349,]
#sdat <- sdat[!sdat$colony_id==349,]



# most
gdat <- subset(gdat, area > log10((pi*(3/2)^2)/10000))
sdat <- subset(sdat, area > log10((pi*(3/2)^2)/10000))
fec <- subset(fec, area > log10((pi*(3/2)^2)/10000))

#gdat <- subset(gdat, area > log10((pi*(1)^2)/10000))
#sdat <- subset(sdat, area > log10((pi*(1)^2)/10000))
#fec <- subset(fec, area > log10((pi*(1)^2)/10000))

#######################################
# When a colony shrinks to under 5cm it dies..? 
#######################################




#######################################
# MODEL PARAMETERS
#######################################
source("R/params.R")
#fig.s1

params[params$spp=="Acy",]
#ggsave("figs/supp/fig.s1.png", fig.s1, width=25, height=13, units="cm", dpi = 300)
params
#write.csv(params, "params.csv")


params[params$spp=="Gre", c("p.int", "p.slp","p.sig")] <- params[params$spp=="Gpe", c("p.int", "p.slp","p.sig")]

params[params$spp=="Gre", c("g.int", "g.slp","g.var")] <- params[params$spp=="Gpe", c("g.int", "g.slp","g.var")]


#######################################
# DEMOGRAPHIC SPACE
#######################################
# size
#ggplot(ss, aes(x=area, col=spp))+
#geom_density()+scale_colour_manual(values=cols)

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
ab2005 <- ggplot(data=d08B, aes(x=N/10, y=Species))+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
ggtitle("2005")+xlab("colonies per m2")
ab2011_14 <- ggplot(data=abun, aes(x=N/10, y=species))+
stat_summary(fun="mean", geom = "bar", width=0.8, col="black", size=0.1)+
stat_summary(fun.data = mean_se, geom = "errorbar", width=0, size=0.2)+
facet_grid(.~year, scales="free_x", space="free_x")+xlab("Colonies per m")

#plot_grid(ab2005, ab2011_14, rel_widths=c(1,1.5))


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

#plot_grid(ggplot(tiles, aes(x=(N_m2_year/eggs), y=Family))+
#geom_boxplot()+scale_x_log10()+ggtitle("Recruitment"),
#ggplot(rsize, aes(x=recsize_1yr, y=Family))+geom_boxplot()+ggtitle("Recruit size"), ncol=1)


# recsize
rec.size.const <- log10(pi*((5/2)/100)^2) # 5cm diameter
params$rec.size <- log10(pi*(params$r.int*(12/12))^2)
agg <- aggregate(rec.size~morphology, params, mean)
params$rec.size <- agg$rec.size[match(params$morphology, agg$morphology)]
r.limit <- log10(pi*((10/100)/2)^2)
params$rec.size <- ifelse(params$rec.size > r.limit,r.limit,params$rec.size)



# recruitment 
rec.const <- 10^-3
#params$rec <- rec.fam$p.rec[match(params$family, rec.fam$Family)]
params$rec <- ifelse(params$family=="Merulinidae", 5.900139e-05, 2*10^-3)


#######################################
# GENERATE IPMS
#######################################

#par(mfcol=c(2, 6))

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
	#image(y, y, t(mod$P)^0.3)   
	#points(sub$area, sub$area_next, cex=0.25)
	#title(sp, line=-1)
	#abline(0, 1, lty=2)
	y.list[[sp]] <- y 
	ipm.k.list[[sp]] <- mod$P
	ipm.p.list[[sp]] <- mod$P
	ipm.r.list[[sp]] <- mod$R
	}

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

params[,c("spp","lam.est", "g.slp")]




#######################################
# REC FIT TO SIZE STRUCTURE
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
	

rec.ll <- function(x) {
  cnt <- size.dist$count[II] # non-recruits
  rec <<- x[1] 
  mod <- bigmatrix()
  eig.vec <- mod$w[II]/sum(mod$w[II])
  return(-sum(cnt * log(eig.vec), na.rm=TRUE)) } # log-likelihood 

rec.ss <- function(x) {
  cnt <- size.dist$count[II] # non-recruits
  cnt <- cnt / sum(cnt)
  rec <<- x[1] 
  mod <- bigmatrix()
  eig.vec <- mod$w[II]/sum(mod$w[II])
  return(sum((cnt - eig.vec)^2)) } # log-likelihood 
  
  
# ss_vec <- 10^seq(-3, -7, -0.01)
# store <- data.frame()
# for (ss in ss_vec) {
	
 # cnt <- size.dist$count[II] # non-recruits
 # cnt <- cnt / sum(cnt)
 # rec <- ss 
 # mod <- bigmatrix()
 # eig.vec <- mod$w[II]/sum(mod$w[II])
 # store <- rbind(store, data.frame(rec=ss, ss=sum((cnt - eig.vec)^2)))

# }  

# plot(ss ~ rec, store)
# store[which.min(store$ss),]

dev.off()  
par(mfcol=c(4, 3), mar=c(3,3,1,1))

n <- 12

params$rec_fit <- NA
params$lam_fit <- NA 

for (sp in spp) {
#sp <- "Asp"
  # Model
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
b <- mesh()$b 
I <- mesh()$I
II <- I
II[1:6] <- FALSE

max.size <- 1
size.dist <- hist(ss$area[ss$spp==sp & ss$area > rec.size & ss$area < max.size], breaks=b, plot=FALSE)

# hist(ss$area[ss$spp==sp & ss$area > rec.size & ss$area < max.size], breaks=b)
#, lower = 0, upper=10000 removed
# x is a startin rec value?
# mle takes a function and finds the most likely value of x ? 

	rec.fit <- optimise(rec.ss, c(0, 1))
	rec <- rec.fit$minimum
	rec
	
	#rec.fit <- mle(rec.ll, start = list(x = 0.001), method="Brent", lower = 0, upper=1 )
	rec.fit <- optimise(rec.ll, c(0, 100))
	rec <- rec.fit$minimum
	rec

#rec <- 0.00001
  mod <- bigmatrix()
	#image(y, y, log(t(mod$K)))
	#title(sp, line=-1)
	#abline(0, 1, lty=2)
mod$lam

	hist(ss$area[ss$spp==sp & ss$area > rec.size & ss$area < max.size], breaks=b, freq=FALSE, main="", ylim=c(0, 2))
	lines(y[II], (mod$w[II]/sum(mod$w[II]))/h/2, col="red")
	abline(v=y[45], lty=2)
	title(sp, line=-1)
text(0,1, round(mod$lam, 3))


  params$lam_fit[params$spp==sp] <- mod$lam
  params$rec_fit[params$spp==sp]<- rec
}

params$rec_fit


ggplot(params, aes(reorder(spp, -lam_fit), lam_fit))+geom_bar(stat="identity", aes(fill=spp))+scale_fill_manual(values=cols)+geom_hline(yintercept=1)



##### CHANGE REC?  #!!!!!!!
params$rec <- params$rec_fit

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
facet_wrap(~morphology, scales="free", ncol=1)+
scale_fill_manual(values=cols)







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
diffs$lam <- params$lam_fit[match(diffs$spp, params$spp)]

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

# Figure S4
#source("figs/supp/fig.s4_ipms.R")	
#fig.s4

#speed <- 0.1 # slow but final
#speed <- 0.5 # fast
#source("figs/fig.2.R")
#fig2
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

comp$lamC <- morph.diff("lam_fit")$x.common
comp$lamR <- morph.diff("lam_fit")$x.rare
comp$lamdiff <- morph.diff("lam_fit")$x.diff
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
	
# method 3 - projecting lamda

#for(sp in spp){
	sp <-"Ahy"
	rec <- params$rec[params$spp==sp] 
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	lam=bigmatrix()$lam
	v=bigmatrix()$v
	w=bigmatrix()$w
	ggplot(NULL, aes(y, v))+geom_point()

pop  <- matrix(nrow=100, ncol=50)
  pop[,1] <- round(v * 1000) # this is the inital population size for the EvoPop
  for(j in 2:50){
    pop[,j] <- pop[,j-1]*lam
  }
head(pop)

head(melt(pop))
ggplot(data=melt(pop), aes(x=Var1, y=value, col=Var2, group=Var2))+geom_line()+scale_y_log10()
psize<-	aggregate(value~Var2, data=melt(pop), sum)
psize$n <- 1000*lam^psize$Var2

	ggplot()+
	geom_point(data=psize, aes(Var2, value))+
	geom_point(data=psize, aes(Var2, n), col="red")+scale_y_log10()
	
	
	
#######################################
# STORMS - WRONG>> currently
#######################################		


params2$lam.pos <- comp$lamdiff[match(params2$morph, comp$morph)]
params2$lam.pos <- ifelse(params2$abundance_pair=="Rare",1.1,  1.1+params2$lam.pos)
#params2$lam.pos[params2$spp=="Gpe"]<-1.01

n1 <- 100
ngen <- 30
pars <- params2[, c("spp","morph", "abundance_pair","lam.pos", "lam.est")]
proj <- merge(pars, data.frame(gen=1:ngen))
proj$n <- n1*proj$lam.pos^proj$gen

proj$mort <- ifelse(proj$morph=="tabular", 0.05,ifelse(proj$morph=="staghorn", 0.35,ifelse(proj$morph=="corymbose", 0.1,ifelse(proj$morph=="corymbose_2",0.1,  0.8))))
proj$n2 <- ifelse(proj$gen>10, proj$n*proj$mort, proj$n)
proj$n3 <- ifelse(proj$gen>20, proj$n2*proj$mort, proj$n2)


ggplot(proj[proj$n3<4000,], aes(gen, n3))+geom_line(aes(group=spp, col=spp))+scale_colour_manual(values=cols)+facet_wrap(~morph, scale="free_y")+
#scale_y_log10()+
geom_segment(data=params2, aes(x=-Inf, xend=Inf, y=abundance_05, yend=abundance_05, col=spp), linetype="dotted")
	

	
	
	
	
	
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
#source("figs/fig.3.R")
#fig3
#ggsave("figs/fig.3.png", fig3, width=14, height=9, units="cm", dpi = 250)



#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
# params_orig <- params # CAREFUL set original


#######################################
# SENSITIVITY ANALYSIS!!  
# Change each param by 10%... 
#######################################

pars <- c("m.int", "m.slp","f.int","f.slp","r.int","p.int","p.slp","p.sig", "s.int", "s.slp") #"s.slp.2"

p.types <- data.frame(pars, type = c("reproduction","reproduction","reproduction","reproduction","growth","growth","growth", "growth", "survival","survival"))
#p.types <- data.frame(pars, type = c("maturity","maturity","fecundity","fecundity","growth","growth","growth", "growth", "survival","survival","survival"))

types <- unique(p.types$type)

sens <- NULL
sens2 <- NULL
for (i in pars) {
#i <- "m.int"
for (sp in spp) {
#sp <- "Ahy"
params <- params_orig 
par.orig <- params[params$spp==sp, i]
#pars2 <- p.types$pars[p.types$type==i]
params[,i] <- params[,i]*1.1
#params[,pars2] <- params[,pars2]*1.1
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
sens <- rbind(sens, data.frame(spp=sp, param=i, lam=bigmatrix()$lam, morph=params$morphology[params$spp==sp], p = "plus10", par.val = params[params$spp==sp, i], par.orig))
params <- params_orig 
params[,i] <- params[,i]*0.9
sens2 <- rbind(sens2, data.frame(spp=sp, param=i, lam=bigmatrix()$lam, morph=params$morphology[params$spp==sp], p ="minus10", par.val = params[params$spp==sp, i], par.orig))
}}
head(sens)

(sens$par.val-sens$par.orig)/sens$par.orig *100

sens$orig <- params$lam.est[match(sens$spp, params$spp)]
sens$change <- (sens$lam-sens$orig)/sens$orig

sens2$orig <- params$lam.est[match(sens2$spp, params$spp)]
sens2$change <- (sens2$lam-sens2$orig)/sens2$orig

sens.df <- data.frame(cbind(sens, sens2))


plot_grid(
ggplot(sens, aes(change, reorder(param, abs(change))))+
geom_bar(stat="identity", aes(fill=spp), position="dodge")+
facet_wrap(~morph, nrow=1)+
scale_fill_manual(values=cols)+geom_vline(xintercept=0)+
labs(x="% change in lambda", y="varying parameter")+
guides(fill="none")
,
ggplot(sens.df, aes(change, change.1))+
geom_abline(slope=-1)+
geom_point(aes(col=spp))+
geom_text_repel(data=sens.df[abs(sens.df$change) >0.2,], aes(label=param), size=2)+
scale_colour_manual(values=cols)+
labs(x="plus 10% change", y="minus 10% change")
, rel_widths=c(1, 0.4))


#######################################
# SENSITIVITY ANALYSIS!!  
# Demo differences affect lam differences..
#######################################

source("R/params_morph.R")
head(p.morph)

rec.morph <- aggregate(rec~morphology, params, mean)
p.morph$rec <- rec.morph$rec[match(p.morph$morph, rec.morph$morph)]

#pars <- c("m.int", "m.slp","f.int","f.slp","r.int","p.int","p.slp","p.sig", "s.int", "s.slp", "s.slp.2")
pars <- c("m.int", "m.slp","f.int","f.slp","g.int","g.slp","g.var", "s.int", "s.slp", "s.slp.2", "rec")

#p.types <- data.frame(pars, type = c("reproduction","reproduction","reproduction","reproduction","growth","growth","growth", "growth", "survival","survival","survival"))
p.types <- data.frame(pars, type = c("maturity","maturity","fecundity","fecundity","growth","growth","growth", "survival","survival","survival", "recruitment"))
#p.types <- data.frame(pars, type = c("intrinsic","intrinsic","intrinsic","intrinsic","intrinsic","extrinsic","extrinsic", "extrinsic", "extrinsic","extrinsic","extrinsic"))

types <- unique(p.types$type)

test <- NULL
test2 <- NULL
for (i in types) {
#i <- "m.int"
#i = "repro"
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
rec <- params$rec[params$spp==sp] * 2
test2<-rbind(test2, data.frame(spp=sp, param=i, lam=bigmatrix()$lam, morph=params$morphology[params$spp==sp], rec="double"))
 } 
}

head(test)
ggplot(rbind(test, test2), aes(lam, param, fill=rec))+geom_bar(stat="Identity", position="dodge")+facet_wrap(~spp)+geom_vline(xintercept=1)

diffs <- NULL
diffs2 <- NULL
for (i in types) {
	#i <- "m.int"
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


plot_grid(
ggplot(rbind(diffs, diffs2), aes(d, reorder(param, -abs(d)), fill=rec))+
geom_bar(stat="Identity", position="dodge")+
facet_wrap(~morph)+
labs(x="Lamda difference (Common - Rare)", y = "Varying parameter")+
geom_vline(xintercept=0),
ggplot(rbind(sums,sums2), aes(lamdiff, d))+geom_point(aes(col=morph, shape=rec))+
scale_y_sqrt()+scale_x_sqrt()+
geom_abline(slope=1)+
labs(x="Real lamda difference", y = "Sum of differences"), 
rel_widths=c(1,0.5))



#######################################
# SENSITIVITY ANALYSIS!!  
# BETWEEN ALL SPECIES !!! 
#######################################

source("R/params_all.R")
head(p.all)

p.all$rec <- mean(params$rec)
p.all$rec.size <- mean(params$rec.size)

pars <- c("m.int", "m.slp","f.int","f.slp","r.int","p.int","p.slp","p.sig", "s.int", "s.slp", "s.slp.2")

p.types <- data.frame(pars, type = c("maturity","maturity","fecundity","fecundity","growth","partial","partial", "partial", "survival","survival","survival"))

pars <- c("m.int", "m.slp","f.int","f.slp","g.int","g.slp","g.var", "s.int", "s.slp", "s.slp.2", "rec", "rec.size")

p.types <- data.frame(pars, type = c("maturity","maturity","fecundity","fecundity","growth","growth","growth", "survival","survival","survival", "recruitment", "recruit size"))


types <- unique(p.types$type)

# general lamda
params <- params_orig
params[,pars] <- p.all[, pars]

	sp <- "Gre"
rec.size <- mean(params$rec.size)
rec <- mean(params$rec) #params$rec.size[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
gen_lam <- bigmatrix()$lam
gen_lam 
 



test <- NULL
for (i in types) {
#i <- "m.int"
#i = "fecundity"
params <- params_orig
params[,pars] <- p.all[, pars]
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
 } 
}

head(test)
ggplot(test, aes(lam, param))+geom_bar(stat="Identity", position="dodge")+facet_wrap(~spp)+geom_vline(xintercept=1)

head(test)
test$lam_orig <- params$lam_fit[match(test$spp, params$spp)]

test$morph <- factor(test$morph, levels=c("massive","digitate","corymbose","staghorn","tabular"))


test$diff <- test$lam - gen_lam #test$lam_orig #gen_lam
ggplot(test, aes(diff, reorder(param, diff)))+geom_bar(stat="Identity", position="dodge", aes(fill=spp))+facet_wrap(~morph, nrow=1)+geom_vline(xintercept=0)+scale_fill_manual(values=cols)+
labs(y="varying parameter", x="change in species lambda when parameter varies")


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
facet_wrap(~morphology, scales="free_x", nrow=2)+
scale_fill_manual(values=cols)

#######################################
# Figure 4
#######################################
source("figs/fig.4.R")
fig4 
#ggsave("figs/fig.4.png", fig4, width=10, height=11, units="cm", dpi = 250)

# long
#ggsave("figs/fig.4.png", fig4, width=10, height=12, units="cm", dpi = 250)


#######################################
# WITHIN VS BETWEEN
#######################################

params <- params_orig
#params <- params[!params$spp %in% c("Gre","Gpe"),]
dems <- c("lam.est", "eR","lam.const", "lam_fit")
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

#######################################
# CORELLATE WITH FITNESS
#######################################



plot_grid(
ggplot(params, aes(r.int, lam.est))+geom_text(aes(label=spp))+
scale_x_log10()+
geom_smooth(method="lm", formula=y~poly(x, 2)),
ggplot(params, aes(f.colony, lam.est))+geom_text(aes(label=spp))+
scale_x_log10()+
geom_smooth(method="lm", formula=y~poly(x, 2)),
ggplot(params, aes(f.int2, lam.est))+geom_text(aes(label=spp))+
geom_smooth(method="lm", formula=y~poly(x, 2)),
ggplot(params, aes(av.surv, lam.est))+geom_text(aes(label=spp))+
geom_smooth(method="lm", formula=y~poly(x, 2))
)

plot_grid(
ggplot(params, aes(r.int, lam.const))+geom_text(aes(label=spp))+
scale_x_log10()+
geom_smooth(method="lm", formula=y~poly(x, 2)),
ggplot(params, aes(f.colony, lam.const))+geom_text(aes(label=spp))+
scale_x_log10()+
geom_smooth(method="lm", formula=y~poly(x, 2)),
ggplot(params, aes(f.int2, lam.const))+geom_text(aes(label=spp))+
geom_smooth(method="lm", formula=y~poly(x, 2)),
ggplot(params, aes(av.surv, lam.const))+geom_text(aes(label=spp))+
geom_smooth(method="lm", formula=y~poly(x, 2))
)

#######################################
# CONSTANT SLOPE MODEL
#######################################

fec.con<-glm.nb(fecundity~area+spp, data=fec[fec$reproductive==1,], link = log)
summary(fec.con)

x2 <- spp # your second predictor (here, sites)
int <-function(mod){ 
	int <- rep(coef(mod)[[1]], length(x2))
	for (i in 2:length(x2)) { int[i] <- int[i] + coef(mod)[grepl(x2[i],
		 names(coef(mod))) & !grepl("area", names(coef(mod)))]}
	int} 
	
slp <-function(mod){ 
	slp <- rep(coef(mod)[[2]], length(x2))
	for (i in 2:length(x2)) { slp[i] <- slp[i] + coef(mod)[grepl(x2[i],
		names(coef(mod))) & grepl("area", names(coef(mod)))]}
	slp} 
# input your model, output in the same order as “x2”
params$f.int2 <- int(fec.con) # intercept = effect of x2
params$f.int


ggplot(params, aes(spp, f.int2))+geom_bar(stat="identity")+facet_wrap(~morphology, scale="free_x")+coord_cartesian(ylim=c(15, 16))

ggplot(params, aes(f.int2, eggC))+geom_text(aes(label=spp))

