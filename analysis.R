
# source("analysis.R")

rm(list = ls())
library("ggrepel")
library("grid")
library("png")

#######################################
# LOAD SCRIPTS
#######################################
source("R/functions.R")
source("R/data_prep.R")

#source("R/goni.R") # removes smalls to fix growth curve
sdat<-sdat[!sdat$colony_id==287,]
gdat<-gdat[!gdat$colony_id==287,]

#gdat <- gdat[!gdat$colony_id==356,]
#gdat <- gdat[!gdat$colony_id==331,]
#gdat <- gdat[!gdat$colony_id==351,]
#gdat <- gdat[!gdat$colony_id==349,]

cutoff_diam<-10
cutoff_a <- log10(a_func(cutoff_diam/2/100))
#gdat <- subset(gdat, area > cutoff_a  )
#gdat <- subset(gdat, area_next > cutoff_a  )

source("R/params.R")
params[,c("spp", "g.slp")]

#######################################
# SADs
#######################################

sad <- read.csv("data/dornelas2008.csv")
plot_grid(
ggplot(params, aes(y=morphology, x=abundance_05))+
geom_line(aes(col=morphology), arrow=arrow(length=unit(1,"mm")))+
xlim(c(-100, max(sad$abundance)+100))+
guides(col="none", fill="none"),
ggplot(sad, aes(x=abundance))+
geom_histogram(bins=30, col="white")+
#scale_y_sqrt(expand=c(0,0), breaks=c(5,25,55,100))+
xlim(c(-100, max(sad$abundance)+100))+
labs(x="Abundance (N)", y="N species"),
ncol=1, align="v", axis="lr", rel_heights=c(0.5,1))

#######################################
# PARAMETER SPACE
#######################################

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
#ggplot(params, aes(f.cm2, f.int))+geom_text(aes(label=spp))

# average size
size.av <- aggregate(area~spp, ss[!is.na(ss$spp),], mean)
params$size.av <- size.av$area[match(params$spp, size.av$spp)]

# egg energy per area
params$en.cm2<-params$f.cm2*params$eggC

# minimum at reproductive maturity
params$min.r<-1/aggregate(area_cm2~spp, fec[fec$reproductive==1,], min)$area_cm2

# PCA of demographic parameters
colnames(params)
rownames(params) <- params$spp
traits <- c("r.int","f.cm2","survcm", "av.surv","p_mort", "min.r", "f.colony")
pca<-prcomp(params[,traits], scale=T, center=T)
biplot(pca)

exp<-round(c(summary(pca)[[1]][1]^2/sum(summary(pca)[[1]]^2),summary(pca)[[1]][2]^2/sum(summary(pca)[[1]]^2)),3)*100
exp

# abundance
ggplot(params, aes(reorder(spp, abundance_05), abundance_05))+geom_bar(stat="identity")+facet_wrap(~morphology, scales="free_x")

# size
ggplot(ss, aes(x=area))+
geom_density()+
facet_wrap(~spp)

# FIGURE 1
source("figs/fig.1.R")
fig.1
#ggsave("figs/fig.1.png", fig.1, width=13.7, height=9, units="cm", dpi = 300)


#######################################
# IPM FUNCTIONS
#######################################

#-------------------------------survival
s.x <- function(x) { 
	u <- params$s.int[params$spp==sp] + 
	  params$s.slp[params$spp==sp] * x + 
	    params$s.slp.2[params$spp==sp] * x^2
  return(inv.logit(u)) 
}

#------------------------------- growth & partial mortality
p.yx <- function(y, x) {
  # x <- -3
  g <- a_func(r_func(10^x) + params$r.int[params$spp==sp] )
  #+ 1.96 * params$r.err[params$spp==sp])
  temp <- 10^y / g  # proportion of max reached. 
  temp[temp > 1] <- 1
  dnorm(logit(1 - temp), params$p.int[params$spp==sp] + x * params$p.slp[params$spp==sp], params$p.sig[params$spp==sp])
}

#------------------------------- reproduction
 r.yx <- function(y, x) {	
 	mat<- inv.logit(params$m.int[params$spp==sp] + 
 	  params$m.slp[params$spp==sp] *x)
 	fec<- exp(params$f.int[params$spp==sp] + 
 	  params$f.slp[params$spp==sp] *x) 
   #siz<- rnorm(y,mean=params$rec.size[params$spp==sp], sd=0.05) 
   out <- (rec* mat * fec)
   out[x < rec.size | y >= rec.size] <- 0 #if x is below recruitment size
   return(out)
   } 

 #------------------------------- kernel
pmort<-T

bigmatrix <- function() {
    if (pmort) {
    G <- h * outer(y, y, p.yx)
  } else {
    G <- h * outer(y, y, g.yx)
  }
  G <- t(t(G) / apply(G, 2, sum))
  S <- s.x(y)
  P <- G 
  for(i in 1:n) P[,i]=G[,i]*S[i]
  R <- h * outer(y, y, r.yx) 
  #R <- h * outer(y, pmin(y, smax), r.yx) #  ceiling
  K <- P + R
  lam <- Re(eigen(K)$values[1])
	w <- abs(Re(eigen(K)$vectors[,1])) 
	v <- abs(Re(eigen(t(K))$vectors[,1]))
	return(list(K=K, lam=lam, w=w, v=v, G=G, S=S, R=R, P=P)) }
		
#######################################
# MESH AND BOUNDARIES
#######################################

max.size <- 1 
n <- 100

mesh <- function(){
	#min.size <- -4
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

# RECRUITMENT MODEL, Assume:
# (A) closed system at trimodal
# (B) 11 species constitute all acroporids/favids

density <- params$abundance_05/2700# 275*1*10m 
area.m2 <- log10(density*(10^(params$size.av)))
params$fec.m2 <- exp(params$f.int+params$f.slp*area.m2) 
fam <- aggregate(fec.m2~family, params, sum)
tiles$eggs<-fam$fec.m2[match(tiles$Family, fam$family)]
tiles$N_m2_year[tiles$N_m2_year==0] <-1
tiles$p.rec <- tiles$N_m2_year/tiles$eggs
rec.fam <- aggregate(p.rec~Family, tiles, median)

ggplot(tiles, aes(x=(N_m2_year/eggs), y=Family))+
geom_boxplot()+scale_x_log10()

# recruit size
ggplot(rsize, aes(x=recsize_1yr, y=Family))+geom_boxplot()
rsize.fam <- aggregate(recsize_1yr~Family, rsize, median)

#######################################
# ESTIMATE REC/LAM
#######################################

# recsize growth
params$rsize.gr <- log10(pi*(params$r.int*(9/12))^2)
r.limit <- log10(pi*((10/100)/2)^2)
params$rsize.gr <- ifelse(params$rsize.gr>r.limit,r.limit,params$rsize.gr)

# recsize family
params$rsize.fam <- log10(pi*((rsize.fam$recsize_1yr[match(params$family, rsize.fam$Family)]/2)/100)^2)

# recsize constant
params$rsize.const <- log10(pi*((5/2)/100)^2)

# recruitment family
params$rec.fam <- rec.fam$p.rec[match(params$family, rec.fam$Family)]

# recruitment energy
#eggCscale <- (scale(params$eggC)/max(abs(scale(params$eggC))))*params$rec.fam
params$recEn <- params$rec.fam * (10^(scale(params$eggC)/max(abs(scale(params$eggC))))) 

# SELECT!!!
params$rec <- params$rec.fam
params$rec.size <- params$rsize.gr

lam.est <- NULL
for (sp in spp) {
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
mod <- bigmatrix()
lam.est<-rbind(lam.est, data.frame(sp, lam=bigmatrix()$lam))}
 lam.est
params$lam.est <- lam.est$lam


#######################################
# PLOT IPMS
#######################################

par(mfcol=c(2, 6))

lam_const <- NULL
ipm.k.list <- list()
ipm.p.list <- list()
ipm.r.list <- list()
y.list <- list()

for (sp in spp) {
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

# Figure S2
source("figs/fig.S2.R")	
#fig.s2


#######################################
# COMPARE MORPHS
#######################################

# Is AL (middle man) classed as a rare or common?
# params$abundance_pair[params$spp=="Asp"] <- "Rare"
params2<-rbind(params, params[6,]) #duplicate AN, 7 (or AM=6)
params2$morph<-as.character(params2$morph)
params2$morph[c(7, 12)]<-c("corymbose_2","corymbose_2")#AN/AM
comp<-dcast(params2, morph~abundance_pair, value.var="spp")
comp$AC<-params2$abundance_05[match(comp$Common, params2$spp)]
comp$AR<-params2$abundance_05[match(comp$Rare, params2$spp)]
comp$diff <- comp$AC/comp$AR
comp

colsC <- cols[names(cols) %in% comp$Common]
names(colsC)<-comp$morph[match(names(colsC),comp$Common)]
colsC

#######################################
# MULTIPLE REC
#######################################

rec.x <- 10^(seq(-6,-2,0.5))

store<-data.frame()
for (sp in spp) {
	for (rec in rec.x){
	  rec.size <- params$rec.size[params$spp==sp]
	  h <- mesh()$h
	  y <- mesh()$y
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  store<-rbind(store, data.frame(sp=sp, rec=rec, lam=bigmatrix()$lam))
   } 
  }
 
 pairs <- NULL 
	for (rec in rec.x){
		df <- store[store$rec==rec,]
		temp <- comp
		temp$c.lam <- df$lam[match(temp$Common, df$sp)]
		temp$r.lam <- df$lam[match(temp$Rare, df$sp)]
		#temp$diff <- temp$c.lam - temp$r.lam
		#temp$ratio <- temp$c.lam / temp$r.lam
		temp$logdiff <- log(temp$c.lam) - log(temp$r.lam)
		temp$rec <- rec
		pairs <- rbind(pairs, temp)
		}

plot_grid(
ggplot(store, aes(x=rec, y=log(lam), col=sp))+
geom_hline(yintercept=0)+geom_line()+
 scale_x_log10()+scale_colour_manual(values=cols)+
 labs(x="probability of recruitment", y=expression(log(lambda))), 
	ggplot(pairs,aes(rec, logdiff, col=morph))+
	geom_hline(yintercept=0)+
		geom_line()+
		scale_x_log10()+
		labs(x="probability of recruitment", y=expression(Delta~log*(lambda)~(common~-~rare))))
		
# Rgrowth = rec required to grow..

Rgrow <- subset(store, lam>1)
Rg <- NULL
for (sp in unique(Rgrow$sp)){
	#sp <- "A"
	sub <- Rgrow[Rgrow$sp==sp,]
	min <- sub[which.min(sub$lam),"rec"]
	Rg <- rbind(Rg, data.frame(sp, min))}
params$Rgrow <- Rg$min[match(params$spp, Rg$sp)]

ggplot(params, aes(y=reorder(spp, -Rgrow), x=Rgrow,fill=spp))+geom_point(shape=21, size=4)+
scale_fill_manual(values=cols)+guides(fill="none")+scale_x_log10()


#######################################
# MULTIPLE REC. SIZE
#######################################

recsize.x <- seq(min(params$rec.size)-0.2, max(params$rec.size)+0.2,0.1)

store2 <- NULL
for (sp in spp) {
	for (rec.size in recsize.x){
	  rec <- params$rec[params$spp==sp]
	  h <- mesh()$h
	  y <- mesh()$y
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  store2<-rbind(store2, data.frame(sp=sp, rec.size=rec.size, lam=bigmatrix()$lam))
   } 
  }
 
pairs2 <- NULL 
	for (rec.size in recsize.x){
		df <- store2[store2$rec.size==rec.size,]
		temp <- comp
		temp$c.lam <- df$lam[match(temp$Common, df$sp)]
		temp$r.lam <- df$lam[match(temp$Rare, df$sp)]
		#temp$diff <- temp$c.lam - temp$r.lam
		#temp$ratio <- temp$c.lam / temp$r.lam
		temp$logdiff <- log(temp$c.lam) - log(temp$r.lam)
		temp$rec.size <- rec.size
		pairs2 <- rbind(pairs2, temp)
		}
			
plot_grid(ggplot(store2, aes(sqrt(((10^rec.size)*10000)/pi)*2, log(lam), col=sp))+
geom_hline(yintercept=0)+
 geom_line(size=1)+
 scale_colour_manual(values=cols)+
 labs(x="rec size (diam)", y=expression(log(lambda))),
 ggplot(pairs2,aes(sqrt(((10^rec.size)*10000)/pi)*2, logdiff, col=morph))+geom_line(size=1)+labs(x="rec.size (diam)", y=expression(Delta~log*(lambda)~(common~-~rare))))
		
# Rgrowth = rec required to grow..

Sgrow <- subset(store2, lam>1)
Sg <- NULL
for (sp in unique(Sgrow$sp)){
	#sp <- "Ahy"
	sub <- Sgrow[Sgrow$sp==sp,]
	min <- sub[which.min(sub$lam),"rec.size"]
	Sg <- rbind(Sg, data.frame(sp, min))}
params$Sgrow <- Sg$min[match(params$spp, Sg$sp)]

ggplot(params, aes(y=reorder(spp, -Sgrow), x=Sgrow,fill=spp))+geom_point(shape=21, size=4)+
scale_fill_manual(values=cols)+guides(fill="none")
		
#######################################
# MULTIPLE REC AND REC.SIZE
#######################################
 
rec.x <- 10^(seq(-6,-2,0.5))

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
  #store3<-rbind(store3, data.frame(sp=sp, rec=rec, lam=mod$lam))
  }
 } 
 store3[[sp]] <- temp  
  }
 
  minlam <- min(log(do.call("rbind", store3)))
  maxlam <- max(log(do.call("rbind", store3)))
  
par(mfrow=c(3,4))
for(sp in spp){
  image(log10(rec.x), recsize.x, log(store3[[sp]]), main=sp, zlim=c(minlam, maxlam))
  abline(h=params[params$spp==sp, "rec.size"], lty=2)
  contour(log10(rec.x), recsize.x, log(store3[[sp]]), add=TRUE, levels=0)
  }
   contour(log10(rec.x), recsize.x, log(store3[[spp[1]]]), levels=0, col=cols[1], lwd=3)
   #abline(h=params[params$spp==spp[1], "rec.size"], lty=2, col=cols[1], lwd=3)
 for(i in 2:length(spp)){
 	sp <- spp[i]
 	#abline(h=params[params$spp==sp, "rec.size"], lty=2, col=cols[i], lwd=3)
  contour(log10(rec.x), recsize.x, log(store3[[sp]]), add=TRUE, levels=0, col=cols[i], lwd=3)
  }


pairs3 <- NULL
for(m in comp$morph){
	#m <- "tabular"
	spC <- comp$Common[comp$morph==m]
	spR <- comp$Rare[comp$morph==m]
	logdiff <- log(store3[[spC]]) - log(store3[[spR]])
	temp <- reshape(data.frame(logdiff), direction="long",
    varying=list(c(1:ncol(logdiff))), ids=rec.x, times=recsize.x,
    new.row.names=c(1:(length(recsize.x)*length(rec.x))))
    #names(temp) <- c("recsize", "lam", "rec")
    temp$scaled <- (temp$X1-mean(temp$X1))/sd(temp$X1)
    pairs3 <- rbind(pairs3, cbind(temp, morph=m))
}
head(pairs3)

pairs3$lam2 <- ifelse(pairs3$X1 < 0, NA, pairs3$X1)


ggplot()+
geom_raster(data=pairs3, aes(x=id, y=time, fill=lam2))+
geom_point(data=params2, aes(rec.fam, rsize.fam ), col="white", shape=3, stroke=1, size=0.3)+
geom_point(data=params2, aes(recEn, rec.size ), col="white", shape=3)+
	scale_y_continuous(expand=c(0,0))+
	scale_x_log10(expand=c(0,0))+
	scale_fill_distiller(palette="Spectral")+
	facet_wrap(~morph, scales="free")+
	labs(x="rec rate", y="rec size")+
	theme(strip.background=element_blank())

#######################################
# MORPH PAIRS
#######################################
	
	
comp$rec <- params2$rec[match(comp$morph, params2$morph)] 
comp$rec.size <- params2$rec.size[match(comp$morph, params2$morph)] 
comp$rec.size.cm <- sqrt(((10^comp$rec.size)*10000)/pi)*2		
comp$lamC <- params$lam.est[match(comp$Common, params$spp)]	
comp$lamR <- params$lam.est[match(comp$Rare, params$spp)]

comp$lamdiff <- comp$lamC - comp$lamR
comp$logdiff <- log(comp$lamC) - log(comp$lamR)
comp$doub.time <- log(2)/log(1+comp$lamdiff)
comp$difftime <- log(comp$diff)/log(1+comp$lamdiff)
		
ggplot(comp[comp$difftime>0,], aes(x=difftime, y=reorder(morph, -difftime)))+geom_bar(stat="Identity")

#######################################
# FIG2
#######################################

# Figure 2
source("figs/fig.2.R")
fig2
#ggsave("figs/fig.2.png", fig2, width=13, height=10.5, units="cm", dpi = 300)



#######################################
# ELASTICITY & IPM MEASURES
#######################################
 # eigen-things can be combined to obtain the sensitivity and elasticity matrices.

s.list <- list()
eK.list <- list()
eR.list <- list()
eP.list <- list()
demovals <- NULL
sizevals <- NULL

for (sp in spp) {
	sp <- "Gre"
	# run IPM
	rec <- params$rec[params$spp==sp]
	#rec<-1*10^-4
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	K <- bigmatrix()$K
	P <- bigmatrix()$P
	R <- bigmatrix()$R
	image(t(K^0.2))
	# population growth
	lam <- Re(eigen(K)$values[1])
	
	# sable size dist /right eigenvec
	w.eigen<-Re(eigen(K)$vectors[,1]) 
	stable.size <- w.eigen/sum(w.eigen)
	
	# reprodutive values /left eigenvec (size contributions) # Gre goes wrong
	v.eigen<-Re(eigen(t(K))$vectors[,1]) 
	repro.val <- v.eigen/v.eigen[1] # relative
	
	# reproductive val * stable size
	v.dot.w<-sum(stable.size*repro.val)*h 
	
	# sensitivity matrix
	sens<-outer(repro.val,stable.size)/v.dot.w  
	
	# elasticity matrices
	K.elas<-matrix(as.vector(sens)*(as.vector(K)/h)/lam, nrow=n)
	
	# survival elasticity
	P.elas<-(P/h)*sens/lam
	eP=sum(P.elas)*h^2

	# reproduction elasticity
	R.elas<-(R/h)*sens/lam
	eR=sum(R.elas)*h^2
	
	# Net reproductive rate from IPMbook monocarp
	N <- solve(diag(n)-P)
	R0 <- abs(eigen(R %*% N)$values[1])
	
	# Generation time
	GT <- log(R0)/log(lam)
	
	demovals <- rbind(demovals, data.frame(spp=sp, eR, eP, R0, GT))
	sizevals <- rbind(sizevals, data.frame(spp=sp, area=y, stable.size, repro.val, v.dot.w))
		s.list[[sp]] <- sens
		eR.list[[sp]] <- R.elas
		eP.list[[sp]] <- P.elas
		eK.list[[sp]] <- K.elas
		}

summary(s.list)

demovals
head(sizevals)
params[,colnames(demovals)]<- demovals[match(demovals$spp, params$spp),]

#sum(analyses()$eK)*h^2
lapply(eK.list, function(x){sum(x)*h^2})
rowSums(demovals[,c("eP","eR")])

plot_grid(
ggplot(sizevals, aes(x=area, y=repro.val, col=spp))+geom_line()+
scale_colour_manual(values=cols)+
scale_y_log10()+guides(col="none"),
ggplot(sizevals, aes(x=area, y=stable.size, col=spp))+geom_line()+
scale_colour_manual(values=cols)+guides(col="none"),
ggplot(params, aes(reorder(spp, -eR), eR, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none"), 
ggplot(params, aes(reorder(spp, -R0), R0, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none")+scale_y_log10(),
ggplot(params, aes(reorder(spp, -GT), GT, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none")+scale_y_log10()
)


#######################################
# R.ELASTICITY across sizes
#######################################

#spp2 <- spp[!spp=="Gre"]
spp2 <- spp
ek.hist <- NULL
for (sp in spp2) {
	#sp <- "Gpe"
	y <- y.list[[sp]]
	e.k <- reshape(data.frame(t(eR.list[[sp]])),
	direction="long", varying=list(c(1:100)), ids=y, times=y,
	new.row.names=c(1:10000))
	e.k.sum <- aggregate(X1~id, e.k, sum)
	ek.hist <- rbind(ek.hist, cbind(e.k.sum, sp))
	}
head(ek.hist)
colnames(ek.hist)[3]<-"spp"


plot_grid(
ggplot(ek.hist, aes(id, X1))+
geom_line(aes(col=sp), size=0.5)+
scale_y_sqrt()+
guides(col="none")+
xlim(c(-4,1))+
geom_vline(xintercept=log10(pi*((15/2)/100)^2))+
scale_colour_manual(values=cols)
,
ggplot(sdat, aes(area))+
geom_density(aes(col=spp))+
guides(col="none")+
xlim(c(-4,1))+
scale_colour_manual(values=cols),
ncol=1)

ggplot()+
geom_density(data=sdat, aes(x=area, fill=spp))+
#scale_y_sqrt()+
geom_line(data=ek.hist, aes(id, X1/10, col=sp), size=0.5)+
guides(col="none")+
xlim(c(-4,1))+
facet_wrap(~spp)+
geom_vline(xintercept=log10(pi*((15/2)/100)^2))+
scale_colour_manual(values=cols)+
scale_fill_manual(values=cols)


#######################################
# R.ELASTICITY vs LAMBDA
#######################################
store[store$rec==max(rec.x),]
#params$lam_4 <- store[store$rec==max(rec.x),"lam"]
#params$lam_1 <- store[store$rec==min(rec.x),"lam"]
params$lam_4 <- store[store$rec==0.01,"lam"]
params$lam_1 <- store[store$rec==0.0001,"lam"]


  ggplot(data=params, aes(x=eR, y=log(lam_4)-log(lam_1)))+
geom_smooth(method="lm", formula=y~x, size=0.1, col="black")+
 geom_point(aes(fill=spp), shape=21, stroke=0.1, size=2)+
 scale_fill_manual(values=cols)+
 guides(fill="none", col="none")
#Pvals <- big$P / h
#P.elas <- Pvals * K.sens / lambda

#######################################
# Generation VS R0
#######################################

ggplot(demovals, aes(GT, R0))+geom_text(aes(label=spp))+scale_x_log10()
# could calculate mean size of flowering plants. *stable size by p_maturity



#source("R/fig.3.R")
#fig.3

ggplot(params, aes(eR, log(lam.est)))+geom_smooth(method="lm")+geom_hline(yintercept=0)+geom_point(aes(size=abundance_05, col=morphology))+geom_text(label=spp)


#######################################
# WITHIN VS BETWEEN
#######################################


params

dems <- c("lam_1","lam_4", "eR", "lam.est")

vars <- NULL
for (t in c(traits, "eggC", dems, "abundance_05")){
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
vars$type <- ifelse(vars$t %in% dems, "demo", ifelse(vars$t=="abundance_05", "abun", "trait"))

ggplot(vars, aes(y=reorder(t, -within), x=within))+geom_point(aes(col=type))+labs(y="% variation within morphologies", x="parameter")+
coord_cartesian(xlim=c(0,55))

params[,c("abundance_pair","spp")]

vars2 <- NULL
for (t in c(traits, "eggC", dems)){
	#t <- "f.cm2"
	abun <- "abundance_pair"
	params$abun <- params[,abun]
	params$t <- scale(params[,t])
mod <- aov(t~abun+morphology, data=params)
summary(mod)
within <- summary(mod)[[1]]["Residuals","Sum Sq"]
between <- summary(mod)[[1]]["abun","Sum Sq"]
morph <- summary(mod)[[1]]["morphology","Sum Sq"]
pval <- summary(mod)[[1]]["abun","Pr(>F)"]
vars2 <- rbind(vars2, data.frame(t, within=100*(within/(within+between+morph)), between=100*(between/(within+between+morph)), pval))
}
vars2$type <- ifelse(vars2$t %in% dems, "demo", ifelse(vars2$t=="abundance_05", "abun", "trait"))
vars2
ggplot(vars2, aes(reorder(t, -between), between))+geom_point(aes(col=type))+coord_flip()+labs(y="% variation between", x="parameter")


vars2$within_morph <- vars$within[match(vars2$t, vars$t)]
ggplot(vars2, aes(within_morph, between))+geom_text(aes(label=t, col=type))


#######################################
# MODEL ABUNDANCE
#######################################


vars3 <- NULL
for (t in c(traits, "eggC", dems)){
	#t <- "f.cm2"
	abun <- "abundance_05"
	params3 <- params
	params3$r.int <- log(params$r.int)
	params3$abun <- log(params3[,abun])
	params3$t <- scale(params3[,t])
	mod <-lm(abun ~ t * r.int, data=params3) 
    summary(mod)
   coefs<-coef(mod)
	se.x <- coef(summary(mod))["t","Std. Error"]
	pval.x <- coef(summary(mod))["t","Pr(>|t|)"]
	vars3 <- rbind(vars3, data.frame(x=coefs[2], se.x, pval.x, sd.resid=sd(residuals(mod, type="response")), t=t, row.names=NULL))
}
vars3$type <- ifelse(vars3$t %in% dems, "demo", ifelse(vars3$t=="abundance_05", "abun", "trait"))
vars3$e.size <- vars3$x
vars3$e.se <- vars3$se.x *1.96 
vars3$sig <- ifelse(vars3$pval.x < 0.05, '*', "")
vars3

ggplot(vars3, aes(e.size/sd.resid, reorder(t, -e.size/sd.resid)))+
geom_vline(xintercept=0)+
geom_point(aes(col=type))+labs(x="effect size", y="parameter")

#######################################
# MODEL ABUNDANCE II
#######################################


vars3 <- NULL
for (t in c(traits, "eggC", dems)){
	#t <- "f.cm2"
	params3 <- params
	#params3$r.int <- log(params$r.int)
	params3$abun <- params3$abundance_pair
	params3$abun <- ifelse(params3$abun=="Common", 1, 0)
	params3$t <- scale(params3[,t])
	mod <-glm(abun ~ t, data=params3, family="binomial") 
    summary(mod)
   coefs<-coef(mod)
	se.x <- coef(summary(mod))["t","Std. Error"]
	pval.x <- coef(summary(mod))["t","Pr(>|z|)"]
	vars3 <- rbind(vars3, data.frame(x=coefs[2], se.x, pval.x, sd.resid=sd(residuals(mod, type="response")), t=t, row.names=NULL))
}
vars3$type <- ifelse(vars3$t %in% dems, "demo", ifelse(vars3$t=="abundance_05", "abun", "trait"))
vars3$e.size <- vars3$x
vars3$e.se <- vars3$se.x *1.96 
vars3$sig <- ifelse(vars3$pval.x < 0.05, '*', "")
vars3

ggplot(vars3, aes(e.size/sd.resid, reorder(t, -e.size/sd.resid)))+
geom_vline(xintercept=0)+
geom_point(aes(col=type))+labs(x="effect size", y="parameter")






#######################################
# PROJECTING LAMDA
#######################################

# Keep the same lams... some pops go to zero
# impossible to get the same diffs 
# but can get the same poportion diffs 
# but some in small poulations

n1 <- 100
ngen <- 50
pars <- params2[, c("spp","morph", "abundance_pair","lam.est")]
pars

proj <- merge(pars, data.frame(gen=1:ngen))
proj$n <- n1*proj$lam.est^proj$gen
proj$id <- paste(proj$morph, proj$gen, sep="_")
proj$abun <- params$abundance_05[match(proj$spp, params$spp)]

ggplot(proj, aes(gen, n))+geom_line(aes(group=spp, col=spp))+scale_colour_manual(values=cols)+facet_wrap(~morph, scales="free_y")

a_diff <- unique(pairs[,c("morph","AC","AR")])
a_diff$diff <- a_diff$AC / a_diff$AR

proj2 <- data.frame(acast(morph+gen~abundance_pair, data=proj, value.var="n"))
proj2[, c("morph","gen")]<- proj[match(rownames(proj2), proj$id), c("morph","gen")]
proj2$diff <- proj2$Common/proj2$Rare
proj2$real <- a_diff$diff[match(proj2$morph, a_diff$morph)]
head(proj2)


ggplot(subset(proj2, gen==ngen), aes( real, diff, col=morph))+geom_point()+scale_y_log10()+scale_x_log10()+guides(col="none")+scale_colour_manual(values=colsC)

proj3 <- proj2
proj3 <- proj3[proj3$diff<proj2$real,]
gens <- aggregate(gen~morph, proj3, max)
proj4 <- proj2
proj4$max <- gens$gen[match(proj4$morph, gens$morph)] +2 #to cross the line
proj4 <- proj4[proj4$gen<proj4$max,]


ggplot(proj4, aes(gen, diff))+geom_line(aes(col=morph))+
geom_segment(aes(x=-Inf, xend=Inf, y=real, yend=real, col=morph), linetype="dotted")+
#scale_y_log10()+
scale_colour_manual(values=colsC)

#######################################
# HIDDEN DIMENSION
#######################################


# size at maximum survival 
max.surv<-do.call(rbind, lapply(split(s.pred, as.factor(s.pred$spp)), function(x) {return(x[which.max(x$pred),])}))
params$safe.size2<-max.surv$area
rec.cm <- 11
params$survcm2<-aggregate(pred~spp, s.pred[s.pred$area<log10(pi*(rec.cm/100/2)^2),], FUN=mean)$pred

# recruit survival
min.S <- aggregate(area~morph, aggregate(area~spp+morph, s.pred, min), max)
s.pred$min.s <- min.S$area[match(s.pred$morph, min.S$morph)]
max.S <- aggregate(area~morph, aggregate(area~spp+morph, s.pred, max), min)
s.pred$max.s <- max.S$area[match(s.pred$morph, min.S$morph)]
quarts <- aggregate(area~morph, s.pred, FUN=function(x){quantile(x, 0.2 )})
s.pred$quarts <- quarts$area[match(s.pred$morph, quarts$morph)]
s.pred$less <-ifelse(s.pred$area < s.pred$quarts & s.pred$area>s.pred$min.s & s.pred$area < s.pred$max.s, "y", "n")
head(s.pred)
ggplot()+
geom_line(data=s.pred, aes(area, pred, group=spp), size=0.1, col="grey")+
geom_line(data=s.pred[s.pred$less=="y",], aes(area, pred, col=spp))+scale_colour_manual(values=cols)+guides(col="none")
rec.surv<- do.call(rbind, lapply(split(s.pred, as.factor(s.pred$spp)), function(x) {
	x[x$less=="y",]
	}))
head(rec.surv)
avrec<-aggregate(pred~spp, rec.surv, mean)
params$rec.surv <- avrec$pred[match(params$spp, avrec$spp)]

# minimum at reproductive maturity
params$min.r2<-aggregate(area_cm2~spp, fec[fec$reproductive==1,], min)$area_cm2

# size 50% maturity (how delayed)
mat50 <- aggregate(area~spp, m.pred[m.pred$pred>0.5,], min)
params$mat50<-mat50$area[match(params$spp, mat50$spp)]

# barplots
bars <- function(t){
	params$x <- params[,t]
ggplot(params, aes(x=reorder(spp, -x), y=x, fill=spp))+guides(fill="none")+
geom_bar(stat="identity")+
ylab(t)+ggtitle(t)+
theme(axis.title.x=element_blank())+
facet_wrap(~morphology, scales="free_x", nrow=1)+
scale_fill_manual(values=cols)}

plot_grid(
bars("eggC"),
bars("f.slp")+coord_cartesian(ylim=c(2.3,2.8)),
ncol=1)

# lms
lm.plot<-function(x, y, n){
	dat<-params
	dat$x<-params[,x]
	dat$y<-params[,y]
    ggplot(dat, aes(x,y))+	
    geom_smooth(method="lm", formula=y~poly(x,n))+
    geom_path(aes(group=morphology), linetype="dotted")+
	geom_point(aes(fill=spp, size=abun), shape=21)+
	geom_text(aes(label=spp), size=1.5)+
	labs(x=x, y=y)+
	scale_fill_manual(values=cols)+guides(fill="none")+	theme(axis.text=element_text(size=6),axis.title=element_text(size=10), plot.title=element_text(size=7))
}

plot_grid(
lm.plot("mat50","f.slp", 1)+
ggtitle("delayed maturity = delayed eggs")+
ylab("higher slope = more delayed"),
lm.plot("f.slp", "rec.surv", 2)+scale_y_log10()+
ggtitle("delayed eggs = juvenile survival"),
lm.plot("eggC", "f.slp", 1),
lm.plot("f.slp", "f.cm2", 1),
lm.plot("r.int", "f.slp", 1)
)


# pca
head(params)
params$safe.sizex<-10^params$safe.size
pca2<-prcomp(log10(params[,c("f.cm2", "f.slp", "eggC", "eR",  "lam.est", "Sgrow")]), scale=T, center=T)
biplot(pca2)
params$Rpca1 <- pca2$x[,1][match(params$spp, rownames(pca2$x))]
params$Rpca2 <- pca2$x[,2][match(params$spp, rownames(pca2$x))]
Rvecs <- data.frame(pca2$rotation, lab=rownames(pca2$rotation))


ggplot()+
geom_point(data=params, aes(Rpca1, Rpca2, col=spp, size=abundance_05))+
geom_segment(data=Rvecs, aes(0,0, xend=PC1, yend=PC2))+
geom_text(data=Rvecs, aes(PC1*1.2, PC2*1.2, label=lab))+
scale_colour_manual(values=cols)+guides(col="none", size="none")

ggplot(params, aes(eggC, Sgrow))+geom_point(aes(col=spp, size=abundance_05))+scale_colour_manual(values=cols)+guides(col="none", size="none")+geom_smooth(method="lm")


gdat$growth <- 10^gdat$area_next - 10^gdat$area
params$av.growth <- aggregate(growth~spp, gdat, mean)$growth
params$av.pmort <- aggregate(p_mort~spp, gdat, mean)$p_mort
gr.plots<-plot_grid(
lm.plot("r.int","min.r",1)+scale_x_log10()+scale_y_log10(),
lm.plot("r.int","p_mort",1)+scale_x_log10(), 
lm.plot("r.int","safe.size",1)+scale_x_log10(),
lm.plot("r.int","survcm",2)+scale_x_log10(),
lm.plot("r.int","av.surv",2)+scale_x_log10(),
lm.plot("r.int","f.cm2",2)+scale_x_log10(),
lm.plot("r.int","eggC",2)+scale_x_log10(),
lm.plot("r.int","av.growth",1)+scale_x_log10(),
lm.plot("r.int","av.pmort",1)+scale_x_log10(),
align="hv")
gr.plots




# escape in size - see Bak and Meesters 1998, Trapon MEPS, Babcock and mundy 96 






  ############ SIZE VS NUMBER 
  
  plot_grid(
ggplot(params, aes(size.av, abundance_05))+geom_text(aes(label=spp))+
geom_path(aes(group=morphology), linetype="dotted")+
geom_point(aes(col=spp))+scale_colour_manual(values=cols)+
scale_y_log10()+guides(col="none")
,
ggplot(params, aes(eggC, f.slp))+geom_text(aes(label=spp))+
geom_path(aes(group=morphology), linetype="dotted")+
geom_point(aes(col=spp))+scale_colour_manual(values=cols)+
scale_y_log10()+guides(col="none")
)


