
rm(list = ls())
library("MASS") 
library("reshape2")
library("ggplot2")
library("cowplot")
library("viridis")
library("png")
theme_set(theme_cowplot())
inv.logit <- function(x) {exp(x)/(1+exp(x))}
#source("data/data_prep.R")

#######################################
# DATA
#######################################
setwd("~/Documents/PostDoc/01_Trimodal/analysis_mm")
params<-read.csv("data/final/00_info.csv")
dat<-read.csv("data/final/01_areas.csv")
gdat<-read.csv("data/final/02_growth.csv")
sdat<-read.csv("data/final/03_surv.csv")
fec<-read.csv("data/final/04_fecundity.csv")
ss<-read.csv("data/final/05_size_structure.csv")
spp<-params$spp[order(params$spp)]
params$abun<-factor(params$abun, levels=c("Rare", "Common", "Dominant"))
cols<-as.character(params$cols)


#######################################
# FECUNDITY - PREP
#######################################
# Aggregate fecundity
fec<-aggregate(.~id+spp_code+morphology, fec[,c("id","area_cm2","eggs","spp_code", "polyps_cm2","eggC","morphology")], mean)
head(fec)
# --------
fec$reproductive <- ifelse(fec$eggs>0, 1, 0) # reproductive status
fec$area <- log10(fec$area_cm2 /10000) 
# --------
fec$f.cm2<-fec$polyps_cm2*fec$eggs
fec$f.colony<-fec$f.cm2*fec$area_cm2 *fec$eggC # original
fec$f.colonyP<-round(fec$f.cm2*fec$area_cm2 * (fec$eggC/max(fec$eggC)) ,0) 
# make energy a survival proportion, and round to 0
fec$fecundity<-fec$f.colonyP # select here



#######################################
# MATURITY 
#######################################
m.int<-NULL
m.slp<-NULL
m.pred<-NULL
for(sp in spp){
	sub<-fec[fec$spp_code==sp,]
	m.mod<-glm(reproductive ~ area, family="binomial", data=sub) 
	m.int<-c(m.int, coef(m.mod)[1])
	m.slp<-c(m.slp, coef(m.mod)[2])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  spp_code=sp, morph=sub$morphology[1])
	new$pred<-predict(m.mod, new, type="response")
	m.pred<-rbind(m.pred, new)
}
params$m.int<-m.int
params$m.slp<-m.slp
# -------- plot
ggplot()+
  geom_jitter(data=fec, 
    aes(x= area, reproductive), 
      height=0.02, shape=21, col="grey")+
  geom_line(data=m.pred, 
    aes(area, pred, col=spp_code), 
      size=1)+
  scale_colour_manual(values=cols)


#######################################
# FECUNDITY  
#######################################
f.int<-NULL
f.slp<-NULL
f.pred<-NULL
for(sp in spp){
	sub<-fec[fec$spp_code==sp & fec$reproductive==1,]
	f.mod<-glm.nb(fecundity ~ area, sub)
	f.int<-c(f.int, coef(f.mod)[1])
	f.slp<-c(f.slp, coef(f.mod)[2])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1),
	spp_code=sp, morph=sub$morphology[1])
	new$pred<-predict(f.mod, new, type="response")
	f.pred<-rbind(f.pred, new)
	}
params$f.int<-f.int
params$f.slp<-f.slp
# -------- plot
ggplot()+ 
  geom_point(data=fec[fec$reproductive==1, ], 
    aes(area, fecundity), 
      shape=21, col="grey")+
  geom_line(data=f.pred, 
    aes(area, pred, col=spp_code), 
      size=1)+
  scale_y_log10()+
  scale_colour_manual(values=cols)


#######################################
# GROWTH 
#######################################
gdat$area <- log10(gdat$area1 / 10000) 
gdat$area_next <- log10(gdat$area2 / 10000) 

g.int<-NULL
g.slp<-NULL
g.var<-NULL
g.pred<-NULL
for(sp in spp){
	sub<-gdat[gdat$spp_code==sp,]
	g.mod<-lm(area_next ~ area, data=sub) 
	g.int<-c(g.int, coef(g.mod)[1])
	g.slp<-c(g.slp, coef(g.mod)[2])
	g.var<-c(g.var, var(residuals(g.mod)))
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1),
	  spp_code=sp, morph=sub$morph[1])
	new$pred<-predict(g.mod, new, type="response")
	g.pred<-rbind(g.pred, new)
}
params$g.int<-g.int
params$g.slp<-g.slp
params$g.var<-g.var

ggplot()+ 
  geom_abline(slope=1, linetype="dotted")+
  geom_point(data=gdat, 
    aes(area, area_next), 
      shape=21, col="grey")+
  geom_line(data=g.pred, 
    aes(area, pred, col=spp_code), 
      size=1)+
  scale_colour_manual(values=cols)+
  facet_wrap(~spp_code)
  
  
 
dorn<-read.csv("data/other/dornelas2017/dornelas2017.csv")
trysp<-"GR"
ggplot()+
  geom_text(data=gdat[gdat$spp_code==trysp ,], 
    aes(x=log(10^area*10000), y=log(10^area_next*10000), label=colony_id))+
  geom_text(data=dorn[dorn$species_code==trysp ,], 
    aes(x=log(area_timestep), y=log(area_nexttimestep), label=colony_id), 
      col="red")+
  geom_smooth(data=dorn[dorn$species_code==trysp ,], 
    aes(x=log(area_timestep), y=log(area_nexttimestep)),
       col="red", se=F, method="lm", formula=y~x)+
  geom_smooth(data=gdat[gdat$spp_code==trysp ,], 
    aes(x=log(10^area*10000), y=log(10^area_next*10000)), 
      se=F, method="lm", formula=y~x)+
  geom_abline(slope=1)




#######################################
# SURVIVAL 
#######################################
sdat$area<-log10(sdat$area_cm2/10000)
sdat$area_sq<-sdat$area^2
#sdat<-sdat[sdat$area>-3,]
sdat<-sdat[!is.na(sdat$surv),]
#sdat<-sdat[!sdat$year==2012,]
#d.cm2(-3)

#sdat<-sdat[sdat$area<quantile(sdat$area, 0.995) & sdat$area>quantile(sdat$area, 0.005),]

s.int<-NULL
s.slp<-NULL
s.slp2<-NULL
s.pred<-NULL
for(sp in spp){
	sub<-sdat[sdat$spp_code==sp,]
	#s.mod<-glm(surv ~ area, family="binomial", data=sub) 
	s.mod2<-glm(surv ~ area + area_sq, family="binomial", data=sub) 
	s.int<-c(s.int, coef(s.mod2)[1])
	s.slp<-c(s.slp, coef(s.mod2)[2])
	s.slp2<-c(s.slp2, coef(s.mod2)[3])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  spp_code=sp, morph=sub$morph[1])
	new$area_sq <- new$area^2
	new$pred<-predict(s.mod2, new, type="response")
	s.pred<-rbind(s.pred, new)
	}
params$s.int<-s.int
params$s.slp<-s.slp
params$s.slp.2<-s.slp2

ggplot()+ 
  geom_jitter(data=sdat, 
    aes(area, surv), 
      shape=21, col="grey",height=0.02)+
  geom_line(data=s.pred, 
    aes(area, pred, col=spp_code), 
      size=0.75)+
      #facet_wrap(~spp_code)+
  scale_colour_manual(values=cols)




# JM Survival:
# lots of 2012 data missing, especially deaths:
# counted F6 (final timepoint) as 1, not NA 
# areas are associated with "Previous year" not current year
# some wrong survivals ("photographed"):
  # -   287 F5 - small AS
  # -   115 F5 - small AI
# some wrong deaths e.g. huge AC 29 in F5
madin<-read.csv("data/other/madin2014/madin2014.csv")
madin$area<-log10(exp(madin$log_area_cm2))
madin$area_sq<-madin$area^2
madin$surv<-ifelse(madin$mortality==1, 0, 1)
mad.p<-NULL
s2.int<-NULL
s2.slp<-NULL
s2.slp2<-NULL
for(sp in spp){
	sub<-madin[madin$spp_code==sp,]
	#s.mod2<-glm(surv ~ area, family="binomial", data=sub) 
	s.mod2<-glm(surv ~ area + area_sq, family="binomial", data=sub) 
	s2.int<-c(s2.int, coef(s.mod2)[1])
	s2.slp<-c(s2.slp, coef(s.mod2)[2])
	s2.slp2<-c(s2.slp2, coef(s.mod2)[3])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), spp_code=sp)
	new$area_sq <- new$area^2
	new$pred<-predict(s.mod2, new, type="response")
	mad.p<-rbind(mad.p, new)
}
 # USE MADIN PARAMS
#params$s.int<-s2.int
#params$s.slp<-s2.slp
#params$s.slp.2<-s2.slp2

ggplot()+
geom_vline(xintercept=-3)+
  geom_jitter(data=madin, aes(area, surv), shape=21, col="grey", height=0.02)+
    #geom_text(data=madin2, aes(area, surv, label=colony_id),  size=2, position=position_jitter(height=0.03))+
  geom_text(data=sdat, aes(area, surv, label=colony_id), 
      col="red",  size=2, position=position_jitter(height=0.03))+
  geom_line(data=mad.p, aes(area, pred))+
  #geom_line(data=mad2.p, aes(area, pred), linetype="dotted")+
  geom_line(data=s.pred, aes(area, pred), col="red")+
  facet_wrap(~spp_code)
  
  
  
spx<-"GR"
ggplot()+
geom_vline(xintercept=-3)+
  geom_jitter(data=madin[madin$spp_code==spx,], aes(area, surv), shape=21, col="grey", height=0.01)+
  geom_text(data=madin[madin$spp_code==spx,], aes(area, surv, label=colony_id),  size=2, position=position_jitter(height=0.03))+
  geom_text(data=sdat[sdat$spp_code==spx,], aes(area, surv, label=colony_id), col="red",  size=2, position=position_jitter(height=0.01))+
    geom_point(data=sdat[sdat$spp_code==spx & sdat$year==2012,], aes(area, surv),  size=2, shape=21)+
  geom_line(data=mad.p[mad.p$spp_code==spx,], aes(area, pred))+
  geom_line(data=s.pred[s.pred$spp_code==spx,], aes(area, pred), col="red")
# numbers: 
nrow(madin[madin$spp_code==spx & madin$surv==1,])
nrow(sdat[sdat$spp_code==spx & sdat$surv==1,])
nrow(madin[madin$spp_code==spx & madin$surv==0,])
nrow(sdat[sdat$spp_code==spx & sdat$surv==0,])



#######################################
# MAXIMUM GROWTH  - look into AD
#######################################
library("quantreg")
gdat$radius1<-sqrt(gdat$area1/pi)/100
gdat$radius2<-sqrt(gdat$area2/pi)/100
gdat$g_radius<-gdat$radius2-gdat$radius1

r.int<-NULL
r.slp<-NULL
r.err<-NULL
r.pred<-NULL
for(sp in spp){
	sub<-gdat[gdat$spp_code==sp,]
	r.mod <-rq(g_radius ~ 1 , data=sub, tau=0.98) # no slope
	#r.mod <-rq(g_radius ~ area , data=sub, tau=0.99)# slope
	r.int<-c(r.int, coef(r.mod)[1])
	#r.slp<-c(r.slp, coef(r.mod)[2])
	r.err <- c(r.err, summary(r.mod, se='boot')$coef[[2]])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  spp_code=sp, morph=sub$morph[1])
	new$pred<-predict(r.mod, new, type="response")
	r.pred<-rbind(r.pred, new)
	}
params$r.int<-r.int
params$r.slp<-r.slp
params$r.err<-r.err


ggplot()+
 geom_point(data=gdat, aes(x=area, y=g_radius), col="grey", shape=21)+
     geom_line(data=r.pred, aes(x=area, y=pred, col=spp_code))+
   scale_colour_manual(values=c(cols))#+facet_wrap(~spp_code)





#######################################
# PARTIAL MORTALITY  - subset??
#######################################
logit <- function(x) { log(x/(1-x)) }
inv.logit <- function(x) { exp(x)/(1+exp(x)) }

gdat$max_g<-params$r.int[match(gdat$spp_code, params$spp)] 
gdat$max_area2<-pi*(gdat$max_g + gdat$radius1)^2 # eq. 2.1 
gdat$p_mort<-  1 - ((gdat$area2/10000)/gdat$max_area2) # eq. 2.2
gdat$p_stasis<-  1 - ((gdat$area1/10000)/gdat$max_area2) 
pdat<-subset(gdat, p_mort>0.001) # less than max growth line
#pdat<-subset(dat, p_mort>0) # zero
pdat$pm_logit <- logit(pdat$p_mort) 

p.int<-NULL
p.slp<-NULL
p.pred<-NULL
p.sig<-NULL
for(sp in spp){
	sub<-pdat[pdat$spp_code==sp,]
	p.mod<-lm(pm_logit ~ area , data=sub)
	p.int <- c(p.int, coef(p.mod)[[1]])
	p.slp <- c(p.slp, coef(p.mod)[[2]])
	p.sig <- c(p.sig, sigma(p.mod))
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  spp_code=sp, morph=sub$morph[1])
	new$pred<-predict(p.mod, new, type="response")
	p.pred<-rbind(p.pred, new)
	}
params$p.int<-p.int
params$p.slp<-p.slp
params$p.sig<-p.sig

ggplot()+
  geom_point(data=pdat, 
    aes(x=area, y=pm_logit), 
       shape=21, col="grey")+
    geom_line(data=pdat, aes(x=area, y=p_stasis), linetype="dotted")+
   geom_line(data=p.pred, 
     aes(x=area, y=pred))+facet_wrap(~spp_code)




#######################################
# MIN/MAX SIZE
#######################################
head(ss)
params$smax<-aggregate(area_cm2~spp_code, dat, max)$area
params$smin<-aggregate(area_cm2~spp_code, dat, min)$area
#params$smax<-apply(aggregate(list(a=dat$area, b=dat$area_next), by=list(dat$spp_code),max)[,c("a","b")], 1, max)
#params$smin<-apply(aggregate(list(a=dat$area, b=dat$area_next), by=list(dat$spp_code),min)[,c("a","b")], 1, min)




##################################################
#### ######  ##  ##   ####
 ##  ##  ## ## ## ##  ##
 ##  #####  ## ## ##    ##
 ##  ##     ## ## ##     ##
#### ##     ## ## ##  ####
##################################################
options(stringsAsFactors=FALSE)

#------------------------------- growth
g.yx <- function(y, x) {
	dnorm(y, mean=params$g.int[params$spp==sp] + 
	  params$g.slp[params$spp==sp]*x,
	    sd=sqrt(params$g.var[params$spp==sp]))
	    }
	
#-------------------------------survival
s.x <- function(x) { 
	u <- params$s.int[params$spp==sp] + 
	  params$s.slp[params$spp==sp] * x + 
	    params$s.slp.2[params$spp==sp] * x^2
  return(inv.logit(u)) 
  }
  
#------------------------------- reproduction
 r.yx <- function(y, x) {	
 	mat<- inv.logit(params$m.int[params$spp==sp] + 
 	  params$m.slp[params$spp==sp] *x)
 	fec<- exp(params$f.int[params$spp==sp] + 
 	  params$f.slp[params$spp==sp] *x) 
 	#siz<- dnorm(y,mean=log10(pi*params$r.int[params$spp==sp]^2), sd=0.05) 
 	#siz<- dnorm(y,mean=rec.size, sd=0.05) 
 	siz <- dnorm(y,mean=log10(0.6*pi*params$r.int[params$spp==sp]^2), sd=params$g.var[params$spp==sp])
   out <- rec* mat * fec * siz 
   return(out)
   } 

# growth functions
a_func <- function(r) {
  pi * r^2
}
r_func <- function(a) {
  sqrt(a / pi)
}
circularity <- function(area, perimeter) {
  (4 * pi * area)/(perimeter^ 2)
}
   
#------------------------------- growth & partial morality
p.yx <- function(y, x) {
  # x <- -4
  g <- a_func(r_func(10^x) + params$r.int[params$spp==sp] )
  #+ 1.96 * params$r.err[params$spp==sp])
  
  temp <- 10^y / g
  temp[temp > 1] <- 1
  dnorm(logit(1 - temp), params$p.int[params$spp==sp] + x * params$p.slp[params$spp==sp], params$p.sig[params$spp==sp])
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
  #R <- h * outer(y, y, r.yx) 
  R <- h * outer(y, pmin(y, smax), r.yx) #  ceiling
  K <- P + R
  lam <- Re(eigen(K)$values[1])
	w <- abs(Re(eigen(K)$vectors[,1])) 
	v <- abs(Re(eigen(t(K))$vectors[,1]))
	return(list(K=K, lam=lam, w=w, v=v, G=G, S=S, R=R, P=P)) }
		
	
#######################################
# MESH AND BOUNDARIES
#######################################

min.size <- -3.5 
max.size <- 1 
rec.size <- -2.21 
n <- 100
b <- min.size + c(0:n) * (max.size - min.size)/n
y <- 0.5 * (b[1:n]+b[2:(n+1)])
h <- y[2] - y[1]
I <- y >= rec.size

#-------------------------------
d.cm2<-function(x){ (sqrt((10^x)*10000)/pi)*2 }
d.cm2(rec.size) # rec diam in cm
d.cm2(log10(pi*params$r.int[params$spp=="AH"]^2))


#######################################
# FECUNDITY VS RECRUITMENT...
#######################################

rmat<-NULL
  for (sp in spp) {
  	rec<-1
  	smax<-params[params$spp==sp, "smax"]
rmat<-rbind(rmat,data.frame(s=y, rep=apply(t(bigmatrix()$R), 1, sum), spp_code=sp))}
head(rmat)
ggplot()+
  geom_line(data=rmat, aes(s, rep, col=spp_code), linetype="dotted")+
  geom_line(data=f.pred,aes(area, pred, col=spp_code))+
  scale_y_log10()+
  facet_wrap(~spp_code)+
  scale_colour_manual(values=cols)


#######################################
# PLOT IPMS
#######################################

par(mfcol=c(2, 6))
lam_const <- NULL
t.rec <- NULL
for (sp in spp) {
	smax<-params[params$spp==sp, "smax"]
    rec <-     1*10^-8
  	sub<-dat[dat$spp_code==sp,]
    mod <- bigmatrix()
	image(y, y, t(mod$K)^0.3)    
	points(sub$area, sub$area_next, cex=0.25)
	title(sp, line=-1)
	abline(0, 1, lty=2)
  lam_const<-c(lam_const,bigmatrix()$lam) }
params$lam_const<-lam_const
params$lam_const


#######################################
# MULTIPLE REC/LAMBDAS
#######################################

store<-data.frame()
recx<-c(1*10^-3, 1*10^-4, 1*10^-5, 1*10^-6, 1*10^-7, 1*10^-8)
for (sp in spp) {
	for (rec in recx){
	smax<-params[params$spp==sp, "smax"]
  	sub<-dat[dat$spp_code==sp,]
  	rec<-rec
   mod <- bigmatrix()
  store<-rbind(store, data.frame(sp=sp, rec=rec, lam=bigmatrix()$lam))
  } }
recdat<- dcast(store, sp~rec, value.var='lam')
params$lam_0<-recdat[,"1e-08"]
params$lam_1<-recdat[,"1e-07"]
params$lam_2<-recdat[,"1e-06"]
params$lam_3<-recdat[,"1e-05"]
params$lam_4<-recdat[,"1e-04"]

#######################################
# EXPLORE LAMBDA
#######################################
lam<-expression(lambda)
lamlong<-melt(params, measure.vars=c("lam_0","lam_1","lam_2","lam_3", "lam_4"), value.name="lambda")
#-------------------------------
params2<-rbind(params, params[6,]) #duplicate millepora(should it be nasuta?) 
params2$morph<-as.character(params2$morph)
params2$morph[c(7, 12)]<-c("corymbose_2","corymbose_2")
comp<-dcast(params2, morph~abunJM, value.var="spp")
comp$sp<-apply(comp[,c(2,3)], 1, paste, collapse = ":")
comp$lam_0<-transform(dcast(params2, morph~abunJM, value.var="lam_0"), x=Common/Rare)$x
comp$lam_1<-transform(dcast(params2, morph~abunJM, value.var="lam_1"), x=Common/Rare)$x
comp$lam_2<-transform(dcast(params2, morph~abunJM, value.var="lam_2"), x=Common/Rare)$x
comp$lam_3<-transform(dcast(params2, morph~abunJM, value.var="lam_3"), x=Common/Rare)$x
comp$lam_4<-transform(dcast(params2, morph~abunJM, value.var="lam_4"), x=Common/Rare)$x
comp<-melt(comp, id.vars=c("morph","sp","Common","Rare"))
#------------------------------- 
comp$morph<-factor(comp$morph, levels=c("tabular","staghorn", "corymbose", "corymbose_2", "digitate","massive"))
ggplot(comp)+
geom_bar(aes(x=morph, y=value),stat="identity", col="black", fill="grey")+
geom_text(aes(x=morph, y=value+0.01, label=sp), size=2, angle=90)+
coord_flip()+scale_y_log10()+
theme(axis.line.y=element_blank(), strip.background=element_blank(), axis.title.y=element_blank())+
ylab("lambda ratio")+
facet_wrap(~variable, nrow=2)


#######################################
# ELASTICITY
#######################################
 # eigen-things can be combined to obtain the sensitivity and elasticity matrices.
 
rec<-1*10^-6

analyses<-function(){
	#sp<-"GR"
	K<-bigmatrix()$K
	lam <- Re(eigen(K)$values[1])
  w.eigen<-Re(eigen(K)$vectors[,1]) # w=right eigen=stable size dist. 
  v.eigen<-Re(eigen(t(K))$vectors[,1]) # v=left, reproductive value=contribution of each size class to pop size. 
stable.size <- w.eigen/sum(w.eigen)
repro.val <- v.eigen/v.eigen[1] 
v.dot.w<-sum(stable.size*repro.val)*h 
sens<-outer(repro.val,stable.size)/v.dot.w  
elas<-matrix(as.vector(sens)*(as.vector(K)/h)/lam, nrow=n)
P.elas<-(bigmatrix()$P/h)*sens/lam
R.elas<-(bigmatrix()$R/h)*sens/lam
return(list(st=stable.size, rv=repro.val, sa=sens, eK=elas, eP=P.elas, eR=R.elas))}

sum(analyses()$eK)*h^2
  
props<-NULL 
for (sp in spp) {
props<-rbind(props, data.frame(eP=sum(analyses()$eP)*h^2, eR=sum(analyses()$eR)*h^2, sp=sp))
	}
props
params$eP<-props$eP
params$eR<-props$eR
rowSums(props[,c(1,2)])

#######################################
# R.ELASTICITY VS RECRUITMENT SENSITIVITY
#######################################

ggplot(data=params, aes(x=eR, y=lam_4-lam_1))+
  geom_text(aes(label=spp))+
  geom_smooth(method="lm", formula=y~x)+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="R.elasticity", y="change in lambda")


#######################################
# R.ELASTICITY vs MORPH/ABUNDANCE
#######################################

params$R<-ifelse(params$morph=="massive", "S", ifelse(params$morph=="digitate", "R", "R"))
mod<-lm(eR~log(abun_05), subset(params, R!="S"))
summary(mod)

plot_grid(
ggplot(data=params, aes(y=eR, x=abun_05))+
geom_smooth(method="lm", formula=y~poly(x,2))+geom_text(aes(label=spp))+
scale_x_log10()+
scale_y_log10()+
labs(y="R.elasticity")#+facet_wrap(~R)
,
ggplot(data=params, aes(x=reorder(morph, -eR), y=eR))+geom_boxplot()+
geom_point()+labs(y="R.elasticity")+
scale_y_log10()+
geom_text(aes(label=spp), nudge_x=0.15, col="red")
)


#Pvals <- big$P / h
#P.elas <- Pvals * K.sens / lambda

ggplot(params, aes(r.int, eR, col=spp))+
geom_line(aes(group=morph))+
geom_point(aes(size=abun_05))+
  scale_x_log10()+
  scale_y_log10()+
scale_colour_manual(values=cols)+scale_x_sqrt()#+facet_wrap(~R)







#######################################
# EXPLORE LAMBDA II
#######################################

ggplot(data=lamlong)+ 
geom_hline(yintercept=1, linetype="dotted")+
geom_line(aes(x=abun, y=lambda, col=spp, group=variable))+
 geom_point(aes(x=abun, y=lambda, fill=spp, size=variable), shape=21)+
 guides(fill="none", col="none")+labs(y=lam)+ scale_fill_manual(values=as.character(params$cols))+
 scale_colour_manual(values=as.character(params$cols))+
 facet_wrap(~morph, scales="free_x")+
 theme(axis.text.x=element_text(angle=90, size=5))


# Lamda abundance... 
pdat<-lamlong[lamlong$variable=="lam_3",]
#pdat<-lamlong[lamlong$variable!="lam_4",]
#pdat<-lamlong
plot_grid(ggplot(pdat,aes(x=abun_05, y=lambda, col=variable, fill=variable))+
geom_smooth(method="lm", formula=y~poly(x,2), size=0.5)+
#geom_hline(yintercept=1, linetype="dotted")+
geom_text(aes(label=spp), size=3)+
labs(x="abundance", y=lam)+theme(legend.title=element_blank(), axis.text=element_text(size=8))+
scale_x_log10(),
ggplot(pdat,aes(x=r.int, y=lambda, col=variable, fill=variable))+
geom_smooth(method="lm", formula=y~poly(x,2), size=0.5)+
#geom_hline(yintercept=1, linetype="dotted")+
geom_text(aes(label=spp), size=3)+
labs(x="growth", y=lam)+theme(legend.title=element_blank(), axis.text=element_text(size=8))+
scale_x_log10(),
nrow=1)



#######################################
#-------------------------------  min.size/max.size-lamda relationship
#######################################

store<-data.frame()
minx<-c(-6,-5,-4,-3,-2,-1)
for (sp in spp) {
	for(min.size2 in minx){
		min.size2<-min.size2
		max.size <- max.size #2
		rec.size <- rec.size
		n <- 100
		b <- min.size2 + c(0:n) * (max.size - min.size2)/n
		y <- 0.5 * (b[1:n]+b[2:(n+1)])
		h <- y[2] - y[1]
		store<-rbind(store, data.frame(sp=sp, min.size=min.size2, lam=bigmatrix()$lam)) }}
ggplot(data=store, aes(x=min.size, y=lam, colour=sp))+
geom_line()+geom_hline(yintercept=1)+
geom_vline(xintercept=min.size, linetype="dotted")+
ggtitle("relative differences")+
scale_colour_manual(values=as.character(params$cols))


store<-data.frame()
maxx <- c( 0,1,2,3, 4,5)
for (sp in spp) {
	for(max.size2 in maxx){
		max.size2<-max.size2
		min.size <- min.size #2
rec.size <- -2.7
n <- 100
b <- min.size + c(0:n) * (max.size2 - min.size)/n
y <- 0.5 * (b[1:n]+b[2:(n+1)])
h <- y[2] - y[1]
		store<-rbind(store, data.frame(sp=sp, max.size=max.size2, lam=bigmatrix()$lam)) }}
ggplot(data=store, aes(x=max.size, y=lam, colour=sp))+
geom_line()+geom_hline(yintercept=1)+
geom_vline(xintercept=max.size, linetype="dotted")+
#lims(x=c(0, 0.00001), y=c(0,4))+
ggtitle("relative differences")+
scale_colour_manual(values=as.character(params$cols))







  


#######################################
# TRADEOFFS 
#######################################
head(params)

# survival
params$p_mort<-inv.logit((params$p.slp*log10(0.01))+params$p.int)
params$av.surv<-aggregate(pred~spp_code, s.pred,mean)$pred
params$survcm<-aggregate(pred~spp_code, s.pred[s.pred$area<log10(pi*(11/100/2)^2),], FUN=mean)$pred
max.surv<-do.call(rbind, lapply(split(s.pred,as.factor(s.pred$spp)), function(x) {return(x[which.max(x$pred),])}))
params$safe.size<-max.surv$area

# growth
gdat$growth<-gdat$area2-gdat$area1	
params$av.growth<-aggregate(growth~spp_code, gdat[gdat$growth>0,], mean)$growth

# reproduction
params$f.cm2<-aggregate(f.cm2~spp_code, fec,mean)$f.cm2
params$en.cm2<-params$f.cm2*params$eggC
params$min.r<-aggregate(area_cm2~spp_code, fec[fec$reproductive==1,], min)$area_cm2

#params2<-params
#params2<-rbind(params2, rep(NA, ncol(params)))
#params2[12,"spp"]<-"CR"
#params2[12,"r.int"]<-0.2/1000
#params2[12,"min.r"]<-pi*((2/10/2)^2) # 1.2-2 mm diam


colnames(params)
pca<-prcomp(params[,c("r.int", "f.cm2","survcm","av.surv","p_mort")], scale=T, center=T)
#pca<-prcomp(params[,c("r.int", "f.cm2","av.surv", "p_mort")], scale=T, center=T)
pcdat<-data.frame(spp=params$spp, x=pca$x[,1], y=pca$x[,2], abun=params$abun_05)
rot<-data.frame(spp=rownames(pca$rotation), x=pca$rotation[,1], y=pca$rotation[,2])
exp<-round(c(summary(pca)[[1]][1]^2/sum(summary(pca)[[1]]^2),summary(pca)[[1]][2]^2/sum(summary(pca)[[1]]^2)),3)*100

ggplot()+
geom_segment(data=rot, aes(x=0,y=0, xend=x, yend=y), arrow=arrow(length=unit(1, "mm")))+
geom_text(data=rot, aes(x*1.2, y*1.2, label=spp))+
geom_point(data=pcdat, aes(x,y, col=spp, size=abun))+
geom_text(data=pcdat, aes(x,y, label=spp), size=3)+
scale_colour_manual(values=c(cols))+
guides(colour="none", size="none")+scale_radius(range=c(5,12))+
labs(x=paste("PC1 (",exp[1],"%)", sep=""),y=paste("PC2 (",exp[2],"%)", sep=""))+
theme_bw()





# linear models..

lm.plot<-function(x, y, n){
	dat<-params
	dat$x<-params[,x]
	dat$y<-params[,y]
    ggplot(dat, aes(x,y))+	
    geom_smooth(method="lm", formula=y~poly(x,n))+
	geom_point(aes(fill=spp), shape=21, size=3.5)+
	geom_text(aes(label=spp), size=1.5)+
	labs(x=x, y=y)+
	scale_fill_manual(values=paste(cols))+guides(fill="none")+
	theme(axis.text=element_text(size=6),axis.title=element_text(size=10))
}

# growth
plot_grid(
lm.plot("r.int","min.r",1)+scale_x_log10()+scale_y_log10(),
lm.plot("r.int","p_mort",1)+scale_x_log10(), 
lm.plot("r.int","safe.size",1)+scale_x_log10(),
lm.plot("r.int","survcm",2)+scale_x_log10(),
lm.plot("r.int","av.surv",2)+scale_x_log10(),
lm.plot("r.int","f.cm2",2)+scale_x_log10(),
lm.plot("r.int","en.cm2",1)+scale_x_log10(),
align="hv")
  



























##################################################
 ####  ###### ##  ## ##### #####
##  ##   ##   ##  ## ##    ##  ##
##  ##   ##   ###### ####  #####
##  ##   ##   ##  ## ##    ##  ##
 ####    ##   ##  ## ##### ##   ##
##################################################
#library("nlme") 
#library("stats4")
#library("boot")
#library("grid")
#library("gridExtra")
  
#######################################
# FECUNDITY II
#######################################
# Data:
# Every row is a POLYP in fec data
# Year > Species > Colony > Branch > Polyp > N eggs
# mostly 4-6 branches per colony (goni=1 branch)
# 6 polyps per branch (consistent)
# Current problems:
# If a polyp has no eggs, we assume the colony has no eggs (wrong)
# Currently model eggs per polyp consistent accross colony
# Eggs per polyp is poorly predicted by size
#######################################
fec$area <- log10(fec$area_cm2 /10000) 
#######################################
# Colony reproductive? 
head(fec)
colmat<-aggregate(.~id+spp_code, fec[,c("id","area","eggs","spp_code")], mean)
colmat$mat<-ifelse(colmat$eggs>0, 1,0)
head(colmat)

m.mod <- glm(mat ~ area * spp_code, family="binomial", colmat)
m.pred<-pred.sp(m.mod, colmat)
m.plot<-plot.sp(colmat, m.pred, "mat",0.02, "Colony maturity")
m.plot
int.sp(m.mod)

#######################################
# Number reproductive polyps?
head(fec)
fec$polypID<-paste(fec$id, fec$branch, fec$polyp, sep=".")
fec$colmat<-colmat$mat[match(fec$id, colmat$id)]
fec2<-subset(fec, colmat==1)
#fec2<-fec
polmat<-aggregate(.~polypID+spp_code, fec2[,c("polypID","area","eggs","spp_code")], mean)
polmat$mat<-ifelse(polmat$eggs>0, 1,0)
head(polmat)
m.mod2 <- glm(mat ~ area * spp_code, family="binomial", polmat)
m.pred2<-pred.sp(m.mod2, polmat)
m.plot2<-plot.sp(polmat, m.pred2, "mat",0.02, "Polyp maturity")
m.plot2
#######################################
# Number eggs per reproductive polyp?
fec2$polmat<-polmat$mat[match(fec2$polypID, polmat$polypID)] 
fec3<-subset(fec2, polmat==1)
polfec<-aggregate(.~id+spp_code, fec3[,c("id","area","eggs","spp_code")], mean)
polfec$eggs<-round(polfec$eggs,0)
head(polfec)
nrow(polfec)
nrow(fec)
f.mod <- glm.nb(eggs ~ area * spp_code, polfec)
summary(f.mod)
f.pred<-pred.sp(f.mod, polfec)
f.plot<-plot.sp(polfec, f.pred, "eggs",0, "Polyp fec")+scale_y_log10()
f.plot



#######################################
# IMAGES
#######################################
tab<-readPNG("data/silhouettes/tabular.png")
tab<-rasterGrob(tab, interpolate=TRUE)
mas<-readPNG("data/silhouettes/massive.png")
mas<-rasterGrob(mas, interpolate=TRUE)
cor<-readPNG("data/silhouettes/corymbose.png")
cor<-rasterGrob(cor, interpolate=TRUE)
cor<-readPNG("data/silhouettes/corymbose.png")
cor<-rasterGrob(cor, interpolate=TRUE)
dig<-readPNG("data/silhouettes/digitate.png")
dig<-rasterGrob(dig, interpolate=TRUE)
brn<-readPNG("data/silhouettes/branching.png")
brn<-rasterGrob(brn, interpolate=TRUE)
#######################################






#######################################
# MAXIMUM GROWTH by MORPHOLOGY
#######################################
dat$year<-as.factor(dat$year)
qreg<-function(q){mor.95<-rq(area_next ~ area + year*morph, data=dat, tau=q)
summary(mor.95)
Yr<-length(unique(dat$year))
Mo<-length(unique(dat$morph))
i95y <-rep(coef(mor.95)[[1]], Yr) +c(0, as.numeric(coef(mor.95)[3:6])) 
i95m <-rep(coef(mor.95)[[1]], Mo) +c(0, as.numeric(coef(mor.95)[7:10])) 
m95<-replicate(Yr, i95m)+t(replicate(Mo, i95y))
m95<-m95+rbind(0,cbind(0, t(matrix(as.numeric(coef(mor.95)[11:26]), ncol=Yr-1))))
rownames(m95)<-sort(unique(dat$morph))
colnames(m95)<-sort(unique(dat$year))
m95}
g.upper<-qreg(0.95)
g.lower<-qreg(0.05)
#------------------------------- plot
ggplot(melt(g.upper), aes(x=reorder(Var1, value), y=10^(value)))+
geom_hline(yintercept=10^(mean(g.upper)), linetype="dotted")+
geom_violin(fill="beige")+geom_point()+labs(y="intercept")+
theme(axis.title.x=element_blank()) #Dornelas 2017




# LIZARD>>>>


tax<-data.frame(rbind(cbind("Acropora.humilis.Group","digitate"),
cbind(c("Acropora.hyacinthus","Acropora.hyacinthus.Group","Acropora.divaricata.Group"),"tabular"), cbind(c("Acropora.aspera.Group","Acropora.florida.Group", "Acropora.muricata.Group","Acropora.robusta.Group","Acropora.lovelli.Group","Acropora.horrida.Group"), "staghorn"), cbind(c("Acropora.echinata.Group", "Acropora.latistella.Group", "Acropora.loripes.Group", "Acropora.nasuta.Group", "Acropora.selago.Group"), "corymbose"), cbind("Goniastrea","massive")))
tax



tax<-data.frame(rbind(
cbind("Acropora.humilis.Group","AS"),
cbind("Acropora.humilis.Group","AD"),
cbind(c("Acropora.hyacinthus","Acropora.hyacinthus.Group"),"AH"), 
cbind(c("Acropora.hyacinthus","Acropora.hyacinthus.Group"),"AC"), 
cbind("Acropora.muricata.Group", "AI"), 
cbind("Acropora.robusta.Group", "AR"), 
cbind("Acropora.nasuta.Group", "AN"), 
cbind(c("Acropora.latistella.Group", "Acropora.loripes.Group", "Acropora.selago.Group"), "AM"), 
cbind(c("Acropora.latistella.Group", "Acropora.loripes.Group", "Acropora.selago.Group"), "AL"), 
cbind("Goniastrea","GR"),
cbind("Goniastrea","GP")))
tax



lizard<-read.csv("~/Documents/PostDoc/01_Trimodal/analysis_mm/data/other/lizard_pratchett/Coral Data Lizard Island 95-19.csv")
unique(lizard$Year)
#ggplot(aggregate(lizard$Total.Coral.Cover, by=list(lizard$Year, lizard$Site, lizard$Zoneno), FUN=mean), aes(x=Group.1, y=x, col=Group.3))+geom_point()+geom_line(aes(group=interaction(Group.2, Group.3)))+facet_wrap(~Group.2)
liz<-subset(lizard, select=-c(X.,Total.Coral.Cover,NCHS,Total.soft.coral))
liz<-melt(liz, id.vars=c("Year","Site","Zoneno","Replicate", "Method"), variable.name="Taxon", value.name="x")
unique(liz$Taxon)
liz<-liz[order(liz$Year),] 
liz$Taxon<- tax$X2[match(liz$Taxon, tax$X1)]    # MORE GROUPS WITH # !!!
liz<-aggregate(x~Year+Site+Zoneno+Replicate+Method+Taxon, liz, sum) # with above
liz<-liz[liz$Taxon %in% unique(c(tax$X1, tax$X2)),]
liz<-liz[liz$Year %in% c(2002,2009,2011, 2017), ] #2002, 1999, 1995, 2017
liz<-liz[liz$Zone %in% c("2. Crest"), ]
#liz<-subset(liz, Site!="2. Washing Machine")
liz$Year<-factor(liz$Year) 

ggplot()+
#geom_path(data=liz, aes(x=Year, y=x, group=paste(Site, Zoneno, Replicate)), size=0.15, col="grey", alpha=0.5)+
geom_rect(data=liz, aes(xmin="2009",xmax="2011",ymin=0,ymax=20), fill="grey90")+
geom_path(data=aggregate(x~Taxon+Year+Zoneno+Site, data=liz, mean), aes(x=Year, y=x, group=paste(Site, Zoneno), col=Site, linetype=Zoneno), size=0.35, arrow=arrow(length=unit(0.1, "cm"), type="closed"))+
geom_path(data=aggregate(x~Taxon+Year, data=liz, mean), aes(x=Year, y=x, group=1), arrow=arrow(length=unit(0.25, "cm")), size=0.5)+
scale_y_sqrt()+
scale_colour_manual(values=c("#8dd3c7","#bebada","#fb8072", "#80b1d3"))+
facet_wrap(~Taxon)+labs(y="mean % cover")+ #scales="free_y"
theme(strip.text=element_text(size=8), axis.text=element_text(size=8), strip.background=element_blank())



# only for morph
aliz<-aggregate(x~Taxon+Year, data=liz[liz$Year %in% c(2002,2011),], mean)
dif<-data.frame(dcast(aliz, Taxon~Year, value.var="x"))
dif$x<-(dif$X2011/dif$X2002)/(11-2)

params$dcov<-dif$x[match(params$spp, dif$Taxon)]
ggplot(params, aes(y=dcov, x=lam_3))+
geom_text(aes(label=spp))+
geom_smooth(method="lm", formula=y~x)+
geom_abline(slope=1)+
scale_y_log10()













#######################################
# IPM + RECRUITMENT -> LAMBDA  *******CHANGING SIZE WINDOW
#######################################

2*sqrt(10^rec.size*10000)/pi  #rec.diam

lam_const <- NULL
par(mfcol=c(2, 6))
for (sp in spp) {
	sub<-dat[dat$spp_code==sp,]
min.size<-0.9*min(c(sub$area, sub$area_next), na.rm=T)
max.size<-1.1*max(c(sub$area, sub$area_next), na.rm=T)
#rec.size <- -2.5 #-2
n <- 100
b <- min.size + c(0:n) * (max.size - min.size)/n
y <- 0.5 * (b[1:n]+b[2:(n+1)]) # midpoints
h <- y[2] - y[1]
I <- y ==y #>= rec.size
	#rec  <-   0.000001
	rec  <- 0.000005 #*50.6 
	#rec<-0.0001
  mod <- bigmatrix()
	image(y, y, t(mod$K)) 
	points(sub$area, sub$area_next, cex=0.25)
	title(sp, line=-1)
	abline(0, 1, lty=2)
lam_const<-c(lam_const,bigmatrix()$lam)
}
#dev.off()
params$lam_const<-lam_const




#######################################
# FIT REC TO SIZE STRUCTURE 
#######################################
head(ss)
#ggplot(ss)+geom_density(aes(x=area, linetype=as.factor(year), col=as.factor(year)))+scale_y_continuous(expand=c(0,0))+theme(legend.title=element_blank())+facet_wrap(~species_code)


#------------------------------- recruitment
rec.ll <- function(x) {
  cnt <- size.dist$count[I] # non-recruits
  rec <<- x[1]
  mod <- bigmatrix()
  eig.vec <- mod$w[I]/sum(mod$w[I])
  return(-sum(cnt * log(eig.vec))) } # log-likelihood 


#-------------------------------
fitrec<-NULL
for (sp in spp) {
	smax<-params[params$spp==sp, "smax"]

	size.dist <- hist(ss$area[ss$spp_code==sp & ss$area > min.size & ss$area < max.size], breaks=b, plot=FALSE)
	rec.fit <- mle(rec.ll, start = list(x = 10), method="Brent", lower = 0, upper=10000)
	rec<-coef(rec.fit)
	#rec<-0.000005
	mod<-bigmatrix()
	fitrec<-rbind(fitrec,data.frame(x=y[I], y=(mod$w[I]/sum(mod$w[I])/h), rec, spp_code=sp, mn=min.size, mx=max.size)) }
params$rec<-unique(fitrec[,c("spp_code","rec")])$rec
#-------------------------------
head(fitrec)

ggplot()+
geom_density(data=ss, aes(x=area, linetype=as.factor(year), col=as.factor(year)))+
#geom_density(data=dat, aes(x=area), col="grey")+
geom_line(data=fitrec, aes(x=x, y=y))+
#geom_vline(xintercept=c(min.size, max.size), linetype="dotted")+
scale_y_continuous(expand=c(0,0))+theme(legend.title=element_blank())+
facet_wrap(~spp_code)

lam_fit<-NULL
for (sp in spp) {
	rec<-params[params$spp==sp,"rec"]
	lam_fit<-c(lam_fit,bigmatrix()$lam)}
  params$lam_fit<-lam_fit

  plot_grid(
  ggplot(params)+geom_point(aes(x=reorder(spp, -rec), y=rec)),
  ggplot(params,aes(y=lam_fit, x=rec))+geom_text(aes(label=spp)))
  
  



d.cm2(-2.5)


#######################################
# RECRUIT SIZE
#######################################
# caribbean agaricia spp grows at ~1cm/month. 12cm/yr = similar to adults?
# plot recruit size dists
  recsizes<-NULL
  for (sp in spp) {
  	siz<- dnorm(y,mean=log10(pi*(params$r.int[params$spp==sp]^2)), sd=0.05)
  	recsizes<-rbind(recsizes, data.frame(size=y, d=siz*h, sp=sp))
  	} 
  	aggregate(d~sp, recsizes, sum)
  ggplot(recsizes, aes(x=size, y=d, col=sp))+geom_line()+
  scale_colour_manual(values=as.character(params$cols))
# mean rec size in cm
  rec_cm2<-NULL
  for (sp in spp) {
rec_cm2<-rbind(rec_cm2, data.frame(d=d.cm2(area(params$r.int[params$spp==sp])), sp=sp))}
ggplot(rec_cm2, aes(x=reorder(sp, -d), y=d, fill=sp))+geom_bar(stat="identity")+
scale_fill_manual(values=as.character(params$cols))+labs(x="sp",y="rec.diam (cm)")



  
  dmin<-min(apply(aggregate(list(a=dat$area, b=dat$area_next), by=list(dat$spp_code),min)[,c("a","b")], 1, min))
dmax<-max(apply(aggregate(list(a=dat$area, b=dat$area_next), by=list(dat$spp_code),max)[,c("a","b")], 1, max))









#######################################
#-------------------------------  recruitment-lamda relationship
#######################################
store<-data.frame()
recx <- c(0.0000001, 0.0000005, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001)
for (sp in spp) {
	for(rec in recx){
		rec<-rec
		store<-rbind(store, data.frame(sp=sp, rec=rec, lam=bigmatrix()$lam))
  }}
  store

ggplot(data=store, aes(x=rec, y=lam, colour=sp))+
geom_line()+
geom_hline(yintercept=1)+
#lims(x=c(0, 0.00001), y=c(0,4))+
ggtitle("relative differences")+
scale_x_log10()+
scale_colour_manual(values=as.character(params$cols))






#######################################
# PLOT ELASTICITY
#######################################

sensdat<-NULL
for (sp in spp) {
	sensdat<-rbind(sensdat, data.frame(melt(analyses()$sa), spp=sp))}
  #ggplot(sensdat, aes(x=Var2, y=Var1, col=value))+ #^0.3 geom_point()+ggtitle("sensitivity analysis")+scale_colour_distiller(palette="Spectral")+facet_wrap(~spp)

# check whether axes are arranged
elasdat<-NULL
for (sp in spp) {
	elasdat<-rbind(elasdat, data.frame(melt(analyses()$eK), spp=sp))}
  #ggplot(elasdat, aes(x=Var2, y=Var1, col=value))+geom_point()+ggtitle("elasticity analysis")+scale_colour_distiller(palette="Spectral")+ facet_wrap(~spp)
  
 aggregate(value~spp, elasdat, sum)








q.pred<-NULL
r.int<-NULL
for(sp in spp){
	sub<-dat[dat$spp_code==sp,]
	q.95a <-rq(g_radius ~ 1, data=sub, tau=0.95)# no slope
	#q.95b <-rq(g_radius ~ area, data=dat, tau=0.95)# slope
	r.int<-c(r.int, coef(q.95a)[1])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), spp_code=sp)
	new$pred<-predict(q.95a, new, type="response")
	q.pred<-rbind(q.pred, new)
	}
params$r.int <- r.int

ggplot(data=q.pred)+
  geom_line(aes(x=area, y=pred, col=spp_code))+
  scale_colour_manual(values=cols)
