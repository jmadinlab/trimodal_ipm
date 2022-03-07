
#######################################
# MATURITY 
#######################################
fec$reproductive <- ifelse(fec$eggs>0, 1, 0) # reproductive status
m.int<-NULL
m.slp<-NULL
m.pred<-NULL
for(sp in spp){
	sub<-fec[fec$spp==sp,]
	m.mod<-glm(reproductive ~ area, family="binomial", data=sub) 
	m.int<-c(m.int, coef(m.mod)[1])
	m.slp<-c(m.slp, coef(m.mod)[2])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  spp=sp, morph=sub$morphology[1])
	new$pred<-predict(m.mod, new, type="response")
	m.pred<-rbind(m.pred, new)
}
params$m.int<-m.int
params$m.slp<-m.slp

m.pred$spp <- factor(m.pred$spp, levels=order)

#######################################
# FECUNDITY  
#######################################

fec$f.cm2<-fec$polyps_cm2*fec$eggs  # eggs per cm2
fec$fecundity<-round(fec$f.cm2*fec$area_cm2) # eggs per colony

f.int<-NULL
f.slp<-NULL
f.pred<-NULL
for(sp in spp){
	#sp <- "Ahy"
	sub<-fec[fec$spp==sp & fec$reproductive==1,]
	f.mod<-glm.nb(fecundity ~ area, sub)
	f.int<-c(f.int, coef(f.mod)[1])
	f.slp<-c(f.slp, coef(f.mod)[2])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1),
	spp=sp, morph=sub$morphology[1])
	new$pred<-predict(f.mod, new, type="response")
	f.pred<-rbind(f.pred, new)
	}
params$f.int<-f.int
params$f.slp<-f.slp

#######################################
# FECUNDITY - CONSTANT SLOPE
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
params$f.int.const <- int(fec.con) # intercept = effect of x2
params$f.int.const

#######################################
# GROWTH 
#######################################

g.int<-NULL
g.slp<-NULL
g.var<-NULL
g.pred<-NULL
for(sp in spp){
	sub<-gdat[gdat$spp==sp,]
	g.mod<-lm(area_next ~ area, data=sub) 
	g.int<-c(g.int, coef(g.mod)[1])
	g.slp<-c(g.slp, coef(g.mod)[2])
	g.var<-c(g.var, var(residuals(g.mod)))
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1),
	  spp=sp, morph=sub$morphology[1])
	new$pred<-predict(g.mod, new, type="response")
	g.pred<-rbind(g.pred, new)
}
params$g.int<-g.int
params$g.slp<-g.slp
params$g.var<-g.var

#######################################
# SURVIVAL 
#######################################
sdat$area_sq<-sdat$area^2
sdat<-sdat[!is.na(sdat$surv),]
#sdat<-sdat[sdat$area>-3,]
s.int<-NULL
s.slp<-NULL
s.slp2<-NULL
s.pred<-NULL
#s.rec <- NULL
for(sp in spp){
	sub<-sdat[sdat$spp==sp,]
	s.mod<-glm(surv ~ area, family="binomial", data=sub) 
	s.mod2<-glm(surv ~ area + area_sq, family="binomial", data=sub) 
	s.mod.c <- s.mod2
	s.int<-c(s.int, coef(s.mod.c)[1])
	s.slp<-c(s.slp, coef(s.mod.c)[2])
	s.slp2<-c(s.slp2, coef(s.mod.c)[3])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  spp=sp, morph=sub$morphology[1])
	new$area_sq <- new$area^2
	new$pred<-predict(s.mod.c, new, type="response")
	#new2<-data.frame(area=params[params$spp==sp,"rec.size"], spp=sp, morph=sub$morphology[1])
	#new2$area_sq <- new2$area^2
	#new2$pred<-predict(s.mod2, new2, type="response")
	#s.rec <- rbind(s.rec, new2)
	s.pred<-rbind(s.pred, new)
	}
params$s.int<-s.int
params$s.slp<-s.slp
params$s.slp.2<-s.slp2

#######################################
# MAXIMUM GROWTH  
#######################################
library("quantreg")
gdat$radius1<-sqrt(gdat$area_cm2/pi)/100
gdat$radius2<-sqrt(gdat$area_cm2_next/pi)/100
gdat$g_radius<-gdat$radius2-gdat$radius1

r.int<-NULL
r.slp<-NULL
r.err<-NULL
r.pred<-NULL
for(sp in spp){
	sub<-gdat[gdat$spp==sp,]
	r.mod <-rq(g_radius ~ 1 , data=sub, tau=0.98) # no slope
	#r.mod <-rq(g_radius ~ area , data=sub, tau=0.99)# slope
	r.int<-c(r.int, coef(r.mod)[1])
	#r.slp<-c(r.slp, coef(r.mod)[2])
	r.err <- c(r.err, summary(r.mod, se='boot')$coef[[2]])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  spp=sp, morph=sub$morph[1])
	new$pred<-predict(r.mod, new, type="response")
	r.pred<-rbind(r.pred, new)
	}
params$r.int<-r.int
params$r.slp<-r.slp
params$r.err<-r.err


#######################################
# PARTIAL MORTALITY  
#######################################
logit <- function(x) { log(x/(1-x)) }
inv.logit <- function(x) { exp(x)/(1+exp(x)) }

gdat$max_g<-params$r.int[match(gdat$spp, params$spp)] 
gdat$max_area2<-pi*(gdat$max_g + gdat$radius1)^2 # eq. 2.1 
gdat$p_mort<-  1 - ((gdat$area_cm2_next/10000)/gdat$max_area2) # eq. 2.2
gdat$p_stasis<-  1 - ((gdat$area_cm2/10000)/gdat$max_area2) 
pdat<-subset(gdat, p_mort>0.001) # less than max growth line
#pdat<-subset(dat, p_mort>0) # zero
pdat$pm_logit <- logit(pdat$p_mort) 

p.int<-NULL
p.slp<-NULL
p.pred<-NULL
p.sig<-NULL
for(sp in spp){
	sub<-pdat[pdat$spp==sp,]
	p.mod<-lm(pm_logit ~ area , data=sub)
	p.int <- c(p.int, coef(p.mod)[[1]])
	p.slp <- c(p.slp, coef(p.mod)[[2]])
	p.sig <- c(p.sig, sigma(p.mod))
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  spp=sp, morph=sub$morph[1])
	new$pred<-predict(p.mod, new, type="response")
	p.pred<-rbind(p.pred, new)
	}
params$p.int<-p.int
params$p.slp<-p.slp
params$p.sig<-p.sig
