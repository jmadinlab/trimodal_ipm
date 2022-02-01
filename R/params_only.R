

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
# eggs per cm2
fec$f.cm2<-fec$polyps_cm2*fec$eggs 
# eggs per colony
 fec$fecundity<-round(fec$f.cm2*fec$area_cm2)
# egg energy per colony
# make energy a survival proportion, and round to 0
#fec$fecundity<-round(fec$f.cm2*fec$area_cm2 * (fec$eggC/max(fec$eggC)) ,0) 

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


