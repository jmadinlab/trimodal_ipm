

morph <- unique(fec$morphology)
p.morph <- data.frame(morph)

#######################################
# MATURITY 
#######################################
fec$reproductive <- ifelse(fec$eggs>0, 1, 0) # reproductive status
m.int<-NULL
m.slp<-NULL
m.pred<-NULL
for(m in morph){
	sub<-fec[fec$morphology==m,]
	m.mod<-glm(reproductive ~ area, family="binomial", data=sub) 
	m.int<-c(m.int, coef(m.mod)[1])
	m.slp<-c(m.slp, coef(m.mod)[2])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), morph=m)
	new$pred<-predict(m.mod, new, type="response")
	m.pred<-rbind(m.pred, new)
}
p.morph$m.int<-m.int
p.morph$m.slp<-m.slp

# -------- plot
#pdf("figs/p1_maturity.pdf" )

p1<-ggplot()+
  geom_jitter(data=fec, aes(x= area, reproductive), height=0.02, shape=21, col="grey")+
  geom_line(data=m.pred, aes(area, pred, col=morph))+
  theme_classic()
p1

#show(p1)
#dev.off()

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
for(m in morph){
	#sp <- "Ahy"
	sub<-fec[fec$morph==m & fec$reproductive==1,]
	f.mod<-glm.nb(fecundity ~ area, sub)
	f.int<-c(f.int, coef(f.mod)[1])
	f.slp<-c(f.slp, coef(f.mod)[2])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), morph=m)
	new$pred<-predict(f.mod, new, type="response")
	f.pred<-rbind(f.pred, new)
	}
p.morph$f.int<-f.int
p.morph$f.slp<-f.slp

p2<-ggplot()+ 
  geom_point(data=fec[fec$reproductive==1, ], aes(area, fecundity), shape=21, col="grey")+
  geom_line(data=f.pred, aes(area, pred, col=morph))+
  scale_y_log10()+
  #facet_wrap(~morph)+
  #scale_colour_manual(values=cols)+
  theme_classic()
p2
#show(p2)
#dev.off()


#######################################
# GROWTH 
#######################################

g.int<-NULL
g.slp<-NULL
g.var<-NULL
g.pred<-NULL
for(m in morph){
	sub<-gdat[gdat$morph==m,]
	g.mod<-lm(area_next ~ area, data=sub) 
	g.int<-c(g.int, coef(g.mod)[1])
	g.slp<-c(g.slp, coef(g.mod)[2])
	g.var<-c(g.var, var(residuals(g.mod)))
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), morph=m)
	new$pred<-predict(g.mod, new, type="response")
	g.pred<-rbind(g.pred, new)
}
p.morph$g.int<-g.int
p.morph$g.slp<-g.slp
p.morph$g.var<-g.var

#pdf("figs/p3_growth.pdf")

p3<-ggplot()+ 
  geom_abline(slope=1, linetype="dotted")+
  geom_point(data=gdat, 
    aes(area, area_next), 
      shape=21, col="grey")+
  geom_line(data=g.pred, 
    aes(area, pred, col=morph))+
  #scale_colour_manual(values=cols)+
  #facet_wrap(~morph)+
 theme_classic()

#show(p3)
#dev.off() 

p.morph[p.morph$morph=="massive", c("p.int", "p.slp","p.sig")] <- params[params$spp=="Gpe", c("p.int", "p.slp","p.sig")]
p.morph[p.morph$morph=="massive", c("g.int", "g.slp","g.var")] <- params[params$spp=="Gpe", c("g.int", "g.slp","g.var")]

  

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
for(m in morph){
	sub<-sdat[sdat$morph==m,]
	s.mod<-glm(surv ~ area, family="binomial", data=sub) 
	s.mod2<-glm(surv ~ area + area_sq, family="binomial", data=sub) 
	s.mod.c <- s.mod2
	s.int<-c(s.int, coef(s.mod.c)[1])
	s.slp<-c(s.slp, coef(s.mod.c)[2])
	s.slp2<-c(s.slp2, coef(s.mod.c)[3])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), morph=m)
	new$area_sq <- new$area^2
	new$pred<-predict(s.mod.c, new, type="response")
	#new2<-data.frame(area=p.morph[p.morph$morph==sp,"rec.size"], morph=sp, morph=sub$morphology[1])
	#new2$area_sq <- new2$area^2
	#new2$pred<-predict(s.mod2, new2, type="response")
	#s.rec <- rbind(s.rec, new2)
	s.pred<-rbind(s.pred, new)
	}
p.morph$s.int<-s.int
p.morph$s.slp<-s.slp
p.morph$s.slp.2<-s.slp2

#pdf("figs/p4_survival.pdf")

p4<-ggplot()+ 
  geom_jitter(data=sdat, 
    aes(area, surv), 
      shape=21, col="grey",height=0.02)+
  geom_line(data=s.pred, 
    aes(area, pred, col=morph))+
      #facet_wrap(~morph)+
  #scale_colour_manual(values=cols)+
 theme_classic()
p4
#show(p4)
#dev.off() 



#######################################
# MAXIMUM GROWTH  - look into AD
#######################################
library("quantreg")
gdat$radius1<-sqrt(gdat$area_cm2/pi)/100
gdat$radius2<-sqrt(gdat$area_cm2_next/pi)/100
gdat$g_radius<-gdat$radius2-gdat$radius1

r.int<-NULL
r.slp<-NULL
r.err<-NULL
r.pred<-NULL
for(m in morph){
	sub<-gdat[gdat$morph==m,]
	r.mod <-rq(g_radius ~ 1 , data=sub, tau=0.98) # no slope
	#r.mod <-rq(g_radius ~ area , data=sub, tau=0.99)# slope
	r.int<-c(r.int, coef(r.mod)[1])
	#r.slp<-c(r.slp, coef(r.mod)[2])
	r.err <- c(r.err, summary(r.mod, se='boot')$coef[[2]])
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  morph=m)
	new$pred<-predict(r.mod, new, type="response")
	r.pred<-rbind(r.pred, new)
	}
p.morph$r.int<-r.int
p.morph$r.slp<-r.slp
p.morph$r.err<-r.err

#pdf("figs/p5_maxgrowth.pdf", )

p5<-ggplot()+
 geom_point(data=gdat, aes(x=area, y=g_radius), col="grey", shape=21)+
     geom_line(data=r.pred, aes(x=area, y=pred, col=morph))+
   #+facet_wrap(~morph)+
    theme_classic()
p5
#show(p5)
#dev.off() 
   


#######################################
# PARTIAL MORTALITY  - subset??
#######################################
logit <- function(x) { log(x/(1-x)) }
inv.logit <- function(x) { exp(x)/(1+exp(x)) }

gdat$max_g<-p.morph$r.int[match(gdat$morph, p.morph$morph)] 
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
for(m in morph){
	sub<-pdat[pdat$morph==m,]
	p.mod<-lm(pm_logit ~ area , data=sub)
	p.int <- c(p.int, coef(p.mod)[[1]])
	p.slp <- c(p.slp, coef(p.mod)[[2]])
	p.sig <- c(p.sig, sigma(p.mod))
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1), 
	  morph=m)
	new$pred<-predict(p.mod, new, type="response")
	p.pred<-rbind(p.pred, new)
	}
p.morph$p.int<-p.int
p.morph$p.slp<-p.slp
p.morph$p.sig<-p.sig

#pdf("figs/p6_partialmort.pdf")

p6<-ggplot()+
  geom_point(data=pdat, 
    aes(x=area, y=pm_logit), 
       shape=21, col="grey")+
    geom_line(data=pdat, aes(x=area, y=logit(p_stasis), group=morphology), linetype="dotted")+
   geom_line(data=p.pred, 
     aes(x=area, y=pred, col=morph))+
     #scale_colour_manual(values=c(cols))+
     #facet_wrap(~morph)+
     theme_classic()
p6

ggplot()+
  geom_point(data=pdat, 
    aes(x=area, y=inv.logit(pm_logit)), 
       shape=21, col="grey")+
    geom_line(data=pdat, aes(x=area, y=p_stasis, group=morphology, col=morphology), linetype="dotted")+
   geom_line(data=p.pred, 
     aes(x=area, y=inv.logit(pred), col=morph))+
     #scale_colour_manual(values=c(cols))+
     #facet_wrap(~morph)+
     theme_classic()

#show(p6)
#dev.off() 

plot_grid(p1, p2, p3, p4, p5, p6)


