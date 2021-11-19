

p.all <- data.frame(spp="all")

#######################################
# MATURITY 
#######################################
fec$reproductive <- ifelse(fec$eggs>0, 1, 0) # reproductive status

	m.mod<-glm(reproductive ~ area, family="binomial", data=fec) 
	m.int<-coef(m.mod)[1]
	m.slp<-coef(m.mod)[2]
	new<-data.frame(area=seq(min(fec$area), max(fec$area), 0.1))
	new$pred<-predict(m.mod, new, type="response")
	m.pred<-new
	
p.all$m.int<-m.int
p.all$m.slp<-m.slp

# -------- plot
#pdf("figs/p1_maturity.pdf" )

p1<-ggplot()+
  geom_jitter(data=fec, aes(x= area, reproductive), height=0.02, shape=21, col="grey")+
  geom_line(data=m.pred, aes(area, pred))+
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

sub <-  fec[ fec$reproductive==1,]
	f.mod<-glm.nb(fecundity ~ area, sub)
	f.int<- coef(f.mod)[1]
	f.slp<-coef(f.mod)[2]
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1))
	new$pred<-predict(f.mod, new, type="response")
	f.pred<-new
	
	
p.all$f.int<-f.int
p.all$f.slp<-f.slp

p2<-ggplot()+ 
  geom_point(data=fec[fec$reproductive==1, ], aes(area, fecundity), shape=21, col="grey")+
  geom_line(data=f.pred, aes(area, pred))+
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


	g.mod<-lm(area_next ~ area, data=gdat) 
	g.int<-coef(g.mod)[1]
	g.slp<-coef(g.mod)[2]
	g.var<-var(residuals(g.mod))
	new<-data.frame(area=seq(min(gdat$area), max(gdat$area), 0.1))
	new$pred<-predict(g.mod, new, type="response")
	g.pred<-new
	
	
p.all$g.int<-g.int
p.all$g.slp<-g.slp
p.all$g.var<-g.var

#pdf("figs/p3_growth.pdf")

p3<-ggplot()+ 
  geom_abline(slope=1, linetype="dotted")+
  geom_point(data=gdat, 
    aes(area, area_next), 
      shape=21, col="grey")+
  geom_line(data=g.pred, 
    aes(area, pred))+
  #scale_colour_manual(values=cols)+
  #facet_wrap(~morph)+
 theme_classic()

#show(p3)
#dev.off() 
  

#######################################
# SURVIVAL 
#######################################
sdat$area_sq<-sdat$area^2
sdat<-sdat[!is.na(sdat$surv),]


	s.mod<-glm(surv ~ area, family="binomial", data=sdat) 
	s.mod2<-glm(surv ~ area + area_sq, family="binomial", data=sdat) 
	s.mod.c <- s.mod2
	s.int<-coef(s.mod.c)[1]
	s.slp<-coef(s.mod.c)[2]
	s.slp2<-coef(s.mod.c)[3]
	new<-data.frame(area=seq(min(sdat$area), max(sdat$area), 0.1))
	new$area_sq <- new$area^2
	new$pred<-predict(s.mod.c, new, type="response")
	s.pred<-new

p.all$s.int<-s.int
p.all$s.slp<-s.slp
p.all$s.slp.2<-s.slp2

#pdf("figs/p4_survival.pdf")

p4<-ggplot()+ 
  geom_jitter(data=sdat, 
    aes(area, surv), 
      shape=21, col="grey",height=0.02)+
  geom_line(data=s.pred, 
    aes(area, pred))+
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

	r.mod <-rq(g_radius ~ 1 , data=gdat, tau=0.98) # no slope
	#r.mod <-rq(g_radius ~ area , data=sub, tau=0.99)# slope
	r.int<-coef(r.mod)[1]
	#r.slp<-c(r.slp, coef(r.mod)[2])
	r.err <- summary(r.mod, se='boot')$coef[[2]]
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1))
	new$pred<-predict(r.mod, new, type="response")
	r.pred<-new
	
p.all$r.int<-r.int
p.all$r.slp<-r.slp
p.all$r.err<-r.err

#pdf("figs/p5_maxgrowth.pdf", )

p5<-ggplot()+
 geom_point(data=gdat, aes(x=area, y=g_radius), col="grey", shape=21)+
     geom_line(data=r.pred, aes(x=area, y=pred))+
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

gdat$max_g<-p.all$r.int
gdat$max_area2<-pi*(gdat$max_g + gdat$radius1)^2 # eq. 2.1 
gdat$p_mort<-  1 - ((gdat$area_cm2_next/10000)/gdat$max_area2) # eq. 2.2
gdat$p_stasis<-  1 - ((gdat$area_cm2/10000)/gdat$max_area2) 
pdat<-subset(gdat, p_mort>0.001) # less than max growth line
#pdat<-subset(dat, p_mort>0) # zero
pdat$pm_logit <- logit(pdat$p_mort) 


	p.mod<-lm(pm_logit ~ area , data=pdat)
	p.int <- coef(p.mod)[[1]]
	p.slp <- coef(p.mod)[[2]]
	p.sig <- sigma(p.mod)
	new<-data.frame(area=seq(min(sub$area), max(sub$area), 0.1))
	new$pred<-predict(p.mod, new, type="response")
	p.pred<-new

p.all$p.int<-p.int
p.all$p.slp<-p.slp
p.all$p.sig<-p.sig

#pdf("figs/p6_partialmort.pdf")

p6<-ggplot()+
  geom_point(data=pdat, 
    aes(x=area, y=pm_logit), 
       shape=21, col="grey")+
    geom_line(data=pdat, aes(x=area, y=logit(p_stasis)), linetype="dotted")+
   geom_line(data=p.pred, 
     aes(x=area, y=pred))+
     #scale_colour_manual(values=c(cols))+
     #facet_wrap(~morph)+
     theme_classic()
p6

ggplot()+
  geom_point(data=pdat, 
    aes(x=area, y=inv.logit(pm_logit)), 
       shape=21, col="grey")+
    geom_line(data=pdat, aes(x=area, y=p_stasis), linetype="dotted")+
   geom_line(data=p.pred, 
     aes(x=area, y=inv.logit(pred)))+
     #scale_colour_manual(values=c(cols))+
     #facet_wrap(~morph)+
     theme_classic()

#show(p6)
#dev.off() 

plot_grid(p1, p2, p3, p4, p5, p6)


