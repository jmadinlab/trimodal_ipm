
#######################################
# LAMBDA BOOTSTRAP
#######################################

#result_list<-lapply(1:10, function(n)source("R/boot.R"))
#lams.x<-do.call(rbind,lapply(result_list,function(x) x$value))
#nrow(lams.x)/11
#ggplot(data=lams.x, aes(x=lam, fill=spp))+geom_density()+scale_fill_manual(values=cols)



logit <- function(x) { log(x/(1-x)) }
inv.logit <- function(x) { exp(x)/(1+exp(x)) }

	 r.yx <- function(y, x) {	
 	mat<- m.int + m.slp*x
 	fec<- exp(f.int + f.slp *x)  
   out <- (rec* mat * fec) 
   out[x < rec.size | y >= rec.size] <- 0 
   return(out)
   } 
   
   s.x <- function(x) { 
	u <- s.int + s.slp * x + s.slp.2 * x^2
  return(inv.logit(u)) 
}

g.yx <- function(y, x) {
	dnorm(y, mean=g.int+ g.slp*x,
	    sd=sqrt(g.var))
	    }
	
	p.yx <- function(y, x) {
  # x <- -3
  g <- a_func(r_func(10^x) + r.int )
  #+ 1.96 * params$r.err[params$spp==sp])
  temp <- 10^y / g  # proportion of max reached. 
  temp[temp > 1] <- 1
  dnorm(logit(1 - temp), p.int + x * p.slp, p.sig)
}



lamboot <- NULL 

n.out <- 5

for(sp in spp){
	#sp <- "AC"
	# sampling
cols.x <- unique(dat$colony_id[dat$spp==sp])
samp <- sample(cols.x, length(cols.x)-n.out, replace=F)
gboot <- gdat[gdat$spp==sp & gdat$colony_id %in% samp,]
sboot <- sdat[sdat$spp==sp & sdat$colony_id %in% samp,]
fcols <- unique(fec$id[fec$spp==sp])
fsamp <- sample(fcols, (length(fcols)-(n.out*2)), replace=F)
fboot <- fec[fec$spp==sp & fec$id %in% fsamp,]

# params
m.mod<-glm(reproductive ~ area, family="binomial", data=fboot) 
	m.int<-coef(m.mod)[1]
	m.slp<-coef(m.mod)[2]
fboot$f.cm2<-fboot$polyps_cm2*fboot$eggs 
fboot$fecundity<-fboot$f.cm2*fboot$area_cm2 
f.mod<-glm.nb(fecundity ~ area, fboot[fboot$reproductive==1,])
	f.int<-coef(f.mod)[1]
	f.slp<-coef(f.mod)[2]
g.mod<-lm(area_next ~ area, data=gboot) 
	g.int<-coef(g.mod)[1]
	g.slp<-coef(g.mod)[2]
	g.var<-var(residuals(g.mod))
s.mod<-glm(surv ~ area, family="binomial", data=sboot) 
s.mod2<-glm(surv ~ area + area_sq, family="binomial", data=sboot) 
	s.int<-coef(s.mod2)[1]
	s.slp<-coef(s.mod2)[2]
	s.slp.2<-coef(s.mod2)[3]	
gboot$radius1<-sqrt(gboot$area_cm2/pi)/100
gboot$radius2<-sqrt(gboot$area_cm2_next/pi)/100
gboot$g_radius<-gboot$radius2-gboot$radius1
r.mod <-rq(g_radius ~ 1 , data=gboot, tau=0.98) # no slope
	#r.mod <-rq(g_radius ~ area , data=sub, tau=0.99)# slope
	r.int<- coef(r.mod)[1]
	#r.slp<-coef(r.mod)[2]
	r.err <- summary(r.mod, se='boot')$coef[[2]]
gboot$max_g<- r.int
gboot$max_area2<-pi*(gboot$max_g + gboot$radius1)^2 # eq. 2.1 
gboot$p_mort<-  1 - ((gboot$area_cm2_next/10000)/gboot$max_area2) # eq. 2.2
gboot$p_stasis<-  1 - ((gboot$area_cm2/10000)/gboot$max_area2) 
pboot<-subset(gboot, p_mort>0.001) # less than max growth line	
	pboot$pm_logit <- logit(pboot$p_mort) 
	p.mod<-lm(pm_logit ~ area , data=pboot)
	p.int <- coef(p.mod)[[1]]
	p.slp <- coef(p.mod)[[2]]
	p.sig <- sigma(p.mod)
# IPM
	rec <- 1*10^-3
	rec.size <- params$rec.size[params$spp==sp]
	#b <- seq(rec.size, max.size, length=n)
	b <- seq(min.size, max.size, length=n)
	h <- b[2] - b[1]
	b <- c(min(b)-(2*h), min(b)-h, b)
	y <- 0.5 * (b[1:n]+b[2:(n+1)])
	I <- y >= rec.size
	
    mod <- bigmatrix()
    mod$lam
    lamboot <- rbind(lamboot, data.frame(lam=mod$lam, spp=sp))
    }
     lamboot

