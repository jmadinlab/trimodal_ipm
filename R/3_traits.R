
#######################################
# DEMOGRAPHIC TRAITS
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

#######################################
# AVERAGE SIZE
#######################################

# average size in demo models
size.dat <- aggregate(area_cm2~spp, dat[!is.na(dat$spp),], mean)
params$size.dat <- size.dat$area_cm2[match(params$spp, size.dat$spp)]
params$size.dat <- log10(params$size.dat/10000)

# average size in size structure data
size.ss <- aggregate(area~spp, ss[!is.na(ss$spp),], mean)
params$size.ss <- size.ss$area[match(params$spp, size.ss$spp)]





# write.csv(params, "data/params.csv")

