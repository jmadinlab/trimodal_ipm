

#######################################
# DATA
#######################################

params<-read.csv("data/info.csv")
dat<- read.csv("data/trimodal.csv")
gdat<-read.csv("data/growth.csv")
sdat<-read.csv("data/survival.csv")
fec<-read.csv("data/fecundity.csv")
ss<-read.csv("data/size_structure.csv")
spp<-params$spp[order(params$spp)]
params$abundance<-ifelse(params$abundance_05<500, "Rare", ifelse(params$abundance_05>1500, "Dominant","Common"))
params$abundance<-factor(params$abundance, levels=c("Rare", "Common", "Dominant"))
cols<-as.character(params$cols)


#################################### 
# ---- remove F9 (hurricane)
unique(gdat$year)
nrow(gdat)
gdat <- gdat[!gdat$year==2014,] #
nrow(gdat)

unique(sdat$year)
nrow(sdat)
sdat <- sdat[!sdat$year==2015,]
nrow(sdat)
sdat$surv <- ifelse(sdat$year==2014, NA, sdat$surv)


# ----------------------------# mean polyp density
pd <- read.csv("data/polyp_density.csv")
pd.mean <- aggregate(list(polyps_cm2=pd$polyps_cm2),list(spp=pd$spp),mean)
pd.mean <- rbind(pd.mean, data.frame(spp="AM", polyps_cm2=pd.mean$polyps_cm2[pd.mean$spp=="AL"]))
params<-merge(params, pd.mean)
# ---------------------------# mean energetics
ec <- read.csv("data/egg_energy.csv")
ec.mean <- aggregate(list(eggC=ec$Carbon_ug_corrected, eggN=ec$Nitrogen_ug) ,list(spp=ec$spp), mean)
params<-merge(params, ec.mean)
# fec
fec$eggC <- params$eggC[match(fec$spp, params$spp)]
fec$polyps_cm2 <- params$polyps_cm2[match(fec$spp, params$spp)]
fec$morphology <- params$morphology[match(fec$spp, params$spp)]


# --- Aggregate fecundity
##################################
fec<-aggregate(.~id+spp+morphology, fec[,c("id","area_cm2","eggs","spp", "polyps_cm2","eggC","morphology")], mean)

# --- All areas in log10 m2
##################################
fec$area <- log10(fec$area_cm2 /10000) 
gdat$area <- log10(gdat$area_cm2 / 10000) 
gdat$area_next <- log10(gdat$area_cm2_next / 10000) 
sdat$area <- log10(sdat$area_cm2 / 10000) 
ss$area <- log10(ss$area_cm2 / 10000) 
gdat$morphology <- params$morphology[match(gdat$spp, params$spp)]
sdat$morphology <- params$morphology[match(sdat$spp, params$spp)]






