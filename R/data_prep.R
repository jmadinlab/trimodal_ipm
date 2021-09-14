
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
names(cols) <- params$spp

#spp names
#################################### 
head(dat)
dat$spp <- params$spp[match(dat$species, params$species)]
gdat$spp <- params$spp[match(gdat$species, params$species)]
sdat$spp <- params$spp[match(sdat$species, params$species)]
fec$spp <- params$spp[match(fec$species, params$species)]
ss$spp <- params$spp[match(ss$species, params$species)]


# ordering
#################################### 
order <- c("Ahy","Acy","Ain","Aro","Ana","Asp","Ami","Adi","Ahu", "Gre","Gpe")
params$spp <- factor(params$spp, levels=order)
labs <- params$species[order(params$spp)]
labs

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
pd$spp <- params$spp[match(pd$species, params$species)] 
pd.mean <- aggregate(list(polyps_cm2=pd$polyps_cm2),list(spp=pd$spp),mean)
pd.mean <- rbind(pd.mean, data.frame(spp="Ami", polyps_cm2=pd.mean$polyps_cm2[pd.mean$spp=="Asp"]))
params<-merge(params, pd.mean)
# ---------------------------# mean energetics
ec <- read.csv("data/egg_energy.csv")
ec$spp <- params$spp[match(ec$species, params$species)] 
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


# --- Lit review for recruit success/size
##################################

tiles <- read.csv("data/recruitment/recruitment.csv")
tiles <- aggregate(N_m2_year~Study+Location+Family,tiles, mean) # sum?? No mean

rsize <- read.csv("data/recruitment/recsize.csv")


#######################################
# ABUNDANCE TRIMODAL
#######################################

ab <- read.csv("data/abundance/lit_merged_20180412.csv")

ab$species[ab$species %in% c("Acropora fat_digitifera", "Acropora difitifera_fat", "Acropora sp. Fat dig")] <- "Acropora cf. digitifera"

ab$species[ab$species=="Acropora robustea"]<-"Acropora robusta"

ab$species[ab$species=="Monitpora crassituberculata"]<-"Montipora crassituberculata"

ab$genus <- substr(ab$species, 1, 5)

# post 1990s
ab <- subset(ab, campaign >1997)

# all transects
unique(ab[,c("site","transect", "campaign")])

# all species
unique(ab$species)

# For TB, Goniastrea sp. == retiformis
ab$species <- ifelse(ab$species=="Goniastrea sp." & ab$observer=="TB", "Goniastrea retiformis", ab$species)

# IDs
ab$ID <- paste(ab$site, ab$campaign, ab$transect, sep="_")

# trimodal only
tri <- subset(ab, site=="Trimodal")
tri <- tri[tri$species %in% params$species,]
tri <- tri[tri$campaign %in% c(2011, 2014),]
head(tri)

# cover and N
tri2 <- aggregate(intercept_cm~., subset(tri, select=-c(site, observer,genus)), length)
tri2$cover <- aggregate(intercept_cm~., subset(tri, select=-c(site, observer,genus)), sum)$intercept_cm/10
head(tri2)

# add zeros
abun<-data.frame(ID=rep(unique(tri2$ID), length(unique(params$species))), species=rep(unique(params$species), each=length(unique(tri2$ID))))
abun$N <- tri2$intercept_cm[match(paste(abun$ID, abun$species), paste(tri2$ID, tri2$species))]
abun$cover <- tri2$cover[match(paste(abun$ID, abun$species),paste(tri2$ID, tri2$species))]
abun$cover[is.na(abun$cover)]<-0
abun$N[is.na(abun$N)]<-0
abun[,c("year")]<-tri2[match(abun$ID, tri2$ID),"campaign"]
abun$tran<-substr(abun$ID, 15,15)
head(abun)



