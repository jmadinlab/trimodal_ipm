
#######################################
# DEMOGRAPHIC DATA
#######################################

# data per colony per year
params<-read.csv("data/info.csv") # basic info
gdat<-read.csv("data/growth.csv") #area(t)/(t+1)
sdat<-read.csv("data/survival.csv") # (1/0) 
fec<-read.csv("data/fecundity.csv") # N per branch
spp<-params$spp[order(params$spp)]

#spp names
#################################### 
gdat$spp <- params$spp[match(gdat$species, params$species)]
sdat$spp <- params$spp[match(sdat$species, params$species)]
fec$spp <- params$spp[match(fec$species, params$species)]

#################################### 
# ---- remove F9 (hurricane/post hurricane)
unique(gdat$year)
gdat <- gdat[!gdat$year==2014,] 
sdat <- sdat[!sdat$year==2015,]
sdat$surv <- ifelse(sdat$year==2014, NA, sdat$surv)
unique(sdat$year)

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
params$eggC.se<-aggregate(Carbon_ug_corrected~spp, ec, sd)$Carbon_ug_corrected

# -- Bad colonies
##################################
sdat<-sdat[!sdat$colony_id==287,]
gdat<-gdat[!gdat$colony_id==287,]

# --- Aggregate fecundity to colony-level
##################################
fec<-aggregate(.~id+spp+morphology, fec[,c("id","area_cm2","eggs","spp", "polyps_cm2","eggC","morphology")], mean)

# --- All areas in log10 m2
##################################
fec$area <- log10(fec$area_cm2 /10000) 
gdat$area <- log10(gdat$area_cm2 / 10000) 
gdat$area_next <- log10(gdat$area_cm2_next / 10000) 
sdat$area <- log10(sdat$area_cm2 / 10000) 
gdat$morphology <- params$morphology[match(gdat$spp, params$spp)]
sdat$morphology <- params$morphology[match(sdat$spp, params$spp)]

# --- remove tiny colonies (<3cm diameter)
##################################
diam <- 3
min.area <- log10(pi*(diam/2/100)^2)
gdat <- subset(gdat, area > min.area)
sdat <- subset(sdat, area > min.area)
fec <- subset(fec, area > min.area)


#######################################
# SIZE STRUCTURE
#######################################

ss<-read.csv("data/size_structure.csv") 
ss$spp <- params$spp[match(ss$species, params$species)]
ss$area <- log10(ss$area_cm2 / 10000) 


#######################################
# ABUNDANCE
#######################################

# BELT TRANSCTS (2005)

a1 <- read.csv("data/BT_counts_2005.csv")
a1$Species[a1$Species == "Acropora fat dig"] <- "Acropora cf. digitifera"
a1A <- a1[a1$Species %in% params$species,]

# add zeros 
a1B<-data.frame(ID=rep(unique(a1A$ID), length(unique(a1A$Species))), Species=rep(unique(a1A$Species), each=length(unique(a1A$ID))))
a1B$N <- a1A$Abundance[match(paste(a1B$ID, a1B$Species), paste(a1A$ID, a1A$Species))]
a1B$N[is.na(a1B$N)]<-0

abun.BT <- a1B
abun.BT[,c("spp","morphology","abundance_05")] <- params[match(abun.BT$Species, params$species), c("spp","morphology","abundance_05")]


# LIT TRANSECTS (2011-14)
#######################################

a2 <- read.csv("data/LIT_counts_1995_2017.csv")
a2$species[a2$species %in% c("Acropora fat_digitifera", "Acropora difitifera_fat", "Acropora sp. Fat dig")] <- "Acropora cf. digitifera"
a2$species[a2$species=="Acropora robustea"]<-"Acropora robusta"
a2$species <- ifelse(a2$species=="Goniastrea sp." & a2$observer=="TB", "Goniastrea retiformis", a2$species)
a2$genus <- substr(a2$species, 1, 5)

# post 1990s
a2 <- subset(a2, campaign >1997)

# IDs
a2$ID <- paste(a2$site, a2$campaign, a2$transect, sep="_")

# trimodal study site only
a2A <- subset(a2, site=="Trimodal")
a2A <- a2A[a2A$species %in% params$species,]
a2A <- a2A[a2A$campaign %in% c(2011, 2014),]
head(a2A)

# cover and N
a2B <- aggregate(intercept_cm~., subset(a2A, select=-c(site, observer,genus)), length)
a2B$cover <- aggregate(intercept_cm~., subset(a2A, select=-c(site, observer,genus)), sum)$intercept_cm/10
head(a2B)

# add zeros
a2C<-data.frame(ID=rep(unique(a2B$ID), length(unique(params$species))), species=rep(unique(params$species), each=length(unique(a2B$ID))))
a2C$N <- a2B$intercept_cm[match(paste(a2C$ID, a2C$species), paste(a2B$ID, a2B$species))]
a2C$cover <- a2B$cover[match(paste(a2C$ID, a2C$species),paste(a2B$ID, a2B$species))]
a2C$cover[is.na(a2C$cover)]<-0
a2C$N[is.na(a2C$N)]<-0
a2C[,c("year")]<-a2B[match(a2C$ID, a2B$ID),"campaign"]
a2C$tran<-substr(a2C$ID, 15,15)
head(a2C)

abun.LIT <- a2C
abun.LIT[,c("spp","morphology","abundance_05")] <- params[match(abun.LIT$species, params$species), c("spp","morphology","abundance_05")]

#write.csv(abun.BT, "data/abun.BT.csv")
#write.csv(abun.LIT, "data/abun.LIT.csv")

