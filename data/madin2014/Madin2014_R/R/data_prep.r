rm(list = ls())
setwd("~/Desktop/Archive")

# Data prep for mortality analysis

# Find most recent data file from trimodal (i.e., tagged_YYYYMMDD.csv)
files <- dir("data/", "tagged*")
current <- sort(files, decreasing = T)[1]
tag <- read.csv(paste("data/", current, sep = ""), as.is = T)

# Only keep F3-F6 data (F6 is F11 in database...)
tag$fieldtrip_id[tag$fieldtrip_id == "F11"] <- "F6"
tag <- tag[tag$fieldtrip_id %in% c("F3", "F4", "F5", "F6", "F7", "F8", "F9"),]
unique(tag$fieldtrip_id)
# Only keep area data where photos and outlines are acceptable
tag <- tag[as.character(tag$outline_acceptable) != "f" & as.character(tag$photo_acceptable) != "f",]

# Only keep colonies that were:
tag <- tag[tag$action == "photographed" | tag$action == "dead" | tag$action == "gone",]
# tag <- tag[tag$action == "photographed" | tag$action == "dead" | tag$action == "gone" | tag$action == "not found",]
# tag <- tag[tag$action == "not found",]

# Remove GP colony under hyacinthus
tag <- tag[tag$colony_id != 329,]


mor <- data.frame(
  species_code=tapply(tag$species_code, list(tag$colony_id), max),
  colony_id=unique(tag$colony_id),
  tapply(log(tag$area_cm2/10000), list(tag$colony_id, tag$fieldtrip_id), mean),
  tapply(tag$action, list(tag$colony_id, tag$fieldtrip_id), max)
)

mor$form[mor$species_code == "AC" | mor$species_code == "AH"] <- "table"
mor$form[mor$species_code == "AD" | mor$species_code == "AS"] <- "digitate"
mor$form[mor$species_code == "AI" | mor$species_code == "AR"] <- "branching"
mor$form[mor$species_code == "GR" | mor$species_code == "GP"] <- "massive"
mor$form[mor$species_code == "AL" | mor$species_code == "AM" | mor$species_code == "AN"] <- "corymbose"

head(mor)

mor$F6[mor$F6.1=="photographed"] <- 1
#mor$F7[mor$F7.1=="photographed"] <- 1 # not present -mike
#mor$F8[mor$F8.1=="photographed"] <- 1
#mor$F9[mor$F9.1=="photographed"] <- 1


# again, not elegant, but concatenate mortality data and add year column
sr34 <- mor[!is.na(mor$F3) & mor$F3.1 == "photographed",]
sr45 <- mor[!is.na(mor$F4) & mor$F4.1 == "photographed",]
sr56 <- mor[!is.na(mor$F5) & mor$F5.1 == "photographed",]
sr67 <- mor[!is.na(mor$F6) & mor$F6.1 == "photographed",]
#sr78 <- mor[!is.na(mor$F7) & mor$F7.1 == "photographed",]
#sr89 <- mor[!is.na(mor$F8) & mor$F8.1 == "photographed",]

sr34$survival[!is.na(sr34$F4.1) & sr34$F4.1 == "photographed"] <- 0
sr34$survival[!is.na(sr34$F4.1) & (sr34$F4.1 == "dead" | sr34$F4.1 == "gone")] <- 1
# sr34$survival[!is.na(sr34$F4.1) & (sr34$F4.1 == "dead" | sr34$F4.1 == "gone" | sr34$F4.1 == "not found")] <- 1
sr45$survival[!is.na(sr45$F5.1) & sr45$F5.1 == "photographed"] <- 0
sr45$survival[!is.na(sr45$F5.1) & (sr45$F5.1 == "dead" | sr45$F5.1 == "gone")] <- 1
# sr45$survival[!is.na(sr45$F5.1) & (sr45$F5.1 == "dead" | sr45$F5.1 == "gone" | sr45$F5.1 == "not found")] <- 1
sr56$survival[!is.na(sr56$F6.1) & sr56$F6.1 == "photographed"] <- 0
sr56$survival[!is.na(sr56$F6.1) & (sr56$F6.1 == "dead" | sr56$F6.1 == "gone")] <- 1
# sr56$survival[!is.na(sr56$F6.1) & (sr56$F6.1 == "dead" | sr56$F6.1 == "gone" | sr56$F6.1 == "not found")] <- 1

#sr67$survival[!is.na(sr67$F7.1) & sr67$F7.1 == "photographed"] <- 0
#sr67$survival[!is.na(sr67$F7.1) & (sr67$F7.1 == "dead" | sr67$F7.1 == "gone")] <- 1

#sr78$survival[!is.na(sr78$F8.1) & sr78$F8.1 == "photographed"] <- 0
#sr78$survival[!is.na(sr78$F8.1) & (sr78$F8.1 == "dead" | sr78$F8.1 == "gone")] <- 1

#sr89$survival[!is.na(sr89$F9.1) & sr89$F9.1 == "photographed"] <- 0
#sr89$survival[!is.na(sr89$F9.1) & (sr89$F9.1 == "dead" | sr89$F9.1 == "gone")] <- 1



# add area squared for quadratic models
sr34$F3_2 <- sr34$F3^2
sr45$F4_2 <- sr45$F4^2
sr56$F5_2 <- sr56$F5^2

head(sr34)

# combine
srt <- data.frame(
  species_code = c(sr34$species_code, sr45$species_code, sr56$species_code),
  form = c(sr34$form, sr45$form, sr56$form),
  action = c(sr34$F4.1, sr45$F5.1, sr56$F6.1),
  colony_id = c(sr34$colony_id, sr45$colony_id, sr56$colony_id),
  year = c(rep("3_4", length(sr34$species_code)), rep("4_5", length(sr45$species_code)), rep("5_6", length(sr56$species_code))),
  area = c(sr34$F3, sr45$F4, sr56$F5),
  area_2 = c(sr34$F3_2, sr45$F4_2, sr56$F5_2),
  survival = c(sr34$survival, sr45$survival, sr56$survival
  )
)

head(srt)

# combine Andrew
dat_andrew <- data.frame(
  species_code = c(sr34$species_code, sr45$species_code, sr56$species_code, sr67$species_code, sr78$species_code, sr89$species_code),
  form = c(sr34$form, sr45$form, sr56$form, sr67$form, sr78$form, sr89$form),
  action = c(sr34$F4.1, sr45$F5.1, sr56$F6.1, sr67$F7.1, sr78$F8.1, sr89$F9.1),
  year = c(rep("3_4", length(sr34$species_code)), rep("4_5", length(sr45$species_code)), rep("5_6", length(sr56$species_code)), rep("6_7", length(sr67$species_code)), rep("7_8", length(sr78$species_code)), rep("8_9", length(sr89$species_code))),
  survival = c(sr34$survival, sr45$survival, sr56$survival, sr67$survival, sr78$survival, sr89$survival)
)

dat <- data.frame(srt)

dat_sav <- data.frame(colony_id=dat$colony_id, species_code=dat$species_code, growth_form=dat$form, log_area_cm2=dat$area, mortality=dat$survival, year=dat$year)
dat_sav <- dat_sav[!is.na(dat_sav$mortality),]
dat_sav <- data.frame(lapply(dat_sav, as.character), stringsAsFactors=FALSE)
dat_sav$growth_form[dat_sav$growth_form == "table"] <- "tabular"
dat_sav$growth_form[dat_sav$growth_form == "branching"] <- "staghorn"
#dat_sav$species[dat_sav$species == "AC"] <- "Acropora cytherea"
#dat_sav$species[dat_sav$species == "AH"] <- "Acropora hyacinthus"
#dat_sav$species[dat_sav$species == "AD"] <- "Acropora cf digitifera"
#dat_sav$species[dat_sav$species == "AS"] <- "Acropora humilis"
#dat_sav$species[dat_sav$species == "AL"] <- "Acropora spathulata"
#dat_sav$species[dat_sav$species == "AI"] <- "Acropora intermedia"
#dat_sav$species[dat_sav$species == "AR"] <- "Acropora robusta"
#dat_sav$species[dat_sav$species == "AM"] <- "Acropora millepora"
#dat_sav$species[dat_sav$species == "AN"] <- "Acropora nasuta"
#dat_sav$species[dat_sav$species == "GR"] <- "Goniastrea retiformis"
#dat_sav$species[dat_sav$species == "GP"] <- "Goniastrea pectinata"

# dat_sav$survival <- abs(as.numeric(dat_sav$survival) - 1)

dat_sav <- dat_sav[order(dat_sav$species),]




write.csv(dat_sav, "output/data_mortality.csv", row.names = FALSE, quote=FALSE)



# Find most recent data file from fractal (i.e., fractal_YYYYMMDD.txt)
csf1 <- read.delim("output/output_take_15.txt", as.is = T)
csf2 <- read.delim("output/output_tps_15.txt", as.is = T)

names(csf2) <- names(csf1)
csf <- rbind(csf1, csf2)

# Group E1 and E2
csf$exp[csf$exp == "E2"] <- "E1"
csf$exp <- factor(csf$exp)

# Remove encrusting and folecious growth forms
csf <- csf[csf$gf != "EN",]
csf <- csf[csf$gf != "FL",]
csf <- csf[csf$gf != "CN",]

csf$gf[csf$gf == "TB"] <- "table"
csf$gf[csf$gf == "DT"] <- "digitate"
csf$gf[csf$gf == "MS"] <- "massive"
csf$gf[csf$gf == "CB"] <- "corymbose"
csf$gf[csf$gf == "BR"] <- "branching"

csf$gf <- factor(csf$gf)

# THESE NEED CHECKING FOR WHY NA's OCCUR!!
csf <- csf[is.na(csf$frac_W) == F,]
csf <- csf[is.na(csf$csf_W) == F,]
csf <- csf[-890,]  # presky outlier -- huge branching with very small base -- LI_E5_016
csf <- csf[!is.na(csf$area_T),]

csf$lcsf <- log(csf$csf_W)
csf$larea <- log(csf$area_T)

csf_temp <- csf[csf$exp=="E4" | csf$exp=="E5",]

csf_sav <- data.frame(growth_form=csf_temp$gf, log_area_cm2=csf_temp$larea, log_csf=csf_temp$lcsf)
csf_sav <- data.frame(lapply(csf_sav, as.character), stringsAsFactors=FALSE)
csf_sav$growth_form[csf_sav$growth_form == "table"] <- "tabular"
csf_sav$growth_form[csf_sav$growth_form == "branching"] <- "arborescent"
csf_sav <- csf_sav[order(csf_sav$growth_form),]
write.csv(csf_sav, "output/data_csf.csv", row.names = FALSE, quote=FALSE)
