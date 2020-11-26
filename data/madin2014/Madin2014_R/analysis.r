library(mgcv)



# mortality analysis
rm(list=ls())
source("R/data_prep.r")
source("R/functions.r")

# Get area/gf to csf conversions for colonies at trimodal site
csf_mod <- glm(lcsf ~ larea * gf, data = csf[csf$exp=="E4" | csf$exp=="E5",])
co <- coef(csf_mod)
summary(csf_mod)

write.table(round(coef(summary(csf_mod)), 3), "output/growthform_csf.csv", sep = ",")

# Comparing species and growth forms
# spp_q = glm(survival ~ area:species_code + area_2:species_code + species_code, family = binomial, data = dat)
# summary(spp_q)
#
# for_q = glm(survival ~ area:form + area_2:form + form, family = binomial, data = dat)
# summary(for_q)
#
# jAICc(spp_q)
# jAICc(for_q)

# Large robusta that was "gone" because whole reef was removed...
# missing large robusta
# http://acropora.bio.mq.edu.au:3000/corals/251
# dat <- dat[c(-803),]

br <- run_test("branching")
br
write.table(round(coef(summary(br$r)), 3), "output/branching.csv", sep = ",")

tb <- run_test("table")
tb
write.table(round(coef(summary(tb$gfMOD)), 3), "output/table.csv", sep = ",")

cb <- run_test("corymbose")
cb
write.table(round(coef(summary(cb$gfMOD)), 3), "output/corymbose.csv", sep = ",")

dt <- run_test("digitate")
dt
write.table(round(coef(summary(dt$gfMOD)), 3), "output/digitate.csv", sep = ",")

ms <- run_test("massive")
ms
write.table(round(coef(summary(ms$gfMOD)), 3), "output/massive.csv", sep = ",")

# Plotting

source("R/make_figs.r")

# growth, longevity

gro <- read.table("data/hya_growth.txt", header= TRUE) / 10000 # m^2
gro <- gro[gro[,2] > 0,][-33,] # colonies that dies and an outlier
gro$linit <- log(gro$Initial)
gro$lnext <- log(gro$Next)

mod_gro <- lm(lnext ~ linit, data=gro)


x <- log(0.005)
p <- c(1)
for (i in 1:10) {
	x <- c(x, predict(mod_gro, list(linit = x[i])))
	p <- c(p, p[i] * predict(tb$gfMOD, list(area = x[i+1], area_2 = x[i+1]^2), type="response"))
}

tb$gfMOD


# For Andrew's mortality paper

species_mort <- function(sp, yr) {
  per <- sum(dat_andrew$survival[dat_andrew$species_code == sp & dat_andrew$year==yr], na.rm = T) / sum(!is.na(dat_andrew$survival[dat_andrew$species_code == sp & dat_andrew$year==yr]))
  num <- sum(!is.na(dat_andrew$survival[dat_andrew$species_code == sp & dat_andrew$year==yr]))
  return(per)
}

spp <- c("AC", "AD", "AH", "AI", "AL", "AM", "AN", "AR", "AS", "GP", "GR")

keep <- data.frame(spp=spp, f3_4=as.vector(sapply(spp, species_mort, yr="3_4")))
keep <- data.frame(keep, f4_5=as.vector(sapply(spp, species_mort, yr="4_5")))
keep <- data.frame(keep, f5_6=as.vector(sapply(spp, species_mort, yr="5_6")))
keep <- data.frame(keep, f6_7=as.vector(sapply(spp, species_mort, yr="6_7")))
keep <- data.frame(keep, f7_8=as.vector(sapply(spp, species_mort, yr="7_8")))
keep <- data.frame(keep, f8_9=as.vector(sapply(spp, species_mort, yr="8_9")))
keep[,2:7] <- round(keep[,2:7], 3) * 100
write.csv(keep, "output/yearly_mort.csv")
