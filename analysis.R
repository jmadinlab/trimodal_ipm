
rm(list = ls())
library("grid")
library("png")
library("MASS") 
library("reshape2")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("stats4")
inv.logit <- function(x) {exp(x)/(1+exp(x))}

#######################################
# LOAD DEMOGRAPHIC DATA 
#######################################

 source("R/1_data_prep.R") # process raw data
 source("R/2_params.R") # get model parameters
 source("R/3_traits.R") # demographic traits

#params <- read.csv("output/params.csv") #skip 1-2-3
spp<-params$spp[order(params$spp)]

# ordering
#################################### 
order <- c("Ahy","Acy","Ain","Aro","Ana", "Asp","Ami","Adi","Ahu", "Gre","Gpe")
params$spp <- factor(params$spp, levels=order)
labs <- params$species[order(params$spp)]
cols<-as.character(params$cols)
names(cols) <- params$spp
cols <- cols[order]

#######################################
# MORPHOLOGICAL PAIRS
#######################################

# Ana declined, do not compare
params2 <- params[!params$spp=="Ana",]
params2$morph <- params2$morphology
comp<-dcast(params2[!params2$spp=="Ana",], morph~abundance_pair, value.var="spp")
comp

colsC <- cols[names(cols) %in% comp$Common]
names(colsC)<-comp$morph[match(names(colsC), comp$Common)]

#######################################
# AVERAGE SIZE
#######################################

ggplot(ss, aes(x=area))+geom_histogram()+ facet_wrap(~species)

# average size in demo models
size.dat <- aggregate(area_cm2~spp, sdat[!is.na(sdat$spp),], mean)
params$size.dat <- size.dat$area_cm2[match(params$spp, size.dat$spp)]
params$size.dat <- log10(params$size.dat/10000)

# average size in size structure data
size.ss <- aggregate(area~spp, ss[!is.na(ss$spp),], mean)
params$size.ss <- size.ss$area[match(params$spp, size.ss$spp)]
10^params$size.ss
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
# DEMOGRAPHIC SPACE (PCA)
#######################################

rownames(params) <- params$spp
traits <- c("r.int","f.cm2","survcm", "av.surv","p_mort", "min.r", "f.colony")

pca<-prcomp(params[,traits], scale=T, center=T)
biplot(pca)
exp<-round(c(summary(pca)[[1]][1]^2/sum(summary(pca)[[1]]^2),summary(pca)[[1]][2]^2/sum(summary(pca)[[1]]^2)),3)*100 # explained variation

library("psych") # PCA matrix 
pairs.panels(log(params[,traits]), scale=T, cex.cor=2)

#######################################
# ABUNDANCE
#######################################

#abun.BT <- read.csv("output/abun.BT.csv")
#abun.LIT <- read.csv("output/abun.LIT.csv")

aggregate(N/10~spp, abun.BT, mean) # N per m2
aggregate(N/10~spp+year, abun.LIT, mean) # N per m

#######################################
# FIGURE 1
#######################################

source("figs/fig.1.R")
fig.1 # ggsave("figs/fig.1.jpeg", fig.1, width=15, height=9.5, units="cm", dpi = 600)

#######################################
# IPMS
#######################################

source("R/3_ipms.R") # IPM functions
max.size <- 0.3 # max(ss$area)
n <- 100

# recruit size: av max growth within morphs
agg <- aggregate(r.int~morphology, params, mean)
params$r.int2  <- agg$r.int[match(params$morphology, agg$morphology)]
params$rec.size <- log10(pi*(params$r.int2*0.85)^2)

# p_rec placeholder (same within families)
params$rec <- ifelse(params$family=="Merulinidae",0.000316, 0.00113)

#######################################
# RECRUITMENT SENSITIVITY I
#######################################	

rec.xDET <- 10^(seq(-5,-1, 0.1))
storeDET<-data.frame()
for (sp in spp) {
	for (rec in rec.xDET){
	  rec.size <- params$rec.size[params$spp==sp]
	  h <- mesh()$h
	  y <- mesh()$y
  storeDET<-rbind(storeDET, data.frame(spp=sp, rec=rec, lam=bigmatrix()$lam))
   } 
  }

ggplot(storeDET, aes(rec, lam))+
geom_line(aes(col=spp))+scale_x_log10()+
scale_colour_manual(values=cols)+
coord_cartesian(ylim=c(0.6,2))

#######################################
# ESTIMATE LAM
#######################################

# Mean-family recruitment based on transects 
source("R/4_recruitment.R")
params$log.rec.tran <- log10(params$rec.tran)
agg <- aggregate(log.rec.tran~family, params, mean)
params$rec <- 10^agg[match(params$family, agg$family),2]

lam.est <- data.frame()
for (sp in spp) {
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
lam.est<-rbind(lam.est, data.frame(spp=sp, lam.est=bigmatrix()$lam, R0=bigmatrix()$R0))
 }

params$lam.est <- lam.est$lam.est
params$R0 <- lam.est$R0
params$GT <- log(params$R0)/log(params$lam.est)
params[,c("spp","lam.est", "R0", "GT")]

#######################################
# BOOTSTRAP LAM
#######################################	

# source("R/bootstap.R")
boot <- read.csv("output/lam.range.csv")
boot$GT <- log(boot$r)/log(boot$lam)

ggplot(boot, aes(GT, log(lam)))+
geom_hline(yintercept=0)+
geom_point(shape=21, aes(col=spp), alpha=0.2)+
coord_cartesian(xlim=c(min(boot$GT),100))+
scale_colour_manual(values=cols)+scale_x_log10()+
geom_point(data=params, aes(GT, log(lam.est)))

#######################################
# FIG 2
#######################################	
 
source("figs/fig.2.R") 
fig.2 #ggsave("figs/fig.2.jpeg", fig.2, width=15, height=12.7, units="cm", dpi = 600)

#######################################
# SUPPLEMENT
#######################################

# Plot abundance / demographic models
source("figs/SUPPLEMENT/supp.fig2.R") # abundance
source("figs/SUPPLEMENT/supp.fig3.R") # models

# Plot IPMs
source("figs/SUPPLEMENT/supp.fig5.R")

#######################################
# COMPARE LAMS
#######################################

params2$lam.est <- params$lam.est[match(params2$spp, params$spp)]

morph.diff <- function(x){
	params2$x.temp <-params[match(params2$spp, params$spp), x]
	comp$x.common <- params2$x.temp[match(comp$Common, params2$spp)]
	comp$x.rare <- params2$x.temp[match(comp$Rare, params2$spp)]
	comp$x.diff <- comp$x.common-comp$x.rare
	comp$x.div <- comp$x.common/comp$x.rare
	comp
	}

comp$AC<-morph.diff("abundance_05")$x.common
comp$AR<-morph.diff("abundance_05")$x.rare
comp$abundiff <- morph.diff("abundance_05")$x.div
comp$lamC <- morph.diff("lam.est")$x.common
comp$lamR <- morph.diff("lam.est")$x.rare
comp$lamdiff <- morph.diff("lam.est")$x.diff

ggplot(comp)+ geom_bar(aes(y=reorder(morph, -lamdiff), x=lamdiff), stat="identity")

#######################################
# TIME NEEDED TO GET DIFFERENCES
#######################################

# projecting lamda from basic size dist.. 

proj.size <- NULL
for(sp in spp){
	rec <- params$rec[params$spp==sp] 
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	lam=bigmatrix()$lam
	v=bigmatrix()$v
	ngen <- 100
	npop <- 100
	v <- matrix(c(1, rep(0,99)))
pop  <- matrix(nrow=npop, ncol=ngen) # sizes/gens.
pop[,1] <- round(v * npop) # inital pop. size
for(j in 2:ngen){ pop[,j] <- bigmatrix()$K %*% as.matrix(pop[,j-1]) }
pop.long <- melt(pop, value.name = "NperSize")
colnames(pop.long)[1]<-"Size"
colnames(pop.long)[2]<-"gen"
psize<-	aggregate(NperSize~gen, data=pop.long, sum)
psize$n <- npop*lam^psize$gen # simpler method
proj.size <- rbind(proj.size, cbind(psize, spp=sp, abun=params$abundance_pair[params$spp==sp], morph=params$morph[params$spp==sp]))
}

proj.size <- proj.size[!proj.size$spp=="Ana",]	

proj <- dcast(morph+gen~abun, data=proj.size, value.var="NperSize")
proj$diff <- proj$Common/proj$Rare
proj$real <- comp$abundiff[match(proj$morph, comp$morph)]
gens <- aggregate(gen~morph, proj[proj$diff<proj$real,], max)
proj$maxgen <- gens$gen[match(proj$morph, gens$morph)] +1 
	
ggplot(proj, aes(gen, diff, col=morph))+
geom_line(linetype="dotted")+
geom_line(data=proj[proj$gen<=proj$maxgen,])+
geom_point(data=proj[proj$gen==proj$maxgen,])+ 
scale_y_sqrt()+scale_x_sqrt()+
coord_cartesian(ylim=c(0,60), xlim=c(0,100))

#######################################
# WITHIN VS BETWEEN
#######################################

Gen5 <- subset(proj.size, gen==5)
params$Gen5<-Gen5$NperSize[match(params$spp, Gen5$spp)]
Gen10 <- subset(proj.size, gen==10)
params$Gen10<-Gen10$NperSize[match(params$spp, Gen10$spp)]

dems <- c("lam.est", "Gen5",  "Gen10", "abundance_05", "rec")
vars <- NULL
for (t in c(traits, "eggC",dems)){
	#t <- "Gen15"
	params$t <- scale(params[,t])
mod <- aov(t~morphology, data=params)
summary(mod)
within <- summary(mod)[[1]]["Residuals","Sum Sq"]
between <- summary(mod)[[1]]["morphology","Sum Sq"]
pval <- summary(mod)[[1]]["morphology","Pr(>F)"]
vars <- rbind(vars, data.frame(t, within=100*(within/(within+between)), between=100*(between/(within+between)), pval))
}

ggplot(vars, aes(within,reorder(t, -within)))+geom_bar(stat="identity")

#######################################
# DRIVERS OF FITNESS DIFFS (within morphs)
#######################################

params_orig <- params # CAREFUL set original 
	
source("R/6_params_morph.R")
# p.morph <- read.csv("output/params_morph.csv") 
head(p.morph)

pars <- c("m.int", "m.slp","f.int","f.slp","g.int","g.slp","g.var", "s.int", "s.slp", "s.slp.2")
p.types <- data.frame(pars, type = c("reproduction","reproduction","reproduction","reproduction","growth","growth","growth", "survival","survival","survival"))
types <- unique(p.types$type)

test <- NULL
for (i in types) {
params <- params_orig 
params[,pars] <- p.morph[match(params$morphology, p.morph$morph), pars]
pars2 <- p.types$pars[p.types$type==i]
params[,pars2] <- params_orig[match(params$spp, params_orig$spp), pars2]
for (sp in spp) {
	#sp <- "Ahy"
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
test<-rbind(test, data.frame(spp=sp, param=i, lam=bigmatrix()$lam, morph=params$morphology[params$spp==sp], rec="normal"))
 } 
}

m.diff <- NULL
for (i in types) {
	#i <- "maturity"
	params <- params_orig 
	sub <- test[test$param==i, ]
	params$x <- sub$lam[match(params$spp, sub$spp)]
	m.diff <- rbind(m.diff, data.frame(morph=comp$morph, d=(morph.diff("x")$x.common - morph.diff("x")$x.rare), param=i, rec="normal"))
	}
head(m.diff)

sums <- aggregate(d~morph+rec, m.diff, sum)
sums$lamdiff <- comp$lamdiff[match(sums$morph, comp$morp)]
sums # compare d with the real lamdiff

ggplot()+
geom_bar(data=m.diff, aes(d, reorder(morph, -d), fill=param), stat="identity", position="stack", col="black")+geom_vline(xintercept=0)

#######################################
# FIGURE 3
#######################################

source("figs/fig.3.R")
fig.3 # ggsave("figs/fig.3.jpeg", fig.3, width=14, height=12, units="cm", dpi = 600)

#######################################
# REPRODUCTIVE TRADEOFF
#######################################

params <- params_orig 

params$fec1cm <- exp(params$f.int +  params$f.slp * log10(1)) /10000
params2$fec1cm <- params$fec1cm[match(params2$spp, params$spp)]

ggplot(params, aes(fec1cm, eggC))+geom_text(aes(label=spp))+geom_smooth(data=params[params$family=="Acroporidae",], method="lm")

source("figs/fig.4.R")
fig.4 # ggsave("figs/fig.4.jpeg", fig.4, width=5.5, height=8, units="cm", dpi = 600)

#######################################
# SUPPLEMENT
#######################################

# sensitivity analysis
source("R/7_sensitivity.R")
source("figs/SUPPLEMENT/supp.fig6.R")
