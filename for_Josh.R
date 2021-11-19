
rm(list = ls())

library("MASS") 
library("reshape2")
library("ggplot2")
library("cowplot")
library("stats4")

source("R/functions.R")
source("R/data_prep.R")

gdat <- subset(gdat, area > log10((pi*(3/2)^2)/10000))
sdat <- subset(sdat, area > log10((pi*(3/2)^2)/10000))
fec <- subset(fec, area > log10((pi*(3/2)^2)/10000))

source("R/params.R")

params[params$spp=="Gre", c("p.int", "p.slp","p.sig")] <- params[params$spp=="Gpe", c("p.int", "p.slp","p.sig")]

params[params$spp=="Gre", c("g.int", "g.slp","g.var")] <- params[params$spp=="Gpe", c("g.int", "g.slp","g.var")]

#######################################
# RECRUITMENT LIT REVIEW
#######################################

#  Assume: (A) closed system, (B) 11 species = all acroporids/favids
# average size
size.av <- aggregate(area~spp, ss[!is.na(ss$spp),], mean)
params$size.av <- size.av$area[match(params$spp, size.av$spp)]
density <- params$abundance_05/2700# 275*1*10m 
area.m2 <- log10(density*(10^(params$size.av)))
params$fec.m2 <- exp(params$f.int+params$f.slp*area.m2) 
fam <- aggregate(fec.m2~family, params, sum)
tiles$eggs<-fam$fec.m2[match(tiles$Family, fam$family)]
tiles$N_m2_year[tiles$N_m2_year==0] <-1
tiles$p.rec <- tiles$N_m2_year/tiles$eggs
rec.fam <- aggregate(p.rec~Family, tiles, median)
rec.fam

# recsize
rec.size.const <- log10(pi*((5/2)/100)^2) # 5cm diameter
params$rec.size <- log10(pi*(params$r.int*(12/12))^2)
agg <- aggregate(rec.size~morphology, params, mean)
params$rec.size <- agg$rec.size[match(params$morphology, agg$morphology)]
r.limit <- log10(pi*((10/100)/2)^2)
params$rec.size <- ifelse(params$rec.size > r.limit,r.limit,params$rec.size)

# recruitment 
rec.const <- 10^-3
params$rec <- rec.fam$p.rec[match(params$family, rec.fam$Family)]
#params$rec <- ifelse(params$family=="Merulinidae", 5.900139e-05, 2*10^-3)


#######################################
# IPM MESH AND BOUNDARIES
#######################################

max.size <- 0.3 # max(ss$area)
n <- 100
mesh <- function(){
	min.size <- rec.size
	b <- seq(min.size, max.size, length=n)
	h <- b[2] - b[1]
	b <- c(min(b)-(2*h), min(b)-h, b)
	y <- 0.5 * (b[1:n]+b[2:(n+1)])
	I <- y >= rec.size
	return(list(b=b,h=h,y=y,I=I, rec.size = rec.size))}
	

#######################################
# ESTIMATE LAM
#######################################

lam.est <- NULL
lam.const <- NULL
for (sp in spp) {
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
lam.est<-rbind(lam.est, data.frame(sp, lam=bigmatrix()$lam))
rec.size <- rec.size.const
rec <- rec.const
h <- mesh()$h
y <- mesh()$y
lam.const<-rbind(lam.const, data.frame(sp, lam=bigmatrix()$lam))
}
params$lam.est <- lam.est$lam
params$lam.const <- lam.const$lam

ggplot(params, aes(reorder(spp, -lam.est), lam.est))+geom_bar(stat="identity", aes(fill=spp))+scale_fill_manual(values=cols)+geom_hline(yintercept=1)

params[,c("spp","lam.est", "g.slp")]


#######################################
# REC FIT TO SIZE STRUCTURE
#######################################

rec.ll <- function(x) {
  cnt <- size.dist$count[II] # non-recruits
  rec <<- x[1] 
  mod <- bigmatrix()
  eig.vec <- mod$w[II]/sum(mod$w[II])
  return(-sum(cnt * log(eig.vec), na.rm=TRUE)) } # log-likelihood 

rec.ss <- function(x) {
  cnt <- size.dist$count[II] # non-recruits
  cnt <- cnt / sum(cnt)
  rec <<- x[1] 
  mod <- bigmatrix()
  eig.vec <- mod$w[II]/sum(mod$w[II])
  return(sum((cnt - eig.vec)^2)) } # log-likelihood 
  
  
# ss_vec <- 10^seq(-3, -7, -0.01)
# store <- data.frame()
# for (ss in ss_vec) {
	
 # cnt <- size.dist$count[II] # non-recruits
 # cnt <- cnt / sum(cnt)
 # rec <- ss 
 # mod <- bigmatrix()
 # eig.vec <- mod$w[II]/sum(mod$w[II])
 # store <- rbind(store, data.frame(rec=ss, ss=sum((cnt - eig.vec)^2)))

# }  

# plot(ss ~ rec, store)
# store[which.min(store$ss),]

dev.off()  
par(mfcol=c(4, 3), mar=c(3,3,1,1))

n <- 12

params$rec_fit <- NA
params$lam_fit <- NA 

for (sp in spp) {
#sp <- "Asp"
  # Model
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
b <- mesh()$b 
I <- mesh()$I
II <- I
II[1:6] <- FALSE

max.size <- 1
size.dist <- hist(ss$area[ss$spp==sp & ss$area > rec.size & ss$area < max.size], breaks=b, plot=FALSE)

# hist(ss$area[ss$spp==sp & ss$area > rec.size & ss$area < max.size], breaks=b)
#, lower = 0, upper=10000 removed
# x is a startin rec value?
# mle takes a function and finds the most likely value of x ? 

	rec.fit <- optimise(rec.ss, c(0, 1))
	rec <- rec.fit$minimum
	rec
	
	#rec.fit <- mle(rec.ll, start = list(x = 0.001), method="Brent", lower = 0, upper=1 )
	rec.fit <- optimise(rec.ll, c(0, 100))
	rec <- rec.fit$minimum
	rec

#rec <- 0.00001
  mod <- bigmatrix()
	#image(y, y, log(t(mod$K)))
	#title(sp, line=-1)
	#abline(0, 1, lty=2)
mod$lam

	hist(ss$area[ss$spp==sp & ss$area > rec.size & ss$area < max.size], breaks=b, freq=FALSE, main="", ylim=c(0, 2))
	lines(y[II], (mod$w[II]/sum(mod$w[II]))/h/2, col="red")
	abline(v=y[45], lty=2)
	title(sp, line=-1)
text(0,1, round(mod$lam, 3))


  params$lam_fit[params$spp==sp] <- mod$lam
  params$rec_fit[params$spp==sp]<- rec
}

params$rec_fit


ggplot(params, aes(reorder(spp, -lam_fit), lam_fit))+geom_bar(stat="identity", aes(fill=spp))+scale_fill_manual(values=cols)+geom_hline(yintercept=1)



