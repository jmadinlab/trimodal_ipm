
rm(list = ls())
library("ggrepel")
library("grid")
library("png")

#######################################
# LOAD SCRIPTS
#######################################
source("R/functions.R")
source("R/data_prep.R")

#source("R/goni.R") # removes smalls to fix growth curve
sdat<-sdat[!sdat$colony_id==287,]
gdat<-gdat[!gdat$colony_id==287,]

cutoff_diam<-6
cutoff_a <- log10(a_func(cutoff_diam/2/100))
gdat <- subset(gdat, area > cutoff_a  )
gdat <- subset(gdat, area_next > cutoff_a  )
dat <- dat[log10(dat$area_cm2/10000)>cutoff_a,]

source("R/params.R")
params[,c("spp", "g.slp")]

#######################################
# GONIASTREA
#######################################

head(dat)
grdat <- subset(dat, spp=="Gre")
grdat <- aggregate(area_cm2~colony_id+fieldtrip_id+year, grdat, mean)
ggplot(grdat, aes(year, area_cm2))+
geom_line(aes(group=colony_id))+
geom_hline(yintercept=40, col="red", linetype="dotted")+
geom_hline(yintercept=10, col="red", linetype="dotted")+
geom_hline(yintercept=1, col="red", linetype="dotted")+
geom_text(data=grdat[grdat$year==2014,], aes(year+0.1, area_cm2, label=colony_id), size=2)+
geom_point(data=grdat[grdat$colony_id==349,], col="red")+
#geom_smooth()+
scale_y_log10()+
theme_classic()

dat[dat$colony_id==356,]


head(gdat)
grgdat <- subset(gdat, spp=="Gre")
grgdat <- grgdat[!grgdat$colony_id==356,]
grgdat <- grgdat[!grgdat$colony_id==331,]
grgdat <- grgdat[!grgdat$colony_id==351,]
grgdat <- grgdat[!grgdat$colony_id==349,]
ggplot(grgdat, aes( area, area_next))+
geom_text(data=gdat[gdat$spp=="Gre",],aes(label=colony_id), col="red")+
geom_text(aes(label=colony_id))+
geom_smooth(method="lm", formula=y~x, se=F)+
#scale_y_log10()+scale_x_log10()+
geom_abline(slope=1)
summary(lm(area_next~area, grgdat))


gdat <- gdat[!gdat$colony_id==356,]
gdat <- gdat[!gdat$colony_id==331,]
gdat <- gdat[!gdat$colony_id==351,]
gdat <- gdat[!gdat$colony_id==349,]

#sdat <- sdat[!sdat$colony_id==356,]
#sdat <- sdat[!sdat$colony_id==331,]
#sdat <- sdat[!sdat$colony_id==351,]
#sdat <- sdat[!sdat$colony_id==349,]

#######################################
# AL
#######################################
# Is AL (middle man) classed as a rare or common?
# i.e. do we compare 2 common with AM (Rare) 
# or 2 rare with AN (common)
# seems as though it's closer to AM and should be Rare
# but morphologically more similar to AM, so good pair

# params$abundance_pair[params$spp=="AL"] <- "Rare"
params2<-rbind(params, params[6,]) #duplicate AN (or AM=6)
params2$morph<-as.character(params2$morph)
params2$morph[c(7, 12)]<-c("corymbose_2","corymbose_2")#AN/AM
comp<-dcast(params2, morph~abundance_pair, value.var="spp")
comp$AC<-params2$abundance_05[match(comp$Common, params2$spp)]
comp$AR<-params2$abundance_05[match(comp$Rare, params2$spp)]
comp$diff <- comp$AC/comp$AR
comp


#######################################
# IPM FUNCTIONS
#######################################

#######################################
# QUADRATIC IN SURVIVAL???
#######################################
#-------------------------------survival
s.x <- function(x) { 
	u <- params$s.int[params$spp==sp] + 
	  params$s.slp[params$spp==sp] * x + 
	    params$s.slp.2[params$spp==sp] * x^2
  return(inv.logit(u)) 
}

#------------------------------- growth & partial mortality
p.yx <- function(y, x) {
  # x <- -3
  g <- a_func(r_func(10^x) + params$r.int[params$spp==sp] )
  #+ 1.96 * params$r.err[params$spp==sp])
  temp <- 10^y / g  # proportion of max reached. 
  temp[temp > 1] <- 1
  dnorm(logit(1 - temp), params$p.int[params$spp==sp] + x * params$p.slp[params$spp==sp], params$p.sig[params$spp==sp])
}

#------------------------------- reproduction
 r.yx <- function(y, x) {	
 	mat<- inv.logit(params$m.int[params$spp==sp] + 
 	  params$m.slp[params$spp==sp] *x)
 	fec<- exp(params$f.int[params$spp==sp] + 
 	  params$f.slp[params$spp==sp] *x) 
   #siz<- rnorm(y,mean=params$rec.size[params$spp==sp], sd=0.05) 
   out <- (rec* mat * fec)
   out[x < rec.size | y >= rec.size] <- 0 #if x is below recruitment size
   return(out)
   } 

 #------------------------------- kernel
pmort<-T

bigmatrix <- function() {
    if (pmort) {
    G <- h * outer(y, y, p.yx)
  } else {
    G <- h * outer(y, y, g.yx)
  }
  G <- t(t(G) / apply(G, 2, sum))
  S <- s.x(y)
  P <- G 
  for(i in 1:n) P[,i]=G[,i]*S[i]
  R <- h * outer(y, y, r.yx) 
  #R <- h * outer(y, pmin(y, smax), r.yx) #  ceiling
  K <- P + R
  lam <- Re(eigen(K)$values[1])
	w <- abs(Re(eigen(K)$vectors[,1])) 
	v <- abs(Re(eigen(t(K))$vectors[,1]))
	return(list(K=K, lam=lam, w=w, v=v, G=G, S=S, R=R, P=P)) }
		
		
#######################################
# MESH AND BOUNDARIES
#######################################

params$rec.size <- log10(pi*(params$r.int*(10/12))^2)
max.size <- 1 
n <- 100

# BEGINNING OF LIFE DEFINED AS FIRST NON-PROPAGULE STAGE - BURNS et al 2010

mesh <- function(){
	#min.size <- -4
	min.size <- rec.size
	b <- seq(min.size, max.size, length=n)
	h <- b[2] - b[1]
	b <- c(min(b)-(2*h), min(b)-h, b)
	y <- 0.5 * (b[1:n]+b[2:(n+1)])
	I <- y >= rec.size
	return(list(b=b,h=h,y=y,I=I, rec.size = rec.size))}

#######################################
# PLOT IPMS
#######################################

par(mfcol=c(2, 6))

lam_const <- NULL
ipm.k.list <- list()
ipm.p.list <- list()
ipm.r.list <- list()
y.list <- list()

for (sp in spp) {
	rec <- 1
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	sub<-gdat[gdat$spp==sp,]
    mod <- bigmatrix()
   # plot
	image(y, y, t(mod$P)^0.3)   
	points(sub$area, sub$area_next, cex=0.25)
	title(sp, line=-1)
	abline(0, 1, lty=2)
	y.list[[sp]] <- y 
	ipm.k.list[[sp]] <- mod$P
	ipm.p.list[[sp]] <- mod$P
	ipm.r.list[[sp]] <- mod$R
	}

# Figure S2
source("R/fig.S2.R")	
#fig.s2

#######################################
# MULTIPLE REC
#######################################

rec.x <- 10^(seq(-6,-1,1))

store<-data.frame()
for (sp in spp) {
	for (rec in rec.x){
	  rec.size <- params$rec.size[params$spp==sp]
	  #rec.size <- max(params$rec.size)
	  h <- mesh()$h
	  y <- mesh()$y
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  store<-rbind(store, data.frame(sp=sp, rec=rec, lam=bigmatrix()$lam))
   } 
  }
 

params$rec.fam <- ifelse(params$morphology=="massive", 10^-4, 10^-3)



#######################################
# ELASTICITY & IPM MEASURES
#######################################
 # eigen-things can be combined to obtain the sensitivity and elasticity matrices.

s.list <- list()
eK.list <- list()
eR.list <- list()
eP.list <- list()
demovals <- NULL
sizevals <- NULL

for (sp in spp) {
	#sp <- "GR"
	# run IPM
	rec <- params$rec.fam[params$spp==sp]
	#rec<-1*10^-4
	rec.size <- params$rec.size[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
	K <- bigmatrix()$K
	P <- bigmatrix()$P
	R <- bigmatrix()$R
	
	# population growth
	lam <- Re(eigen(K)$values[1])
	
	# sable size dist /right eigenvec
	w.eigen<-Re(eigen(K)$vectors[,1]) 
	stable.size <- w.eigen/sum(w.eigen)
	
	# reprodutive values /left eigenvec (size contributions)
	v.eigen<-Re(eigen(t(K))$vectors[,1]) 
	repro.val <- v.eigen/v.eigen[1] # relative
	
	# reproductive val * stable size
	v.dot.w<-sum(stable.size*repro.val)*h 
	
	# sensitivity matrix
	sens<-outer(repro.val,stable.size)/v.dot.w  
	
	# elasticity matrices
	K.elas<-matrix(as.vector(sens)*(as.vector(K)/h)/lam, nrow=n)
	
	# survival elasticity
	P.elas<-(P/h)*sens/lam
	eP=sum(P.elas)*h^2

	# reproduction elasticity
	R.elas<-(R/h)*sens/lam
	eR=sum(R.elas)*h^2
	
	# Net reproductive rate from IPMbook monocarp
	N <- solve(diag(n)-P)
	R0 <- abs(eigen(R %*% N)$values[1])
	
	# Generation time
	GT <- log(R0)/log(lam)
	
	demovals <- rbind(demovals, data.frame(spp=sp, eR, eP, R0, GT))
	sizevals <- rbind(sizevals, data.frame(spp=sp, area=y, stable.size, repro.val, v.dot.w))
		s.list[[sp]] <- sens
		eR.list[[sp]] <- R.elas
		eP.list[[sp]] <- P.elas
		eK.list[[sp]] <- K.elas
		}

demovals
head(sizevals)
params[,colnames(demovals)]<- demovals[match(demovals$spp, params$spp),]

#sum(analyses()$eK)*h^2
lapply(eK.list, function(x){sum(x)*h^2})
rowSums(demovals[,c("eP","eR")])

plot_grid(
ggplot(sizevals, aes(x=area, y=repro.val, col=spp))+geom_line()+
scale_colour_manual(values=cols)+
scale_y_log10()+guides(col="none"),
ggplot(sizevals, aes(x=area, y=stable.size, col=spp))+geom_line()+
scale_colour_manual(values=cols)+guides(col="none"),
ggplot(params, aes(reorder(spp, -eR), eR, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none"), 
ggplot(params, aes(reorder(spp, -R0), R0, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none")+scale_y_log10(),
ggplot(params, aes(reorder(spp, -GT), GT, fill=spp))+geom_bar(stat="identity")+scale_fill_manual(values=cols)+guides(fill="none")+scale_y_log10()
)

#######################################
# R.ELASTICITY vs LAMBDA
#######################################
store[store$rec==max(rec.x),]
#params$lam_4 <- store[store$rec==max(rec.x),"lam"]
#params$lam_1 <- store[store$rec==min(rec.x),"lam"]
params$lam_4 <- store[store$rec==0.01,"lam"]
params$lam_1 <- store[store$rec==0.0001,"lam"]


  ggplot(data=params, aes(x=eR, y=log(lam_4)-log(lam_1)))+
geom_smooth(method="lm", formula=y~x, size=0.1, col="black")+
 geom_point(aes(fill=spp), shape=21, stroke=0.1, size=2)+
 scale_fill_manual(values=cols)+
 guides(fill="none", col="none")
#Pvals <- big$P / h
#P.elas <- Pvals * K.sens / lambda

#######################################
# Generation VS R0
#######################################

ggplot(demovals, aes(GT, R0))+geom_text(aes(label=spp))+scale_x_log10()
# could calculate mean size of flowering plants. *stable size by p_maturity

source("R/fig.3.R")
fig.3

#######################################
# min.size/max.size-lamda relationship
#######################################
  
  store<-data.frame()
  minx<-c(-6,-5,-4,-3,-2,-1)
  for (sp in spp) {
    for(min.size2 in minx){
      min.size2<-min.size2
      max.size <- max.size #2
      #rec.size <- rec.size
      n <- 100
      b <- min.size2 + c(0:n) * (max.size - min.size2)/n
      y <- 0.5 * (b[1:n]+b[2:(n+1)])
      h <- y[2] - y[1]
      store<-rbind(store, data.frame(sp=sp, min.size=min.size2, lam=bigmatrix()$lam)) }}
      
      head(store)
      
  ggplot(data=store, aes(x=min.size, y=lam, colour=sp))+
    geom_line()+geom_hline(yintercept=1)+
    geom_vline(xintercept=min.size2, linetype="dotted")+
    ggtitle("relative differences")+
        geom_segment(data=params, aes(x=rec.size, xend=rec.size, y=Inf, yend=Inf), linetype="dotted")+
    scale_colour_manual(values=as.character(params$cols))
  
  
  store<-data.frame()
  maxx <- c( 0,1,2,3, 4,5)
  for (sp in spp) {
    for(max.size2 in maxx){
      max.size2<-max.size2
      min.size <- min.size #2
      rec.size <- -2.7
      n <- 100
      b <- min.size + c(0:n) * (max.size2 - min.size)/n
      y <- 0.5 * (b[1:n]+b[2:(n+1)])
      h <- y[2] - y[1]
      store<-rbind(store, data.frame(sp=sp, max.size=max.size2, lam=bigmatrix()$lam)) }}
  ggplot(data=store, aes(x=max.size, y=lam, colour=sp))+
    geom_line()+geom_hline(yintercept=1)+
    geom_vline(xintercept=max.size, linetype="dotted")+
    #lims(x=c(0, 0.00001), y=c(0,4))+
    ggtitle("relative differences")+
    scale_colour_manual(values=as.character(params$cols))
  
  



#######################################
# FIT REC TO SIZE STRUCTURE
#######################################
fitrec<-NULL

par(mfrow=c(4,3))
for (sp in spp) {
#sp <- "AC"
  rec.size <- params$rec.size[params$spp==sp]
  #b <- seq(rec.size, max.size, length=n)
  b <- seq(-4, max.size, length=n)
  h <- b[2] - b[1]
  b <- c(min(b)-h, b)
  y <- 0.5 * (b[1:n]+b[2:(n+1)])
  I <- y >= rec.size
    
    size.dist <- hist(ss$area[ss$spp==sp & ss$area > rec.size & ss$area < max.size], breaks=b)
    rec.fit <- mle(rec.ll, start = list(x = 0.1), method="Brent", lower = 0, upper=1)
    optimise(rec.ll, c(0,1))
    rec<-coef(rec.fit)
    mod<-bigmatrix()
    lines(y, mod$w/h, col="red")
    
    
    fitrec<-rbind(fitrec,data.frame(x=y[I], y=(mod$w[I]/sum(mod$w[I])/h), rec, spp=sp, mn=min.size, mx=max.size)) }
  
  params$rec<-unique(fitrec[,c("spp","rec")])$rec
params$rec




