
rm(list = ls())
library("grid")
library("png")

#######################################
# LOAD SCRIPTS
#######################################
source("R/functions.R")
source("R/data_prep.R")

sdat<-sdat[!sdat$colony_id==287,]
gdat<-gdat[!gdat$colony_id==287,]

source("R/params.R")
params[,c("spp", "g.slp")]

#######################################
# IPM FUNCTIONS
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

max.size <- 1 
n <- 100

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
# RECRUITMENT LIT REVIEW
#######################################

# RECRUITMENT MODEL, Assume:
# (A) closed system at trimodal
# (B) 11 species constitute all acroporids/favids

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

ggplot(tiles, aes(x=(N_m2_year/eggs), y=Family))+
geom_boxplot()+scale_x_log10()

# recruit size
ggplot(rsize, aes(x=recsize_1yr, y=Family))+geom_boxplot()
rsize.fam <- aggregate(recsize_1yr~Family, rsize, median)

#######################################
# ESTIMATE REC/LAM
#######################################

# recsize growth
params$rsize.gr <- log10(pi*(params$r.int*(9/12))^2)
r.limit <- log10(pi*((10/100)/2)^2)
params$rsize.gr <- ifelse(params$rsize.gr>r.limit,r.limit,params$rsize.gr)

# recruitment family
params$rec.fam <- rec.fam$p.rec[match(params$family, rec.fam$Family)]

# SELECT!!!
params$rec <- params$rec.fam
params$rec.size <- params$rsize.gr

lam.est <- NULL
for (sp in spp) {
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
mod <- bigmatrix()
lam.est<-rbind(lam.est, data.frame(sp, lam=bigmatrix()$lam))}
 lam.est
params$lam.est <- lam.est$lam


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
#source("figs/fig.S2.R")	
#fig.s2


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
	#sp <- "Gre"
	# run IPM
	rec <- params$rec[params$spp==sp]
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
	
	# reprodutive values /left eigenvec (size contributions) # Gre goes wrong
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

summary(s.list)

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
# COMPARE MORPHS
#######################################

# Is AL (middle man) classed as a rare or common?
# params$abundance_pair[params$spp=="Asp"] <- "Rare"
params2<-rbind(params, params[6,]) #duplicate AN, 7 (or AM=6)
params2$morph<-as.character(params2$morph)
params2$morph[c(7, 12)]<-c("corymbose_2","corymbose_2")#AN/AM
comp<-dcast(params2, morph~abundance_pair, value.var="spp")
comp$AC<-params2$abundance_05[match(comp$Common, params2$spp)]
comp$AR<-params2$abundance_05[match(comp$Rare, params2$spp)]
comp$diff <- comp$AC/comp$AR
comp

colsC <- cols[names(cols) %in% comp$Common]
names(colsC)<-comp$morph[match(names(colsC),comp$Common)]
colsC

#######################################
# MORPH PAIRS
#######################################
	
	
comp$rec <- params2$rec[match(comp$morph, params2$morph)] 
comp$rec.size <- params2$rec.size[match(comp$morph, params2$morph)] 
comp$rec.size.cm <- sqrt(((10^comp$rec.size)*10000)/pi)*2		
comp$lamC <- params$lam.est[match(comp$Common, params$spp)]	
comp$lamR <- params$lam.est[match(comp$Rare, params$spp)]
comp$lamdiff <- comp$lamC - comp$lamR
comp$logdiff <- log(comp$lamC) - log(comp$lamR)
comp$doub.time <- log(2)/log(1+comp$lamdiff)
comp$difftime <- log(comp$diff)/log(1+comp$lamdiff)
		
ggplot(comp[comp$difftime>0,], aes(x=difftime, y=reorder(morph, -difftime)))+geom_bar(stat="Identity")



