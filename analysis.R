
rm(list = ls())

#######################################
# LOAD SCRIPTS
#######################################
source("R/functions.R")
source("R/data_prep.R")
sdat<-sdat[!sdat$colony_id==287,]
gdat<-gdat[!gdat$colony_id==287,]
source("R/params.R")


#######################################
# IPM FUNCTIONS
#######################################

#------------------------------- growth
g.yx <- function(y, x) {
	dnorm(y, mean=params$g.int[params$spp==sp] + 
	  params$g.slp[params$spp==sp]*x,
	    sd=sqrt(params$g.var[params$spp==sp]))
	    }
	
#-------------------------------survival
s.x <- function(x) { 
	u <- params$s.int[params$spp==sp] + 
	  params$s.slp[params$spp==sp] * x + 
	    params$s.slp.2[params$spp==sp] * x^2
  return(inv.logit(u)) 
}

#------------------------------- growth & partial morality
p.yx <- function(y, x) {
  # x <- -4
  g <- a_func(r_func(10^x) + params$r.int[params$spp==sp] )
  #+ 1.96 * params$r.err[params$spp==sp])
  temp <- 10^y / g
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
   out <- rec* mat * fec 
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

params$rec.size <- log10(pi*params$r.int^2)
max.size <- 1 
n <- 100
rec <- 0.001
#rec <- 1*10^-8
#min.size <- -3.5 
#rec.size <- -2.21 
# Maximum size windows, constant or vary by species? 
#params$smax<-aggregate(area_cm2~spp, dat, max)$area
#params$smin<-aggregate(area_cm2~spp, dat, min)$area


#######################################
# FIT REC TO SIZE STRUCTURE
#######################################
fitrec<-NULL

#for (sp in spp) {
sp <- "GR"
  rec.size <- params$rec.size[params$spp==sp]
  b <- seq(rec.size, max.size, length=n)
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
    
    
   # fitrec<-rbind(fitrec,data.frame(x=y[I], y=(mod$w[I]/sum(mod$w[I])/h), rec, spp=sp, mn=min.size, mx=max.size)) }
  
  params$rec<-unique(fitrec[,c("spp","rec")])$rec



  lam_const<-c(lam_const,bigmatrix()$lam) 
  
  #}


#######################################
# PLOT IPMS
#######################################

par(mfcol=c(2, 6))

lam_const <- NULL

for (sp in spp) {
	#smax<-params[params$spp==sp, "smax"]
	rec.size <- params$rec.size[params$spp==sp]
	b <- seq(rec.size, max.size, length=n)
	h <- b[2] - b[1]
	b <- c(min(b)-h, b)
	y <- 0.5 * (b[1:n]+b[2:(n+1)])
	I <- y >= rec.size
  	sub<-gdat[gdat$spp==sp,]
    mod <- bigmatrix()
	image(y, y, t(mod$P)^0.3)    
	points(sub$area, sub$area_next, cex=0.25)
	title(sp, line=-1)
	abline(0, 1, lty=2)
  lam_const<-c(lam_const,bigmatrix()$lam) }

params$lam_const<-lam_const
params$lam_const

#######################################
# MULTIPLE REC
#######################################

store<-data.frame()
recx<-c(1, 0.1, 0.01,0.001, 1*10^-4, 1*10^-5, 1*10^-6)
for (sp in spp) {
	for (rec in recx){
	  #smax<-params[params$spp==sp, "smax"]
	  rec.size <- params$rec.size[params$spp==sp]
	  b <- seq(rec.size, max.size, length=n)
	  h <- b[2] - b[1]
	  b <- c(min(b)-h, b)
	  y <- 0.5 * (b[1:n]+b[2:(n+1)])
	  I <- y >= rec.size
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  store<-rbind(store, data.frame(sp=sp, rec=rec, lam=bigmatrix()$lam))
  } }
recdat<- dcast(store, sp~rec, value.var='lam')
params$lam_0<-recdat[,"0.1"]
params$lam_1<-recdat[,"0.01"]
params$lam_2<-recdat[,"0.001"]
params$lam_3<-recdat[,"1e-04"]
params$lam_4<-recdat[,"1e-05"]


#######################################
# EXPLORE LAMBDA
#######################################
lam<-expression(lambda)
lamlong<-melt(params, measure.vars=c("lam_0","lam_1","lam_2","lam_3", "lam_4"), value.name="lambda")
ggplot(data=lamlong)+ 
  geom_hline(yintercept=1, linetype="dotted")+
  geom_line(aes(x=abundance, y=lambda, col=spp, group=variable))+
  geom_point(aes(x=abundance, y=lambda, fill=spp, size=variable), shape=21)+
  guides(fill="none", col="none")+labs(y=lam)+ scale_fill_manual(values=as.character(params$cols))+
  scale_colour_manual(values=as.character(params$cols))+
  facet_wrap(~morphology, scales="free_x")+
  theme(axis.text.x=element_text(angle=90, size=5))





#######################################
# COMPARE SPECIES PAIRS
#######################################
params2<-rbind(params, params[6,]) #duplicate millepora(should it be nasuta?) 
params2$morph<-as.character(params2$morph)
params2$morph[c(7, 12)]<-c("corymbose_2","corymbose_2")
comp<-dcast(params2, morph~abundance_pair, value.var="spp")
comp$sp<-apply(comp[,c(2,3)], 1, paste, collapse = ":")
comp$lam_0<-transform(dcast(params2, morph~abundance_pair, value.var="lam_0"), x=Common/Rare)$x
comp$lam_1<-transform(dcast(params2, morph~abundance_pair, value.var="lam_1"), x=Common/Rare)$x
comp$lam_2<-transform(dcast(params2, morph~abundance_pair, value.var="lam_2"), x=Common/Rare)$x
comp$lam_3<-transform(dcast(params2, morph~abundance_pair, value.var="lam_3"), x=Common/Rare)$x
comp$lam_4<-transform(dcast(params2, morph~abundance_pair, value.var="lam_4"), x=Common/Rare)$x
comp<-melt(comp, id.vars=c("morph","sp","Common","Rare"))
#------------------------------- 
comp$morph<-factor(comp$morph, levels=c("tabular","staghorn", "corymbose", "corymbose_2", "digitate","massive"))
ggplot(comp)+
geom_bar(aes(x=morph, y=value),stat="identity", col="black", fill="grey")+
geom_text(aes(x=morph, y=value+0.01, label=sp), size=2, angle=90)+
coord_flip()+scale_y_log10()+
theme(axis.line.y=element_blank(), strip.background=element_blank(), axis.title.y=element_blank())+
ylab("lambda ratio")+
facet_wrap(~variable, nrow=2)


#######################################
# ELASTICITY
#######################################
 # eigen-things can be combined to obtain the sensitivity and elasticity matrices.
 
rec<-1*10^-6

analyses<-function(){
	#sp<-"GR"
	K<-bigmatrix()$K
	lam <- Re(eigen(K)$values[1])
  w.eigen<-Re(eigen(K)$vectors[,1]) # w=right eigen=stable size dist. 
  v.eigen<-Re(eigen(t(K))$vectors[,1]) # v=left, reproductive value=contribution of each size class to pop size. 
stable.size <- w.eigen/sum(w.eigen)
repro.val <- v.eigen/v.eigen[1] 
v.dot.w<-sum(stable.size*repro.val)*h 
sens<-outer(repro.val,stable.size)/v.dot.w  
elas<-matrix(as.vector(sens)*(as.vector(K)/h)/lam, nrow=n)
P.elas<-(bigmatrix()$P/h)*sens/lam
R.elas<-(bigmatrix()$R/h)*sens/lam
return(list(st=stable.size, rv=repro.val, sa=sens, eK=elas, eP=P.elas, eR=R.elas))}

sum(analyses()$eK)*h^2
  
props<-NULL 
for (sp in spp) {
props<-rbind(props, data.frame(eP=sum(analyses()$eP)*h^2, eR=sum(analyses()$eR)*h^2, sp=sp))
	}
props
params$eP<-props$eP
params$eR<-props$eR
rowSums(props[,c(1,2)])


#######################################
# R.ELASTICITY VS RECRUITMENT SENSITIVITY
#######################################

ggplot(data=params, aes(x=eR, y=lam_4-lam_1))+
  geom_text(aes(label=spp))+
  geom_smooth(method="lm", formula=y~x)+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="R.elasticity", y="change in lambda")


#######################################
# R.ELASTICITY vs MORPH/ABUNDANCE
#######################################

params$R<-ifelse(params$morphology=="massive", "S", ifelse(params$morphology=="digitate", "R", "R"))
mod<-lm(eR~log(abundance_05), subset(params, R!="S"))
summary(mod)

plot_grid(
ggplot(data=params, aes(y=eR, x=abundance_05))+
geom_smooth(method="lm", formula=y~poly(x,2))+geom_text(aes(label=spp))+
scale_x_log10()+
scale_y_log10()+
labs(y="R.elasticity")#+facet_wrap(~R)
,
ggplot(data=params, aes(x=reorder(morphology, -eR), y=eR))+geom_boxplot()+
geom_point()+labs(y="R.elasticity")+
scale_y_log10()+
geom_text(aes(label=spp), nudge_x=0.15, col="red")
)


#Pvals <- big$P / h
#P.elas <- Pvals * K.sens / lambda



  
  #######################################
  #-------------------------------  min.size/max.size-lamda relationship
  #######################################
  
  store<-data.frame()
  minx<-c(-6,-5,-4,-3,-2,-1)
  for (sp in spp) {
    for(min.size2 in minx){
      min.size2<-min.size2
      max.size <- max.size #2
      rec.size <- rec.size
      n <- 100
      b <- min.size2 + c(0:n) * (max.size - min.size2)/n
      y <- 0.5 * (b[1:n]+b[2:(n+1)])
      h <- y[2] - y[1]
      store<-rbind(store, data.frame(sp=sp, min.size=min.size2, lam=bigmatrix()$lam)) }}
  ggplot(data=store, aes(x=min.size, y=lam, colour=sp))+
    geom_line()+geom_hline(yintercept=1)+
    geom_vline(xintercept=min.size, linetype="dotted")+
    ggtitle("relative differences")+
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
  
  
  
  
  










