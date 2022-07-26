
gdat_orig <- gdat
sdat_orig <-sdat
fec_orig <-fec

##################################
##################################
##################################

# STEP 1: sample colonies with replacement from the data until the original sample size (n) is reached, then repeat x times

NREP <- 1000

g.cols <- NULL
s.cols <- NULL
f.cols <- NULL
#sp <- "Ahy"
#i <- 1
fec$time <- substr(fec$id, 1,2)
yrs <- unique(fec$time)
for(sp in spp){
for(i in 1:NREP){
g.pool <- unique(gdat[gdat$spp==sp,"colony_id"])
s.pool <- unique(sdat[sdat$spp==sp,"colony_id"])
n1 <- length(c(unique(g.pool),unique(s.pool)))
g.samp <- sample(g.pool, n1, replace=T)
g.cols <- rbind(g.cols, data.frame(cols=g.samp, samp=i, spp=sp))
s.samp <- sample(s.pool, n1, replace=T)
s.cols <- rbind(s.cols, data.frame(cols=s.samp, samp=i, spp=sp))
f.pool <- unique(fec[fec$spp==sp,"id"]) 
n2 <- length(unique(f.pool))
f.samp <- sample(f.pool, n2, replace=T)
f.cols <- rbind(f.cols, data.frame(cols=f.samp, samp=i, spp=sp))
}}

# STEP 2: Re-fit demographic models for each sample and re-calculate lambda

p.range <- NULL
lam.range <- NULL
#i <- 2
for(i in 1:NREP){
	#i <- 18
gdat <- gdat_orig
sdat <-sdat_orig 
fec <-fec_orig 
gdat <- gdat[gdat$colony_id %in% g.cols$cols[g.cols$samp==i], ]
sdat <- sdat[sdat$colony_id %in% s.cols$cols[s.cols$samp==i], ]
fec <- fec[fec$id %in% f.cols$cols[f.cols$samp==i], ]

source("R/params.R")
params[params$spp=="Gre", c("g.int", "g.slp","g.var")] <- params[params$spp=="Gpe", c("g.int", "g.slp","g.var")]

for(sp in spp){
	tryCatch({
		#sp <-"Gpe"
	rec.size <- params$rec.size[params$spp==sp]
	rec <- params$rec[params$spp==sp]
	params[params$spp==sp,]
	h <- mesh()$h
	y <- mesh()$y
	lam.boot <- bigmatrix()$lam
	P <- bigmatrix()$P
	R <- bigmatrix()$R
	N <- solve(diag(n)-P)
	R0 <- abs(eigen(R %*% N)$values[1])
lam.range<-rbind(lam.range, data.frame(spp=sp, samp=i, lam=lam.boot, r=R0))	
},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
}

head(lam.range)

lam.range$morphology <- params$morphology[match(lam.range$spp, params$spp)]

ggplot(lam.range)+
geom_density(aes(x=log(lam), col=spp, fill=spp), alpha=0.5)+
geom_vline(xintercept=0, linetype="dotted")+
facet_wrap(~morphology, scale="free_y", ncol=1)+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
theme_classic()

gdat <- gdat_orig
sdat <-sdat_orig 
fec <-fec_orig 

#write.csv(lam.range, "data/lam.range.csv")
