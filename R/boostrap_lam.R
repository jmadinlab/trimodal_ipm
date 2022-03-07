

gdat <- gdat_orig
sdat <-sdat_orig 
fec <-fec_orig 

#samp_sizes <- aggregate(colony_id~spp, gdat, function(x){length(unique(x))})
#samp_sizes$surv <- aggregate(colony_id~spp, sdat, function(x){length(unique(x))})$colony_id
#samp_sizes$fec <- aggregate(id~spp, fec, function(x){length(unique(x))})$id
#samp_sizes



REMOVE <- 4
F.REMOVE <- 15
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


# Adi 1, Gre Loads

p.range <- NULL
lam.range <- NULL
#i <- 2
for(i in 1:NREP){
	#i <- 18
gdat <- gdat_orig
sdat <-sdat_orig 
fec <-fec_orig 

length(unique(gdat$colony_id))
gdat <- gdat[gdat$colony_id %in% g.cols$cols[g.cols$samp==i], ]
length(unique(gdat$colony_id))
sdat <- sdat[sdat$colony_id %in% s.cols$cols[s.cols$samp==i], ]
fec <- fec[fec$id %in% f.cols$cols[f.cols$samp==i], ]

# remove tiny colonies (<3cm diameter)
diam <- 3
gdat <- subset(gdat, area > log10((pi*(diam/2)^2)/10000))
sdat <- subset(sdat, area > log10((pi*(diam/2)^2)/10000))
fec <- subset(fec, area > log10((pi*(diam/2)^2)/10000))

sdat<-sdat[!sdat$colony_id==287,]
gdat<-gdat[!gdat$colony_id==287,]

source("R/params.R")

params[params$spp=="Gre", c("g.int", "g.slp","g.var")] <- params[params$spp=="Gpe", c("g.int", "g.slp","g.var")]

p.range <- rbind(p.range, cbind(samp=i, params[,c("g.slp","g.int","g.var", "f.slp","f.int", "s.slp", "s.int","s.slp.2")]))
for(sp in spp){
	#sp <- "Ahu"

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
#lam.range

lam.range$morphology <- params$morphology[match(lam.range$spp, params$spp)]

ggplot(lam.range)+
geom_density(aes(x=lam, col=spp, fill=spp), alpha=0.5)+
geom_vline(xintercept=1, linetype="dotted")+
facet_wrap(~morphology, scale="free_y", ncol=1)+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
xlim(c(0.5, 1.8))+
#scale_x_log10(limits=c(0.5, 2))+
theme_classic()

table(lam.range$spp)

ggplot(lam.range, aes(y=reorder(morphology, -lam), x=lam, fill=spp))+geom_boxplot()+geom_vline(xintercept=1)+scale_fill_manual(values=cols)+guides(fill="none")

gdat <- gdat_orig
sdat <-sdat_orig 
fec <-fec_orig 

#write.csv(lam.range, "data/lam.range.csv")

# aggregate(eggs~spp, fec, function(x){ })
