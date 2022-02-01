

gdat_orig <- gdat
sdat_orig <- sdat
fec_orig <- fec

#samp_sizes <- aggregate(colony_id~spp, gdat, function(x){length(unique(x))})
#samp_sizes$surv <- aggregate(colony_id~spp, sdat, function(x){length(unique(x))})$colony_id
#samp_sizes

REMOVE <- 4
F.REMOVE <- 15
NREP <- 100

g.cols <- NULL
s.cols <- NULL
f.cols <- NULL
#sp <- "Ahy"
#i <- 1
for(sp in spp){
for(i in 1:NREP){
g.pool <- unique(gdat[gdat$spp==sp,"colony_id"])
g.samp <- sample(g.pool, length(g.pool)-REMOVE)
g.cols <- rbind(g.cols, data.frame(cols=g.samp, samp=i, spp=sp))
s.pool <- unique(sdat[sdat$spp==sp,"colony_id"])
s.samp <- sample(s.pool, length(s.pool)-REMOVE)
s.cols <- rbind(s.cols, data.frame(cols=s.samp, samp=i, spp=sp))
f.pool <- unique(fec[fec$spp==sp,"id"]) 
f.samp <- sample(f.pool, length(f.pool)-F.REMOVE)
f.cols <- rbind(f.cols, data.frame(cols=f.samp, samp=i, spp=sp))
}}

p.range <- NULL
lam.range <- NULL
#i <- 2
for(i in 1:NREP){
	
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

source("R/params_only.R")

params[params$spp=="Gre", c("g.int", "g.slp","g.var")] <- params[params$spp=="Gpe", c("g.int", "g.slp","g.var")]

p.range <- rbind(p.range, cbind(samp=i, params[,c("g.slp","g.int","g.var", "f.slp","f.int", "s.slp", "s.int","s.slp.2")]))
for(sp in spp){
	rec.size <- params$rec.size[params$spp==sp]
	rec <- params$rec.mean[params$spp==sp]
	h <- mesh()$h
	y <- mesh()$y
lam.range<-rbind(lam.range, data.frame(spp=sp, samp=i, lam=bigmatrix()$lam))	
}}

head(lam.range)

lam.range$morphology <- params$morphology[match(lam.range$spp, params$spp)]

ggplot(lam.range)+
geom_density(aes(x=lam, col=spp, fill=spp), alpha=0.5)+
scale_y_sqrt(expand=c(0,0))+
#scale_x_log10()+
#facet_wrap(~morphology, scale="free")+
scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+
guides(fill="none", col="none")+
theme_classic()


ggplot(lam.range, aes(reorder(spp, -lam), lam))+geom_boxplot()+geom_hline(yintercept=1)

gdat <- gdat_orig
sdat <-sdat_orig 
fec <-fec_orig 

#write.csv(lam.range, "data/lam.range.csv")

# aggregate(eggs~spp, fec, function(x){ })
