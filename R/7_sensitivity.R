
#######################################
# RECRUITMENT SENSITIVITY II
#######################################	


# lambda at different rec/recsize combinations
rec.x <- 10^(seq(-6,-2,0.5))
recsize.x <- seq(min(params$rec.size)-0.1, max(params$rec.size)+0.1,0.1)
3
store3 <- list()
for (sp in spp) {
	temp<-matrix(NA,length(rec.x), length(recsize.x))
	for (i in 1:length(rec.x)){
		rec <- rec.x[i]
		for(j in 1:length(recsize.x)){
			rec.size <- recsize.x[j]
		h <- mesh()$h
		y <- mesh()$y	
   mod <- bigmatrix()
   temp[i,j] <- mod$lam
  }
 } 
 store3[[sp]] <- temp  
  }

minlam <- min(log(do.call("rbind", store3)))
maxlam <- max(log(do.call("rbind", store3)))

par(mfrow=c(3,4))
sp.conts <- NULL
for(sp in spp){
	#sp<-"Ahy"
  image(log10(rec.x), recsize.x, log(store3[[sp]]), main=sp, zlim=c(minlam, maxlam))
  abline(h=params[params$spp==sp, "rec.size"], lty=2)
  contour(log10(rec.x), recsize.x, log(store3[[sp]]), add=TRUE, levels=0)
  }
   contour(log10(rec.x), recsize.x, log(store3[[spp[1]]]), levels=0, col=cols[1], lwd=3)
 for(i in 1:length(spp)){
 	sp <- spp[i]
  contour(log10(rec.x), recsize.x, log(store3[[sp]]), add=TRUE, levels=0, col=cols[i], lwd=3)
  conts <- melt(store3[[spp[i]]])
  sp.conts <- rbind(sp.conts, cbind(conts,rec = rep(rec.x, length(recsize.x)), recsize = rep(recsize.x, each=length(rec.x)), spp=sp))
  }
head(sp.conts)

##### PLOT

params2$rec <- params$rec[match(params2$spp, params$spp)]
params2$rec.size <- params$rec.size[match(params2$spp, params$spp)]

comp$rec <- params2$rec[match(comp$morph, params2$morph)] 
comp$rec.size <- params2$rec.size[match(comp$morph, params2$morph)] 
comp$rec.size.cm <- sqrt(((10^comp$rec.size)*10000)/pi)*2	

# difference in lambda across rec/recsize combinations
pairs3 <- NULL
for(m in comp$morph){
	#m <- "tabular"
	spC <- comp$Common[comp$morph==m]
	spR <- comp$Rare[comp$morph==m]
	logdiff <- log(store3[[spC]]) - log(store3[[spR]])
	temp <- melt(logdiff)
	temp$rec <- rep(rec.x, length(recsize.x))
	temp$recsize <- rep(recsize.x, each=length(rec.x))
    pairs3 <- rbind(pairs3, cbind(temp, morph=m))
}
head(pairs3)

pairs3$lam2 <- ifelse(pairs3$value < 0, NA, pairs3$value)

ggplot()+
geom_raster(data=pairs3, aes(x=rec, y=recsize, fill=lam2))+
geom_point(data=params2, aes(rec, rec.size ), col="white", shape=3)+
	scale_y_continuous(expand=c(0,0))+
	scale_x_log10(expand=c(0,0))+
	scale_fill_distiller(palette="Spectral")+
	facet_wrap(~morph, scales="free")+
	labs(x="rec rate", y="rec size")+
	theme(strip.background=element_blank())
