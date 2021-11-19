 
 
#######################################
# MULTIPLE REC
#######################################

rec.x <- 10^(seq(-6,-2,0.5))

store<-data.frame()
for (sp in spp) {
	for (rec in rec.x){
	  rec.size <- params$rec.size[params$spp==sp]
	  h <- mesh()$h
	  y <- mesh()$y
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  store<-rbind(store, data.frame(sp=sp, rec=rec, lam=bigmatrix()$lam))
   } 
  }
 
 pairs <- NULL 
	for (rec in rec.x){
		df <- store[store$rec==rec,]
		temp <- comp
		temp$c.lam <- df$lam[match(temp$Common, df$sp)]
		temp$r.lam <- df$lam[match(temp$Rare, df$sp)]
		#temp$diff <- temp$c.lam - temp$r.lam
		#temp$ratio <- temp$c.lam / temp$r.lam
		temp$logdiff <- log(temp$c.lam) - log(temp$r.lam)
		temp$rec <- rec
		pairs <- rbind(pairs, temp)
		}


#######################################
# MULTIPLE REC. SIZE
#######################################

recsize.x <- seq(min(params$rsize.gr)-0.05, max(params$rsize.gr)+0.05,0.1)

store2 <- NULL
for (sp in spp) {
	for (rec.size in recsize.x){
	  rec <- params$rec[params$spp==sp]
	  h <- mesh()$h
	  y <- mesh()$y
  sub<-dat[dat$spp==sp,]
   mod <- bigmatrix()
  store2<-rbind(store2, data.frame(sp=sp, rec.size=rec.size, lam=bigmatrix()$lam))
   } 
  }
 
pairs2 <- NULL 
	for (rec.size in recsize.x){
		df <- store2[store2$rec.size==rec.size,]
		temp <- comp
		temp$c.lam <- df$lam[match(temp$Common, df$sp)]
		temp$r.lam <- df$lam[match(temp$Rare, df$sp)]
		#temp$diff <- temp$c.lam - temp$r.lam
		#temp$ratio <- temp$c.lam / temp$r.lam
		temp$logdiff <- log(temp$c.lam) - log(temp$r.lam)
		temp$rec.size <- rec.size
		pairs2 <- rbind(pairs2, temp)
		}
			
 
rec_sens <- plot_grid(
ggplot(store, aes(x=rec, y=log(lam), col=sp))+
geom_hline(yintercept=0)+geom_line(size=1)+
 scale_x_log10(breaks=c(10^-6, 10^-4,10^-2), labels=c(
		expression(10^-6),
		expression(10^-4),
		expression(10^-2)))+scale_colour_manual(values=cols)+
 theme_classic()+
 theme(legend.title=element_blank(), legend.key.height=unit(1,"mm"))+
 labs(x="probability of recruitment", y=expression(log(lambda))), 
	ggplot(pairs,aes(rec, logdiff, col=morph))+
	geom_hline(yintercept=0)+
		geom_line(size=1)+
		scale_x_log10(breaks=c(10^-6, 10^-4,10^-2), labels=c(
		expression(10^-6),
		expression(10^-4),
		expression(10^-2)))+
		scale_colour_manual(values=colsC)+
		labs(x="probability of recruitment", y=expression(Delta~log*(lambda)~(common~-~rare)))+
		theme_classic()+theme(legend.title=element_blank(), legend.key.height=unit(1,"mm")),
		ncol=1, align="v", axis="lr")
	rec_sens	
		
					
rsize_sens<-plot_grid(ggplot(store2, aes(sqrt(((10^rec.size)*10000)/pi)*2, log(lam), col=sp))+
geom_hline(yintercept=0)+
 geom_line(size=1)+
 scale_colour_manual(values=cols)+
 theme_classic()+
 theme(legend.title=element_blank(), legend.key.height=unit(1,"mm"))+
 labs(x="recruit diameter (cm)", y=expression(log(lambda))),
 ggplot(pairs2,aes(sqrt(((10^rec.size)*10000)/pi)*2, logdiff, col=morph))+geom_line(size=1)+
  theme_classic()+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=colsC)+
 theme(legend.title=element_blank(), legend.key.height=unit(1,"mm"))+
 labs(x="recruit diameter (cm)", y=expression(Delta~log*(lambda)~(common~-~rare))),
 ncol=1, align="v", axis="lr")
rsize_sens



plot_grid(rec_sens, rsize_sens, nrow=1)




#######################################
# RECS NEEDED TO GROW
#######################################

tiles.a <- subset(tiles, Family=="Acroporidae")
rsize.a <- subset(rsize, Family=="Acroporidae")
tiles.m <- subset(tiles, Family=="Merulinidae")
rsize.m <- subset(rsize, Family=="Merulinidae")
log10(pi*((5/2)/100)^2)

params$rec0 <- c(8*10^-4, 1*10^-1, 1*10^-1, 3.3*10^-4, 1.2*10^-3, 1*10^-1, 1*10^-1, 1.7*10^-3, 1.3*10^-3, 6.3*10^-4, 5.9*10^-3)

sp.conts2 <- sp.conts
for(sp in spp){
sp.conts2 <- sp.conts2[!(sp.conts2$spp==sp & sp.conts2$recsize > params$rsize.gr[params$spp==sp]), ]}


ggplot()+  
#geom_rect(data=NULL, aes(xmin=min(tiles.a$N_m2_year/tiles.a$eggs), xmax=max(tiles.a$N_m2_year/tiles.a$eggs), ymin=-Inf, ymax=Inf), fill="grey90", col="grey90")+
geom_rect(data=NULL, aes(xmin=quantile(tiles.a$N_m2_year/tiles.a$eggs, 0.25), xmax=quantile(tiles.a$N_m2_year/tiles.a$eggs, 0.75), ymin=params$rsize.gr2[params$spp=="Ahu"], ymax=Inf), fill="grey90", col='grey90')+
geom_rect(data=NULL, aes(xmin=quantile(tiles.m$N_m2_year/tiles.m$eggs, 0.25), xmax=quantile(tiles.m$N_m2_year/tiles.m$eggs, 0.75), ymin=-Inf, ymax=params$rsize.gr2[params$spp=="Gre"]), fill="grey90", col=NA)+
stat_contour(data=sp.conts, aes(x=rec, y=recsize, z = log(value), col=spp), breaks=c(0), size=0.35, linetype="dashed")+
stat_contour(data=sp.conts2, aes(x=rec, y=recsize, z = log(value), col=spp), breaks=c(0), size=1)+
scale_x_log10( breaks=c( 10^-5, 10^-4, 10^-3, 10^-2))+
#geom_segment(data=params, inherit.aes=F, aes(x=10^-5, xend=rec0, y=rsize.gr,yend=rsize.gr, col=spp), linetype="dotted")+
labs(x="Probability of larval survival", y="Recruit size (cm)")+
geom_text(data=NULL,  aes(9*10^-5, -3, label="Population \nshrinking"), size=2.5)+
geom_text(data=NULL,  aes(5*10^-3, -2.2, label="Population \ngrowing"), size=2.5)+
geom_point(data=params,inherit.aes=F, aes(y=rsize.gr, x=rec0, col=spp), size=2)+
scale_fill_manual(values=cols)+
scale_colour_manual(values=cols)+
scale_y_continuous(breaks=log10(pi*((c(2,4,6,8)/2)/100)^2), labels=c(2,4,6,8))+
coord_cartesian(xlim=c(1*10^-5, 10^-2))+
#ggtitle("loglam = 0")+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
axis.text.x=element_text(angle=30, vjust=1, hjust=1))+
guides(colour="none", fill="none")







