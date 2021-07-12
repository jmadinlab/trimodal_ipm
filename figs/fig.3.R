
colsC <- params$cols[params$abundance_pair=="Common"]
names(colsC) <- params$morphology[params$abundance_pair=="Common"]


demovals$eR[demovals$spp=="GR"]<-5e-3
demovals$lam_4 <- store[store$rec==0.01,"lam"]
demovals$lam_1 <- store[store$rec==0.0001,"lam"]
demovals$ldiff <- log(demovals$lam_4)-log(demovals$lam_1)
demovals$morphology <- params$morphology
demovals$eRarr <- demovals$eR
demovals$Larr <- demovals$ldiff
demovals$Larr <- ifelse(demovals$morphology=="staghorn", demovals$ldiff+0.1, demovals$Larr) 
demovals$Larr <- ifelse(demovals$morphology=="tabular", demovals$ldiff+0.1, demovals$Larr) 
demovals$Larr <- ifelse(demovals$morphology=="digitate", demovals$ldiff-0.1, demovals$Larr) 
demovals$Larr <- ifelse(demovals$morphology=="corymbose", demovals$ldiff-0.1, demovals$Larr) 
demovals$eRarr <- ifelse(demovals$morphology=="massive", demovals$eR-0.015, demovals$eRarr) 

demovals$spp <- factor(demovals$spp, levels=c(order))
demovals$order <- c(2,3,1,5,7,8,9,6,4,11,10)

eRbar <- ggplot(demovals[order(demovals$order),], aes(x=reorder(morphology, eR), y=eR, fill=spp))+
geom_bar(stat="identity", position = position_dodge2(preserve = "single"), width=0.8, col="black", size=0.05)+
theme_classic()+
#scale_y_sqrt(expand=c(0,0))+
scale_y_continuous(expand=c(0,0))+
scale_colour_manual(values=colsC)+
coord_flip()+
scale_fill_manual(values=cols)+guides(fill="none", col="none")+
theme(axis.title.x=element_blank(), axis.line.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=6), axis.title.y=element_blank(), axis.ticks.x=element_blank(),plot.background=element_blank(),panel.background=element_blank())
eRbar

eRarrow <- ggplot(demovals[order(demovals$order),], aes(x=reorder(morphology, eR), y=eR, fill=spp))+
geom_path(aes(group=morphology, col=morphology),arrow=arrow(length=unit(0.75,"mm"), ends="first"))+
coord_flip()+
scale_colour_manual(values=colsC)+guides(fill="none", col="none")+theme_void()
eRarrow

Relasplot<- ggplot(data=demovals[order(demovals$order),], aes(x=eR, y=ldiff))+
geom_smooth(method="lm", formula=y~x, size=0.1, col="black")+
 geom_point(aes(fill=spp), shape=21, stroke=0.1, size=2)+
# geom_path(aes(x=eRarr, y=Larr, col=morphology), arrow=arrow(length=unit(1,"mm")))+
  #geom_text(aes(label=spp))+
  scale_colour_manual(values=colsC)+
  scale_fill_manual(values=cols)+
  #scale_x_log10()+
  #scale_y_log10()+
  #annotation_custom(ggplotGrob(eRbar), 0.1, 0.35, -0.05, 0.5)+
  #annotation_custom(tab, 0.2, 0.35, 0, 0.5)+
  guides(fill="none", col="none")+
  labs(x="Elasticity to fecundity", y=expression(Delta~log(lambda)~with~recruitment))+
  theme_classic()+
  theme(axis.text=element_text(size=7), 
  axis.title.y=element_text(size=7, vjust=-10),
  axis.title.x=element_text(size=7), plot.background=element_blank(), panel.background=element_blank())

elas2 <- plot_grid(eRbar, NULL,eRarrow, NULL,Relasplot, ncol=1, rel_heights=c(0.3,0,0.1,-0.15,1), align="v", labels=c("A","","B",""), label_size=9)
elas2

params$pc1 <- pca$x[,1]
params$pc2 <- pca$x[,2]


library("plot3D")
scatter3Dx <- function(x, y, z,...)
  { panelfirst <- function(pmat) {
   	# flat against the bottom
      XY <- trans3D(x, y, z = rep(Rz[1], length(z)), pmat = pmat)
      scatter2D(XY$x, XY$y,  pch = 19, col="grey", 
              cex = 0.25, add = TRUE, colkey = FALSE)}
  scatter3D(x, y, z, ..., panel.first=panelfirst) 
}

psize1 <- 1.9
psize2 <- 1.7
labsiz <- 2

png("data/hid2b.png", width = 5, height = 5.5, units = 'in', res = 300)
params$x<-params$pc1
params$y<-params$pc2
params$z<-params$R0
Rx<-c(min(params$x)-0.5,max(params$x)+0.5)
Ry<-c(min(params$y)-0.5,max(params$y)+0.5)
Rz<-c(min(params$z)-0.5,max(params$z)+0.5)
Lx<-"PC1"
Ly<-"PC2"
Lz<-"Ro"
p_large <- subset(params, z>min(params$z)+1.1)
scatter3Dx(p_large$x,p_large$y,p_large$z-0.5,colvar = NULL, col = "black", pch = 19, cex = 0.01, bty = "u", col.panel ="white",col.grid = "grey80", 
 theta=135, phi=40, xlab=Lx, zlab=Lz, ylab=Ly, ticktype = "detailed", xlim=Rx, zlim=Rz, ylim=Ry, expand=1.5, nticks=4, cex.axis=0.001,  cex.lab=labsiz, xaxt="n", type="h") # theta=right
 scatter3Dx(params$x,params$y,params$z,col="black", add=TRUE,colvar = NULL, cex=psize1,pch = 19)
scatter3Dx(params$x,params$y,params$z,col=cols, add=TRUE,colvar = NULL, cex=psize2,pch = 19)
dev.off()

png("data/hid1b.png", width = 5, height = 5.5, units = 'in', res = 300)
params$z<-log10(params$GT)
Rx<-c(min(params$x)-0.5,max(params$x)+0.5)
Ry<-c(min(params$y)-0.5,max(params$y)+0.5)
Rz<-c(min(params$z)-0.2,max(params$z))
Lx<-"PC1"
Ly<-"PC2"
Lz<-"T"
p_large <- subset(params, z>min(params$z)+0)
scatter3Dx(p_large$x,p_large$y,p_large$z-0.1,colvar = NULL, col = "black", pch = 19, cex = 0, bty = "u", col.panel ="white",col.grid = "grey80", 
 theta=135, phi=40, xlab=Lx, zlab=Lz, ylab=Ly, ticktype = "detailed", xlim=Rx, zlim=Rz, ylim=Ry, expand=1.5, nticks=4, cex.axis=0.001,  cex.lab=labsiz, xaxt="n", type="h") # theta=right
 scatter3Dx(params$x,params$y,params$z,col="black", add=TRUE,colvar = NULL, cex=psize1,pch = 19)
scatter3Dx(params$x,params$y,params$z,col=cols, add=TRUE,colvar = NULL, cex=psize2,pch = 19)
dev.off()

hid1 <-readPNG("data/hid1.png")
hid1<-rasterGrob(hid1, interpolate=TRUE)
hid2 <-readPNG("data/hid2.png")
hid2<-rasterGrob(hid2, interpolate=TRUE)

imageread<-function(image){
qplot(1:10, 1:10, geom="blank") +
  annotation_custom(image, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+theme_void()+
   theme(axis.line.x=element_blank(),
   axis.line.y=element_blank(),  
   axis.text.x=element_blank(), 
   axis.text.y=element_blank(), 
   axis.title.x=element_blank(), 
   axis.title.y=element_blank(), 
   plot.title=element_text(hjust=0.5, vjust=0,face="bold", size=5),
   #plot.title=element_blank(),
   panel.grid.major = element_blank(), 
   axis.ticks=element_blank(), 
   panel.grid.minor = element_blank(),
   panel.border=element_blank())}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

threeDs <- plot_grid(
imageread(hid1)+ggtitle("Generation time"),
NULL,
imageread(hid2)+ggtitle("Net reproductive rate"),#+ggtitle(expression(bold(R[0]))),
ncol=1, 
rel_heights=c(1,0.05,1))
threeDs 

# FIGURE 3


demovals$Rarr <- demovals$R0
demovals$GTarr <- demovals$GT
demovals$Rarr <- ifelse(demovals$morphology=="digitate", demovals$Rarr+0.6, demovals$Rarr)
demovals$Rarr <- ifelse(demovals$morphology=="staghorn", demovals$Rarr+0.6, demovals$Rarr)
demovals$GTarr <- ifelse(demovals$spp=="GR", demovals$GT*0.8, demovals$GTarr)
demovals$GTarr <- ifelse(demovals$spp=="GP", demovals$GT*1.3, demovals$GTarr)
demovals$Rarr <- ifelse(demovals$spp=="AC", demovals$R0+0.5, demovals$Rarr)
demovals$Rarr <- ifelse(demovals$spp=="AH", demovals$R0-0.5, demovals$Rarr)
demovals$Rarr <- ifelse(demovals$spp=="AL", demovals$R0-0.4, demovals$Rarr)
demovals$Rarr <- ifelse(demovals$spp=="AM", demovals$R0+0.4, demovals$Rarr)
demovals$Rarr <- ifelse(demovals$spp=="AN", demovals$R0+0.4, demovals$Rarr)




t_off <- ggplot(demovals[order(demovals$order),], aes(GT, R0))+
#geom_text(aes(label=spp))+
geom_point(aes(fill=spp), shape=21, size=2.5, stroke=0.2)+
geom_path(aes(x=GTarr, y=Rarr, 
col=morphology), arrow=arrow(length=unit(1,"mm"),ends="first"))+
scale_x_log10()+
scale_colour_manual(values=colsC)+
scale_fill_manual(values=cols)+
guides(col="none", fill="none")+
theme_classic()+
labs(x="Generation time (T)", y=expression(Net~Reproductive~Rate~(R[0])))+
theme(axis.title=element_text(size=7), axis.text=element_text(size=7))
t_off



RvsGT<-plot_grid(t_off,threeDs, nrow=1, rel_widths=c(1, 0.4), align="h", labels=c("C","D"), label_size=9)
RvsGT



fig.3 <- plot_grid(
  elas2,
RvsGT,
rel_widths=c(0.75,1), labels=c("",""),  label_size=9)
fig.3

#ggsave("figs/fig.3.png", fig.3, width=13, height=5.5, units="cm", dpi = 300)

