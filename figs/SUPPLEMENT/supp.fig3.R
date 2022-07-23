
# first load "params.R"

# -------- plot
#pdf("figs/p1_maturity.pdf" )

p1<-ggplot()+
  geom_jitter(data=fec, 
    aes(x= area, reproductive), 
      height=0.02, shape=21, col="grey")+
  geom_line(data=m.pred, 
    aes(area, pred, col=spp))+
  scale_colour_manual(values=cols, labels=labs)+
  theme_classic()
p1

#show(p1)
#dev.off()

# -------- plot
#pdf("figs/p2_fecundity.pdf")

p2<-ggplot()+ 
  geom_point(data=fec[fec$reproductive==1, ], 
    aes(area, fecundity), 
      shape=21, col="grey")+
  geom_line(data=f.pred, aes(area, pred, col=spp))+
     #geom_line(data=f.pred2, aes(area, pred, col=spp))+
  scale_y_log10()+
  #facet_wrap(~spp)+
  scale_colour_manual(values=cols)+
  theme_classic()
p2
#show(p2)
#dev.off()


#pdf("figs/p3_growth.pdf")

p3<-ggplot()+ 
  geom_abline(slope=1, linetype="dotted")+
  geom_point(data=gdat, 
    aes(area, area_next), 
      shape=21, col="grey")+
  geom_line(data=g.pred, 
    aes(area, pred, col=spp))+
  scale_colour_manual(values=cols)+
  #facet_wrap(~spp)+
 theme_classic()

#show(p3)
#dev.off() 
  

#pdf("figs/p4_survival.pdf")

p4<-ggplot()+ 
  geom_jitter(data=sdat, 
    aes(area, surv), 
      shape=21, col="grey",height=0.02)+
  geom_line(data=s.pred, 
    aes(area, pred, col=spp))+
      #facet_wrap(~spp)+
  scale_colour_manual(values=cols)+
 theme_classic()
p4
#show(p4)
#dev.off() 


#pdf("figs/p5_maxgrowth.pdf", )

p5<-ggplot()+
 geom_point(data=gdat, aes(x=area, y=g_radius), col="grey", shape=21)+
     geom_line(data=r.pred, aes(x=area, y=pred, col=spp))+
   scale_colour_manual(values=c(cols))+
   #+facet_wrap(~spp)+
    theme_classic()

#show(p5)
#dev.off() 
   

#pdf("figs/p6_partialmort.pdf")

p6<-ggplot()+
  geom_point(data=pdat, 
    aes(x=area, y=pm_logit), 
       shape=21, col="grey")+
    geom_line(data=pdat, aes(x=area, y=logit(p_stasis), group=spp), linetype="dotted")+
   geom_line(data=p.pred, 
     aes(x=area, y=pred, col=spp))+
     scale_colour_manual(values=c(cols))+
     #facet_wrap(~spp)+
     theme_classic()


ggplot()+
  geom_point(data=pdat, 
    aes(x=area, y=inv.logit(pm_logit)), 
       shape=21, col="grey")+
    geom_line(data=pdat, aes(x=area, y=p_stasis, group=spp, col=spp), linetype="dotted")+
   geom_line(data=p.pred, 
     aes(x=area, y=inv.logit(pred), col=spp))+
     scale_colour_manual(values=c(cols))+
     facet_wrap(~spp)+
     theme_classic()

#show(p6)
#dev.off() 


#######################################
# PLOT PARAMETERS
#######################################

partheme <- theme(axis.text=element_text(size=7), axis.title=element_text(size=8),
plot.title=element_text(size=8, hjust=0.5, face="bold"), legend.title=element_blank(), legend.text=element_text(face="italic"))

xlab <- expression(log[10]*(area[~t]))


fig.s3 <- plot_grid(plot_grid(
p1+guides(col="none")+partheme+ggtitle("reproductive maturity")+
labs(x=xlab,y="probability of maturity"), 
p2+guides(col="none")+partheme+ggtitle("fecundity")+
labs(x=xlab,y=expression("eggs per colony")),
p3+guides(col="none")+partheme+ggtitle("growth rate")+
labs(x=xlab, y=expression(log[10]*(area[~t~+1]))),  
p4+guides(col="none")+partheme+ggtitle("survival rate")+
labs(x=xlab, y="annual survival rate"), 
p5+guides(col="none")+partheme+ggtitle("maximum growth")+
labs(x=xlab, y=expression(log[10]*(radial~growth))),
p6+guides(col="none")+partheme+ggtitle("partial mortality")+
labs(x=xlab, y="logit(proportion area lost)"), 
nrow=2, align="hv", labels="AUTO", label_size=8),
get_legend(p1+partheme), rel_widths=c(1,0.2))
fig.s3




