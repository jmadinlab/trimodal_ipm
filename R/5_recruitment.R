
#######################################
# FIRST LAM ESTIMATE
#######################################

lam.est <- NULL
for (sp in spp) {
rec.size <- params$rec.size[params$spp==sp]
rec <- params$rec[params$spp==sp]
h <- mesh()$h
y <- mesh()$y
lam.est<-rbind(lam.est, data.frame(sp, lam=bigmatrix()$lam))
}
params$lam.est <- lam.est$lam


#######################################
# CHANGE IN ABUNDANCE (TRANSECTS)
#######################################

head(abun.LIT)
diffs <- data.frame(dcast(abun.LIT, species+tran~year, value.var="N"))
diffs <- aggregate(.~species, subset(diffs, select=-c(tran)), mean)
diffs[diffs==0]<-0.05 # absent species = rare
diffs$diff <- diffs$X2014 - diffs$X2011
diffs$lam.tran <- (diffs$X2014/diffs$X2011)^(1/3)
diffs

params$lam.tran <- diffs$lam.tran[match(params$species, diffs$species)]

ggplot(params, aes(lam.est, lam.tran))+geom_text(aes(label=spp))+geom_abline(slope=1)+geom_smooth(method="lm", se=F, linetype="dotted")+scale_y_log10()+scale_x_log10()

summary(lm(log(lam.tran)~log(lam.est), params))

#######################################
# REC AT TRANSECT LAMS
#######################################	

rec.transect <- NULL
for(sp in spp){
	#sp <- "Gre"
sub <- storeDET[storeDET$spp==sp,]
r.sp <- params$lam.tran[params$spp==sp]
n.rec <- which(abs(sub$lam-r.sp)==min(abs(sub$lam-r.sp)))
if(n.rec==1){
newlam<-	min(sub$lam)+0.01
n.rec2 <- which(abs(sub$lam-newlam)==min(abs(sub$lam-newlam)))
new <- sub[n.rec2,]
} else { new <- sub[n.rec,] }
rec.transect <- rbind(rec.transect, cbind(new, lam.r=r.sp))
}

rec.transect

params$rec.tran <- rec.transect$rec[match(params$spp, rec.transect$spp)]

# rec at 1... 
rec.one <- NULL
for(sp in spp){
sub <- storeDET[storeDET$spp==sp,]
n.rec <- which(abs(sub$lam-1)==min(abs(sub$lam-1)))
rec.one <- rbind(rec.one, cbind(sub[n.rec,], lam.r=r.sp))
}
params$rec.one <- rec.one$rec[match(params$spp, rec.one$spp)]
