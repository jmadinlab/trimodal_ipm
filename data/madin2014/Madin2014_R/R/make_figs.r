# SPECIES

pdf(file = "figs/fig_bestfit.pdf", height=5, width=3.5)
  op <- par(mfrow=c(3, 2), mar=c(1, 1, 2, 1), oma=c(4, 3, 1, 0), ann=FALSE)
  
  sp <- "massive"
  dead = sum(dat$survival[dat$form == sp], na.rm = T) / sum(!is.na(dat$survival[dat$form == sp])) 
  ss = seq(min(dat$area[dat$form == sp]), max(dat$area[dat$form == sp]), 0.1)
  plot(jitter(survival, amount=0) ~ area, data = dat[dat$form == sp,], xlim = c(-8, 1), col = "grey", axes=FALSE)
  spp_g = gam(survival ~ s(area, fx=FALSE, k=-0, bs = "cr"), family = binomial, data = dat[dat$form==sp,])
  lines(ss, bf(predict(spp_g, list(area = ss, species_code = factor(rep(sp, length(ss)))))), lty = 2)
  lines(ss, bf(predict(ms$gfMOD, list(area = ss, area_2 = ss^2, species_code = factor(rep(sp, length(ss)))))), lty = 1)
  text(-1, 0.8, bquote(paste(italic("M"), " = ", .(round(dead, 3)))), cex=0.9)
  mtext("(a) Massive", side=3, adj=0, cex=0.6, line=0.2)
  axis(1, labels=FALSE)
  axis(2, at=c(0,1), las=2, cex.axis=0.9)

  z <- readPNG("figs/silh_massive.png") 
  rasterImage(z, -3.5, 0.45, -1, 0.7) 

  sp <- "digitate"
  dead = sum(dat$survival[dat$form == sp], na.rm = T) / sum(!is.na(dat$survival[dat$form == sp])) 
  ss = seq(min(dat$area[dat$form == sp]), max(dat$area[dat$form == sp]), 0.1)
  plot(jitter(survival, amount=0) ~ area, data = dat[dat$form == sp,], xlim = c(-8, 1), col = "grey", axes=FALSE)
  spp_g = gam(survival ~ s(area, fx=FALSE, k=-0, bs = "cr"), family = binomial, data = dat[dat$form==sp,])
  lines(ss, bf(predict(spp_g, list(area = ss, species_code = factor(rep(sp, length(ss)))))), lty = 2)
  lines(ss, bf(predict(dt$gfMOD, list(area = ss, area_2 = ss^2, species_code = factor(rep(sp, length(ss)))))), lty = 1)
  text(-1, 0.8, bquote(paste(italic("M"), " = ", .(round(dead, 3)))), cex=0.9)
  mtext("(b) Digitate", side=3, adj=0, cex=0.6, line=0.2)
  axis(1, labels=FALSE)
  axis(2, at=c(0,1), labels=FALSE, las=2)

  z <- readPNG("figs/silh_digitate.png") 
  rasterImage(z, -3.5, 0.45, -1, 0.7) 

  sp <- "corymbose"
  dead = sum(dat$survival[dat$form == sp], na.rm = T) / sum(!is.na(dat$survival[dat$form == sp]))	
  ss = seq(min(dat$area[dat$form == sp]), max(dat$area[dat$form == sp]), 0.1)
  plot(jitter(survival, amount=0) ~ area, data = dat[dat$form == sp,], xlim = c(-8, 1), col = "grey", axes=FALSE)
  spp_g = gam(survival ~ s(area, fx=FALSE, k=-0, bs = "cr"), family = binomial, data = dat[dat$form==sp,])
  lines(ss, bf(predict(spp_g, list(area = ss, species_code = factor(rep(sp, length(ss)))))), lty = 2)
  lines(ss, bf(predict(cb$gfMOD, list(area = ss, area_2 = ss^2, species_code = factor(rep(sp, length(ss)))))), lty = 1)
  text(-1, 0.8, bquote(paste(italic("M"), " = ", .(round(dead, 3)))), cex=0.9)
  mtext("(c) Corymbose", side=3, adj=0, cex=0.6, line=0.2)
  axis(1, labels=FALSE)
  axis(2, at=c(0,1), las=2, cex.axis=0.9)

  z <- readPNG("figs/silh_corymbose.png") 
  rasterImage(z, -3, 0.45, -1, 0.7) 

  sp <- "table"
  dead = sum(dat$survival[dat$form == sp], na.rm = T) / sum(!is.na(dat$survival[dat$form == sp])) 
  ss = seq(min(dat$area[dat$form == sp]), max(dat$area[dat$form == sp]), 0.1)
  plot(jitter(survival, amount=0) ~ area, data = dat[dat$form == sp,], xlim = c(-8, 1), col = "grey", pch=21, axes=FALSE)
  spp_g = gam(survival ~ s(area, fx=FALSE, k=-0, bs = "cr"), family = binomial, data = dat[dat$form==sp,])
  lines(ss, bf(predict(spp_g, list(area = ss, species_code = factor(rep(sp, length(ss)))))), lty = 2)
  lines(ss, bf(predict(tb$gfMOD, list(area = ss, area_2 = ss^2, species_code = factor(rep(sp, length(ss)))))), lty = 1)
  text(-1, 0.8, bquote(paste(italic("M"), " = ", .(round(dead, 3)))), cex=0.9)
  mtext("(d) Tabular", side=3, adj=0, cex=0.6, line=0.2)
  axis(1, labels=FALSE)
  axis(2, at=c(0,1), las=2, labels=FALSE)

  z <- readPNG("figs/silh_tabular.png") 
  rasterImage(z, -4, 0.5, 0, 0.7) 

  # Branching 
  sp <- "AI"
  dead = sum(dat$survival[dat$species_code == sp], na.rm = T) / sum(!is.na(dat$survival[dat$species_code == sp])) 
  ss = seq(min(dat$area[dat$species_code == sp]), max(dat$area[dat$species_code == sp]), 0.1)
  plot(jitter(survival, amount=0) ~ area, data = dat[dat$species_code == sp,], xlim = c(-8, 1), col = "grey", axes=FALSE)
  spp_g = gam(survival ~ s(area, fx=FALSE, k=-0.5, bs = "cr"), family = binomial, data = dat[dat$species_code==sp,])
  lines(ss, bf(predict(spp_g, list(area = ss, species_code = factor(rep(sp, length(ss)))))), lty = 2)
  lines(ss, bf(predict(br$spMOD, list(area = ss, area_2 = ss^2, species_code = factor(rep(sp, length(ss)))))), lty = 1)
  text(-1, 0.8, bquote(paste(italic("M"), " = ", .(round(dead, 3)))), cex=0.9)
  mtext(expression(paste("(e) ", italic("Acropora intermedia"))), side=3, adj=0, cex=0.6, line=0.2)
  axis(1, cex.axis=0.9)
  axis(2, at=c(0,1), las=2, cex.axis=0.9)
  
  z <- readPNG("figs/silh_branching.png") 
  rasterImage(z, 0.5, 0.25, 2.5, 0.75, xpd = TRUE) 

  sp <- "AR"
  dead = sum(dat$survival[dat$species_code == sp], na.rm = T) / sum(!is.na(dat$survival[dat$species_code == sp])) 
  ss = seq(min(dat$area[dat$species_code == sp]), max(dat$area[dat$species_code == sp]), 0.1)
  plot(jitter(survival, amount=0) ~ area, data = dat[dat$species_code == sp,], xlim = c(-8, 1), col = "grey", axes=FALSE)
  spp_g = gam(survival ~ s(area, fx=FALSE, k=-0.5, bs = "cr"), family = binomial, data = dat[dat$species_code==sp,])
  lines(ss, bf(predict(spp_g, list(area = ss, species_code = factor(rep(sp, length(ss)))))), lty = 2)
  lines(ss, bf(predict(br$spMOD, list(area = ss, area_2 = ss^2, species_code = factor(rep(sp, length(ss)))))), lty = 1)
  text(-1, 0.8, bquote(paste(italic("M"), " = ", .(round(dead, 3)))), cex=0.9)
  mtext(expression(paste("(f) ", italic("Acropora robusta"))), side=3, adj=0, cex=0.6, line=0.2)
  axis(1, cex.axis=0.9)
  axis(2, at=c(0,1), las=2, labels=FALSE)

  mtext(expression(paste("Colony size ",m^2, " (ln)", sep="")), side=1, outer=TRUE, line=2, cex=0.7)
  mtext("Yearly mortality", side=2, outer=TRUE, line=1.5, cex=0.7)

dev.off()


# CSF vs MORTALITY
# pdf(file = "figs/fig_csf.pdf")

#   par(mfrow=c(2, 3))
#   store = c()
#   for (pp in seq(0.5, 0.01, -0.1)) {
#     ff = c(); csf_mn = c(); csf_se = c()
#     mor_mn = c(); mor_se = c(); mor_lo = c()

#     for (sp in c("AH", "AC", "AD", "AS", "AI", "AR", "GR", "GP", "AL", "AM", "AN")) {

#       temp1 = dat[dat$species_code == sp & !is.na(dat$survival),]
#       tt = round(dim(temp1)[1] * pp) # largest 20% colonies for species
#       tp = temp1[order(temp1$area, decreasing = T),][1:tt,]
#       ff = c(ff, sp)
#       csf_mn = c(csf_mn, csf[which(rownames(csf)==as.vector(tp$form[1])),"mean"])
#       csf_se = c(csf_se, csf[which(rownames(csf)==as.vector(tp$form[1])),1]-csf[which(rownames(csf)==as.vector(tp$form[1])),2])
#       mor_mn = c(mor_mn, qbinom(0.5, tt, mean(tp$survival))/(tt))
#       mor_se = c(mor_se, (qbinom(0.84, tt, mean(tp$survival))/(tt) - qbinom(0.5, tt, mean(tp$survival))/(tt))/sqrt(tt))
#       #mor_lo = c(mor_lo, qbinom(0.16, tt, mean(tp$survival))/(tt))
#     }

#     plot(csf_mn, mor_mn, ylim = c(0, 0.4), xlim = c(-1, 1.5), pch=20, cex=1.2, axes=FALSE, xlab="Slope of size/storm vulnerability relationship ±95% CIs", ylab=paste("Yearly mortality of largest ", pp*100,"% of colonies ±95% CIs", sep=""), main=pp)
#     text(csf_mn+0.1, mor_mn+0.01, ff, cex=0.8)
#     arrows(csf_mn + csf_se/1.96, mor_mn, csf_mn - csf_se/1.96, mor_mn, code=3, angle=90, length=0.05)
#     arrows(csf_mn, mor_mn + mor_se, csf_mn, mor_mn - mor_se, code=3, angle=90, length=0.05)
#     axis(1, las=1)
#     axis(2, at=c(0, 0.2, 0.4), las=2)
#     x = csf_mn[ff != "AI" & ff != "AR"]
#     y = mor_mn[ff != "AI" & ff != "AR"]
    
#     mod = cor.test(csf_mn, mor_mn)
#     temp <- round(cbind(mod$stat, mod$para, mod$p.val, mod$est), 3)
#     print(summary(mod))

#     mod = cor.test(x, y)
#     temp <- cbind(temp, round(cbind(mod$stat, mod$para, mod$p.val, mod$est), 3))
#     store <- rbind(store, temp)
#     print(summary(mod))

#     abline(mod, lty=2)

#     # readline()
#   }
#   write.table(store, "output/csf_mort.csv", sep = ",")

# dev.off()

# CSF 

pdf(file = "figs/fig_csf_relationships.pdf", width=5.5, height=5.5)

  cc <- c()
  forms <- c("branching", "table", "corymbose", "digitate", "massive")
  col_vec <- brewer.pal(length(forms), "Set1")

  par(mar=c(5, 5, 1, 2))
  plot(0,0, xlim=c(-8, 2), ylim=c(-2, 9), xlab=expression(paste("Colony size ",m^2, " (ln)", sep="")), ylab="Colony Shape Factor (ln)", axes=FALSE, type="n")
  
  for (i in 1:length(forms)) {
    ss <- csf[(csf$exp=="E4" | csf$exp=="E5") & csf$gf==forms[i],]
    colt <- col_vec[i]
    points(lcsf ~ larea, data=ss, pch=20, col=colt)
    lines(ss$larea, predict(csf_mod, list(larea=ss$larea, gf=factor(rep(forms[i], nrow(ss)))), se=F), col=col_vec[i])
    cc <- c(cc, nrow(ss))
  }
  
  z <- readPNG("figs/silh_branching.png") 
  rasterImage(z, -1, 7, 0, 9, xpd = TRUE) 
  z <- readPNG("figs/silh_tabular.png") 
  rasterImage(z, 0, 4.5, 2, 5.3, xpd = TRUE) 
  z <- readPNG("figs/silh_corymbose.png") 
  rasterImage(z, -0.9, 2, 0.2, 3, xpd = TRUE) 
  z <- readPNG("figs/silh_digitate.png") 
  rasterImage(z, -1.3, 0.5, 0, 1.5, xpd = TRUE) 
  z <- readPNG("figs/silh_massive.png") 
  rasterImage(z, 1, -2.2, 2.2, -1.3, xpd = TRUE) 

  axis(1)
  axis(2, las=2)

dev.off()



### CSF / MORTALITY

pdf(file = "figs/fig_csf_mort.pdf", width=5.5, height=5.5)
  par(mar=c(5, 5, 1, 2))
  pp <- 0.2
  ff = c(); csf_mn = c(); csf_se = c()
  mor_mn = c(); mor_se = c(); mor_lo = c()

  for (sp in c("AH", "AC", "AD", "AS", "AI", "AR", "GR", "GP", "AL", "AM", "AN")) {

    temp1 = dat[dat$species_code == sp & !is.na(dat$survival),]
    tt = round(dim(temp1)[1] * pp) # largest 20% colonies for species
    tp = temp1[order(temp1$area, decreasing = T),][1:tt,]
    ff = c(ff, sp)

    csf_temp <- mean(predict(csf_mod, list(larea = tp$area, gf = factor(rep(as.vector(tp$form[1]), length(tp$area)))), se = F))
    csf_mn <- c(csf_mn, csf_temp)
    
    mor_mn_temp <- qbinom(0.5, tt, mean(tp$survival))/(tt)
    mor_mn_temp <- sum(tp$survival)/nrow(tp)
    mor_mn = c(mor_mn, mor_mn_temp)
    mor_se_temp <- (qbinom(0.84, tt, mean(tp$survival))/(tt) - qbinom(0.5, tt, mean(tp$survival))/(tt))/sqrt(tt)
    mor_se_temp <- mor_mn_temp*(1-mor_mn_temp)/(nrow(tp)-1)
    mor_se = c(mor_se, mor_se_temp)  
  }

  plot(csf_mn, mor_mn, ylim = c(0, 0.4), xlim = c(-2, 10), pch=20, cex=1.2, axes=FALSE, xlab=expression("Mechanical vulnerability of largest colonies (ln"~italic(CSF)~")"), ylab=paste("Yearly mortality of largest colonies", sep=""))

  axis(1, las=1)
  axis(2, at=c(0, 0.2, 0.4), las=2)
  x = csf_mn[ff != "AI" & ff != "AR"]
  y = mor_mn[ff != "AI" & ff != "AR"]
  
  mod = lm(y ~ x)
  lines(seq(-1, 5, 0.2), predict(mod, list(x=seq(-1, 5, 0.2))), lty=2)

  text(8, 0.07, "Aro", cex=0.8)
  text(8, 0.011, "Ain", cex=0.8)
  z <- readPNG("figs/silh_branching.png") 
  rasterImage(z, 8.5, 0, 9.8, 0.09, xpd = TRUE) 

  text(4.5, 0.26, "Ahy", cex=0.8)
  text(4, 0.3, "Acy", cex=0.8)
  z <- readPNG("figs/silh_tabular.png") 
  rasterImage(z, 4.5, 0.26, 7.6, 0.292, xpd = TRUE) 

  text(1.4, 0.215, "Ami", cex=0.8)
  text(1.5, 0.185, "Asp", cex=0.8)
  text(2.5, 0.145, "Ana", cex=0.8)
  z <- readPNG("figs/silh_corymbose.png") 
  rasterImage(z, 2.6, 0.155, 4.4, 0.2, xpd = TRUE) 

  text(1.2, 0.145, "Ahu", cex=0.8)
  text(1.5, 0.07, "Adi", cex=0.8)
  z <- readPNG("figs/silh_digitate.png") 
  rasterImage(z, -0.7, 0.08, 1, 0.125, xpd = TRUE) 

  text(-0.75, 0.015, "Gre", cex=0.8)
  text(0.5, 0.013, "Gpe", cex=0.8)
  z <- readPNG("figs/silh_massive.png") 
  rasterImage(z, -2.4, -0.01, -0.7, 0.025, xpd = TRUE) 

dev.off()

#####

par(mfrow=c(2, 5))
store = c()
sean = c()

for (pp in seq(0.5, 0.01, -0.1)) {
  ff = c(); csf_mn = c(); csf_se = c()
  mor_mn = c(); mor_se = c(); mor_lo = c()

  for (sp in c("AH", "AC", "AD", "AS", "AI", "AR", "GR", "GP", "AL", "AM", "AN")) {

    temp1 = dat[dat$species_code == sp & !is.na(dat$survival),]
    tt = round(dim(temp1)[1] * pp) # largest 20% colonies for species
    tp = temp1[order(temp1$area, decreasing = T),][1:tt,]
    ff = c(ff, sp)

	csf_temp <- mean(predict(csf_mod, list(larea = tp$area, gf = factor(rep(as.vector(tp$form[1]), length(tp$area)))), se = F))
    csf_mn <- c(csf_mn, csf_temp)
    
	mor_mn_temp <- qbinom(0.5, tt, mean(tp$survival))/(tt)
	mor_mn_temp <- sum(tp$survival)/nrow(tp)
    mor_mn = c(mor_mn, mor_mn_temp)
	mor_se_temp <- (qbinom(0.84, tt, mean(tp$survival))/(tt) - qbinom(0.5, tt, mean(tp$survival))/(tt))/sqrt(tt)
	mor_se_temp <- mor_mn_temp*(1-mor_mn_temp)/(nrow(tp)-1)
    mor_se = c(mor_se, mor_se_temp)
	
	sean <- rbind(sean, data.frame(group="largest", proportion=pp, species=sp, csf_mean=csf_temp, mort_est=mor_mn_temp, mort_var=mor_se_temp, n=tt))
	
  }

  # plot(csf_mn, mor_mn, ylim = c(0, 0.4), xlim = c(-2, 10), pch=20, cex=1.2, axes=FALSE, xlab="Mean mechanical vulnerability of largest colonies", ylab=paste("Yearly mortality of largest ", pp*100,"% of colonies ±95% CIs", sep=""), main=paste("Largest", pp))
  # text(csf_mn+0.1, mor_mn+0.01, ff, cex=0.8)

  # arrows(csf_mn, mor_mn + mor_se, csf_mn, mor_mn - mor_se, code=3, angle=90, length=0.05)
  # axis(1, las=1)
  # axis(2, at=c(0, 0.2, 0.4), las=2)
  x = csf_mn[ff != "AI" & ff != "AR"]
  y = mor_mn[ff != "AI" & ff != "AR"]
  
  mod = cor.test(csf_mn, mor_mn)
  temp <- round(cbind(mod$stat, mod$para, mod$p.val, mod$est), 3)
  print(summary(mod))

  mod = cor.test(x, y)
  temp <- cbind(temp, round(cbind(mod$stat, mod$para, mod$p.val, mod$est), 3))
  store <- rbind(store, temp)
  print(summary(mod))

  abline(mod, lty=2)

  # readline()
}

for (pp in seq(0.5, 0.01, -0.1)) {
  ff = c(); csf_mn = c(); csf_se = c()
  mor_mn = c(); mor_se = c(); mor_lo = c()

  for (sp in c("AH", "AC", "AD", "AS", "AI", "AR", "GR", "GP", "AL", "AM", "AN")) {

    temp1 = dat[dat$species_code == sp & !is.na(dat$survival),]
    tt = round(dim(temp1)[1] * pp) # largest 20% colonies for species
    tp = temp1[order(temp1$area, decreasing = F),][1:tt,]
    ff = c(ff, sp)

	csf_temp <- mean(predict(csf_mod, list(larea = tp$area, gf = factor(rep(as.vector(tp$form[1]), length(tp$area)))), se = F))
    csf_mn <- c(csf_mn, csf_temp)

	mor_mn_temp <- qbinom(0.5, tt, mean(tp$survival))/(tt)
	mor_mn_temp <- sum(tp$survival)/nrow(tp)
    mor_mn = c(mor_mn, mor_mn_temp)
	mor_se_temp <- (qbinom(0.84, tt, mean(tp$survival))/(tt) - qbinom(0.5, tt, mean(tp$survival))/(tt))/sqrt(tt)
	mor_se_temp <- mor_mn_temp*(1-mor_mn_temp)/(nrow(tp)-1)
    mor_se = c(mor_se, mor_se_temp)
	
	sean <- rbind(sean, data.frame(group="smallest", proportion=pp, species=sp, csf_mean=csf_temp, mort_est=mor_mn_temp, mort_var=mor_se_temp, n=tt))
  }

  # plot(csf_mn, mor_mn, ylim = c(0, 0.5), xlim = c(-2, 10), pch=20, cex=1.2, axes=FALSE, xlab="Mean mechanical vulnerability of largest colonies", ylab=paste("Yearly mortality of largest ", pp*100,"% of colonies ±95% CIs", sep=""), main=paste("Smallest", pp))
  # text(csf_mn+0.1, mor_mn+0.01, ff, cex=0.8)

  # arrows(csf_mn, mor_mn + mor_se, csf_mn, mor_mn - mor_se, code=3, angle=90, length=0.05)
  # axis(1, las=1)
  # axis(2, at=c(0, 0.2, 0.4), las=2)
  x = csf_mn[ff != "AI" & ff != "AR"]
  y = mor_mn[ff != "AI" & ff != "AR"]
  
  mod = cor.test(csf_mn, mor_mn)
  temp <- round(cbind(mod$stat, mod$para, mod$p.val, mod$est), 3)
  print(summary(mod))

  mod = cor.test(x, y)
  temp <- cbind(temp, round(cbind(mod$stat, mod$para, mod$p.val, mod$est), 3))
  store <- rbind(store, temp)
  print(summary(mod))

  abline(mod, lty=2)

  # readline()
}

save(sean, file="output/sean.RData")
write.table(store, "output/csf_mort_new.csv", sep = ",")

# dev.off()

