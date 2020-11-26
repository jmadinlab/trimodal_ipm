if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
if (!require("png")) {
  install.packages("png")
  library(png)
}


# binomial function
bf = function(x) {exp(x) / (1 + exp(x))}

# second-order AIC, where LL is the logLik, K is the the number of parameters in the fitted model, and n is the number of observations
# Was tested against AICc in library(AICcmodavg)
jAIC <- function(mod) {-2*logLik(mod)[1] + 2*length(coef(mod))}
jAICc <- function(mod, n) {
  K <- length(coef(mod))
  n <- length(residuals(mod))
  jAIC(mod) + 2*K*(K+1)/(n-K-1) 
}

# Comparing species and groth form model fits using stepwise model selection followed by AICc
run_test <- function(gf) {
  sp_q = glm(survival ~ area:species_code + area_2:species_code + species_code, family = binomial, data = dat[dat$form == gf,])
  spMOD <- step(sp_q, test="Chi")

  fo_q = glm(survival ~ area + area_2, family = binomial, data = dat[dat$form == gf,])
  gfMOD <- step(fo_q, test="Chi")

  spAICc <- jAICc(spMOD)
  gfAICc <- jAICc(gfMOD)
  return(list(gf=gf, spMOD=spMOD, gfMOD=gfMOD, spAICc=spAICc, gfAICc=gfAICc))
}
