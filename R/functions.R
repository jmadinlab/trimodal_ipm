# functions 

options(stringsAsFactors=FALSE)

#######################################
# IPM FUNCTIONS
#######################################

#-------------------------------survival
s.x <- function(x) { 
	u <- params$s.int[params$spp==sp] + 
	  params$s.slp[params$spp==sp] * x + 
	    params$s.slp.2[params$spp==sp] * x^2
  return(inv.logit(u)) 
}

#------------------------------- growth
g.yx <- function(y, x) {
	dnorm(y, mean=params$g.int[params$spp==sp] + 
	  params$g.slp[params$spp==sp]*x,
	    sd=sqrt(params$g.var[params$spp==sp]))
	    }

#------------------------------- growth & partial mortality
p.yx <- function(y, x) {
  # x <- -3
  g <- a_func(r_func(10^x) + params$r.int[params$spp==sp] )
  #+ 1.96 * params$r.err[params$spp==sp])
  temp <- 10^y / g  # proportion of max reached. 
  temp[temp > 1] <- 1
  dnorm(logit(1 - temp), params$p.int[params$spp==sp] + x * params$p.slp[params$spp==sp], params$p.sig[params$spp==sp])
}

#------------------------------- reproduction
 r.yx <- function(y, x) {	
 	mat<- inv.logit(params$m.int[params$spp==sp] + 
 	  params$m.slp[params$spp==sp] *x)
 	fec<- exp(params$f.int[params$spp==sp] + 
 	  params$f.slp[params$spp==sp] *x) 
   #siz<- rnorm(y,mean=params$rec.size[params$spp==sp], sd=0.05) 
   out <- (rec* mat * fec)
   out[x < rec.size | y >= rec.size] <- 0 #if x is below recruitment size
   return(out)
   } 

 #------------------------------- kernel
pmort<-T

bigmatrix <- function() {
    if (pmort) {
    G <- h * outer(y, y, p.yx)
  } else {
    G <- h * outer(y, y, g.yx)
  }
  G <- t(t(G) / apply(G, 2, sum)) # ???? 
  S <- s.x(y)
  P <- G  # placeholder, redifine P on next line
  for(i in 1:n) P[,i]=G[,i]*S[i]
  R <- h * outer(y, y, r.yx) 
  #R <- h * outer(y, pmin(y, smax), r.yx) #  ceiling
  K <- P + R
  lam <- Re(eigen(K)$values[1])
	w <- abs(Re(eigen(K)$vectors[,1])) 
	v <- abs(Re(eigen(t(K))$vectors[,1]))
	return(list(K=K, lam=lam, w=w, v=v, G=G, S=S, R=R, P=P)) }
	

#######################################
# OTHER FUNCTIONS
#######################################

inv.logit <- function(x) {exp(x)/(1+exp(x))}


circularity <- function(area, perimeter) {
  (4 * pi * area)/(perimeter^ 2)
}

a_func <- function(r) {  
  pi * r^2  
  }

r_func <- function(a) {
  sqrt(a / pi)
}


rec.ll <- function(x) {
  cnt <- size.dist$count[I] # non-recruits
  rec <<- x[1]
  mod <- bigmatrix()
  eig.vec <- mod$w[I]/sum(mod$w[I])
  return(-sum(cnt * log(eig.vec), na.rm=TRUE)) } # log-likelihood 
  
  