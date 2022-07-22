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
   out[x < rec.size | y > rec.size | y < rec.size-h] <- 0 #if x is below recruitment size/ y is above recruit size. 

   return(out)
   } 
    
#------------------------------- mesh
mesh <- function(){
	min.size <- rec.size # 3.5
	b <- seq(min.size, max.size, length=n)
	h <- b[2] - b[1]
	b <- c(min(b)-(2*h), min(b)-h, b)
	y <- 0.5 * (b[1:n]+b[2:(n+1)])
	I <- y >= rec.size
	return(list(b=b,h=h,y=y,I=I, rec.size = rec.size))}

 #------------------------------- kernel
pmort<-F

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
  N <- solve(diag(n)-P)
  R0 <- abs(eigen(R %*% N)$values[1])
	return(list(K=K, lam=lam, w=w, v=v, G=G, S=S, R=R, P=P, R0=R0)) }
	
