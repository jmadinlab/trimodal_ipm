# functions 

options(stringsAsFactors=FALSE)

library("MASS") 
library("reshape2")
library("ggplot2")
library("cowplot")
library("stats4")

inv.logit <- function(x) {exp(x)/(1+exp(x))}




#------------------------------- growth functions


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

