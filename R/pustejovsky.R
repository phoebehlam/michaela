#' cohen's d from dichotomizing/extreme group design to r
#'
#' cohen's d from dichotomizing/extreme group design to correlation coefficient\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns extreme group r (referred to as r (subscript eg) in pustejovsky paper)
#'
#' @param d cohen's d
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#' +1 for positive associations (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
#'
#' @examples
#' dicho.d.r(.80, 1/3, 1/3, 15, 15, -1)
#'
#' @export
dicho.d.r <- function (d, p1, p2, n1, n2, dir) {
  #auxilliary stats
  f = n1/(n1+n2)
  c1 = qnorm(p1)
  c2 = qnorm(1-p2)
  v1 = dnorm(c1)/p1
  v2 = dnorm(c2)/p2
  a = ((v1+v2)^2) / (f*v1*(v1+c1)+(1-f)*v2*(v2-c2))
  b = 1/sqrt(f*v1*(v1+c1)+(1-f)*v2*(v2-c2))
  
  r.eg = (b*d)/sqrt(d^2 + a) * dir
  
  return(r.eg)
}

#' cohen's d from dichotomizing/extreme group design to z
#'
#' cohen's d from dichotomizing/extreme group design to fisher's z\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns z based on taylor series approximation
#'
#' @param d cohen's d
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param k kth term for taylor series approximation (enter 1 to 5)
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#' +1 for positive associations (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
#'
#' @examples
#' dicho.d.z(.80, 1/3, 1/3, 15, 15, 5, -1)
#'
#' @export
dicho.d.z <- function (d, p1, p2, n1, n2, k, dir) {
  
  #auxilliary stats
  f = n1/(n1+n2)
  c1 = qnorm(p1)
  c2 = qnorm(1-p2)
  v1 = dnorm(c1)/p1
  v2 = dnorm(c2)/p2
  a = ((v1+v2)^2) / (f*v1*(v1+c1)+(1-f)*v2*(v2-c2))
  b = 1/sqrt(f*v1*(v1+c1)+(1-f)*v2*(v2-c2))
  
  r.eg = (b*d)/sqrt(d^2 + a) * dir
  r.pbs = d/sqrt(d^2+a) * dir
  
  #taylor series approximation, up to 5
  k0 = ((log(1+r.pbs) - log(1-r.pbs))/2)*(((b-1)^0)*(r.pbs^0)/factorial(0))
  k1 = (1/(1-r.pbs^2)) * (((b-1)^1)*(r.pbs^1)/factorial(1))
  k2 = (2*r.pbs)/(1-r.pbs^2)^2 * (((b-1)^2)*(r.pbs^2)/factorial(2))
  k3 = (2+6*r.pbs^2)/(1-r.pbs^2)^3 * (((b-1)^3)*(r.pbs^3)/factorial(3))
  k4 = (24*r.pbs + 24*r.pbs^3)/(1-r.pbs^2)^4 * (((b-1)^4)*(r.pbs^4)/factorial(4))
  k5 = (24 + 240*r.pbs^2 + 120*r.pbs^4)/(1-r.pbs^2)^5 * (((b-1)^5)*(r.pbs^5)/factorial(5))
  
  if (k==1) {
    z1 = k0+k1
    return(z1)
  } else if (k==2) {
    z2 = k0+k1+k2
    return(z2)
  } else if (k==3) {
    z3 = k0+k1+k2+k3
    return (z3)
  } else if (k==4) {
    z4 = k0+k1+k2+k3+k4
    return (z4)
  } else if (k==5) {
    z5 = k0+k1+k2+k3+k4+k5
    return (z5)
  }
  
}

#' variance of z converted from dichotomizing/extreme group design d
#'
#' compute the variance of z from a cohen's d derived from dichotomizing/extreme group design\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#'
#' @param d cohen's d
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param type select one (read \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} for recommendations\cr
#' there are situations where the "controlled experiment" formulas are recommended even for stats derived from dichotomize/extreme group designs\cr
#' so please refer to his paper)\itemize{
#' \item"vd.ce" compute variance of d using the controlled experiment formulas
#' \item"vd.eg" compute variance of d using the extreme group formulas}
#'
#' @examples
#' dichod.zvar(.80, 1/3, 1/3, 15, 15, "vd.ce")
#'
#' @export
dichod.zvar <- function (d, p1, p2, n1, n2, type = c("vd.ce", "vd.eg")) {
  
  #auxilliary stats
  f = n1/(n1+n2)
  c1 = qnorm(p1)
  c2 = qnorm(1-p2)
  v1 = dnorm(c1)/p1
  v2 = dnorm(c2)/p2
  a = ((v1+v2)^2) / (f*v1*(v1+c1)+(1-f)*v2*(v2-c2))
  b = 1/sqrt(f*v1*(v1+c1)+(1-f)*v2*(v2-c2))
  
  r.eg = (b*d)/sqrt(d^2 + a)
  r.pbs = d/sqrt(d^2+a)
  
  #vd.eg appendix B7 equation and #20
  k1.1 = -v1*r.eg
  k1.2 = -v2*r.eg
  k1 = k1.1 + k1.2
  
  k2.1 = 1-r.eg^2*v1*(v1+c1)
  k2.2 = 1- r.eg^2*v2*(v2-c2)
  k2 = f*k2.1+(1-f)*k2.2
  
  k3.1 = -(r.eg^3)*v1*(v1*(v1+c1)+(v1+c1)^2-1)
  k3.2 = -(r.eg^3)*v2*(v2*(v2-c2)+(v2-c2)^2-1)
  k3 = k3.1 + k3.2
  
  k4.1 = -(r.eg^4)*v1*((v1+c1)^3+4*v1*(v1+c1)^2+(v1^2-3)*(v1+c1)-v1)
  k4.2 = -(r.eg^4)*v2*((v2-c2)^3+4*v2*(v2-c2)^2+(v2^2-3)*(v2-c2)-v2)
  k4 = f*k4.1+(1-f)*k4.2
  
  vdeg = (((k2.1/f+k2.2/(1-f))/k2)-(k1*k3/k2^2)+(k1^2*(f*k2.1^2+(1-f)*k2.2^2)/(2*k2^3))+(k1^2*k4/(4*k2^3)))/(n1+n2)
  vzeg = (a^2*b^2)/((d^2 + a)*(d^2*(1-b^2)+a)^2) * vdeg
  
  #vd.ce, equation #14 and then #20
  vdce = (n1+n2)/(n1*n2) + d^2/(2*(n1+n2))
  vzce = (a^2*b^2)/((d^2 + a)*(d^2*(1-b^2)+a)^2) * vdce
  
  if (type == "vd.ce") {
    return (vzce)
  } else if (type == "vd.eg") {
    return(vzeg)
  }
  
}

#' means and sds from dichotomizing/extreme group design to r
#'
#' means and sds from dichotomizing/extreme group design to correlation coefficient\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns extreme group r (referred to as r (subscript eg) in pustejovsky paper)
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#' +1 for positive associations (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
#'
#' @examples
#' dicho.meansd.r(15, 13, 1.4, 0.9, 15, 15, 1/3, 1/3)
#'
#' @export
dicho.meansd.r <- function (m1, m2, sd1, sd2, n1, n2, p1, p2) {
  d = (m1-m2)/sqrt ((sd1^2+sd2^2)/2)
  r = dicho.d.r(d, p1, p2, n1, n2, 1) #direction always one, because the d already comes with sign
  return (r)
}

#' means and sd from dichotomizing/extreme group design to z
#'
#' means and sd from dichotomizing/extreme group design to fisher's z\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns extreme group r (referred to as r (subscript eg) in pustejovsky paper) \cr
#' for taylor series approximation, defaults 5th term (use dicho.d.z() if you'd like to adjust to lower term)
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#' +1 for positive associations (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
#'
#' @examples
#' dicho.meansd.z(15, 13, 1.4, 0.9, 15, 15, 1/3, 1/3)
#'
#' @export
dicho.meansd.z <- function (m1, m2, sd1, sd2, n1, n2, p1, p2) {
  d = (m1-m2)/sqrt ((sd1^2+sd2^2)/2)
  z = dicho.d.z(d, p1, p2, n1, n2, 5, 1) #direction always one, because the d already comes with sign
  return (z)
}

#' variance of z converted from dichotomizing/extreme group design mean and sd
#'
#' compute the variance of z from means and sds derived from dichotomizing/extreme group design d\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param type select one (read \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} for recommendations\cr
#' there are situations where the "controlled experiment" formulas are recommended even for stats derived from dichotomize/extreme group designs\cr
#' so please refer to his paper)\itemize{
#' \item"vd.ce" compute variance of d using the controlled experiment formulas
#' \item"vd.eg" compute variance of d using the extreme group formulas}
#'
#' @examples
#' dicho.meansd.zvar(15, 13, 1.4, 0.9, 15, 15, 1/3, 1/3, "vd.eg")
#'
#' @export
dicho.meansd.zvar <- function (m1, m2, sd1, sd2, n1, n2, p1, p2, type = c("vd.ce", "vd.eg")) {
  d = (m1-m2)/sqrt ((sd1^2+sd2^2)/2)
  if (type == "vd.ce") {
    zvar = dichod.zvar(d, p1, p2, n1, n2, type = "vd.ce")
  } else if (type == "vd.eg") {
    zvar = dichod.zvar(d, p1, p2, n1, n2, type = "vd.eg")
  }
  return (zvar)
}

#' t from dichotomizing/extreme group design to r
#'
#' t-statistics from dichotomizing/extreme group design to correlation coefficient\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns extreme group r (referred to as r (subscript eg) in pustejovsky paper) \cr
#' for unadjusted statistics only
#'
#' @param t t-statistics
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#' +1 for positive associations (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
#'
#' @examples
#' dicho.t.r(2.05, 1/3, 1/3, 15, 15, -1)
#'
#' @export
dicho.t.r <- function (t, totn, n1, n2, p1, p2, dir) {
  
  d = (2*t)/sqrt(totn-2)
  
  if(missing(dir)) {
    r = NA_real_
  } else {
    r = dicho.d.r(d, p1, p2, n1, n2, dir)
  }
  
  return (r)
}

#' t from dichotomizing/extreme group design to z
#'
#' t-statistics from dichotomizing/extreme group design to fisher's z\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns z based on taylor series approximation, defaults 5th term (use dicho.d.z() if you'd like to adjust to lower term)\cr
#' for unadjusted statistics only
#'
#' @param t t-statistics
#' @param totn the total n
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param dir if t-statistics is in absolute value, provide the empirical direction:\cr
#' +1 for positive associations (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
#'
#' @examples
#' dicho.t.z(2.05, 15+15, 15, 15, 1/3, 1/3, -1)
#'
#' @export
dicho.t.z <- function (t, totn, n1, n2, p1, p2, dir) {
  d = (2*t)/sqrt(totn-2)
  
  if(missing(dir)) {
    z = NA_real_
  } else {
    z = dicho.d.z(d, p1, p2, n1, n2, 5, dir)
  }
  
  return (z)
}

#' variance of z converted from dichotomizing/extreme group design t
#'
#' compute the variance of z from t-statistics derived from dichotomizing/extreme group design t\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#'
#' @param t t-statistics
#' @param totn the total n
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param type select one from below, please read \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.}for recommendations\cr
#' there are situations where the "controlled experiment" formulas are recommended even for stats derived from dichotomize/extreme group designs\cr
#' so please refer to his paper\itemize{
#' \item"vd.ce" compute variance of d using the controlled experiment formulas
#' \item"vd.eg" compute variance of d using the extreme group formulas}
#'
#' @examples
#' dicho.t.zvar(2.05, 15+15, 15, 15, 1/3, 1/3, "vd.eg")
#'
#' @export
dicho.t.zvar <- function (t, totn, n1, n2, p1, p2, type = c("vd.ce", "vd.eg")) {
  d = (2*t)/sqrt(totn-2)
  
  if (type == "vd.ce") {
    zvar = dichod.zvar(d, p1, p2, n1, n2, type = "vd.ce")
  } else if (type == "vd.eg") {
    zvar = dichod.zvar(d, p1, p2, n1, n2, type = "vd.eg")
  }
  
  return (zvar)
}

#' means and ci from dichotomizing/extreme group design to r
#'
#' means and ci from dichotomizing/extreme group design to correlation coefficient\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns extreme group r (referred to as r (subscript eg) in pustejovsky paper)
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param cil1 lower bound of 95 confidence interval for group 1 mean
#' @param ciu1 upper bound of 95 confidence interval for group 1 mean
#' @param cil2 lower bound of 95 confidence interval for group 2 mean
#' @param ciu2 upper bound of 95 confidence interval for group 2 mean
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param k number of predictors
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param result choose one:\itemize{
#' \item"cil" estimates se from lower bound of ci
#' \item"ciu" estimates se from upper bound of ci
#' \item"sediff" computes the difference in estimated se computed from upper vs. lower bound ci
#' }
#'
#' @examples
#' dicho.meanci.r(5, 3, 3, 6, 1, 4, 30, 31, 0, .33, .33)
#'
#' @export
dicho.meanci.r <- function (m1, m2, cil1, ciu1, cil2, ciu2, n1, n2, k, p1, p2, result = c("sediff", "cil", "ciu")){
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors); for raw means without covariates, enter k = 0')
  }#k = 0 would give df = n - 1 for each mean
  
  se1_ciu = suppressWarnings(ci.se(m1, cil1, ciu1, n1, k, result = "ciu"))
  se1_cil = suppressWarnings(ci.se(m1, cil1, ciu1, n1, k, result = "cil"))
  
  se2_ciu = suppressWarnings(ci.se(m2, cil2, ciu2, n2, k, result = "ciu"))
  se2_cil = suppressWarnings(ci.se(m2, cil2, ciu2, n2, k, result = "cil"))
  
  if (missing (result)) {
    
    r = dicho.meanse.r(m1, m2, se1_ciu, se2_ciu, n1, n2, p1, p2)
    
    warning ('hi! please use result = "sediff" with this same function to check se discreapncies')
    return (r)
    
  } else if (result == "ciu") {
    r = dicho.meanse.r(lm1, 2, se1_ciu, se2_ciu, n1, n2, p1, p2)
    
    warning ('hi! please use result = "sediff" with this same function to check se discreapncies')
    return (r)
  } else if (result == "cil") {
    r = dicho.meanse.r(lm1, m2, se1_cil, se2_cil ,n1, n2, p1, p2)
    
    warning ('hi! please use result = "sediff" with this same function to check se discreapncies')
    return (r)
  } else if (result == "sediff") {
    diff1 = abs(se1_ciu - se1_cil)
    diff2 = abs(se2_ciu - se2_cil)
    
    diff1.df <- as.data.frame (table (diff1))
    diff2.df <- as.data.frame (table (diff2))
    
    list <- list ("se 1 diff" = diff1.df, "se 2 diff" = diff2.df)
    
    return (list)
  }
}

#' means and ci from dichotomizing/extreme group design to z
#'
#' means and confidence interval from dichotomizing/extreme group design to fisher's z\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns extreme group r (referred to as r (subscript eg) in pustejovsky paper) \cr
#' for taylor series approximation, defaults 5th term (use dicho.d.z() if you'd like to adjust to lower term)
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param cil1 lower bound of 95 confidence interval for group 1 mean
#' @param ciu1 upper bound of 95 confidence interval for group 1 mean
#' @param cil2 lower bound of 95 confidence interval for group 2 mean
#' @param ciu2 upper bound of 95 confidence interval for group 2 mean
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param k number of predictors
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#'
#' @examples
#' dicho.meanci.z(5, 3, 3, 6, 1, 4, 30, 31, 0, .33, .33)
#'
#' @export
dicho.meanci.z <- function (m1, m2, cil1, ciu1, cil2, ciu2, n1, n2, k, p1, p2, result = c("sediff", "cil", "ciu")){
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors); for raw means without covariates, enter k = 0')
  }#k = 0 would give df = n - 1 for each mean
  
  se1_ciu = suppressWarnings(ci.se(m1, cil1, ciu1, n1, k, result = "ciu"))
  se1_cil = suppressWarnings(ci.se(m1, cil1, ciu1, n1, k, result = "cil"))
  
  se2_ciu = suppressWarnings(ci.se(m2, cil2, ciu2, n2, k, result = "ciu"))
  se2_cil = suppressWarnings(ci.se(m2, cil2, ciu2, n2, k, result = "cil"))
  
  if (missing (result)) {
    
    z = dicho.meanse.z(m1, m2, se1_ciu, se2_ciu, n1, n2, p1, p2)
    
    warning ('hi! please use result = "sediff" with this same function to check se discreapncies')
    return (z)
    
  } else if (result == "ciu") {
    z = dicho.meanse.z(m1, m2, se1_ciu, se2_ciu, n1, n2, p1, p2)
    
    warning ('hi! please use result = "sediff" with this same function to check se discreapncies')
    return (z)
  } else if (result == "cil") {
    z = dicho.meanse.z(m1, m2, se1_cil, se2_cil ,n1, n2, p1, p2)
    
    warning ('hi! please use result = "sediff" with this same function to check se discreapncies')
    return (z)
  } else if (result == "sediff") {
    diff1 = abs(se1_ciu - se1_cil)
    diff2 = abs(se2_ciu - se2_cil)
    
    diff1.df <- as.data.frame (table (diff1))
    diff2.df <- as.data.frame (table (diff2))
    
    list <- list ("se 1 diff" = diff1.df, "se 2 diff" = diff2.df)
    
    return (list)
  }
}

#' variance of z converted from dichotomizing/extreme group design mean and ci
#'
#' compute the variance of z from t-statistics derived from dichotomizing/extreme group design mean and ci\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param cil1 lower bound of 95 confidence interval for group 1 mean
#' @param ciu1 upper bound of 95 confidence interval for group 1 mean
#' @param cil2 lower bound of 95 confidence interval for group 2 mean
#' @param ciu2 upper bound of 95 confidence interval for group 2 mean
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param k number of predictors
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param type select one from below, please read \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.}for recommendations\cr
#' there are situations where the "controlled experiment" formulas are recommended even for stats derived from dichotomize/extreme group designs\cr
#' so please refer to his paper\itemize{
#' \item"vd.ce" compute variance of d using the controlled experiment formulas
#' \item"vd.eg" compute variance of d using the extreme group formulas}
#'
#' @examples
#' dicho.meanci.zvar(5, 3, 3, 6, 1, 4, 30, 31, 0, .33, .33, "vd.eg")
#'
#' @export
dicho.meanci.zvar <- function (m1, m2, cil1, ciu1, cil2, ciu2, n1, n2, k, p1, p2, type = c("vd.ce", "vd.eg")){
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors); for raw means without covariates, enter k = 0')
  }#k = 0 would give df = n - 1 for each mean
  
  se1_ciu = suppressWarnings(ci.se(m1, cil1, ciu1, n1, k, result = "ciu"))
  se2_ciu = suppressWarnings(ci.se(m2, cil2, ciu2, n2, k, result = "ciu"))
  
  if (type == "vd.ce") {
    zvar = dicho.meanse.zvar(m1, m2, se1_ciu, se2_ciu, n1, n2, p1, p2, type = "vd.ce")
  } else if (type == "vd.eg") {
    zvar = dicho.meanse.zvar(m1, m2, se1_ciu, se2_ciu, n1, n2, p1, p2, type = "vd.eg")
  }
  
  return(zvar)
  
}

#' means and se from dichotomizing/extreme group design to r
#'
#' means and standard errors from dichotomizing/extreme group design to correlation coefficient\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns extreme group r (referred to as r subscript eg in pustejovsky paper)
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param se1 standard error of group 1
#' @param se2 standard error of group 2
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' 
#' @examples
#' dicho.meanse.r(1.1, 0.5, 0.6, 0.2, 10, 18, .25, .25)
#' 
#' @export
dicho.meanse.r <- function (m1, m2, se1, se2, n1, n2, p1, p2) {
  sd1 = se.sd(se1, n1)
  sd2 = se.sd(se2, n2)
  r = dicho.meansd.r (m1, m2, sd1, sd2, n1, n2, p1, p2)
  
  return (r)
  
}

#' means and se from dichotomizing/extreme group design to z
#'
#' means and se from dichotomizing/extreme group design to fisher's z\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#' returns extreme group r (referred to as r (subscript eg) in pustejovsky paper) \cr
#' for taylor series approximation, defaults 5th term (use dicho.d.z() if you'd like to adjust to lower term)
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param se1 standard error of group 1
#' @param se2 standard error of group 2
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#'
#' @examples
#' dicho.meanse.z(1.1, 0.5, 0.6, 0.2, 10, 18, .25, .25)
#'
#' @export
dicho.meanse.z <- function (m1, m2, se1, se2, n1, n2, p1, p2) {
  sd1 = se.sd(se1, n1)
  sd2 = se.sd(se2, n2)
  z = dicho.meansd.z (m1, m2, sd1, sd2, n1, n2, p1, p2)
  
  return (z)
  
}

#' variance of z converted from dichotomizing/extreme group design mean and se
#'
#' compute the variance of z from t-statistics derived from dichotomizing/extreme group design mean and se\cr
#' see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param se1 standard error of group 1
#' @param se2 standard error of group 2
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param type select one from below, please read \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.}for recommendations\cr
#' there are situations where the "controlled experiment" formulas are recommended even for stats derived from dichotomize/extreme group designs\cr
#' so please refer to his paper\itemize{
#' \item"vd.ce" compute variance of d using the controlled experiment formulas
#' \item"vd.eg" compute variance of d using the extreme group formulas}
#' 
#' @examples
#' dicho.meanse.z(1.1, 0.5, 0.6, 0.2, 10, 18, .25, .25, "vd.ce")
#'
#'@export
dicho.meanse.zvar <- function (m1, m2, se1, se2, n1, n2, p1, p2, type = c("vd.ce", "vd.eg")) {
  sd1 = se.sd(se1, n1)
  sd2 = se.sd(se2, n2)
  
  if (type == "vd.ce") {
    zvar = dicho.meansd.zvar(m1, m2, sd1, sd2, n1, n2, p1, p2, type = "vd.ce")
  } else if (type == "vd.eg") {
    zvar = dicho.meansd.zvar(m1, m2, sd1, sd2, n1, n2, p1, p2, type = "vd.eg")
  }
  
  return (zvar)
}









