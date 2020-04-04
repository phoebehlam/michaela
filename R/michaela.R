.onAttach <- function(libname,pkgname) {
  packageStartupMessage(' m i c h a e l a  |  beta. please report any bugs.')
}

#' r to z
#'
#' convert correlation coefficient to fisher's z\cr
#' correlation coefficient can be either bivariate pearson's r or partial correlation coefficient (hedges & olkin, 1985)
#'
#' @param r correlation coefficient
#'
#' @examples
#' single value
#' r.z(.3)
#'
#' dplyr example: dataset named "dat" with a column of correlation coefficient values named "corr", create a new var called "z"
#' dat %>% mutate (z = r.z(corr)) -> dat
#'
#' @export
r.z <- function (r) {
  z = 0.5 * log((1 + r)/(1 - r))
  return (z)
}

#' z to r
#'
#' convert fisher's z to correlation coefficient
#'
#' @param z fisher's z
#'
#' @examples
#' single value
#' z.r(.3)
#'
#' dplyr example: dataset named "dat" with a column of correlation coefficient values named "fishz", create a new var called "r"
#' dat %>% mutate (r = z.r(fishz)) -> dat
#'
#' @export
#z to r
z.r <- function (z) {

  r = (exp(2*z) - 1)/(exp(2*z) + 1)

  return (r)
}

#' variance of fisher's z
#'
#' compute the sampling variance of fisher's z
#'
#' @param n sample size from which the fisher's z is derived
#'
#' @examples
#' zvar(100)
#'
#' @export
zvar <- function (n) {
  zvar = 1/(n - 3)
  return (zvar)
}

#' r squared to r
#'
#' compute the semi-partial correlation coefficient from delta r squared
#'
#' @param rsq delta r squared from regression
#' @param dir r squared is in absolute value, so please supply the empirical direction of the link (+1 or -1)
#'
#' @examples
#' # a paper reports delta r squared as .04 and in text describes a negative association
#' rsq.r(.04, -1)
#'
#' @export
rsq.r <- function (rsq, dir) {
  r = dir*sqrt (rsq)
  return (r)
}


#' or to r
#'
#' odd ratio to correlation coefficient: please only use this conversion for "true" odd ratios calculated by
#' or = ad/bc (a = # of exposed cases; b = # of exposed non-cases, c = # of unexposed cases, d = # of unexposed non-cases)
#' for estiamted or derived from exponentiating the coefficient of a logistic regression, use expb.r() instead
#'
#' @param or "true" odd ratio calculated from or = ad/bc (see description above)
#' @param n1 cell size of group 1, if unknown, leave blank-- defaults correction factor a = 4 in estimation of d (assumes n1 = n2)
#' @param n2 cell size of group 2, if unknown, leave blank-- defaults correction factor a = 4 in estimation of d (assumes n1 = n2)
#'
#' @examples
#' or.r(1.25)
#' or.r(1.25, 25, 34)
#'
#' @export
or.r <- function (or, n1, n2) {
  d = log(or)*sqrt(3)/pi

  if (missing(n1)) {
    r = d/sqrt(d^2 + 4)
  } else if (missing(n2)) {
    r = d/sqrt(d^2 + 4)
  } else {
    a = ((n1 + n2)^2)/(n1*n2)
    r = d/sqrt(d^2 + a)
  }

 return (r)
}

#' ci to se
#'
#' compute standard error from confidence interval\cr
#' we recommend checking the computations by using result = "diff". see below.\cr
#' computations are based on t-distribution, not z\cr
#'
#' @param est the point estimate (e.g., a regression coefficient, a mean)
#' @param cil the lower bound of the confidence interval
#' @param ciu the upper bound of the confidence interval
#' @param n the sample size
#' @param k the total number of predictors. for converting the ci associated with a raw mean, enter k = 0.
#' @param result select one of: \itemize{
#' \item"ciu" will compute se based on upper bound of ci (this is the default if unspecified)
#' \item"cil" will compute se based on lower bound of ci
#' \item"avg" will compute se from both ciu and cil and take the average
#' \item"diff" will give the difference between the se computed using ciu and using cil. for transformed point estimate, see sediff() function.
#' }
#'
#' @examples
#' # e.g., regression b=.156, 95CI [.10, .225], n=166, 3 total predictors
#' ci.se(.156, .10, .225, 166, 3)
#'
#' @export
ci.se <- function (est, cil, ciu, n, k, result = c("ciu", "cil", "avg", "diff")) {

  if (missing(k)) {
    stop ("hi, please specify k (the number of predictors); if determining se for raw means, enter k = 0.")
  } #raw means would enter k = 0 that way df would be n - 1

  df = n - k - 1
  qt <- qt (.025, df, lower.tail = FALSE)
  se_ciu = (ciu - est)/qt
  se_cil = (cil - est)/-qt
  diff = abs(se_ciu - se_cil)

  warning ("hi! please also use sediff() or result = 'diff' in this function to check discrepancies using se computed from upper vs. lower ci")
  if (missing (result)) {
    return (se_ciu)
  } else if (result == "ciu"){
    return (se_ciu)
  } else if (result == "cil") {
    return (se_cil)
  } else if (result == "avg"){
    avg = (se_ciu + se_cil)/2
    return (avg)
  } else if (result == "diff") {
    return(diff)
  }
}

#' differences in computed se
#'
#' the differences in the computed standard error using the upper bound vs. lower bound of 95 confidence interval\cr
#' computations are based on t-distribution, not z\cr
#' returns and prints a value/vector of difference(s)
#'
#' @param est the point estimate (e.g., a regression coefficient, a mean)
#' @param cil the lower bound of the 95 confidence interval
#' @param ciu the upper bound of the 95 confidence interval
#' @param n the sample size
#' @param k the total number of predictors. for converting the ci associated with a raw mean, enter k = 0.
#' @param type select one of:\itemize{
#' \item"reg" for regular coefficient/mean \cr
#' \item"exp b" for exponentiated regression coefficients (e.g., estimated odd ratios from exponentiated logistic regression coefficients) and cis \cr
#' \item"geo" for geometric means (exponentiated arithmetic means of logged raw units) and cis \cr
#' \item"percent" for percent changes (exponentiated regression coefficients expressed in percent changes) and cis
#'}
#'
#' @examples
#' # e.g., exponentiated logistic regression coefficient (estimated odd ratios) and cis
#' se.diff(4.05, 1.56, 10.51, 141, 8, "exp b")
#'
#' @export
sediff <- function (est, cil, ciu, n, k, type = c("reg", "exp b", "geo", "percent")) {
  if (missing(k)) {
    stop ("hi, please specify k (the number of predictors); if determining se for raw means enter 0.")
  }

  if (missing(type)) {
    stop ('hi, please specify type-- options are "reg" (regression), "exp b" (exponentiated regression b), "geo" (geometric means), or "percent"(percent change from transformed exponentiated b)')
  }

  df = n - k - 1
  qt <- qt (.025, df, lower.tail = FALSE)

  if (type == "reg") {
    se_ciu = (ciu - est)/qt
    se_cil = (cil - est)/-qt
  }
  else if (type == "exp b") {
    lnb = log (est)
    lncil = log (cil)
    lnciu = log (ciu)

    se_ciu = (lnciu - lnb)/qt
    se_cil = (lncil - lnb)/-qt
  }
  else if (type == "geo") {
    lnb = log (est)
    lncil = log (cil)
    lnciu = log (ciu)

    se_ciu = (lnciu - lnb)/qt
    se_cil = (lncil - lnb)/-qt
  }
  else if (type == "percent") {
    b = log (est/100 + 1)
    cil2 = log (cil/100 + 1)
    ciu2 = log (ciu/100 + 1)

    se_ciu = (ciu2 - b)/qt
    se_cil = (cil2 - b)/-qt
  }
  diff = abs(se_ciu - se_cil)
  return (diff)
  print (diff)
}

#' b and se to t
#'
#' regression coefficient and standard error to t-statistics
#'
#' @param b the regression coefficient (non-transformed only, see other functions for transformed regression coefficients)
#' @param se the standard error
#'
#' @examples
#' bse.t(.25, .09)
#' dat %>% mutate (t_from_b_se = bse.t(reg_coef, reg_se)) -> dat
#'
#' @export
bse.t <- function (b, se) {
  t = b/se
  return (t)
}

#' t to r
#'
#' t-statistics to correlation coefficient\cr
#' unadjusted t-statistics estimate bivariate correlation coefficient (enter k = 1)\cr
#' adjusted t-statistics (e.g., from a regression) estimate partial correlation coefficient \href{https://www.journals.uchicago.edu/doi/abs/10.5243/jsswr.2013.24}{(aloe & thompson, 2013)}.
#'
#' @param t the t-statistics
#' @param n the sample size
#' @param k the total number of predictors, enter k = 1 if bivariate
#' @param dir empirical direction: +1 or -1. articles sometimes report t-statistics in absolute value, may be prudent to check when extracting stats\cr
#' if the reported t-statistics is positive, read text to ensure it reflects a positive association, enter the t as is and enter dir = 1 if so\cr
#' if the reported t-statistics is positive, but text suggests a negative association, enter the t as is and enter dir = -1\cr
#' if the reported t-statistics is negative, you can simply enter the negative t-value and enter dir = 1\cr
#' note: if your t-statistics is computed from taking the ratio of b and se, e.g., using bse.t(), then it should already have empirical direction "built-in" since b is directional itself, so simply enter dir = 1
#'
#' @examples
#' t.r(2.14, 100, -1)
#' dat %>% mutate (r_from_t = t.r(t_stat, analytical_n, direction)) -> dat
#'
#' @export
t.r <- function (t, n, k, dir) {
  if (missing(k)) {
    stop ("hi, please specify k (the number of predictors); if t-test of two groups with no covariates, enter k = 1")
  } #k = 1 here instead of 0 for no covariates because df for 2-sampled t-test is (n1 - 1) + (n2 - 1) = total n - 2, so k - 1 will give n - 1 - 1
  if (missing (dir)) {
    r = sign(t)*sqrt(t^2/(t^2 + (n - k - 1)))
    } else {
      r = sign(dir)*sqrt(t^2/(t^2 + (n - k - 1)))
      }
  return (r)
}

#' b and se to r
#'
#' regression coefficient and standard error to correlation coefficient\cr
#' unadjusted b and se estimate bivariate correlation coefficient (enter k = 1)\cr
#' adjusted b and se estimate partial correlation coefficient \href{https://www.journals.uchicago.edu/doi/abs/10.5243/jsswr.2013.24}{(aloe & thompson, 2013)}.
#'
#' @param b the regression coefficient (untransformed only, see other functions e.g., expb.r() for transformed coefficients)
#' @param se the standard error
#' @param n the sample size
#' @param k the total number of predictors, enter k = 1 if bivariate
#'
#' @examples
#' bse.r(.25, .09, 100, 3)
#' dat %>% mutate (r_from_bse = bse.r(reg_coef, reg_se, reg_n, reg_predictors)) -> dat
#'
#' @export
bse.r <- function (b, se, n, k) {
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors)')
  }
  t = b/se
  r = sign(b)*sqrt(t^2/(t^2 + (n - k - 1)))
  return (r)
}

#' b and ci to r
#'
#' regression coefficient and confidence interval to correlation coefficient: b and ci -> t -> r\cr
#' unadjusted t-statistics estimate bivariate correlation coefficient (enter k = 1)\cr
#' adjusted t-statistics (e.g., from a regression) estimate partial correlation coefficient \href{https://www.journals.uchicago.edu/doi/abs/10.5243/jsswr.2013.24}{(aloe & thompson, 2013)}.
#'
#' @param b the regression coefficient (untransformed only, see other functions e.g., expb.r() for transformed coefficients)
#' @param cil the lower bound of the 95 confidence interval
#' @param ciu the upper bound of the 95 confidence interval
#' @param n the sample size
#' @param k the total number of predictors, enter k = 1 if bivariate
#'
#' @examples
#' bci.r(.156, .10, .225, 166, 3)
#' dat %>% mutate (r_from_bci = bci.r(reg_coef, reg_lowerci, reg_upperci, reg_n, reg_predictors)) -> dat
#'
#' @export
bci.r <- function (b, cil, ciu, n, k, result = c("cil", "ciu", "avg")) {
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors)')
  }

  if (missing(result)) {
    se = ci.se (b, cil, ciu, n, k)
  } else if (result == "ciu") {
    se = ci.se (b, cil, ciu, n, k, result = "ciu")
  } else if (result == "cil") {
    se = ci.se (b, cil, ciu, n, k, result= "cil")
  } else if (result == "avg")  {
    se = ci.se (b, cil, ciu, n, k, result= "avg")
  }

  r = bse.r (b = b, se = se, n = n, k = k)
  return (r)
}

#' exponentiated b and ci to r
#'
#' exponentiated logistic regression coefficient and confidence interval to correlation coefficient\cr
#' expnentiated b and ci -> t -> r\cr
#' unadjusted t-statistics estimate bivariate correlation coefficient (enter k = 1)\cr
#' adjusted t-statistics (e.g., from a regression) estimate partial correlation coefficient \href{https://www.journals.uchicago.edu/doi/abs/10.5243/jsswr.2013.24}{(aloe & thompson, 2013)}.
#'
#' @param b exponentiated logistic regression coefficient
#' @param cil the lower bound of the 95 confidence interval
#' @param ciu the upper bound of the 95 confidence interval
#' @param n the sample size
#' @param k the total number of predictors, enter k = 1 if bivariate
#'
#' @examples
#' expb.r(4.05, 1.56, 10.51, 141, 8)
#' dat %>% mutate (r_from_expbci = expb.r(expb_coef, expb_lowerci, expb_upperci, expb_n, expb_predictors)) -> dat
#'
#' @export
expb.r <- function (b, cil, ciu, n, k) {
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors)')
  }

  lnb = log (b)
  lncil = log (cil)
  lnciu = log (ciu)

  se = ci.se (lnb, lncil, lnciu, n, k)
  r = bse.r (lnb, se, n, k)

  return (r)
}

#' risk ratio to r
#'
#' risk ratio to correlation coefficient
#' "true" risk ratios calculated from (a/(a+b))/(c/(c+d))\cr
#' a = # of exposed cases\cr
#' b = # of exposed non-cases\cr
#' c = # of non-exposed cases\cr
#' d = # of non-exposed non-cases
#'
#' @param rr the risk ratio
#' @param prev the prevalence of the outcome in the reference group: c/(c+d)\cr
#' c = # of non-exposed cases\cr
#' d = # of non-exposed non-cases
#' @param n1 the cell size of group 1 (if unknown, leave blank, assumes correcting factor a = 4 where n1 = n2)
#' @param n2 the cell size of group 2 (if unknown, leave blank, assumes correcting factor a = 4 where n1 = n2)
#'
#' @examples
#' rr.r(4.05, 1.56, 10.51, 141, 8)
#' dat %>% mutate (r_from_rr = rr.r(rr, rr_lowerci, rr_upperci, rr_n, rr_predictors)) -> dat
#'
#' @export
rr.r <- function (rr, prev, n1, n2) {
  or = (rr*(1-prev))/(1-rr*prev)
  r = or.r(or, n1, n2)
  return (r)
}


#' percent change to r
#'
#' percent change derived from transformed exponentiated regression coefficients and confidence interval to correlation coefficient\cr
#' percent change and ci -> expnentiated b and ci -> t -> r\cr
#' unadjusted t-statistics estimate bivariate correlation coefficient (enter k = 1)\cr
#' adjusted t-statistics (e.g., from a regression) estimate partial correlation coefficient \href{https://www.journals.uchicago.edu/doi/abs/10.5243/jsswr.2013.24}{(aloe & thompson, 2013)}.
#'
#' @param percent the percent change calculated from exponentiated logistic regression coefficient \cr
#' specifically computed such that percent change = 100*(exp(b) - 1); where b is the regression coefficient \cr
#' if percent change expressed in decimals, please *100 first \cr
#' or if percent change not computed with above equation, please compute back to exponentiated b and then use expb.r()
#' @param cil the lower bound of the confidence interval around the percent change
#' @param ciu the upper bound of the confidence interval around the percent change
#' @param n the sample size
#' @param k the total number of predictors, enter k = 1 if bivariate
#' @param result one of the following:\itemize{
#' \item"ciu" uses the upper bound of ci to compute se (defaults if missing input)
#' \item"cil" uses the lower bound of ci to compute se
#' \item"avg" uses both upper bound and lower bound and then take the average of the resulting se's
#' }
#'
#' @examples
#' percent.r(7, 3, 11, 3134, 6)
#' dat %>% mutate (r_from_perc = percent.r(percent, percent_lowerci, percent_upperci, percent_n, percent_predictors)) -> dat
#'
#' @export
percent.r <- function (percent, cil, ciu, n, k, result = c("cil", "ciu", "avg")) {
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors)')
  }

  b = log (percent/100 + 1)
  cil2 = log (cil/100 + 1)
  ciu2 = log (ciu/100 + 1)
  se_ciu = suppressWarnings(ci.se (b, cil2, ciu2, n, k, result = "ciu"))
  se_cil = suppressWarnings(ci.se (b, cil2, ciu2, n, k, result = "cil"))
  se_avg = suppressWarnings(ci.se (b, cil2, ciu2, n, k, result = "avg"))

  if(missing(result)) {
    r = bse.r (b, se_ciu, n, k)
  } else if(result == "ciu") {
    r = bse.r (b, se_ciu, n, k)
  } else if(result == "cil") {
    r = bse.r (b, se_cil, n, k)
  } else if(result == "avg") {
    r = bse.r (b, se_avg, n, k)
  }
  warning ("hi! please also use sediff() to check discrepancies for se computed from upper vs. lower ci")

  return (r)
}

#' f to r
#'
#' f-statistics (df1=1) to correlation coefficient\cr
#' unadjusted statistics only
#'
#' @param f the f-statistics
#' @param df the degrees of freedom2
#' @param dir empirical direction: +1 or -1. articles sometimes report f-statistics in absolute value, may be prudent to check when extracting stats\cr
#' if the reported f-statistics is positive, read text to ensure it reflects a positive association, enter the f as is and enter dir = 1 if so\cr
#' if the reported f-statistics is positive, but text suggests a negative association, enter the f as is and enter dir = -1\cr
#' if the reported f-statistics is negative, you can simply enter the negative f-value and enter dir = 1\cr
#'
#' @examples
#' f.r(3.21, 97, -1)
#' dat %>% mutate (r_from_f = f.r(reg_f, reg_df, reg_emp_direction)) -> dat
#'
#' @export
f.r <- function (f, df, dir) {
  if (missing (dir)) {
    r = NA_real_
  } else {
    r = sign(dir)*sqrt(f/(f + df))
  }
  return (r)
}

# code this in one day
## dichotomized groups: ancova f to r (for adjusted)
# looking into compute.es, a.fes function
# requires covariate outcome correlation or multiple correlation
# then, d <- sqrt(f * (n.1 + n.2)/(n.1 * n.2)) * sqrt(1 - R^2)

#' p to r
#'
#' regression p-values (or other t distribution-based p-values) to correlation coefficient\cr
#' estimates t-statistics and then correlation coefficients\cr
#'
#' @param p the p-value
#' @param n the sample size
#' @param k the number of predictors
#' @param dir empirical direction: +1 for positive associaitons and -1 for negative associations
#'
#' @examples
#' regp.r(.023, 100, 3, 1)
#' dat %>% mutate (r_from_pval = regp.r(reg_pval, reg_n, reg_predictors, reg_emp_direction)) -> dat
#'
#' @export
regp.r <- function (p, n, k, dir) {
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors), for 2-sampled t-test without covariates, enter k=1')
  }
  if (missing(dir)) {
    r = NA_real_
  }
  else {
    t = qt (p/2, n - k - 1, lower.tail= FALSE)
    r = t.r (t, n, k, dir)
  }
  return (r)
}

#' mean and sd to d
#'
#' means and standard deviations of two groups to cohen's d
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#'
#' @examples
#' meansd.d(15.1, 13.2, 2.3, 1.1, 25, 45)
#' dat %>% mutate (d_from_meansd = meansd.d(mean1, mean2, sd1, sd2, n1, n2)) -> dat
#'
#' @export
meansd.d <- function(m1, m2, sd1, sd2, n1, n2) {
  d = (m1-m2)/sqrt ((sd1^2+sd2^2)/2)
  return(d)
}

#' d to t
#'
#' cohen's d to t-statistics
#'
#' @param d cohen's d
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param k total number of predictors
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#' +1 for positive associaitons and -1 for negative associations
#'
#' @examples
#' d.t(.80, 15, 20, 1, -1)
#' dat %>% mutate (t_from_d = d.t(cohend, stress_n, control_n, numpred, emp_dir)) -> dat
#'
#' @export
d.t <- function (d, n1, n2, k, dir) {
  if(missing(dir)){
    t = NA_real_
  } else {
    t = dir*d*sqrt((n1*n2*(n1+n2-k-1)))/(n1+n2)
  }
  return (t)
}

#' se to sd
#'
#' standard error of raw mean to standard deviation
#'
#' @param se standard error
#' @param n cell size
#'
#' @examples
#' se.sd(0.04, 15)
#'
#' @export
se.sd <- function (se, n) {
  sd = se*sqrt(n)

  return (sd)
}

#' mean and se to d
#'
#' means and standard errors of two groups to cohen's d
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param se1 standard deviation of group 1
#' @param se2 standard deviation of group 2
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#'
#' @examples
#' meanse.d(15.1, 13.2, .64, .164, 25, 45)
#'
#' @export
meanse.d <- function (m1, m2, se1, se2, n1, n2) {
  sd1 = se.sd(se1, n1)
  sd2 = se.sd(se2, n2)
  d = meansd.d(m1, m2, sd1, sd2, n1, n2)
  return (d)
}

#' d to r
#'
#' cohen's d to correlation coefficient
#'
#' @param d cohen's d
#' @param n1 cell size of group 1 (if missing, defaults correction factor a = 4, assumes n1=n2)
#' @param n2 cell size of group 2 (if missing, defaults correction factor a = 4, assumes n1=n2)
#' @param k total number of predictors
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#'  +1 for positive associaitons (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
#'
#' @examples
#' d.t(.80, 15, 20, 1, -1)
#' dat %>% mutate (t_from_d = d.t(cohend, stress_n, control_n, numpred, emp_dir)) -> dat
#'
#' @export
d.r <- function (d, n1, n2, dir) {
  if (missing (n1)) {
    a = 4
  } else if (missing (dir)){
    r = NA_real_
  } else {
    a = ((n1 + n2)^2)/(n1*n2)
  }

  r = sign(dir)*d/sqrt (d^2 + a)

  return  (r)
}

#' mean and sd to r
#'
#' means and standard deviations of two groups to correlation coefficient
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param sd1 standard deviation of group 1
#' @param sd2 standard deviation of group 2
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#'
#' @examples
#' meansd.r(15.1, 13.2, 2.3, 1.1, 25, 45)
#' dat %>% mutate (r_from_meansd = meansd.r(mean1, mean2, sd1, sd2, n1, n2)) -> dat
#'
#' @export
meansd.r <- function (m1, m2, sd1, sd2, n1, n2) {
  d = meansd.d(m1, m2, sd1, sd2, n1, n2)
  t = d.t (d, n1, n2, 1, 1) #direction always one, because the d already comes with sign
  r = t.r(t, n1+n2, 1, 1) #direction always one, because the d already comes with sign
  return (r)
}

#' mean and se to r
#'
#' means and standard error of two groups to correlation coefficient
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param se1 standard error of group 1
#' @param se2 standard error of group 2
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#'
#' @examples
#' meanse.r(15.1, 13.2, .64, .164, 25, 45)
#'
#' @export
meanse.r <- function (m1, m2, se1, se2, n1, n2, k) {
  sd1 = se.sd(se1, n1)
  sd2 = se.sd(se2, n2)
  r = meansd.r (m1, m2, sd1, sd2, n1, n2)

  return (r)
}

#' mean and ci to d
#'
#' means and confidence interval of two groups to cohen's d
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param cil lower bound of 95 confidence interval
#' @param ciu upper bound of 95 confidence interval
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param k total number of predictors (k=0 for raw)
#'
#'
#' @export
meanci.d <- function (m1, m2, cil1, ciu1, cil2, ciu2, n1, n2, k) {
  se1_ciu = suppressWarnings(ci.se(m1, cil1, ciu1, n1, k, result = "ciu"))
  se2_ciu = suppressWarnings(ci.se(m2, cil2, ciu2, n2, k, result = "ciu"))

  d = meansd.d(m1, m2, se.sd(se1_ciu, n1), se.sd(se2_ciu, n2))

  return (d)

}

#' mean and ci to r
#'
#' means and confidence interval of two groups to correlation coefficient
#'
#' @param m1 mean of group 1
#' @param m2 mean of group 2
#' @param cil lower bound of 95 confidence interval
#' @param ciu upper bound of 95 confidence interval
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param k total number of predictors (k=0 for raw)
#'
#'
#' @export
meanci.r <- function (m1, m2, cil1, ciu1, cil2, ciu2, n1, n2, k, result = c("sediff", "cil", "ciu")){
  if (missing(k)) {
    stop ('hi, please specify k (the number of predictors); for raw means without covariates, enter k = 0')
  }#k = 0 would give df = n - 1 for each mean

  se1_ciu = suppressWarnings(ci.se(m1, cil1, ciu1, n1, k, result = "ciu"))
  se1_cil = suppressWarnings(ci.se(m1, cil1, ciu1, n1, k, result = "cil"))

  se2_ciu = suppressWarnings(ci.se(m2, cil2, ciu2, n2, k, result = "ciu"))
  se2_cil = suppressWarnings(ci.se(m2, cil2, ciu2, n2, k, result = "cil"))

  if (missing (result)) {

    r = meanse.r(m1, m2, se1_ciu, se2_ciu, n1, n2)

    warning ('hi! please use result = "sediff" with this same function to check se discreapncies')
    return (r)

  } else if (result == "ciu") {
    r = meanse.r(m1, m2, se1_ciu, se2_ciu, n1, n2)

    warning ('hi! please use result = "sediff" with this same function to check se discreapncies')
    return (r)
  } else if (result == "cil") {
    r = meanse.r(m1, m2, se1_cil, se2_cil ,n1, n2)

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



#' median, minimum, and max to mean
#'
#' estimate mean from median, minimum and max. see \href{https://onlinelibrary.wiley.com/doi/full/10.4073/cmpn.2016.3}{polanin & snilstveit, 2016}
#'
#' @param median the median
#' @param min the minimum
#' @param max the maximum
#'
#' @export
median.mean <- function (median, min, max) {
  mean = (min + (2*median) + max)/4
  return(mean)
}

#' median, minimum, and max to sd
#'
#' estimate sd from median, minimum and max. see \href{https://onlinelibrary.wiley.com/doi/full/10.4073/cmpn.2016.3}{polanin & snilstveit, 2016}
#'
#' @param median the median
#' @param min the minimum
#' @param max the maximum
#'
#' @export
median.sd <- function (median, min, max) {

  sd = sqrt(1/12*((min - (2*median) + max)^2)/4 + (max- min)^2)

  return(sd)
}

#' median, minimum, and max to r
#'
#' estimate correlation coefficient from median, minimum and max. see \href{https://onlinelibrary.wiley.com/doi/full/10.4073/cmpn.2016.3}{polanin & snilstveit, 2016}
#'
#' @param median1 the median for group 1
#' @param median2 the median for group 2
#' @param min1 the minimum for group 1
#' @param min2 the minimum for group 2
#' @param max1 the maximum for group 1
#' @param max2 the maximum for group 2
#' @param n1 the cell size for group 1
#' @param n2 the cell size for group 2
#'
#' @export
median.r <- function (median1, median2, min1, max1, min2, max2, n1, n2) {

  mean1 = median.mean(median1, min1, max1)
  mean2 = median.mean(median2, min2, max2)
  sd1 = median.sd(median1, min1, max1)
  sd2 = median.sd(median2, min2, max2)

  r = meansd.r(mean1, mean2, sd1, sd2, n1, n2)

  return(r)
}

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
#' +1 for positive associaitons (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
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
#' @param k kth term for taylor series approximation (1 to 5)
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#' +1 for positive associaitons (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
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

  #taylor series approximation, term = 5
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
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param dir if cohen's d is in absolute value, provide the empirical direction:\cr
#' +1 for positive associaitons (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
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
#' for taylor series approximation, defaults 5th term
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
#' +1 for positive associaitons (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
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
#' +1 for positive associaitons (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
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
#' returns z based on taylor series approximation, defaults 5th term \cr
#' for unadjusted statistics only
#'
#' @param t t-statistics
#' @param totn the total n
#' @param n1 cell size of group 1
#' @param n2 cell size of group 2
#' @param p1 cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param p2 cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)
#' @param dir if t-statistics is in absolute value, provide the empirical direction:\cr
#' +1 for positive associaitons (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)
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
#' for taylor series approximation, defaults 5th term
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




