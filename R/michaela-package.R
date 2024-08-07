#' m  i  c  h  a  e  l  a
#'
#' @section description:
#' michaela is beta. please report any bugs\cr
#' a package for converting various stats metrics to correlation coefficient\cr
#' simple functions aimed to be user-friendly for folks without r knowledge\cr
#' all functions are in the format of "from.to" \cr
#' e.g., \emph{from} regression b and se \emph{to} r: bse.r()\cr
#' e.g., \emph{from} odd ratio \emph{to} r: or.r()\cr
#' \cr
#' all functions are documented, use ?function to see usage, arguments, and examples (e.g., ?bse.r())\cr
#'
#' \emph{\strong{a note on converting unadjusted vs. adjusted metrics:}}\cr
#' given bivariate metrics, all michaela functions can estimate bivariate pearson's correlation\cr
#' given adjusted metrics, \emph{some} michaela functions can estimate the partial correlation coefficient\cr
#' please use with caution: if you extracted adjusted metrics, make sure the function can actually estimate partial correlation\cr
#' for your convenience, equations that can take adjusted metrics are marked with *. \cr
#' \cr
#' see references section below for a list of aloe and becker's work on synthesis of partial effect size and inaccuracies in replacing bivariate corr with regression results
#'
#' @section correlational designs:
#' *r.z(): r to fisher's z \cr
#' *z.r(): fisher's z to r \cr
#' *zvar(): variance of z \cr
#' *rsq.r(): r squared to r (estimates semi-partial correlation instead of partial correlation)\cr
#' *bse.t(): regression coefficient and standard error to t-statistics \cr
#' *t.r(): t-statistics to r \cr
#' *bse.r(): regression coefficient and standard error to r \cr
#' *bci.r(): regression coefficient and confidence interval to r \cr
#' *expb.r(): exponentiated logistic regression coefficient and confidence interval to r \cr
#' *percent.r(): percent change derived from transformed exponentiated regression coefficients and confidence interval to r \cr
#' f.r(): f-statistics (df1=1) to r \cr
#' *regp.r(): regression p-values (or other t distribution-based p-values) to r
#'
#' @section group difference designs:
#' or.r(): odd ratio to r \cr
#' rr.r(): risk ratio to r \cr
#' meansd.d(): means and standard deviations of two groups to cohen's d\cr
#' meanse.d(): means and standard errors of two groups to cohen's d\cr
#' meanci.d(): means and confidence intervals of two groups to cohen's d\cr
#' meansd.r(): means and standard deviations of two groups to r\cr
#' meanse.r(): means and standard errors of two groups to r\cr
#' meanci.r(): means and confidence intervals of two groups to r\cr
#' geo.d(): geometric means and cis to cohen's d\cr
#' geo.r(): geometric means and cis to r\cr
#' d.t(): cohen's d to t-statistics\cr
#' d.r(): cohen's d to r \cr
#' median.r(): medians, minimums, and maximums to r
#'
#' @section dichotomized/extreme group designs:
#' use with caution, please refer and follow recommendations by \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky (2014)}.\cr
#' dicho.d.r(): cohen's d from dichotomizing/extreme group design to r\cr
#' dicho.d.z(): cohen's d from dichotomizing/extreme group design to fisher's z\cr
#' dichod.zvar(): variance of z from a cohen's d derived from dichotomizing/extreme group design\cr
#' dicho.meansd.r(): means and sds from dichotomizing/extreme group design to r\cr
#' dicho.meansd.z(): means and sd from dichotomizing/extreme group design to fisher's z\cr
#' dicho.meansd.zvar(): variance of z converted from dichotomizing/extreme group design mean and sd\cr
#' dicho.t.r(): t-statistics from dichotomizing/extreme group design to r\cr
#' dicho.t.z(): t-statistics from dichotomizing/extreme group design to fisher's z\cr
#' dicho.t.zvar(): variance of z derived from dichotomizing/extreme group design t-statistics\cr
#' dicho.se.r(): mean and se from dichotomizing/extreme group design to r\cr
#' dicho.se.z(): mean and se from dichotomizing/extreme group design to fisher's z\cr
#' dicho.se.zvar(): variance of z derived from dichotomizing/extreme group design mean and se\cr
#' dicho.meanci.r(): mean and ci from dichotomizing/extreme group design to r\cr
#' dicho.meanci.z(): mean and ci from dichotomizing/extreme group design to fisher's z\cr
#' dicho.meanci.zvar(): variance of z derived from dichotomizing/extreme group design mean and ci\cr
#'
#' @section auxillary conversions:
#' *ci.se(): confidence interval to standard error \cr
#' *sediff(): difference between standard errors computed using upper ci and lower ci \cr
#' *geoci.se(): confidence interval around geometric mean to standard error of arithmetic mean of logged raw units\cr
#' se.sd(): standard deviation to standard error\cr
#' median.mean(): median, minimum, and maximum to estimated mean\cr
#' median.sd(): median, minimum, and maximum to estimated standard deviation
#'
#' @section references:
#' re conversion of adjusted effect sizes and combining unadjusted and adjusted effect sizes, please refer and follow recommendations below\enumerate{
#'\item\href{https://www.journals.uchicago.edu/doi/abs/10.5243/jsswr.2013.24}{aloe, a. m., & thompson, c. g. (2013). the synthesis of partial effect sizes. journal of the society for social work and research, 4(4), 390-405.}
#'\item\href{https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1126}{aloe, a. m. (2015). inaccuracy of regression results in replacing bivariate correlations. research synthesis methods, 6(1), 21-27.}
#'\item\href{https://www.semanticscholar.org/paper/Synthesizing-bivariate-and-partial-effect-sizes-Aloe-Tanner-Smith/aee589e75229f314cad5ba8647c9756814ca9e97}{aloe, a. m., tanner‐smith, e. e., becker, b. j., & wilson, d. b. (2016). synthesizing bivariate and partial effect sizes. campbell systematic reviews, 12(1), 1-9.}
#'\item\href{https://projecteuclid.org/euclid.ss/1199285041}{becker, b. j., & wu, m. j. (2007). the synthesis of regression slopes in meta-analysis. statistical science, 22(3), 414-429.}
#'}
#' @docType package
#' @name michaela
NULL
