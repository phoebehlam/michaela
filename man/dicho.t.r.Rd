% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pustejovsky.R
\name{dicho.t.r}
\alias{dicho.t.r}
\title{t from dichotomizing/extreme group design to r}
\usage{
dicho.t.r(t, totn, n1, n2, p1, p2, dir)
}
\arguments{
\item{t}{t-statistics}

\item{n1}{cell size of group 1}

\item{n2}{cell size of group 2}

\item{p1}{cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)}

\item{p2}{cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)}

\item{dir}{if cohen's d is in absolute value, provide the empirical direction:\cr
+1 for positive associations (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)}
}
\description{
t-statistics from dichotomizing/extreme group design to correlation coefficient\cr
see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
returns extreme group r (referred to as r (subscript eg) in pustejovsky paper) \cr
for unadjusted statistics only
}
\examples{
dicho.t.r(2.05, 1/3, 1/3, 15, 15, -1)

}
