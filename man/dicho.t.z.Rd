% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pustejovsky.R
\name{dicho.t.z}
\alias{dicho.t.z}
\title{t from dichotomizing/extreme group design to z}
\usage{
dicho.t.z(t, totn, n1, n2, p1, p2, dir)
}
\arguments{
\item{t}{t-statistics}

\item{totn}{the total n}

\item{n1}{cell size of group 1}

\item{n2}{cell size of group 2}

\item{p1}{cutoff percentile for group 1 (sample-based, please see other formulas for population-based in pustejovsky, 2014)}

\item{p2}{cutoff percentile for group 2 (sample-based, please see other formulas for population-based in pustejovsky, 2014)}

\item{dir}{if t-statistics is in absolute value, provide the empirical direction:\cr
+1 for positive associations (e.g., group 1 - group 2 >= 0) and -1 for negative associations (e.g., group 1 - group 2 < 0)}
}
\description{
t-statistics from dichotomizing/extreme group design to fisher's z\cr
see \href{https://psycnet.apa.org/buy/2013-34335-001}{pustejovsky, 2014. psychological methods.} \cr
returns z based on taylor series approximation, defaults 5th term (use dicho.d.z() if you'd like to adjust to lower term)\cr
for unadjusted statistics only
}
\examples{
dicho.t.z(2.05, 15+15, 15, 15, 1/3, 1/3, -1)

}
