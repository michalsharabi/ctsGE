% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PreparingTheIndexes.R
\name{index}
\alias{index}
\title{Indexing function}
\usage{
index(x, cutoff)
}
\arguments{
\item{x}{A numeric}

\item{cutoff}{A numeric, dermine the threshold for indexing}
}
\value{
Gene expression index
}
\description{
Takes a numeric vector and return an expression index
(i.e., a sequence of 1,-1, and 0)
}
\details{
The function defines limits around the center (median or mean),
+/- cutoff value in median absolute deviation (MAD) or standard deviation
(SD) units respectively.The user defines a parameter cutoff that determines
the limits around the gene-expression center. Then the function calculates
the index value at each time point according to:

\enumerate{
     \item \bold{0:}  standardized value is within the limits (+/- cutoff)
     \item \bold{1:}  standardized value exceeds the upper limit (+ cutoff)
     \item \bold{-1:} standardized value exceeds the lower limit (- cutoff)}
}
\examples{
rawCounts <-
c(103.5, 75.1, 97.3, 27.12, 34.83, 35.53, 40.59, 30.84, 16.39, 29.29)

(sCounts <- scale(rawCounts)[,1])# standardized mean-base scaling

cutoff <- seq(0.2,2,0.1)    # different cutoff produce different indexes

for(i in cutoff){print(index(sCounts,i))}
}
\seealso{
\code{\link{PreparingTheIndexes}}
}

