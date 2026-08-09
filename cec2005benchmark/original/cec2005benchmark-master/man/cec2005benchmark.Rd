\name{cec2005benchmark}
\alias{cec2005benchmark}
\title{Benchmark Functions}
\description{
Common interface to all benchmark functions
}
\usage{
cec2005benchmark(i, x)
}

\arguments{
  \item{i}{Function number between 1 and 25}
  \item{x}{Either a vector with 2, 10, 30 or 50 elements or a matrix with
  the same number of columns and one row for each vector to be evaluated}
}
\value{
A vector with the evaluation of the function for each row of \code{x}
}
\references{
P. N. Suganthan, N. Hansen, J. J. Liang, K. Deb, Y.-P. Chen, A. Auger,
S. Tiwari (2005). Problem Definitions and Evaluation Criteria for the
CEC 2005 Special Session on Real-Parameter Optimization, Technical
report, Nanyang Technological University, Singapore, May 2005,
\url{http://www.ntu.edu.sg/home/EPNSugan/index_files/CEC-05/CEC05.htm}.
}
