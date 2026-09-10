\name{TOthreshda1}
\alias{TOthreshda1}
\title{Data analytic wavelet thresholding routine}
\usage{
TOthreshda1(ywd, alpha = 0.05, verbose = FALSE, return.threshold = FALSE)
}
\arguments{
\item{ywd}{The \code{\link{wd.object}} that you wish to threshold.}
\item{alpha}{The smoothing parameter which is a p-value }
\item{verbose}{Whether messages get printed}
\item{return.threshold}{If TRUE then the threshold value gets returned
	rather than the actual thresholded object}
}
\description{
This function might be better called using the regular
\code{\link{threshold}} function using the \code{op1} policy.

Corresponds to the wavelet thresholding routine developed by
Ogden and Parzen (1994) Data dependent wavelet thresholding in nonparametric
regression with change-point applications. \emph{Tech Rep 176},
University of South Carolina, Department of Statistics.
}
\details{
The TOthreshda1 method operates by testing the max of each set of squared       
   wavelet coefficients to see if it behaves as the nth order statistic of   
   a set of independent chi^2(1) r.v.'s.  If not, it is removed, and the    
   max of the remaining subset is tested, continuing in this fashion until   
   the max of the subset is judged not to be significant.  

   In this situation, the level of the hypothesis tests, alpha, has default  
   value 0.05.  Note that the choice of alpha controls the smoothness of     
   the resulting wavelet estimator -- in general, a relatively large alpha   
   makes it easier to include coefficients, resulting in a more wiggly       
   estimate; a smaller alpha will make it more difficult to include          
   coefficients, yielding smoother estimates.         
}
\value{
Returns the threshold value if \code{return.threshold==TRUE} otherwise
returns the shrunk set of wavelet coefficients.
}
\seealso{\code{\link{threshold}},\code{\link{TOthreshda2}}, \code{\link{wd}}}
\author{Todd Ogden}
\keyword{smooth}
