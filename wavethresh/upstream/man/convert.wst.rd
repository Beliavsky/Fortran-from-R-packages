\name{convert.wst}
\alias{convert.wst}
\title{Convert a non-decimated wst object into a wd object. }
\description{
Convert a packed-ordered non-decimated wavelet transform object into a time-ordered non-decimated wavelet transform object.}
\usage{
\method{convert}{wst}(wst, \dots) 
}
\arguments{
\item{wst}{The \code{\link{wst}} class object that you wish to convert.}
\item{\dots}{any other arguments}
}
\details{
In WaveThresh3 a non-decimated wavelet transform can be ordered in two different ways: as a time-ordered or packet-ordered representation. The coefficients in the two objects are \emph{exactly the same} it is just their internal representation and ordering which is different. The two different representations are useful in different situations. The packet-ordering is useful for curve estimation applications and the time-ordering is useful for time series applications. 

See Nason, Sapatinas and Sawczenko, 1998 for further details on ordering and weaving. 

Note that the input object must be of the non-decimated type. In other words the type component of the input object must be "\code{station}". 
Once the input object has been converted the output can be used with any of the functions suitable for the \code{\link{wd.object}}.
 
The actual weaving permutation for shuffling coefficients from one representation to another is achieved by the \code{\link{getarrvec}} function.  
}
\value{
An object of class \code{\link{wd}} containing exactly the same information as the input object but ordered differently as a packet-ordered object. 
}
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1997 }
\seealso{
\code{\link{convert}}, \code{\link{getarrvec}}, \code{\link{levarr}}, \code{\link{wd}}, \code{\link{wd.object}}, \code{\link{wst}}, \code{\link{wst.object}}. 
}
\examples{
#
# Generate a sequence of 32 random normals (say) and take their
# \code{packed-ordered non-decimated wavelet transform}
#
myrand <- wst(rnorm(32))
#
# Print out the result (to verify the class and type of the object)
#
#myrand
#Class 'wst' : Stationary Wavelet Transform Object:
#       ~~~  : List with 8 components with names
#              wp Carray nlevelsWT filter date 
#
#$WP and $Carray are the coefficient matrices

#
#Created on : Tue Sep 29 12:29:45 1998 
#
#summary(.):
#----------
#Levels:  5 
#Length of original:  32 
#Filter was:  Daub cmpct on least asymm N=10 
#Boundary handling:  periodic 

#Date:  Tue Sep 29 12:29:45 1998 
#
# Yep, the myrand object is of class: \code{\link{wst.object}}.
#
# Now let's convert it to class \code{\link{wd}}. The object
# gets returned and, as usual in S, is printed.
#
convert(myrand)
#Class 'wd' : Discrete Wavelet Transform Object:
#       ~~  : List with 8 components with names
#              C D nlevelsWT fl.dbase filter type bc date 
#
#$ C and $ D are LONG coefficient vectors !
#
#Created on : Tue Sep 29 12:29:45 1998 
#Type of decomposition:  station 
#
#summary(.):
#----------
#Levels:  5 
#Length of original:  32 
#Filter was:  Daub cmpct on least asymm N=10 
#Boundary handling:  periodic 
#Transform type:  station 
#Date:  Tue Sep 29 12:29:45 1998 
#
# The returned object is of class \code{\link{wd}} with a
# type of "station".
# I.e. it has been converted successfully.
}
\keyword{manip}
\author{G P Nason}
