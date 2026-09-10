\name{wd}
\alias{wd}
\title{Wavelet transform (decomposition).}
\description{
This function can perform two types of discrete wavelet transform (DWT). The standard DWT computes the DWT according to Mallat's pyramidal algorithm (Mallat, 1989) (it also has the ability to compute the \emph{wavelets on the interval} transform of Cohen, Daubechies and Vial, 1993).
 
The non-decimated DWT (NDWT) contains all possible shifted versions of the DWT. The order of computation of the DWT is O(n), and it is O(n log n) for the NDWT if n is the number of data points. 
}
\usage{
wd(data, filter.number=10, family="DaubLeAsymm", type="wavelet",
bc="periodic", verbose=FALSE, min.scale=0, precond=TRUE)
}
\arguments{
\item{data}{A vector containing the data you wish to decompose. The length of this vector must be a power of 2.}
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments. 

For the ``wavelets on the interval'' (\code{bc="interval"}) transform the filter number ranges from 1 to 8. See the table of filter coefficients indexed after the reference to Cohen, Daubechies and Vial, 1993.}
\item{family}{specifies the family of wavelets that you want to use. Two popular options are "DaubExPhase" and "DaubLeAsymm" but see the help for \code{\link{filter.select}} for more possibilities.
 
This argument is ignored for the ``wavelets on the interval'' transform (\code{bc="interval"}).

Note that, as of version 4.6.1 you can use the Lina-Mayrand complex-valued
wavelets.
} 
\item{type}{specifies the type of wavelet transform. This can be "wavelet" (default) in which case the standard DWT is performed (as in previous releases of WaveThresh). If type is "station" then the non-decimated DWT is performed. At present, only periodic boundary conditions can be used with the non-decimated wavelet transform.}
\item{bc}{specifies the boundary handling. If \code{bc="periodic"} the default, then the function you decompose is assumed to be periodic on it's interval of definition, if \code{bc="symmetric"} then the function beyond its boundaries is assumed to be a symmetric reflection of the function in the boundary. The symmetric option was the implicit default in releases prior to 2.2. If \code{bc=="interval"} then the ``wavelets on the interval algorithm'' due to Cohen, Daubechies and Vial is used. (The \code{WaveThresh} implementation of the ``wavelets on the interval transform'' was coded by Piotr Fryzlewicz, Department of Mathematics, Wroclaw University of Technology, Poland; this code was largely based on code written by Markus Monnerjahn, RHRK, Universitat Kaiserslautern; integration into \code{WaveThresh} by \code{GPN}. See the nice project report by Piotr on this piece of code). }
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
\item{min.scale}{Only used for the ``wavelets on the interval transform''. The wavelet algorithm starts with fine scale data and iteratively coarsens it. This argument controls how many times this iterative procedure is applied by specifying at which scale level to stop decomposiing. }
\item{precond }{Only used for the ``wavelets on the interval transform''. This argument specifies whether preconditioning is applied (called prefiltering in Cohen, Daubechies and Vial, 1993.) Preconditioning ensures that sequences like 1,1,1,1 or 1,2,3,4 map to zero high pass coefficients. }
}
\details{
If type=="wavelet" then the code implements Mallat's pyramid algorithm (Mallat 1989). For more details of this implementation see Nason and Silverman, 1994. Essentially it works like this: you start off with some data cm, which is a real vector of length \eqn{2^m}, say. 

Then from this you obtain two vectors of length \eqn{2^(m-1)}.
One of these is a set of smoothed data, c(m-1), say. This looks like a smoothed version of cm. The other is a vector, d(m-1), say. This corresponds to the detail removed in smoothing cm to c(m-1). More precisely, they are the coefficients of the wavelet expansion corresponding to the highest resolution wavelets in the expansion. Similarly, c(m-2) and d(m-2) are obtained from c(m-1), etc. until you reach c0 and d0. 

All levels of smoothed data are stacked into a single vector for memory efficiency and ease of transport across the SPlus-C interface. 

The smoothing is performed directly by convolution with the wavelet filter
(\code{filter.select(n)$H}, essentially low- pass filtering), and then dyadic decimation (selecting every other datum, see Vaidyanathan (1990)). The detail extraction is performed by the mirror filter of H, which we call G and is a bandpass filter. G and H are also known quadrature mirror filters. 

There are now two methods of handling "boundary problems". If you know that your function is periodic (on it's interval) then use the bc="periodic" option, if you think that the function is symmetric reflection about each boundary then use bc="symmetric". You might also consider using the "wavelets on the interval" transform which is suitable for data arising from a function that is known to be defined on some compact interval, see Cohen, Daubechies, and Vial, 1993. If you don't know then it is wise to experiment with both methods, in any case, if you don't have very much data don't infer too much about your decomposition! If you have loads of data then don't infer too much about the boundaries. It can be easier to interpret the wavelet coefficients from a bc="periodic" decomposition, so that is now the default. Numerical Recipes implements some of the wavelets code, in particular we have compared our code to "wt1" and "daub4" on page 595. We are pleased to announce that our code gives the same answers! The only difference that you might notice is that one of the coefficients, at the beginning or end of the decomposition, always appears in the "wrong" place. This is not so, when you assume periodic boundaries you can imagine the function defined on a circle and you can basically place the coefficient at the beginning or the end (because there is no beginning or end, as it were). 

The non-deciated DWT contains all circular shifts of the standard DWT. Naively imagine that you do the standard DWT on some data using the Haar wavelets. Coefficients 1 and 2 are added and difference, and also coefficients 3 and 4; 5 and 6 etc. If there is a discontinuity between 1 and 2 then you will pick it up within the transform. If it is between 2 and 3 you will loose it. So it would be nice to do the standard DWT using 2 and 3; 4 and 5 etc. In other words, pick up the data and rotate it by one position and you get another transform. You can do this in one transform that also does more shifts at lower resolution levels. There are a number of points to note about this transform. 

Note that a time-ordered non-decimated wavelet transform object may be converted into a \code{packet-ordered non-decimated wavelet transform} object (and vice versa) by using the \code{\link{convert}} function. 

The NDWT is translation equivariant. The DWT is neither translation invariant or equivariant. The standard DWT is orthogonal, the non-decimated transform is most definitely not. This has the added disadvantage that non-decimated wavelet coefficients, even if you supply independent normal noise. This is unlike the standard DWT where the coefficients are independent (normal noise). 

You might like to consider growing wavelet syntheses using the
\code{\link{wavegrow}} function.
}
\value{
An object of class \code{\link{wd}}.

For boundary conditions apart from \code{bc="interval"} this object is a list with the following components. 

\item{C}{Vector of sets of successively smoothed data. The pyramid structure of Mallat is stacked so that it fits into a vector. The function \code{\link{accessC}} should be used to extract a set for a particular level.}
\item{D}{Vector of sets of wavelet coefficients at different resolution levels. Again, Mallat's pyramid structure is stacked into a vector. The function \code{\link{accessD}} should be used to extract the coefficients for a particular resolution level.} 
\item{nlevelsWT}{The number of resolution levels. This depends on the length of the data vector. If \code{length(data)=2^m}, then there will be m resolution levels. This means there will be m levels of wavelet coefficients (indexed 0,1,2,...,(m-1)), and m+1 levels of smoothed data (indexed 0,1,2,...,m). }
\item{fl.dbase}{There is more information stored in the C and D than is described above. In the decomposition ``extra'' coefficients are generated that help take care of the boundary effects, this database lists where these start and finish, so the "true" data can be extracted.}
\item{filter}{A list containing information about the filter type: Contains the string "wavelet" or "station" depending on which type of transform was performed. }
\item{date}{The date the transform was performed.}
\item{bc}{How the boundaries were handled.}

If the ``wavelets on the interval'' transform is used (i.e. \code{bc="interval"}) then the internal structure of the wd object is changed as follows. 

\itemize{
\item{The coefficient vectors C and D have been replaced by a single vector \code{transformed.vector}. The new single vector contains just the transformed coefficients: i.e. the wavelet coefficients down to a particular scale (determined by \code{min.scale} above). The scaling function coefficients are stored first in the array (there will be \code{2^min.scale} of them. Then the wavelet coefficients are stored as consecutive vectors coarsest to finest of length \code{2^min.scale}, \code{2^(min.scale+1)} up to a vector which is half of the length of the original data.)
 
In any case the user is recommended to use the functions \code{\link{accessC}}, \code{\link{accessD}}, \code{\link{putC}} and \code{\link{putD}} to access coefficients from the \code{\link{wd}} object.} 

\item{The extra component \code{current.scale} records to which level the transform has been done (usually this is \code{min.scale} as specified in the arguments).}
\item{The extra component \code{filters.used} is a vector of integers that record which filter index was used as each level of the decomposition. At coarser scales sometimes a wavelet with shorter support is needed. }
\item{The extra logical component \code{preconditioned} specifies whether preconditioning was turned on or off.}
\item{The component \code{fl.dbase} is still present but only contains data corresponding to the storage of the coefficients that are present in \code{transformed.vector}. In particular, since only one scale of the father wavelet coefficients is stored the component \code{first.last.c} of \code{fl.dbase} is now a three-vector containing the indices of the first and last entries of the father wavelet coefficients and the offset of where they are stored in \code{transformed.vector}. Likewise, the component \code{first.last.d} of \code{fl.dbase} is still a matrix but there are now only rows for each scale level in the \code{transformed.vector} (something like \code{nlevelsWT(wd)-wd$current.scale}). }
\item{The \code{filter} coefficient is also slightly different as the filter coefficients are no longer stored here (since they are hard coded into the wavelets on the interval transform.)} 
}
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 Integration of ``wavelets on the interval'' code by Piotr Fryzlewicz and Markus Monnerjahn was at Version 3.9.6, 1999. }
\seealso{
\code{\link{wd.int}}, \code{\link{wr}}, \code{\link{wr.int}}, \code{\link{wr.wd}}, \code{\link{accessC}}, \code{\link{accessD}}, \code{\link{putD}}, \code{\link{putC}}, \code{\link{filter.select}}, \code{\link{plot.wd}}, \code{\link{threshold}}, \code{\link{wavegrow}}
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
\dontrun{ts.plot(test.data)}
#
# Decompose test.data and plot the wavelet coefficients
#
wds <- wd(test.data)
\dontrun{plot(wds)}
#
# Now do the time-ordered non-decimated wavelet transform of the same thing
#
wdS <- wd(test.data, type="station")
\dontrun{plot(wdS)}
#
# Next examples
# ------------
# The chirp signal is also another good examples to use.
#
# Generate some test data
#
test.chirp <- simchirp()$y
\dontrun{ts.plot(test.chirp, main="Simulated chirp signal")}
#
# Now let's do the time-ordered non-decimated wavelet transform.
# For a change let's use Daubechies least-asymmetric phase wavelet with 8
# vanishing moments (a totally arbitrary choice, please don't read
# anything into it).
#
chirpwdS <- wd(test.chirp, filter.number=8, family="DaubLeAsymm", type="station")
\dontrun{plot(chirpwdS, main="TOND WT of Chirp signal")}
#
# Note that the coefficients in this plot are exactly the same as those
# generated by the packet-ordered non-decimated wavelet transform
# except that they are in a different order on each resolution level.
# See Nason, Sapatinas and Sawczenko, 1998
# for further information.
}
\keyword{math}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
