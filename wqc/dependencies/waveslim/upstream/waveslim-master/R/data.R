#' Autocovariance and Autocorrelation Sequences for a Seasonal Persistent
#' Process
#' 
#' The autocovariance and autocorrelation sequences from the time series model
#' in Figures 8, 9, 10, and 11 of Andel (1986).  They were obtained through
#' numeric integration of the spectral density function.
#' 
#' @usage data(acvs.andel8)
#' @usage data(acvs.andel9)
#' @usage data(acvs.andel10)
#' @usage data(acvs.andel11)
#' @name acvs.andel8
#' @docType data
#' @aliases acvs.andel9 acvs.andel10 acvs.andel11
#' @format A data frame with 4096 rows and three columns: lag, autocovariance
#' sequence, autocorrelation sequence.
#' @references Andel, J. (1986) Long memory time series models,
#' \emph{Kypernetika}, \bold{22}, No. 2, 105-123.
#' @keywords datasets
NULL

#' Simulated AR(1) Series
#' 
#' Simulated AR(1) series used in Gencay, Selcuk and Whitcher (2001).
#' 
#' @usage data(ar1)
#' @name ar1
#' @docType data
#' @format A vector containing 200 observations.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' @keywords datasets
NULL

#' Barbara Test Image
#' 
#' The Barbara image comes from Allen Gersho's lab at the University of
#' California, Santa Barbara.
#' 
#' @usage data(barbara)
#' @name barbara
#' @docType data
#' @format A 256 \eqn{\times}{x} 256 matrix.
#' @source Internet.
#' @keywords datasets
NULL

#' A Piecewise-Constant Function
#' 
#' \deqn{blocks(x) = \sum_{j=1}^{11}(1 + {\rm sign}(x-p_j)) h_j / 2}{%
#' blocks(x) = sum[j=1,11] (1 + sign(x - p_j)) h_j/2}
#' 
#' @usage data(blocks)
#' @name blocks
#' @docType data
#' @format A vector containing 512 observations.
#' @references Bruce, A., and H.-Y. Gao (1996) \emph{Applied Wavelet Analysis
#' with S-PLUS}, Springer: New York.
#' @source S+WAVELETS.
#' @keywords datasets
NULL

#' U.S. Consumer Price Index
#' 
#' Monthly U.S. consumer price index from 1948:1 to 1999:12.
#' 
#' @usage data(cpi)
#' @name cpi
#' @docType data
#' @format A vector containing 624 observations.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' @source Unknown.
#' @keywords datasets
NULL

#' Digital Photograph of Ingrid Daubechies
#' 
#' A digital photograph of Ingrid Daubechies taken at the 1993 AMS winter
#' meetings in San Antonio, Texas.  The photograph was taken by David Donoho
#' with a Canon XapShot video still frame camera.
#' 
#' @usage data(dau)
#' @name dau
#' @docType data
#' @format A 256 \eqn{\times}{x} 256 matrix.
#' @references Bruce, A., and H.-Y. Gao (1996) \emph{Applied Wavelet Analysis
#' with S-PLUS}, Springer: New York.
#' @source S+WAVELETS.
#' @keywords datasets
NULL

#' Sinusoid with Changing Amplitude and Frequency
#' 
#' \deqn{doppler(x) = \sqrt{x(1 - x)} }{% doppler(x) = sqrt{x(1-x)}
#' sin[(2.1*pi)/(x+0.05)]}\deqn{ \sin\left(\frac{2.1\pi}{x+0.05}\right)}{%
#' doppler(x) = sqrt{x(1-x)} sin[(2.1*pi)/(x+0.05)]}
#' 
#' @usage data(doppler)
#' @name doppler
#' @docType data
#' @format A vector containing 512 observations.
#' @references Bruce, A., and H.-Y. Gao (1996) \emph{Applied Wavelet Analysis
#' with S-PLUS}, Springer: New York.
#' @source S+WAVELETS.
#' @keywords datasets
NULL

#' Exchange Rates Between the Deutsche Mark, Japanese Yen and U.S. Dollar
#' 
#' Monthly foreign exchange rates for the Deutsche Mark - U.S. Dollar (DEM-USD)
#' and Japanese Yen - U.S. Dollar (JPY-USD) starting in 1970.
#' 
#' @usage data(exchange)
#' @name exchange
#' @docType data
#' @format A bivariate time series containing 348 observations.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' @source Unknown.
#' @keywords datasets
NULL

#' Sine with Jumps at 0.3 and 0.72
#' 
#' \deqn{heavisine(x) = 4\sin(4{\pi}x) - \mathrm{sign}(x-0.3) - }{%
#' heavisine(x) = 4*sin(4*pi*x) - sign(x-0.3) - sign(0.72-x)}\deqn{
#' \mathrm{sign}(0.72-x)}{% heavisine(x) = 4*sin(4*pi*x) - sign(x-0.3) -
#' sign(0.72-x)}
#' 
#' @usage data(heavisine)
#' @name heavisine
#' @docType data
#' @format A vector containing 512 observations.
#' @references Bruce, A., and H.-Y. Gao (1996) \emph{Applied Wavelet Analysis
#' with S-PLUS}, Springer: New York.
#' @source S+WAVELETS.
#' @keywords datasets
NULL

#' Daily IBM Stock Prices
#' 
#' Daily IBM stock prices spanning May~17, 1961 to November~2, 1962.
#' 
#' @usage data(ibm)
#' @name ibm
#' @docType data
#' @format A vector containing 369 observations.
#' @source Box, G. E.P. and Jenkins, G. M. (1976) \emph{Time Series Analysis:
#' Forecasting and Control}, Holden Day, San Francisco, 2nd edition.
#' @keywords datasets
NULL

#' Japanese Gross National Product
#' 
#' Quarterly Japanese gross national product from 1955:1 to 1996:4.
#' 
#' @usage data(japan)
#' @name japan
#' @docType data
#' @format A vector containing 169 observations.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' 
#' Hecq, A. (1998) Does seasonal adjustment induce common cycles?,
#' \emph{Empirical Economics}, \bold{59}, 289-297.
#' @source Unknown.
#' @keywords datasets
NULL

#' Sine with Jumps at 0.625 and 0.875
#' 
#' \deqn{jumpsine(x) = 10\left( \sin(4{\pi}x) +
#' I_{[0.625 < x \leq 0.875]}\right)}{%
#' jumpsine(x) = 10*(sin(4*pi*x) + I_[0.625 < x <= 0.875])}
#' 
#' @usage data(jumpsine)
#' @name jumpsine
#' @docType data
#' @format A vector containing 512 observations.
#' @references Bruce, A., and H.-Y. Gao (1996) \emph{Applied Wavelet Analysis
#' with S-PLUS}, Springer: New York.
#' @source S+WAVELETS.
#' @keywords datasets
NULL

#' 1995 Kobe Earthquake Data
#' 
#' Seismograph (vertical acceleration, nm/sq.sec) of the Kobe earthquake,
#' recorded at Tasmania University, HobarTRUE, Australia on 16 January 1995
#' beginning at 20:56:51 (GMTRUE) and continuing for 51 minutes at 1 second
#' intervals.
#' 
#' @usage data(kobe)
#' @name kobe
#' @docType data
#' @format A vector containing 3048 observations.
#' @source Data management centre, Washington University.
#' @keywords datasets
NULL

#' Linear Chirp
#' 
#' \deqn{linchirp(x) = \sin(0.125 \pi n x^2)}{% 
#' linchirp(x) = sin(0.125*pi*n*x^2)}
#' 
#' @usage data(linchirp)
#' @name linchirp
#' @docType data
#' @format A vector containing 512 observations.
#' @references Bruce, A., and H.-Y. Gao (1996) \emph{Applied Wavelet Analysis
#' with S-PLUS}, Springer: New York.
#' @source S+WAVELETS.
#' @keywords datasets
NULL

#' Mexican Money Supply
#' 
#' Percentage changes in monthly Mexican money supply.
#' 
#' @usage data(mexm)
#' @name mexm
#' @docType data
#' @format A vector containing 516 observations.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' @source Unknown.
#' @keywords datasets
NULL

#' Nile River Minima
#' 
#' Yearly minimal water levels of the Nile river for the years 622 to 1281,
#' measured at the Roda gauge near Cairo (Tousson, 1925, p. 366-385). The data
#' are listed in chronological sequence by row.
#' 
#' The original Nile river data supplied by Beran only contained only 500
#' observations (622 to 1121).  However, the book claimed to have 660
#' observations (622 to 1281).  The remaining observations from the book were
#' added, by hand, but the series still only contained 653 observations (622 to
#' 1264).
#' 
#' Note, now the data consists of 663 observations (spanning the years
#' 622-1284) as in original source (Toussoun, 1925).
#' 
#' @usage data(nile)
#' @name nile
#' @docType data
#' @format A length 663 vector.
#' @references Beran, J. (1994) \emph{Statistics for Long-Memory Processes},
#' Chapman Hall: Englewood, NJ.
#' @source Toussoun, O. (1925) M\'emoire sur l'Histoire du Nil, Volume 18 in
#' \emph{M\'emoires a l'Institut d'Egypte}, pp. 366-404.
#' @keywords datasets
NULL

#' U.S. Tourism
#' 
#' Quarterly U.S. tourism figures from 1960:1 to 1999:4.
#' 
#' @usage data(tourism)
#' @name tourism
#' @docType data
#' @format A vector containing 160 observations.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' @source Unknown.
#' @keywords datasets
NULL

#' U.S. Unemployment
#' 
#' Monthly U.S. unemployment figures from 1948:1 to 1999:12.
#' 
#' @usage data(unemploy)
#' @name unemploy
#' @docType data
#' @format A vector containing 624 observations.
#' @references Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An
#' Introduction to Wavelets and Other Filtering Methods in Finance and
#' Economics}, Academic Press.
#' @source Unknown.
#' @keywords datasets
NULL

#' Image with Box and X
#' 
#' \deqn{xbox(i,j) = I_{[i=n/4,\;3n/4,\;j;~ n/4 \leq j \leq 3n/4]} + }{%
#' xbox(i,j) = I_[i = n/4, 3n/4, j; n/4 \leq j \leq 3n/4] + I_[n/4 \leq i \leq
#' 3n/4; j = n/4, 3n/4, i]}\deqn{ I_{[n/4 \leq i \leq 3n/4;~
#' j=n/4,\;3n/4,\;i]}}{% xbox(i,j) = I_[i = n/4, 3n/4, j; n/4 \leq j \leq 3n/4]
#' + I_[n/4 \leq i \leq 3n/4; j = n/4, 3n/4, i]}
#' 
#' @usage data(xbox)
#' @name xbox
#' @docType data
#' @format A 128 \eqn{\times}{x} 128 matrix.
#' @references Bruce, A., and H.-Y. Gao (1996) \emph{Applied Wavelet Analysis
#' with S-PLUS}, Springer: New York.
#' @source S+WAVELETS.
#' @keywords datasets
NULL
