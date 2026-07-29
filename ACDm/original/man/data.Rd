\name{DataFiles}

\alias{DataFiles}

\alias{transData}
\alias{durData}
\alias{adjDurData}
\alias{defaultSplineObj}

\title{Time Series Data Sets} 

\description{
    The data file \code{transData} is the base data used in all of the examples. It is a data.frame with rows representing a single transaction and has the columns 'time', 'price', giving the trade price, and 'volume', giving the number of shares traded for the transaction. The data set is based on real transactions but has been obfuscated by transforming the dates, price and volume, for proprietary reasons. It covers two weeks of nearly 100 000 transactions, recorded with 1 second precision.    

The \code{durData} data.frame is simply the trade durations formed from \code{transData} using the function \code{durData <- computeDurations(transData)} 

The \code{adjDurData} data object is in turn created by \code{adjDurData <- diurnalAdj(durData, aggregation = "all")} to add diurnally adjusted durations.

\code{defaultSplineObj} is an estimated cubic spline of the diurnal component using the sample data. It is used when simulating from \code{sim_ACD()} with the argument \code{diurnalFactor} set to \code{TRUE}, when no user \code{splineObj} is provided.
}
