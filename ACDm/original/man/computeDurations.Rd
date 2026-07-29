\name{computeDurations}
\alias{computeDurations}

\title{
Durations computation
}
\description{
Computes durations from a data.frame containing the time stamps of transactions. Trade durations, price durations and volume durations can be computed (if the appropriate data columns are given).
}
\usage{
computeDurations(transactions, open = "10:00:00", close = "18:25:00", 
rm0dur = TRUE, type = "trade", priceDiff = .1, cumVol = 10000)
}
\arguments{
  \item{transactions}{
a data.frame with, at least, transaction time in a column named 'time' (see Details)
}
  \item{open}{
the opening time of the exchange. Transactions done outside the trading hours will be ignored. 
}
  \item{close}{
the closing time of the exchange.
}
  \item{rm0dur}{
if \code{TRUE} zero-durations will be removed and transactions done on the same second will be aggregated, e.g. price will then be the volume weighted avrage price of the aggregated transactions. 
}
  \item{type}{
the type of durations to be computed. Either "trade", "price", or "volume".
}
  \item{priceDiff}{
only if \code{type = "price"}. Price durtions are (here) defind as the duration until the price has changed by at least 'priceDiff' in absolute value.
}
  \item{cumVol}{
only if \code{type = "cumVol"}. Volume durtions are (here) defind as the duration until the cumulative traded volume since the last duration has surpassed 'cumVol'.
}
}
\details{
The data.frame must include a column named 'time' with the time of each transaction, in a time format recognizable by \link{POSIXlt} or strings in format "yyyy-mm-dd hh:mm:ss". If the column 'price' or 'volume' is included its also possible to compute price- and volume durations (see arguments \code{priceDiff} and \code{cumVol}) 
}
\value{
a data.frame with columns:

    \item{time}{
        the calander time of the start of each duration spell.
    }
    
    \item{price}{
        the volume weighted avrage price of the shares traded during the spell of the duration.
    }
    
    \item{volume}{
        the volume (total shares traded) during the duration spell.
    }
    
    \item{Ntrans}{
        number of transactions done during the spell.
    }
    
    \item{durations}{
        the computed durations.
    }
}
\author{
Markus Belfrage
}
\examples{
\dontrun{
#only the first 3 days of data:
data(transData)
durDataShort <- computeDurations(transData[1:56700, ]) 
head(durDataShort)}
}