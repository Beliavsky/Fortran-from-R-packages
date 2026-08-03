\name{getEER}
\alias{getEER}

\encoding{latin1}

\title{Download effective exchange rates data frame from Bank of International Settlement}

\description{
It downloads effective exchange rates data frame from BIS.
}
\usage{
getEER(Areas=c("US","JP"),Freq="Monthly", Type="Real",Basket="Broad")
}

\arguments{
  \item{Areas}{
The currency country of effective exchange rates (EER) of BIS, it allows many countries. The defaults are "US" and "JP".
}

  \item{Freq}{The type of EER, it has two types: "Daily" and "Monthly". The default is "Monthly". "Daily" works only for "Nominal" of NEER.

}

  \item{Type}{The type of EER, it has two types: "Nominal" and "Real". The default is "Real".

}

  \item{Basket}{The basket of EER, it has two types: "Broad" and "Narrow". The default is "Broad".

}

}


\details{

This function connects with <"https://www.bis.org/statistics/eer/"> and downloads the specified data. The rownames of downloaded data embeds the timestamp already, which can be directly transformed into \code{timeSeries} via, as.timeSeries.
}


\value{

  Retrieved BIS effective exchange rates time series, without any time series format.

}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.

}

\examples{
#To save time, the example below is commented.
output=getEER(Areas=c("US","JP"),Freq="Monthly", Type="Real",Basket="Broad")
head(output)

}
