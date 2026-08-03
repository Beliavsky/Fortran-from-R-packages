\name{getFrench.Portfolios}
\alias{getFrench.Portfolios}

\encoding{latin1}

\title{Download 24 asset pricing factors data from the data library of Dr. French}

\description{
It downloads 24 factors data used for asset pricing analysis from the data library of Dr. Kenneth R. French at Dartmouth College.
}
\usage{
getFrench.Portfolios(filename="Portfolios_Formed_on_ME")
}

\arguments{
  \item{filename}{
The name of portfolio data file as listed in (http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html), the default is "Portfolios_Formed_on_ME". So far, this function supports retrieving 24 portfolio data files: \cr
"Portfolios_Formed_on_ME",#Portfolios Formed on Size\cr
"Portfolios_Formed_on_BE-ME",#On Book-to-Market\cr
"Portfolios_Formed_on_OP", #On Operating Profitability\cr
"Portfolios_Formed_on_INV",#On Investment\cr
"6_Portfolios_2x3",        #6 Ports on Size and Book-to-Market\cr
"25_Portfolios_5x5",       #25 Ports on Size and Book-to-Market\cr
"100_Portfolios_10x10",    #100 Ports on Size and Book-to-Market\cr
"6_Portfolios_ME_INV_2x3", #6 Ports on Size and Investment\cr
"25_Portfolios_ME_INV_5x5",#25 Ports on Size and Investment\cr
"100_Portfolios_ME_INV_10x10",#100 Ports on Size and Investment\cr
"25_Portfolios_BEME_OP_5x5",#25 Ports on Book-to-Market and Operating Profitability\cr
"25_Portfolios_BEME_INV_5x5",#25 Ports on Book-to-Market and Investment\cr
"25_Portfolios_OP_INV_5x5",#25 Ports on Operating Profitability and Investment\cr
"32_Portfolios_ME_BEME_OP_2x4x4",#32 Ports on Size,Book-to-Market and Operating Profitability\cr
"32_Portfolios_ME_BEME_INV_2x4x4",#32 Ports on Size,Book-to-Market and Investment\cr
"32_Portfolios_ME_OP_INV_2x4x4",#32 Ports on Size,Operating Profitability and Investment\cr
"5_Industry_Portfolios", # Industry Portfolios by c(5,10,12,17,30,38,48,49)

}


}


\details{

This function connects with <"http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html"> and downloads the specified portfolio data constructed by factors. Currently, we support on retriving and arranging 24 portfolio datasets. Sometimes, the datafile contains multiple data tables, hence the code returns a list. Since the csv spreedsheet on the web is not structured data table, it does not only contain irregular headings, but also 7, or more, data tablesr; and the trivial portfolio data tables labelled by "Portfolio Formed by ..." are omitted. Check the "table.names" of output object.
}


\value{

  \item{data}{The data retrieved and arranged.}
  \item{table.names}{The names of data table.}
  \item{file.name}{The file name of portfolio data.}

}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.

}

\examples{
#To save time, the example below is commented.
output=getFrench.Portfolios()
output$file.name
dim(output$data[[1]])
output$data
output$table.names

}
