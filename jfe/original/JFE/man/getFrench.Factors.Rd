\name{getFrench.Factors}
\alias{getFrench.Factors}

\encoding{latin1}

\title{Download seven asset pricing factors data from the data library of Dr. French}

\description{
It downloads seven factors data used for asset pricing analysis from the data library of Dr. Kenneth R. French at Dartmouth College.
}
\usage{
getFrench.Factors(filename="F-F_Research_Data_Factors")
}

\arguments{

  \item{filename}{
The name of data file as listed in <"http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html">, which is an important database for asset pricing literature, this function supports seven factor files: \cr
"F-F_Research_Data_Factors",# Fama/French 3 Factors  \cr
"F-F_Research_Data_5_Factors_2x3",# Fama/French 5 Factors \cr
"F-F_Momentum_Factor", # Fama/French Momentum Factors \cr
"F-F_Momentum_Factor_daily", # Fama/French Momentum Factors(daily) \cr
"F-F_Research_Data_Factors_weekly",# Fama/French 3 Factors(weekly) \cr
"F-F_Research_Data_Factors_daily", # Fama/French 3 Factors(daily) \cr
"F-F_Research_Data_5_Factors_2x3_daily" # Fama/French 5 Factors (daily)

}


}


\details{

This function connects with <"http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html"> and downloads the specified factors data. Sometimes, the datafile contains multiple data tables, hence the code returns a list. The website may occur somthing unexpected, for example, "F-F_Research_Data_Factors_weekly" is ill-formated hence the .csv data table cannot be retrieved. However, it was OK before. For this problem, we will not fix the code.
}


\value{

  \item{ff.factor}{The data retrieved and arranged.}

}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.

}

\examples{

getFrench.Factors()


}
