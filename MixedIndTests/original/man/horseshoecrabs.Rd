    \name{horseshoecrabs}
     \docType{data}
     \alias{horseshoecrabs}
     \title{Horseshoecrabs dataset}
     \description{Horseshoe Crab Data from Table 3.2 of Agresti(2007). This data set consists of five variables, three of which are categorical, measured on 173 female crabs, each having a male attached in her nest.}
     \usage{data(horseshoecrabs)}
     \format{Data frame with 173 rows and 5 variables:
\itemize{
  \item X1: Color of the female (1: light medium, 2:  medium, 3: dark medium, 4: dark)
 \item X2: Spine condition (1: both good. 2: one worn or broken, 3:  both worn or broken)
 \item X3: Carapace width (cm)
 \item X4: Number of satellites, i.e., other males around the female
 \item X5: Weight (kg)}
 }
\keyword{datasets}
\references{  Agresti, A. (2007). An Introduction to Categorical data analysis,
 John Wiley & Sons, Wiley Series in Probability and Statistics, 2nd edition.}
\examples{
data(horseshoecrabs)
x =data.matrix(horseshoecrabs)
out = TestIndCopula(x,trunc.level=5,graph=TRUE)
}

