\name{FinTS-package}
\alias{FinTS-package}
\alias{FinTS}
\docType{package}
\title{
  \Sexpr[stage=build,results=hide]{pd <- packageDescription("FinTS")}
  \Sexpr[stage=build,results=hide]{lb <- library(help="FinTS", character.only=TRUE)}
  \Sexpr[stage=build]{pd$Title}
}
\description{
  \Sexpr[stage=build]{pd$Description}
}
\details{

  % 2023-08-20 Georgi: not really needed, commenting out.
  %
  % \tabular{ll}{
  % Package     \tab  \Sexpr[stage=build]{pd$Package   } \cr 
  % Version     \tab  \Sexpr[stage=build]{pd$Version   } \cr
  % Date        \tab  \Sexpr[stage=build]{pd$Date      } \cr    
  % Author      \tab  \Sexpr[stage=build]{pd$Author    } \cr     
  % Maintainer  \tab  \Sexpr[stage=build]{pd$Maintainer} \cr    
  % License     \tab  \Sexpr[stage=build]{pd$License   } \cr    
  % URL         \tab  \Sexpr[stage=build]{pd$URL       } \cr    
  % Depends     \tab  \Sexpr[stage=build]{pd$Depends   } \cr    
  % Suggests    \tab  \Sexpr[stage=build]{pd$Suggests  } \cr    
  % Built       \tab  \Sexpr[stage=build]{pd$Built     } \cr
  % }

   % this will print the index but this doesn't seem useful.
   % 
   Index:
   \Sexpr[stage=build,results=rd]{paste("\\\\preformatted{", paste(if(!is.null(lb$info[[2]])) lb$info[[2]] else "No entries",collapse="\n"), "}", sep="")}


    See the \code{scripts} subdirectory of the \code{FinTS} installation
    directory =  system.file(package='FinTS').

    Corrections to the script files provided and contributions to script
    files for other chapters will be graciously accepted.
}
\author{
Spencer Graves

Maintainer: % Spencer Graves <spencer.graves@prodsyse.com>
   \Sexpr[stage=build]{pd$Maintainer}
}
\references{
    Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley)
}
\keyword{ package }
\keyword{ ts }
\seealso{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\examples{
# Where is the 'FinTS' directory?
system.file(package='FinTS')

# View the script file 'ch01.R', which is in the 'scripts'
# subdirectory of the system.file(package='FinTS') directory:
runscript(1, 'view')

# SP statistics in Table 1.2 of Tsay
data(d.ibmvwewsp6203)
FinTS.stats(100*d.ibmvwewsp6203[, "SP"])
}
