\name{rice}
\alias{rice}
\docType{data}
\title{ Rice taste data}
\description{
The rice taste data consists of five inputs and a single output whose values 
are associated with subjective evaluations as follows: xl: flavor, 
x2: appearance, x3: taste, x4: stickiness, x5: toughness, y: overall evaluation.
Sensory test data have been obtained by such subjective evaluations for 105 kinds of
rice (e.g., Sasanishiki, Akita-Komachi, etc.). The data set was used by 
Nozaki et al. (1997) to demonstrate the high performance of a proposed 
for automatically generating fuzzy if-then rules from numerical data.    
}
\usage{data(rice)}
\format{
    A data frame with 105 observations on the following 6 variables: 
    \describe{
        \item{\code{Favor}}{compactness}
        \item{\code{Appearance}}{circularity}
        \item{\code{Taste}}{distance circularity}
        \item{\code{Stickiness}}{radius ratio}
        \item{\code{Toughness}}{principal axis aspect ratio}
        \item{\code{Overall_evaluation}}{maximum length aspect ratio}
    }
}
%%\details{}
\source{
Nozaki, K., Ishibuchi, H, and Tanaka, H. (1997) 
A simple but powerful heuristic method for generating fuzzy rules from numerical data
\emph{Fuzzy Sets and Systems} \bold{86} 3 p. 251--270. 
}
%\references{}
%\examples{}
\keyword{datasets}
