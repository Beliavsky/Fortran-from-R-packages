\name{Pen}
\alias{Pen}
\title{Pen's Parade}
\description{
 plots Pen's Parade of a vector \code{x}}
\usage{
Pen(x, n = rep(1, length(x)), group = NULL,
  scaled = TRUE, abline = TRUE, add = FALSE, segments = NULL,  
  main = "Pen's Parade", ylab = NULL, xlab = NULL, 
  col = NULL, lwd = NULL, las = 1, fill = NULL, \dots)  
}
\arguments{
 \item{x}{a vector containing non-negative elements.}
 \item{n}{a vector of frequencies or weights, must be same length as \code{x}.}
 \item{group}{a factor coding different groups, must be same length as \code{x}. See also details.}
 \item{scaled}{logical. Should Pen's parade be divided by \code{mean(x)}?}
 \item{abline}{logical. Should a horizontal line for the mean be drawn?}
 \item{add}{logical. Should the plot be added to an existing plot?}
 \item{segments}{logical. Should histogram-like segments be drawn?}
 \item{col}{a (vector of) color(s) for drawing the curve.}
 \item{fill}{a (vector of) color(s) for filling the area under the curve.}
 \item{xlab,ylab}{axis labels. Suitable defaults depending on
   \code{scaled} and \code{n} are chosen.}
 \item{main, lwd, las, \dots}{further high-level \code{\link{plot}} parameters.}
}
\details{
Pen's Parade is basically the inverse distribution function
(standardized by \code{mean(x)}).

\code{Pen} allows for fine control of the layout---the graphical parameters \code{col}
and \code{fill} can be vectorized if histogram-like segments are drawn
(\code{segments = TRUE})---but implements several heuristics in choosing its
default plotting parameters. If a grouping factor \code{group} is given,
the default is to draw segments with a grey-shaded filling. If no fill color
is used, the default is to draw a thick blue curve. But as all of these are just
defaults, they can of course easily be changed. See also the examples.
}

\references{
F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam,

F A Cowell: Measuring Inequality, 1995 Prentice Hall/Harvester Wheatshef,

J Pen: Income Distribution, 1971, Harmondsworth: Allen Lane.}

\seealso{\code{\link{Lc}}, \code{\link{plot.Lc}}}

\examples{
# load and attach Philippine income data
data(Ilocos)
attach(Ilocos)
# plot Pen's Parade of income
Pen(income)
Pen(income, fill = hsv(0.1, 0.3, 1))

# income distribution of the USA in 1968 (in 10 classes)
# x vector of class means, n vector of class frequencies
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
n <- c(482, 825, 722, 690, 661, 760, 745, 2140, 1911, 1024)
Pen(x, n = n)
# create artificial grouping variable
myfac <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
Pen(x, n = n, group = myfac)
}
\keyword{hplot}

