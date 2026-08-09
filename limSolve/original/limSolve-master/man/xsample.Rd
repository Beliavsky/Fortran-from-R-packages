\name{xsample}
\alias{xsample}
\title{
  Randomly samples an underdetermined problem with linear equality and
  inequality constraints
}
\description{
  Random sampling of inverse linear problems with linear equality and
  inequality constraints.
  Uses either a "hit and run" algorithm (random or coordinate directions)
  or a mirroring technique for sampling.
 
 The Markov Chain Monte Carlo method produces a sample solution for
 \deqn{Ex=f}
 \deqn{Ax\simeq B}{Ax~=B}
 \deqn{Gx>=h}
 where \eqn{Ex=F} have to be met exactly, and x is distributed
 according to \eqn{ p(\mathbf{x})\propto
  e^{-\frac{1}{2}(\mathbf{Ax-b})^T\mathbf{W}^2(\mathbf{Ax-b})}
  }{p(x)~e^(-.5||W(Ax-B)||^2)}
}

\usage{
xsample(A = NULL, B = NULL, E = NULL, F =NULL, 
        G = NULL, H = NULL, sdB = NULL, W = 1, 
        iter = 3000, outputlength = iter, burninlength = NULL, 
        type = "mirror", jmp = NULL, tol = sqrt(.Machine$double.eps), 
        x0 = NULL, fulloutput = FALSE, test = TRUE, ispos = FALSE, 
        verbose = TRUE, lower = NULL, upper = NULL)
}

\arguments{
  \item{A }{numeric matrix containing the coefficients of the
    (approximate) equality constraints, \eqn{Ax\simeq B}{Ax~=B}.
  }
  \item{B }{numeric vector containing the right-hand side of the
    (approximate) equality constraints.
  }
  \item{E }{numeric matrix containing the coefficients of the (exact)
    equality constraints, \eqn{Ex=F}.
  }
  \item{F }{numeric vector containing the right-hand side of the
    (exact) equality constraints.
  }
  \item{G }{numeric matrix containing the coefficients of the inequality
    constraints, \eqn{Gx>=H}.
  }
  \item{H }{numeric vector containing the right-hand side of the inequality
    constraints.
  }
  \item{sdB }{vector with standard deviation on B. Defaults to \code{NULL}.
  }
  \item{W }{weighting for \eqn{Ax\simeq B}{Ax~=B}. Only used if
    \code{sdB=NULL} and the problem is
    overdetermined. In that case, the error of B around the model Ax is
    estimated based on the residuals of \eqn{Ax\simeq B}{Ax~=B}. This
    error is made proportional to 1/W. If sdB is not NULL, \eqn{W=diag(sdB^-1)}.

  }
  \item{iter }{integer determining the number of iterations.
  }
  \item{outputlength}{ number of iterations kept in the output; at most
    equal to \code{iter}.
  }
  \item{burninlength}{a number of extra iterations, performed at first, to
    "warm up" the algorithm.
  }
  \item{type }{type of algorithm: one of: "mirror", (mirroring algorithm),
    "rda" (random directions algorithm) or "cda" (coordinates directions
    algorithm).
  }
  \item{jmp }{jump length of the transformed variables q: \eqn{x=x0+Zq}
    (only if \code{type}=="mirror"); if jmp is \code{NULL}, a reasonable
    value is determined by xsample, depending on the size of the NULL space.
  }
  \item{tol }{tolerance for equality and inequality constraints; numbers
    whose absolute value is smaller than \code{tol} are set to zero.
  }
  \item{x0 }{initial (particular) solution.
  }
  \item{fulloutput }{if \code{TRUE}, also outputs the transformed variables q.
  }
  \item{test }{if \code{TRUE}, xsample will test for hidden equalities (see
    details). This may be necessary for large problems, but slows down
    execution a bit.
  }
  \item{ispos }{if \code{TRUE}, it is imposed that unknowns are positive
    quantities. Used when generating \code{x0} by linear programming, if x0 is not given.
  }
  \item{verbose }{logical to print warnings and messages.
  }
  \item{upper, lower }{vector containing upper and lower bounds 
    on the unknowns. If one value, it is assumed to apply to all unknowns.
    If a vector, it should have a length equal to the number of unknowns; this
    vector can contain NA for unbounded variables. 
    The upper and lower bounds are added to the inequality conditions G*x>=H.
  }

}
\value{
  a list containing:
  \item{X }{matrix whose rows contain the sampled values of x.
  }
  \item{acceptedratio }{ratio of acceptance (i.e. the ratio of the accepted
    runs / total iterations).
  }
  \item{Q }{only returned if \code{fulloutput} is \code{TRUE}: the
    transformed samples Q.
  }
  \item{p }{only returned if \code{fulloutput} is \code{TRUE}: probability
    vector for all samples (e.g. one value for each row of \code{X}).
  }
  \item{jmp }{the jump length used for the random walk. Can be used to
    check the automated jump length.
  }
}

\author{
  Karel Van den Meersche

  Karline Soetaert <karline.soetaert@nioz.nl>
}

\examples{
#-------------------------------------------------------------------------------
# A simple problem
#-------------------------------------------------------------------------------
# Sample the probability density function of x1,...x4
# subject to:
# x1 + x2       + x4 = 3
#      x2  -x3  + x4 = -1
# xi   > 0

E <- matrix(nrow = 2, byrow = TRUE, data = c(1, 1, 0,  1,
                                             0, 1, -1, 1))
F   <- c(3, -1)

xs  <- xsample(E = E, F = F, lower = 0)
pairs(xs)

#-------------------------------------------------------------------------------
# Sample the underdetermined Mink diet problem
#-------------------------------------------------------------------------------
E <- rbind(Minkdiet$Prey, rep(1, 7))
F <- c(Minkdiet$Mink, 1)

# Here the Requirement x > 0 is been inposed in G and H.
pairs(xsample(E = E, F = F, G = diag(7), H = rep(0, 7), iter = 5000,
      output = 1000, type = "cda")$X,
      main = "Minkdiet 1000 solutions, - cda")
}
\references{
  Van den Meersche K, Soetaert K, Van Oevelen D (2009). xsample(): An R
  Function for Sampling Linear Inverse Problems.
  Journal of Statistical Software, Code Snippets, 30(1), 1-15.

}
\details{The algorithm proceeds in two steps.

  \enumerate{
	  \item  the equality constraints \eqn{Ex=F} are eliminated, and the
      system \eqn{Ex=f}, \eqn{Gx>=h} is rewritten as \eqn{G(p+Zq)>= h},
      i.e. containing only inequality constraints and where Z is a basis for
      the null space of E.
  	\item  the distribution of \eqn{q} is sampled numerically
      using a random walk (based on the Metropolis algorithm).
  }

  There are three algorithms for selecting new samples: \code{rda},
   \code{cda} (two hit-and-run algorithms) and a novel \code{mirror} algorithm.

  \itemize{
	  \item In the \code{rda} algorithm first a random direction is selected,
      and the new sample obtained by uniformly sampling the line
      connecting the old sample and the intersection with the planes defined
      by the inequality constraints.
    \item the \code{cda} algorithm is similar, except that the direction is
      chosen along one of the coordinate axes.
    \item the \code{mirror} algorithm is yet unpublished; it uses the
      inequality constraints as "reflecting planes" along which jumps are
      reflected.
      In contrast to \code{cda} and \code{rda}, this algorithm also works
      with unbounded problems (i.e. for which some of the unknowns can attain
      Inf).
  }
  For more information, see the package vignette \code{vignette(xsample)} or
  the file xsample.pdf in the packages \file{docs} subdirectory.
  
  Raftery and Lewis (1996) suggest a minimum of 3000 iterations to reach
  the extremes.

  If provided, then \code{x0} should be a valid particular solution (i.e.
  \eqn{E*x0=b} and \eqn{G*x0>=h}), else the algorithm will fail.

  For larger problems, a central solution may be necessary as a starting
  point for the \code{rda} and \code{cda} algorithms.  A good starting
  value is provided by the "central" value when running the function
  \code{\link{xranges}} with option \code{central} equal to \code{TRUE}.

  If the particular solution (\code{x0}) is not provided, then the
  parsimonious solution is sought, see \code{\link{ldei}}.

  This may however not be the most efficient way to start the algorithm. The
  parsimonious solution is usually located near the edges, and the
  \code{rda} and \code{cda} algorithms may not get out of this corner.
  The \code{mirror} algorithm is insensitive to that. Here it may be even
  better to start in a corner (as this position will always never be
  reached by random sampling).

  The algorithm will fail if there are hidden equalities. For instance, two
  inequalities may together impose an equality on an unknown,
  or, inequalities may impose equalities on a linear combination of two or
  more unknowns.

  In this case, the basis of the null space Z will be deficient. Therefore,
  \code{xsample} starts by checking if such hidden equalities exist.
  If it is suspected that this is NOT the case, set \code{test} to
  \code{FALSE}. This will speed up execution slightly.

  It is our experience that for small problems either the \code{rda} and
  \code{cda} algorithms are often more efficient.
  For really large problems, the \code{mirror} algorithm is usually much more
  efficient; select a jump length (\code{jmp}) that ensures good random
  coverage, while still keeping the number of reflections reasonable.
  If unsure about the size of jmp, the default will do.

  See \code{\link{E_coli}} for an example where a relatively large problem
  is sampled.
}
\seealso{
  \code{\link{Minkdiet}}, for a description of the Mink diet example.

  \code{\link{ldei}}, to find the least distance solution

  \code{\link{lsei}}, to find the least squares solution

  \code{\link{varsample}}, to randomly sample variables of an lsei problem.

  \code{\link{varranges}}, to estimate ranges of inverse variables.
}

\keyword{algebra}
\keyword{array}
\keyword{optimize}

