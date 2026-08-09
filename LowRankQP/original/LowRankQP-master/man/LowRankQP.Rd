\name{LowRankQP}
\alias{LowRankQP}
\title{Solve Low Rank Quadratic Programming Problems}
\description{
  This routine implements a primal-dual interior point method
  solving quadratic programming problems of the form 
  \tabular{ll}{
    min       \tab \eqn{d^T alpha + 1/2 alpha^T H alpha} \cr
    such that \tab \eqn{A alpha = b} \cr
              \tab \eqn{0 <= alpha <= u} \cr
  }
  with dual
  \tabular{ll}{
    min       \tab \eqn{1/2 alpha^T H alpha + beta^T b + xi^T u} \cr
    such that \tab \eqn{H alpha + c + A^T beta - zeta + xi = 0} \cr
              \tab \eqn{xi, zeta >= 0} \cr
  }
  where \eqn{H=V} if \eqn{V} is square and \eqn{H=VV^T} otherwise.
}
\usage{
LowRankQP(Vmat,dvec,Amat,bvec,uvec,method="PFCF",verbose=FALSE,niter=200,epsterm=1.0E-8)
}
\arguments{
  \item{Vmat}{
    matrix appearing in the quadratic function to be minimized.
  }
  \item{dvec}{
    vector appearing in the quadratic function to be minimized.
  }
  \item{Amat}{
    matrix defining the constraints under which we want to minimize the
    quadratic function.
  }
  \item{bvec}{
    vector holding the values of \eqn{b} (defaults to zero).
  }
  \item{uvec}{
    vector holding the values of \eqn{u}.
  }
  \item{method}{
    Method used for inverting H+D where D is full rank diagonal. 
    If \eqn{V} is square:
    \itemize{
    \item 'LU':   Use LU factorization. (More stable)
    \item 'CHOL': Use Cholesky factorization. (Faster)
    }
    If \eqn{V} is not square:
    \itemize{    
    \item 'SMW':  Use Sherman-Morrison-Woodbury (Faster)
    \item 'PFCF': Use Product Form Cholesky Factorization (More stable) 
    }
  }
  \item{verbose}{
    Display iterations and termination statistics.
  }
  \item{niter}{
    Number of iteration to perform.
  }
  \item{epsterm}{
    Termination tolerance. See equation (12) of Ormerod et al (2008).
  }  
}
\value{
a list with the following components:
\item{alpha}{
vector containing the solution of the quadratic programming problem.
}
\item{beta}{
vector containing the solution of the dual of quadratic programming problem.
}
\item{xi}{
vector containing the solution of the dual quadratic programming problem.
}
\item{zeta}{
vector containing the solution of the dual quadratic programming problem.
}
}

\references{
Ormerod, J.T., Wand, M.P. and Koch, I. (2008). Penalised spline support vector 
classifiers: computational issues, Computational Statistics, 23, 623-641.

Boyd, S. and Vandenberghe, L. (2004). Convex Optimization. Cambridge University Press. 
      
Ferris, M. C. and Munson, T. S. (2003). Interior point methods for massive support 
vector machines. SIAM Journal on Optimization, 13, 783-804.
       
Fine, S. and Scheinberg, K. (2001). Efficient SVM training using low-rank kernel 
representations. Journal of Machine Learning Research, 2, 243-264.

B. Sch\"olkopf and A. J. Smola. (2002). Learning with Kernels. The MIT Press, 
Cambridge, Massachusetts.
}

\examples{
library(LowRankQP)

# Assume we want to minimize: (0 -5 0 0 0 0) \%*\% alpha + 1/2 alpha[1:3]^T alpha[1:3]
# under the constraints:      A^T alpha = b
# with b = (-8,  2,  0 )^T
# and      (-4   2   0 ) 
#      A = (-3   1  -2 )
#          ( 0   0   1 )
#          (-1   0   0 )
#          ( 0  -1   0 )
#          ( 0   0  -1 )
#  alpha >= 0
#
# (Same example as used in quadprog)
#
# we can use LowRankQP as follows:

Vmat          <- matrix(0,6,6)
diag(Vmat)    <- c(1, 1,1,0,0,0)
dvec          <- c(0,-5,0,0,0,0)
Amat          <- matrix(c(-4,-3,0,-1,0,0,2,1,0,0,-1,0,0,-2,1,0,0,-1),6,3)
bvec          <- c(-8,2,0)
uvec          <- c(100,100,100,100,100,100)
LowRankQP(Vmat,dvec,t(Amat),bvec,uvec,method="CHOL")

# Now solve the same problem except use low-rank V

Vmat          <- matrix(c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0),6,3)
dvec          <- c(0,-5,0,0,0,0)
Amat          <- matrix(c(-4,-3,0,-1,0,0,2,1,0,0,-1,0,0,-2,1,0,0,-1),6,3)
bvec          <- c(-8,2,0)
uvec          <- c(100,100,100,100,100,100)
LowRankQP(Vmat,dvec,t(Amat),bvec,uvec,method="SMW")

}

\keyword{optimize}


