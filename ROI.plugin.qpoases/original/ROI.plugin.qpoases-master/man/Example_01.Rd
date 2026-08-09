\name{Example-1}
\alias{Example-1}
\title{Quadratic Problem 1}


\description{
\deqn{maximize \ \ x_1^2 + x_2^2 + x_3^2 - 5 x_2}
\deqn{subject \ to:}
\deqn{-4 x_1 - 3 x_2 +     \geq -8}
\deqn{ 2 x_1 +   x_2 +     \geq 2}
\deqn{       - 2 x_2 + x_3 \geq 0}
\deqn{x_1, x_2, x_3 \geq 0}
}

\examples{
Sys.setenv(ROI_LOAD_PLUGINS = "FALSE")
library("ROI")
library("ROI.plugin.qpoases")

A <- cbind(c(-4, -3, 0), 
           c( 2,  1, 0), 
           c( 0, -2, 1))
x <- OP(Q_objective(diag(3), L =  c(0, -5, 0)),
        L_constraint(L = t(A),
                     dir = rep(">=", 3),
                     rhs = c(-8, 2, 0)))

opt <- ROI_solve(x, solver="qpoases")
opt
## Optimal solution found.
## The objective value is: -2.380952e+00
solution(opt)
## [1]  0.4761905 1.0476190 2.0952381

}

