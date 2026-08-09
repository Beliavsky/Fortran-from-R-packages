\name{Blending}
\docType{data}
\alias{Blending}
\title{
  A linear inverse blending problem
}

\description{
  A manufacturer produces a feeding mix for pet animals.

  The feed mix contains two nutritive ingredients and one ingredient (filler)
  to provide bulk.

  One kg of feed mix must contain a minimum quantity of each of four nutrients
  as below:

  \tabular{rlllll}{
     \tab Nutrient  \tab A  \tab B  \tab C  \tab D \cr
     \tab gram      \tab 80 \tab 50 \tab 25 \tab 5 \cr
  }

  The ingredients have the following nutrient values and cost

  \tabular{rllllll}{
    \tab  (gram/kg)    \tab A   \tab   B \tab C  \tab D    \tab Cost/kg \cr
    \tab Ingredient 1  \tab 100 \tab  50 \tab 40 \tab 10   \tab 40  \cr
    \tab Ingredient 2  \tab 200 \tab 150 \tab 10 \tab  -   \tab 60  \cr
    \tab Filler        \tab -   \tab -   \tab -  \tab  -   \tab 0   \cr
  }


  The problem is to find the composition of the feeding mix that minimises
  the production costs subject to the constraints above.

  Stated otherwise: what is the optimal amount of ingredients in one kg of
  feeding mix?


  Mathematically this can be estimated by solving a linear programming problem:
  \deqn{\min(\sum {Cost_i*x_i})} subject to
  \deqn{x_i>=0}
  \deqn{Ex=f}
  \deqn{Gx>=h}

  Where the \code{Cost} (to be minimised) is given by:
  \deqn{x_1*40+x_2*60}

  The \code{equality} ensures that the sum of the three fractions equals 1:
  \deqn{1 = x_1+x_2+x_3}

  And the \code{inequalities} enforce the nutritional constraints:
  \deqn{100*x_1+200*x_2>80}
  \deqn{50*x_1+150*x_2>50}
  and so on
 
  The solution is Ingredient1 (x1) = 0.5909, Ingredient2 (x2)=0.1364 and
  Filler (x3)=0.2727.
}

\usage{Blending}

\format{
  A list with matrix \code{G} and vector \code{H} that contain the inequality
  conditions and with vector \code{Cost}, defining the cost function.

  Columnnames of \code{G} or names of \code{Cost} are the names of the
  ingredients, rownames of \code{G} and names of \code{H} are the nutrients.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>.
}

\examples{
# Generate the equality condition (sum of ingredients = 1)
E <- rep(1, 3)
F <- 1

G <- Blending$G
H <- Blending$H

# add positivity requirement
G <- rbind(G, diag(3))
H <- c(H, rep(0, 3))

# 1. Solve the model with linear programming
res <- linp(E = t(E), F = F, G = G, H = H, Cost = Blending$Cost)


# show results
print(c(res$X, Cost = res$solutionNorm))

dotchart(x = as.vector(res$X), labels = colnames(G),
         main = "Optimal blending with ranges",
         sub = "using linp and xranges", pch = 16, 
         xlim = c(0, 1))

# 2. Possible ranges of the three ingredients
(xr <- xranges(E, F, G, H))
segments(xr[,1], 1:ncol(G), xr[,2], 1:ncol(G))
legend ("topright", pch = c(16, NA), lty = c(NA, 1),
        legend = c("Minimal cost", "range"))

# 3. Random sample of the three ingredients
# The inequality that all x > 0 has to be added!
xs <- xsample(E = E, F = F, G = G, H = H)$X

pairs(xs, main = "Blending, 3000 solutions with xsample")

# Cost associated to these random samples
Costs <- as.vector(varsample(xs, EqA = Blending$Cost))
hist(Costs)
legend("topright", c("Optimal solution",
       format(res$solutionNorm, digits = 3)))
}

\seealso{
  \code{\link{linp}} to solve a linear programming problem.
}
\keyword{datasets}

