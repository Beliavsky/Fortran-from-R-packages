\name{glht}
\alias{glht}
\alias{glht.matrix}
\alias{glht.character}
\alias{glht.expression}
\alias{glht.mcp}
\alias{glht.mlf}
\alias{mcp}
\title{ General Linear Hypotheses }
\description{
  General linear hypotheses and multiple comparisons 
  for parametric models, including generalized linear models, 
  linear mixed effects models, and survival models.
}
\usage{
glht(model, linfct, ...)
\method{glht}{matrix}(model, linfct, 
    alternative = c("two.sided", "less", "greater"), 
    rhs = 0, ...)
\method{glht}{character}(model, linfct, ...)
\method{glht}{expression}(model, linfct, ...)
\method{glht}{mcp}(model, linfct, ...)
\method{glht}{mlf}(model, linfct, ...)
mcp(..., interaction_average = FALSE, covariate_average = FALSE)
}
\arguments{
  \item{model}{ a fitted model, 
               for example an object returned by \code{\link{lm}}, 
               \code{\link{glm}}, or \code{\link{aov}} etc. It is 
               assumed that \code{\link{coef}} and 
               \code{\link{vcov}} methods are available for 
               \code{model}. For multiple comparisons of
               means, methods \code{\link{model.matrix}}, 
               \code{\link{model.frame}} and \code{\link{terms}}
               are expected to be available for \code{model} as well.}
  \item{linfct}{ a specification of the linear hypotheses to be tested. 
            Linear functions can be specified by either the matrix
            of coefficients or by symbolic descriptions of 
            one or more linear hypotheses. Multiple comparisons
            in AN(C)OVA models are specified by objects returned from
            function \code{mcp}.}.
  \item{alternative}{ a character string specifying the alternative 
      hypothesis, must be one of '"two.sided"' (default), '"greater"' or
      '"less"'.  You can specify just the initial letter.}
  \item{rhs}{ an optional numeric vector specifying the 
              right hand side of the hypothesis.}
  \item{interaction_average}{logical indicating if comparisons are
      averaging over interaction terms. Experimental!}
  \item{covariate_average}{logical indicating if comparisons are
      averaging over additional covariates. Experimental!}
  \item{\dots}{ additional arguments to function \code{\link{modelparm}} in all
                \code{glht} methods. For function \code{mcp},
                multiple comparisons are defined by 
                matrices or symbolic descriptions specifying contrasts
                of factor levels where the arguments correspond
                to factor names.}
}
\details{

  A general linear hypothesis refers to null hypotheses of the form
  \eqn{H_0: K \theta = m} for some parametric model
  \code{model} with parameter estimates \code{coef(model)}. 

  The null hypothesis is specified by a linear function \eqn{K \theta}, 
  the direction of the alternative and the right hand side \eqn{m}.
  Here, \code{alternative} equal to \code{"two.sided"} refers to 
  a null hypothesis \eqn{H_0: K \theta = m}, whereas
  \code{"less"} corresponds to \eqn{H_0: K \theta \ge m} and  
  \code{"greater"} refers to 
  \eqn{H_0: K \theta \le m}. The right hand side vector \eqn{m} can be defined
  via the \code{rhs} argument.

  The generic method \code{glht} dispatches on its second argument
  (\code{linfct}). There are three ways, and thus methods, 
  to specify linear functions to be tested:

  1) The \code{matrix} of coefficients \eqn{K} can be specified directly
  via the \code{linfct} argument. In this case,
  the number of columns of this matrix needs to correspond to the number of
  parameters estimated by \code{model}. It is assumed that
  appropriate \code{coef} and \code{vcov} methods are available
  for \code{model} (\code{\link{modelparm}} deals with some exceptions). 

  2) A symbolic description,
  either a \code{character} or \code{expression} vector passed to \code{glht}
  via its \code{linfct} argument, can be used to define
  the null hypothesis. A symbolic description must be interpretable as a valid
  R expression consisting of both the left and right hand side 
  of a linear hypothesis.
  Only the names of \code{coef(model)} must be used as variable
  names. The alternative is given by
  the direction under the null hypothesis (\code{=} or \code{==}
  refer to \code{"two.sided"}, \code{<=} means
  \code{"greater"} and \code{>=} indicates 
  \code{"less"}). Numeric vectors of length one
  are valid values for the right hand side.
  
  3) Multiple comparisons of means are defined by objects
  of class \code{mcp} as returned by the \code{mcp} function.
  For each factor, which is included in \code{model} 
  as independent variable,
  a contrast matrix or a symbolic description of the contrasts
  can be specified as arguments to \code{mcp}. A symbolic
  description may be a \code{character} or \code{expression} 
  where the factor levels
  are only used as variables names. In addition,
  the \code{type} argument to the contrast generating function
  \code{\link{contrMat}} may serve as a symbolic description of 
  contrasts as well.

  4) The \code{lsm} function in package \code{lsmeans} offers a symbolic
  interface for the definition of least-squares means for factor combinations
  which is very helpful when more complex contrasts are of special interest.

  The \code{mcp} function must be used with care when defining parameters
  of interest in two-way ANOVA or ANCOVA models. Here, the definition
  of treatment differences (such as Tukey's all-pair comparisons or Dunnett's
  comparison with a control) might be problem specific. 
  Because it is impossible to determine the parameters of interest
  automatically in this case, \code{mcp} in multcomp
  version 1.0-0 and higher generates comparisons for the main effects
  only, ignoring covariates and interactions (older versions
  automatically averaged over interaction terms). A warning is given. We refer to
  \bibcitet{multcomp::hsu1996}, Chapter 7, and
  \bibcitet{multcomp::searle1971}, Chapter 7.3, for further discussions and examples on this
  issue.

  \code{glht} extracts the number of degrees of freedom
  for models of class \code{lm} (via \code{\link{modelparm}}) and the
  exact multivariate t distribution is evaluated. For all other
  models, results rely on the normal approximation. Alternatively, the 
  degrees of freedom to be used for the evaluation of multivariate
  t distributions can be given by the additional \code{df} argument to
  \code{\link{modelparm}} specified via \code{\dots}.

  \code{glht} methods return a specification of the null hypothesis
  \eqn{H_0: K \theta = m}. The value of the linear function
  \eqn{K \theta} can be extracted using the \code{\link{coef}} method and
  the corresponding covariance matrix is available from the 
  \code{\link{vcov}} method. Various simultaneous and univariate tests and 
  confidence intervals are available from \code{\link{summary.glht}}
  and \code{\link{confint.glht}} methods, respectively.

  A more detailed description of the underlying methodology is
  available from \bibcitet{hothorn2008b} and 
  \bibcitet{Bretz+Hothorn+Westfall_2010}.

}
\value{

  An object of class \code{glht}, more specifically a list with elements
  \item{model}{a fitted model, used in the call to \code{glht}}
  \item{linfct}{ the matrix of linear functions}
  \item{rhs}{ the vector of right hand side values \eqn{m}}
  \item{coef}{ the values of the linear functions}
  \item{vcov}{ the covariance matrix of the values of the linear functions}
  \item{df}{ optionally, the degrees of freedom when the exact t 
             distribution is used for inference}
  \item{alternative}{ a character string specifying the alternative
      hypothesis}
  \item{type}{ optionally, a character string giving the name of the specific procedure}
  with \code{print}, \code{\link{summary}}, 
  \code{\link{confint}}, \code{\link{coef}} and \code{\link{vcov}} 
  methods being available. When called with \code{linfct} being an
  \code{mcp} object, an additional element \code{focus} is available
  storing the names of the factors under test.

}
\references{\bibshow{*}}
\examples{

  ### multiple linear model, swiss data
  lmod <- lm(Fertility ~ ., data = swiss)

  ### test of H_0: all regression coefficients are zero 
  ### (ignore intercept)

  ### define coefficients of linear function directly
  K <- diag(length(coef(lmod)))[-1,]
  rownames(K) <- names(coef(lmod))[-1]
  K

  ### set up general linear hypothesis
  glht(lmod, linfct = K)

  ### alternatively, use a symbolic description 
  ### instead of a matrix 
  glht(lmod, linfct = c("Agriculture = 0",
                        "Examination = 0",
                        "Education = 0",
                        "Catholic = 0",
                        "Infant.Mortality = 0"))


  ### multiple comparison procedures
  ### set up a one-way ANOVA
  amod <- aov(breaks ~ tension, data = warpbreaks)

  ### set up all-pair comparisons for factor `tension'
  ### using a symbolic description (`type' argument 
  ### to `contrMat()')
  glht(amod, linfct = mcp(tension = "Tukey"))

  ### alternatively, describe differences symbolically
  glht(amod, linfct = mcp(tension = c("M - L = 0", 
                                      "H - L = 0",
                                      "H - M = 0")))

  ### alternatively, define contrast matrix directly
  contr <- rbind("M - L" = c(-1, 1, 0),
                 "H - L" = c(-1, 0, 1), 
                 "H - M" = c(0, -1, 1))
  glht(amod, linfct = mcp(tension = contr))

  ### alternatively, define linear function for coef(amod)
  ### instead of contrasts for `tension'
  ### (take model contrasts and intercept into account)
  glht(amod, linfct = cbind(0, contr \%*\% contr.treatment(3)))


  ### mix of one- and two-sided alternatives
  warpbreaks.aov <- aov(breaks ~ wool + tension,
                      data = warpbreaks)

  ### contrasts for `tension'
  K <- rbind("L - M" = c( 1, -1,  0),     
             "M - L" = c(-1,  1,  0),       
             "L - H" = c( 1,  0, -1),     
             "M - H" = c( 0,  1, -1))

  warpbreaks.mc <- glht(warpbreaks.aov, 
                        linfct = mcp(tension = K),
                        alternative = "less")

  ### correlation of first two tests is -1
  cov2cor(vcov(warpbreaks.mc))

  ### use smallest of the two one-sided
  ### p-value as two-sided p-value -> 0.0232
  summary(warpbreaks.mc)

  ### more complex models: Continuous outcome logistic
  ### regression; parameters are log-odds ratios
  if (require("tram", quietly = TRUE, warn.conflicts = FALSE)) {
      confint(glht(Colr(breaks ~ wool + tension, 
                        data = warpbreaks), 
                   linfct = mcp("tension" = "Tukey")))
  }
}
\keyword{htest}
