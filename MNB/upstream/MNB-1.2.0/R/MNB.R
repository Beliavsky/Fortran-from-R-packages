#' Diagnostic tools for a multivariate negative binomial model
#'
#' Diagnostic tools as residual analysis, global, local and total-local influence
#'      for the multivariate model from the random intercept Poisson-GlG mode. Including
#'      also, the estimation process by maximum likelihood and generating multivariate
#'      negative binomial data.
#'
#' @section MNB package functions:
#' \itemize{
#' \item \code{\link[MNB]{rMNB}}
#' \item \code{\link[MNB]{fit.MNB}}
#' \item \code{\link[MNB]{qMNB}}
#' \item \code{\link[MNB]{re.MNB}}
#' \item \code{\link[MNB]{envelope.MNB}}
#' \item \code{\link[MNB]{global.MNB}}
#' \item \code{\link[MNB]{local.MNB}}
#'}
#'
#' @author Jalmar M F Carrasco <carrascojalmar@gmail.com>,
#' Cristian M Villegas Lobos <master.villegas@gmail.com> and Lizandra C Fabio <lizandrafabio@gmail.com>
#' @references
#' \itemize{
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023).
#' Diagnostic tools for a multivariate negative binomial
#' model for fitting correlated data with overdispersion. Communications in
#' Statistics - Theory and Methods, 52, 1833–1853.
#' \item Fabio, L. C., Villegas, C., Mamun, A. S., and Carrasco, J. M. F. (2025).
#' Residual analysis for discrete correlated data in the multivariate approach.
#' Brazilian Journal of Biometrics, 43, e43728.
#' }
#'
#'@aliases package
#'@name MNB
#
NULL
#' Seizures data
#'
#' @description The data set described in Diggle et.al (2013) refers to an experiment in which 59 epileptic patients
#' were randomly assigned to one of two treatment groups: treatment (progabide drug) and placebo groups.
#' The number of seizures experienced by each patient during the baseline period (week eight) and the four consecutive
#' periods (every two weeks) was recorded. The main goal of this application is to analyze the drug effect with
#' respect to the placebo. Two dummies covariates are considered in this study; Group which assumes values equal
#' to 1 if the patient belongs to treatment group and 0 otherwise, and Period which assumes values equal to 1 if the
#' number of seizures are recorded during the treatment and 0 if are measured in the baseline period. It is taking
#' into account the Time covariate which represents the number of weeks required for the counting of seizures in each
#' patient of the placebo and treatment groups.
#'
#' @docType data
#' @keywords dataSets
#' @name seizures
#' @usage data(seizures)
#' @format This data frame contains the following columns:
#' \itemize{
#' \item Y: The number epileptic seizure.
#' \item trt: Treatment: binary indicators for the prograbide and placebo groups.
#' \item period: binary indicator for the baseline period.
#' \item week: number od weeks
#' \item ind: Indicator on the ith patient.
#' }
#' @references
#' \itemize{
#' \item Diggle, P. J., Liang, K. Y., and Zeger, S. L. (2013). Analysis of Longitudinal Data. Oxford
#' University Press, N.Y., 2 edition.
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023).
#' Diagnostic tools for a multivariate negative binomial
#' model for fitting correlated data with overdispersion. Communications in
#' Statistics - Theory and Methods, 52, 1833–1853.
#' \item Fabio, L. C., Villegas, C., Mamun, A. S., and Carrasco, J. M. F. (2025).
#' Residual analysis for discrete correlated data in the multivariate approach.
#' Brazilian Journal of Biometrics, 43, e43728.
#' }
#' @examples
#'
#' data(seizures)
#' head(seizures)
#'
NULL
#' Alzheimer data
#'
#' @description The Alzheimer’s data is presented in Hand and Taylor (1987) and Hand and
#' Crowder (1996) to assess deterioration aspects of intellect, self-care and
#' personality in senile patients with Alzheimer’s disease. Two groups of patients
#' were compared, one of which received a placebo and the other treatment with
#' lecithin. In the data, each of the subjects, 26 in the placebo group and 22
#' in the lecithin group, were measured on five occasions (initially, 1st, 2nd, 4th and 6th).
#' The measurements were the number of words that the patients could recalled from lists of words.
#'
#' @docType data
#' @keywords dataSets
#' @name alzheimer
#' @usage data(alzheimer)
#' @format This data frame contains the following columns:
#' \itemize{
#' \item Y: The number of words that the patients could recalled from lists of words.
#' \item trt: Placebo ano lecithin groups.
#' \item ind: Indicator on the ith patient.
#' \item time:  initially, 1st, 2nd, 4th and 6th visit.
#' }
#' @references
#' \itemize{
#' \item Hand, D. J. and Crowder, M. (1996). Practical Longitudinal Data Analysis. London: Chapman and Hall.
#' \item Hand, D. J. and Taylor, C. C. (1987). Analysis of Variance and Repeated Measures. London: Chapman and Hall.
#' \item Fabio, L. C., Villegas, C., Carrasco, J. M. F., and de Castro, M. (2023). Diagnostic tools for a multivariate negative binomial model for fitting correlated data with overdispersion. Communications in Statistics - Theory and Methods, 52, 1833–1853.
#' \item Fabio, L. C., Villegas, C., Mamun, A. S., and Carrasco, J. M. F. (2025). Residual analysis for discrete correlated data in the multivariate approach. Brazilian Journal of Biometrics, 43, e43728.
#' }
#' @examples
#'
#' data(alzheimer)
#' head(alzheimer)
#'
NULL
