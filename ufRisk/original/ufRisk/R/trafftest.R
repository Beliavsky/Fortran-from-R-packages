#' Backtesting of Value-at-Risk and Expected Shortfall via Traffic Light Tests
#'
#' Backtesting methods, most importantly traffic light tests, are applied to
#' previously calculated Value-at-Risk and Expected Shortfall series.
#'
#'
#' @param obj A list returned by the \code{varcast} function, that contains
#' different estimated Value-at-Risk and Expected Shortfall series; any
#' other list that follows the name conventions of the \code{varcast} function
#' can be used as well.
#'
#' @importFrom stats 'pnorm' 'pbinom' 'pt'
#'
#' @export
#'
#' @details
#' The Traffic Light Test for backtesting the Value-at-Risk (VaR) was proposed
#' by the Basel Committee on Banking Supervision (1996). A formal mathematical
#' description was given by Constanzino and Curran (2018). Following
#' Constanzino and Curran (2018), define the Value-at-Risk breach indicator by
#'
#' \deqn{X_{VaR}^{(i)}(\alpha) = 1_{ \{L_i \geq VaR_i(\alpha)\} },}
#'
#' where \eqn{i} defines the corresponding trading day, \eqn{L_i} is the loss
#' (denoted as a positive value) on day \eqn{i} and \eqn{\alpha} is the confidence
#' level of the VaR (e.g. if \eqn{\alpha = 0.95}, the 95\%-VaR is considered). The
#' total number of breaches over all trading days \eqn{i = 1, 2, ..., N} is then
#' given by
#'
#' \deqn{X_{VaR}^{N}(\alpha) = \sum_{i=1}^{N} 1_{\{L_i \geq VaR_i(\alpha)\}}.}
#'
#' Following a Binomial Distribution, the cumulative probabilities of observing
#' a specific number of breaches or less can be computed. Under the hypothesis
#' that the selected volatility model is true, the cumulative probability of
#' observing \eqn{X_{VaR}^N (\alpha)} breaches is therefore easily obtainable. The
#' Basel Committee on Banking Supervision (1996) defined three zones. Depending
#' on the zone the calculated cumulative probability can be sorted into, the
#' suitability of the selected model can be assessed. Models with calculated
#' cumulative probabilities < 95\% belong to the green zone and are considered
#' appropriate. Furthermore, if the probabilities are greater or equal to 95\%
#' but smaller than 99.99\%, the corresponding models are categorized into the
#' yellow zone. The red zone is for models with cumulative probabilities greater
#' or equal to 99.99\%. If the test results in a yellow zone classification, the
#' respective VaR values require additional monitoring. Moreover, the Basel
#' Committee recommended to consider additional capital requirements of a bank,
#' if its model used is in the yellow zone. Models in the red zone are
#' considered to be heavily flawed.
#'
#' Based on the same three-zone approach with the same zone boundaries,
#' Constanzino and Curran (2018) proposed a traffic light test for the Expected
#' Shortfall (ES). The total severity of breaches is given by
#'
#' \deqn{X_{ES}^N(\alpha) = \sum_{i=1}^N(1 - (1 - F(L_i))/(1 - \alpha))
#'  * 1_{\{L_i \geq VaR_i(\alpha)\}},}
#'
#' with \eqn{F(L_i)} being the cumulative distribution of the loss at day \eqn{i}. As stated
#' by Constanzino and Curran (2018), \eqn{X_{ES}^N(\alpha)} is approximately normally
#' distributed \eqn{\mathcal{N}(\mu_{ES}}, \eqn{N \sigma_{ES}^2)} for large samples, where
#' \eqn{\mu_{ES} = 0.5(1 - \alpha)N} and
#' \eqn{\sigma_{ES}^2 = (1 - \alpha)(4 - 3(1 - \alpha)) / 12}, from which cumulative
#' probabilities for the observed breaches \eqn{X_{ES}^N} can be easily obtained.
#'
#' For semiparametric models, the backtesting of the VaR is analogous to the
#' described approach. Backtesting the ES, however, requires minor adjustments.
#' Given that the model's underlying innovations follow a standardized
#' t-distribution with degrees of freedom \eqn{\nu}, define by \eqn{r_t} the demeaned
#' returns and by \eqn{\hat{s}_t} the estimated total volatility.
#'
#' \deqn{\hat{\epsilon}_t^* = -r_t / \hat{s}_t \sqrt{\nu / (\nu - 2)}}
#'
#' are now suitable to calculate the total severity of breaches under the
#' assumption that \eqn{\epsilon_t^*} are identically and independently
#' distributed t-distributed random variables.
#'
#' This function uses an object returned by the \code{varcast} function
#' of the \code{ufRisk} package as an input for the
#' function argument \code{obj}. A list with different elements, such as
#' the cumulative probabilities for the VaR and ES series within \code{obj},
#' is returned. Instead of the list, only the traffic light backtesting results
#' are printed to the R console.
#'
#' NOTE:
#'
#' More information on VaR and ES can be found in the documentation of the
#' \code{varcast} function of the \code{ufRisk} package
#' \code{\link{varcast}}.
#'
#' @return A list of class \code{ufRisk} is returned with the following elements.
#' \describe{
#' \item{model}{selected model for estimation}
#' \item{p_VaR.e}{cumulative probability of observing the number of
#' breaches or fewer for (1 - \code{a.e})100\%-VaR}
#' \item{p_VaR.v}{cumulative probability of observing the number of
#' breaches or fewer for (1 - \code{a.v})100\%-VaR}
#' \item{p_ES}{cumulative probability of observing the number of
#' breaches or fewer for (1 - \code{a.e})100\%-ES}
#' \item{pot_VaR.e}{number of exceedances for (1 - \code{a.e})100\%-VaR}
#' \item{pot_VaR.v}{number of exceedances for (1 - \code{a.v})100\%-VaR}
#' \item{potES}{number of exceedances for (1 - \code{a.e})100\%-ES}
#' \item{br.sum}{sum of breaches for (1 - \code{a.e})100\%-ES}
#' \item{WAD}{weighted absolute deviations - a model selection criterion}
#' \item{a.v}{coverage level for the (1-\code{a.v})100\% VaR}
#' \item{a.e}{coverage level for (1-\code{a.e})100\% VaR}
#' }
#'
#' @author
#'\itemize{
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'}
#'
#' @references
#' Basel Committee on Banking Supervision (1996). Supervisory Framework For The
#' Use of Back-Testing in Conjunction With The Internal Models Approach to
#' Market Risk Capital Requirements.
#' Available online: \url{https://www.bis.org/publ/bcbs22.htm} (accessed on 23 June
#' 2020).
#'
#' Constanzino, N., and Curran, M. (2018). A Simple Traffic Light Approach to
#' Backtesting Expected Shortfall. In: Risks 6.1.2.
#'
#' @examples
#'
#'# Example for Walmart Inc. (WMT)
#' prices = WMT$price.close
#' output = varcast(prices)
#' trafftest(output)
#'


trafftest <- function(obj) {

    ret.out <- obj$ret.out
    loss <- -ret.out
    sig.fc <- obj$sig.fc
    VaR.e <- obj$VaR.e
    VaR.v <- obj$VaR.v
    ES <- obj$ES
    dfree <- obj$dfree
    if (is.na(dfree)) sdev <- 1
    if (!is.na(dfree)) sdev <- sqrt(dfree / (dfree - 2))
    a.v <- obj$a.v
    a.e <- obj$a.e
    model <- obj$model
    n.out <- length(ret.out)
    mean.ret <- obj$mean
    if (is.null(mean.ret)) mean.ret <- 0

    br_ind <- which(VaR.e < loss)
    loss_br <- -(ret.out - mean.ret) / sig.fc
    loss_br <- loss_br[br_ind] * sdev
    if (is.na(dfree)) br <- 1 - (1 - pnorm(loss_br))/a.e
    if (!is.na(dfree)) br <- 1 - (1 - pt(loss_br, dfree))/a.e
    br.sum <- sum(br)
    pot_VaR.e <- length(VaR.e[VaR.e < loss])
    pot_VaR.v <- length(VaR.v[VaR.v < loss])
    potES <- length(ES[ES < loss])
    p_VaR.e <- pbinom(pot_VaR.e, n.out, a.e)
    p_VaR.v <- pbinom(pot_VaR.v, n.out, a.v)

    mean.e <- n.out * a.e
    mean.v <- n.out *a.v
    mean.br <- 0.5 * a.e * n.out
    p_ES <- pnorm(br.sum, mean = mean.br, sd = sqrt(a.e * (4 - 3 * a.e)/12
                                                    * n.out))

    WAD <- abs(pot_VaR.e - mean.e)/mean.e + abs(pot_VaR.v - mean.v)/mean.v +
            abs(br.sum - mean.br)/mean.br

    results <- list(model = model, p_VaR.v = p_VaR.v, p_VaR.e = p_VaR.e,
                    p_ES = p_ES, pot_VaR.e = pot_VaR.e, pot_VaR.v = pot_VaR.v,
                    potES = potES, br.sum = br.sum, WAD = WAD, a.v = a.v,
                    a.e = a.e, garchOrder = obj[["garchOrder"]])
    class(results) <- "ufRisk"
    attr(results, "function") <- "trafftest"
    results
}



