#'Function to compute gross returns.
#'
#'Calculate gross returns from prices.
#'
#'@param x data object containing ordered price observations
#'
#'@return An object of the same class as x with the gross returns
#'
#'@export
g_ret<- function(x){
  return(sort(exp(diff(log(x)))))
}
#'
#' Function to find the most-left distribution set.
#'
#' An adaptive clustering algorithm identifies sub-groups of gross returns at each iteration by approximating their distribution with a
#' sequence of two-component log-normal mixtures.
#'
#' @param y vector or data frame
#' @param shift double
#' @param n_it integer
#' @param plot option
#'
#' @return data object
#'
#' @export
#'
#'
tail_mixture<-function(y,shift, n_it, plot){
  vec=which(y>shift)
  y=y[vec]-shift
  GMModel=normalmixEM(log(y), lambda = NULL, mu = NULL, sigma = NULL, k = 2,  mean.constr = NULL, sd.constr = NULL, epsilon = 1e-08, maxit = 1000, maxrestarts=500,  verb = FALSE, fast=FALSE, ECM = FALSE, arbmean = TRUE, arbvar = TRUE)
  zsol=numeric()
  zsol[1]=GMModel$mu[1]
  zsol[3]=GMModel$mu[2]
  zsol[2]=GMModel$sigma[1]
  zsol[4]=GMModel$sigma[2]
  pa=GMModel$lambda[1]
  h1=pa*dlnorm(min(log(y)),zsol[1],zsol[2])
  h2=(1-pa)*dlnorm(min(log(y)), zsol[3],zsol[4])
  indec=as.numeric(which.max(c(h1,h2)))
  if(indec==2){
    mu1=zsol[3]
    sigma1=zsol[4]
    ppa=(1-pa)
    mu2=zsol[1]
    sigma2=zsol[2]
  }
  else {mu1=zsol[1]
  sigma1=zsol[2]
  ppa=pa
  mu2=zsol[3]
  sigma2=zsol[4]}
  delta=(mu2*sigma1^2-mu1*sigma2^2)^2+(mu2^2*sigma1^2-mu1^2*sigma2^2+2*sigma1^2*sigma2^2*log(sigma2*ppa/(sigma1*(1-ppa))))*(sigma2^2-sigma1^2)
  if(delta<=0) {
    if(n_it==1){
      flag=1
      a=0
      inters=1
      inters1=1
      ppa=1
      stop("no subpopulations", "\n")
    }
    else{
      flag=1
      a=0
      inters=1
      inters1=1
      ppa=1
      return("only one subpopulation exists", "\n")
    }}
  else {
    aa=min(exp((-(mu2*sigma1^2-mu1*sigma2^2)-sqrt(delta))/(sigma2^2-sigma1^2)),
           exp((-(mu2*sigma1^2-mu1*sigma2^2)+sqrt(delta))/(sigma2^2-sigma1^2)))
  }
  flag=0
  inters=plnorm(aa,mu2,sigma2)
  inters1=(1-plnorm(aa,mu1,sigma1))
  a=aa+shift
  x1=seq(from=min(y), to=max(y), by=(max(y)-min(y))/1000)
  mix1=dlnorm(x1,mu1, sigma1)
  mix2=dlnorm(x1, mu2,sigma2)
  kern<-density(y, n=1001)
  f=kern$y
  ff1=ppa*mix1+(1-ppa)*mix2
  if (plot == "T") {
    plot_x=plot(x1, f)
    hist(y, freq = FALSE)
    lines(x1,ff1, col='red')
    lines(x1, mix1*ppa)
    lines(x1,mix2*(1-ppa))
    points(aa, ppa*dlnorm(aa,mu1,sigma1), col="red")
    Sys.sleep(5)
    return (c(a,flag,mu1,sigma1,mu2,sigma2,ppa,inters,inters1, plot_x))}
  if (plot == "F") {
    return (c(a,flag,mu1,sigma1,mu2,sigma2,ppa,inters,inters1))}
}
#' Procedure to find the most-left distribution set.
#'
#' Estimation of the vector of unknown parameters for the density functions associated with the
#' two mixture components.
#' @param y object of class "g_ret"
#' @param plot_cp option
#' @return An object of class "infoset" is a list containing the following components for the firse two iterations (k=2):
#'  \describe{
#'  \item{change.points}{a vector of change points.}
#'  \item{prior.probability}{the a priori probabilities.}
#'  \item{first.type.errors}{the cumulative distribution functions associated with the leftmost component of the mixture.}
#'  \item{second.type.errors}{the cumulative distribution functions associated with the rightmost component of the mixture.}
#'  \item{mean}{the parameters (drift) of the left-hand component of the log-normal mixture.}
#'  \item{sd}{the parameters (volatility) of the left-hand component of the log-normal mixture.}
#' }
#'
#'@references
#'Mariani, F., Polinesi, G., Recchioni, M. C. (2022). A tail-revisited Markowitz mean-variance approach and a portfolio network centrality. Computational Management Science, 19(3), 425-455.
#'
#'Mariani, F., Ciommi, M., Chelli, F. M., Recchioni, M. C. (2020). An iterative approach to stratification: Poverty at regional level in Italy. Social Indicators Research, 1-31.
#'
#'@examples
#'
#' \donttest{gross.ret<-as.data.frame(lapply(sample.data, g_ret))
#' infoset(gross.ret$ETF_1, plot_cp = "T")}
#'
#'############################################################
#' ## EXAMPLE 1: Clustering ETFs
#'############################################################
#'
#'\donttest{
#' gross.ret<-as.data.frame(lapply(sample.data, g_ret))
#' result<-NULL
#' for(i in 1:ncol(gross.ret)){
#'  result[[i]]<-infoset(gross.ret[,i], plot_cp = "F")
#' }
#' output<-matrix(unlist(result),12,ncol=ncol(gross.ret)) # output contains the information set
#' output<-t(output)
#' rownames(output)<-colnames(gross.ret)
#' colnames(output)<-c("ch_1","ch_2","priori_1","priori_2","first_1",
#'                    "first_2","second_1","second_2","mean_1","mean_2","dev_1", "dev_2")
#' output<- as.data.frame(output)
#' group_label <- as.factor(asset.label$label)
#' d <- dist(output, method = 'euclidean')
#' hc_SIMS <- hclust(d, method = 'complete')
#' library(dendextend)
#' library(colorspace)
#' dend_SIMS <- as.dendrogram(hc_SIMS)
#' dend_SIMS <- color_branches(dend_SIMS, k = 4, col = c(1:4))
#' labels_colors(dend_SIMS) <-
#'        rainbow_hcl(5)[sort_levels_values(as.numeric(group_label)[order.dendrogram(dend_SIMS)])]
#' labels(dend_SIMS) <- paste(as.character(group_label)[order.dendrogram(dend_SIMS)],
#'        '(', labels(dend_SIMS), ')', sep = '')
#' dend_SIMS <- hang.dendrogram(dend_SIMS, hang_height = 0.001)
#' dend_SIMS <- assign_values_to_leaves_nodePar(dend_SIMS, 0.5, 'lab.cex')
#' dev.new()
#' old_par <- par(no.readonly = TRUE)
#' on.exit(par(old_par))
#' par(mar = c(1.8, 1.8, 1.8, 1))
#' plot(dend_SIMS, main = 'Complete linkage (the labels give the true ETF class)',
#'      horiz = TRUE, nodePar = list(cex = 0.007))
#' legend('topleft', legend = c('emerging equity Asia', 'emerging equity America',
#'                            'corporate bond', 'commodities', 'aggregate bond'),
#'              fill = c('#BDAB66', '#65BC8C', '#C29DDE', '#E495A5', '#55B8D0'), border = 'white')
#'}
#'
#' @export
#'
#'
infoset <- function(y, plot_cp){
  shift = NULL
  shift[1] = 0
  k = 1
  flag = 0
  meant = NULL
  dev = NULL
  meanr = NULL
  devr = NULL
  pa = NULL
  cum = NULL
  cum1 = NULL
  out = NULL
  shift = NULL
  shiftt = 0
  while(flag == 0 & shiftt < median(y) & k<3){
    if(plot_cp == "T"){
      out <- tail_mixture(y, shiftt, k, plot = "T")
    }
    if(plot_cp == "F"){
      out <- tail_mixture(y, shiftt, k, plot = "F")
    }
    k = k+1
    shiftt = out[1]
    if(shiftt>=median(y)){
      message("Not valid change point: time series too short")
    }
    shift[(k-1)] <- out[1]
    flag <- out[2]
    meant[(k-1)] <- out[3]
    dev[(k-1)] <- out[4]
    meanr[(k-1)] <- out[5]
    devr[(k-1)] <- out[6]
    pa[(k-1)] <- out[7]
    cum[(k-1)] <- out[8]
    cum1[(k-1)] <- out[9]
  }
  if(flag == 1){
    shift = shift[1:(length(shift) - 1)]
    cum = cum[1:(length(cum) - 1)]
    cum1 = cum1[1:(length(cum1) - 1)]
    meant = meant[1:(length(meant) - 1)]
    dev = dev[1:(length(dev) - 1)]
    pa = pa[1:(length(pa) - 1)]
  }
  list('change point' = shift, 'prior probability' = pa, 'first type error' = cum,
       'second type error' = cum1, 'mean' = meant, 'sd' = dev)
}
#'
#'Function to create overlapping windows.
#'
#'
#' @param data vector or data frame
#' @param FT Window size. By default set to 1290 training days (five years).
#' @param ov Number of different days for two consecutive time windows.. By default set to 125 training days (six months).
#'
#' @return a list containing the rolling windows
#'
#' @export
#'
#'
create_overlapping_windows <- function(data, FT = 1290, ov = 125) {

  # Check if FT and ov are positive numbers
  if (FT <= 0 || ov <= 0) {
    stop("Both FT (window size) and ov (overlap) must be positive numbers.")
  }

  # Check if the data has enough rows for at least one window
  if (nrow(data) < FT) {
    stop("The dataset is too small for the given window size.")
  }

  # List to store windows
  TW <- list()

  # Calculate the number of windows
  T <- floor((nrow(data) - FT) / ov)

  # Loop to create the overlapping windows
  for (t in 0:T) {
    TW[[(t + 1)]] <- data[(1 + t * ov):(FT + t * ov), ]
  }

  # Return the list of windows
  return(TW)
}
#'
#' Function to compute Left risk measure.
#'
#'
#' @param data A (T x N) matrix or data.frame containing the N time series over period T
#' @param FT Window size.
#' @param ov umber of different days for two consecutive time windows.
#'
#' @examples
#'
#' \donttest{
#' LR <- LR_cp(sample.data, FT= 1290, ov = 125)
#' df <- as.data.frame(matrix(unlist(LR), nrow = length(LR), ncol = ncol(sample.data)))
#' colnames(df) <- c(paste("tw", rep(1:16)))
#' plot(df[,1], pch=19, col=asset.label$label, ylab="LR_cp", xlab="ETFs")
#' }
#'
#' @return A (N x T) data.frame containing the LR_cp measure for the N time series over time windows
#'
#' @export
#'
#'
LR_cp<-function(data, FT, ov){
  if(is.data.frame(data)==F){
    gross.ret<-as.data.frame(g_ret(data))
    data<-as.data.frame(data)
  }
  if(is.data.frame(data)==T){
    gross.ret<-as.data.frame(lapply(data, g_ret))
  }
  result<-NULL
  for(i in 1:ncol(gross.ret)){
    result[[i]]<-infoset(gross.ret[,i], plot_cp = "F")
  }
  output<-matrix(unlist(result),12,ncol=ncol(gross.ret)) # output contains the information set
  output<-t(output)
  rownames(output)<-colnames(gross.ret)
  colnames(output)<-c("ch_1","ch_2","priori_1","priori_2","first_1",
                      "first_2","second_1","second_2","mean_1","mean_2","dev_1", "dev_2")
  output<- as.data.frame(output)
  FT = 1290
  ov= 125
  W <- create_overlapping_windows(data, FT = FT, ov = ov)
  ret <- list()
  y <- list()
  for(i in 1:ncol(data)){
    for (t in 1:length(W)){
      ret[[t]] <- matrix(0,nrow = nrow(as.data.frame(W[[1]])) - 1, ncol = ncol(data))
      y[[t]] <- matrix(0, nrow = nrow(as.data.frame(W[[1]])) - 1, ncol = ncol(data))
    }
  }
  ch <- log(output$ch_1)
  h <- list()
  diff <- list()
  n <- list()
  a_int <- list()
  hh <- list()
  LR_cp_measure <- list()
  hh_p <- list()
  p_loss <- list()
  expected_loss <- list()
  i <- ncol(data)
  for(t in 1:length(W)){
    h[[t]] <- vector('list', i)
    diff[[t]] <- vector('list', i)
    n[[t]] <- vector('list', i)
    a_int[[t]] <- vector('list', i)
    hh[[t]] <- vector('list', i)
    LR_cp_measure[[t]]<- vector('list', i)
    hh_p[[t]] <- vector('list', i)
    p_loss[[t]] <- vector('list', i)
    expected_loss[[t]] <- vector('list', i)
  }

  for(i in 1:ncol(data)){
    for (t in 1:length(W)){
      W[[t]]<-as.data.frame(W[[t]])
      ret[[t]][,i] <- diff(log(W[[t]][, i]))
      y[[t]][,i] <- exp(ret[[t]][, i])
      y[[t]][,i] <- sort(y[[t]][, i])
      h[[t]][[i]] <- hist(log(y[[t]][, i]))
      diff[[t]][[i]] <- h[[t]][[i]]$breaks-ch[[i]]
      n[[t]][[i]] <- length(diff[[t]][[i]][diff[[t]][[i]] < 0])
      a_int[[t]][[i]] <- h[[t]][[i]]$breaks[2] - h[[t]][[i]]$breaks[1]
      hh[[t]][[i]] <- h[[t]][[i]]$mids*h[[t]][[i]]$density*a_int[[t]][[i]]
      hh_p[[t]][[i]] <- h[[t]][[i]]$density*a_int[[t]][[i]]
      expected_loss[[t]][[i]] <- sum(hh[[t]][[i]][1:(n[[t]][[i]] - 1)]) +
        h[[t]][[i]]$density[n[[t]][[i]]]*(ch[[i]] - h[[t]][[i]]$breaks[n[[t]][[i]]])*0.5*
        (ch[[i]] + h[[t]][[i]]$breaks[n[[t]][[i]]])
      p_loss[[t]][[i]] <- sum(hh_p[[t]][[i]][1:(n[[t]][[i]] - 1)]) +
        h[[t]][[i]]$density[n[[t]][[i]]]*(ch[[i]] - h[[t]][[i]]$breaks[n[[t]][[i]]])
      LR_cp_measure[[t]][[i]] <- -(expected_loss[[t]][[i]]/p_loss[[t]][[i]])
    }
  }
  return( LR_cp_measure)
}
#'
#' Plot methods for a LR_cp object
#'
#' @param LR_cp_measure object of class LR_cp
#' @param asset_label vector containing asset label
#'
#' @return plot of LR_cp measures by asset classes
#'
plot_LR_cp <- function(LR_cp_measure, asset_label) {
  tw<-length(LR_cp_measure)
  ncol<-length(LR_cp_measure[[1]])
  # Convert the list to a data frame
  df <- as.data.frame(matrix(unlist(LR_cp_measure),nrow= tw, ncol=ncol))
  # Check that asset_label has a column named 'label'
  if (!"label" %in% colnames(asset_label)) {
    stop("Error: asset_label must contain a column named 'label'")
  }

  # Create the plot
  plot(
    df[,1],
    pch = 19,
    col = asset_label$label,
    ylab = "LR_cp",
    xlab = "ETFs"
  )
}
#'
#' Function to compute portfolio values
#'
#' @param data A (T x N) matrix or data.frame containing the N time series over period T
#' @param LR_cp_measure object of class LR_cp (only for "C_M" and "C_EDC" asset allocation strategies)
#' @param ptf Type of portfolio to be computed. Asset allocation strategies available are:
#'            "M" is the Markowitz portfolio, "C_M" is the combined Markowitz portfolio,
#'            "EDC" uses the extreme downside correlation and "C_EDC" is the combined  extreme downside correlation portfolio
#' @param FT Window size.
#' @param ov Overlap.
#'
#' @return An object of class "ptf_construction" is a list containing the following components for all the time windows considered:
#'  \describe{
#'   \item{ptf oos value}{a vector of out of sample returns.}
#'   \item{weigths}{portfolio weights.}
#'   }
#'
#' @export
ptf_construction<-function (data, FT, ov, LR_cp_measure, ptf = c("M", "C_M", "EDC",
                                                                 "C_EDC"))
{
  W <- create_overlapping_windows(data, FT, ov)
  ret <- list()
  y <- list()
  for (i in 1:ncol(data)) {
    for (t in 1:length(W)) {
      ret[[t]] <- matrix(0, nrow = nrow(as.data.frame(W[[1]])) -
                           1, ncol = ncol(data))
      y[[t]] <- matrix(0, nrow = nrow(as.data.frame(W[[1]])) -
                         1, ncol = ncol(data))
    }
  }
  i <- ncol(data)
  diff <- list()
  for (t in 1:length(W)) {
    diff[[t]] <- vector("list", i)
  }
  for (i in 1:ncol(data)) {
    for (t in 1:length(W)) {
      ret[[t]][, i] <- diff(log(W[[t]][, i]))
    }
  }
  r <- list()
  meanret <- list()
  COV_ret <- list()
  for (t in 1:length(W)) {
    r[[t]] <- matrix(colMeans(ret[[t]]))
    meanret[[t]] <- sum(r[[t]])/ncol(data)
    COV_ret[[t]] <- cov(ret[[t]][((FT - 1) - ((ov + 1) *
                                                2)):(FT - 1), ])
    COV_ret[[t]] <- nearPD(COV_ret[[t]])$mat
  }
  tw <- length(W) - 1
  n <- ncol(data)
  B <- list()
  f <- list()
  sol <- list()
  w <- list()
  Pvalue_ret <- NULL
  Pvalue_M <- list()
  oos.values <- list()
  if (ptf == "M") {
    for (t in 1:tw) {
      B[[t]] <- matrix(1, 1, n)
      B[[t]] <- rbind(B[[t]], t(r[[t]]), diag(n), -diag(n))
      f[[t]] <- c(1, meanret[[t]], rep(0, n), rep(-1,
                                                  n))
      sol[[t]] <- solve.QP(Dmat = COV_ret[[t]], dvec = 0 *
                             r[[t]], Amat = t(B[[t]]), bvec = f[[t]], meq = 2)
      w[[t]] <- round(sol[[t]]$solution, 6)
      Pvalue_ret[t] <- sum(w[[t]] * W[[t]][(FT - 1), ])
      Pvalue_M[[t]] <- (t(w[[t]]) %*% t(W[[t + 1]][((FT -
                                                       1) - ov):(FT - 1), ]) - Pvalue_ret[t])/Pvalue_ret[t]
      oos.values[[t]] <- Pvalue_M[[t]]
    }
  }
  if (ptf == "C_M") {
    for (t in 1:length(W)) {
      LR_cp_measure[[t]] <- matrix(unlist(LR_cp_measure[[t]]),
                                   nrow = 1, ncol = ncol(data))
    }
    B <- list()
    f <- list()
    sol <- list()
    w <- list()
    Pvalue_ret <- NULL
    Pvalue_C <- list()
    oos.values <- list()
    lambda <- 1e-04
    for (t in 1:tw) {
      B[[t]] <- matrix(1, 1, n)
      B[[t]] <- rbind(B[[t]], t(r[[t]]), diag(n), -diag(n))
      f[[t]] <- c(1, meanret[[t]], rep(0, n), rep(-1,
                                                  n))
      sol[[t]] <- solve.QP(Dmat = COV_ret[[t]], dvec = -lambda *
                             LR_cp_measure[[t]], Amat = t(B[[t]]), bvec = f[[t]],
                           meq = 2)
      w[[t]] <- round(sol[[t]]$solution, 6)
      Pvalue_ret[t] <- sum(w[[t]] * W[[t]][(FT - 1), ])
      Pvalue_C[[t]] <- (t(w[[t]]) %*% t(W[[t + 1]][((FT -
                                                       1) - ov):(FT - 1), ]) - Pvalue_ret[t])/Pvalue_ret[t]
      oos.values[[t]] <- Pvalue_C[[t]]
    }
  }
  if (ptf == "EDC") {
    ret <- list()
    y <- list()
    for (i in 1:ncol(data)) {
      for (t in 1:length(W)) {
        ret[[t]] <- matrix(0, nrow = nrow(W[[1]]) -
                             1, ncol = ncol(data))
        y[[t]] <- matrix(0, nrow = nrow(W[[1]]) - 1,
                         ncol = ncol(data))
      }
    }
    for (i in 1:ncol(data)) {
      for (t in 1:length(W)) {
        ret[[t]][, i] <- diff(log(W[[t]][, i]))
        y[[t]][, i] <- exp(ret[[t]][, i])
      }
    }
    W_out <- list()
    for (t in 1:length(W)) {
      W_out[[(t)]] = ret[[t]][((FT - 1) - ((ov + 1) *
                                             2)):(FT - 1), ]
    }
    quant <- matrix(0, nrow = length(W), ncol = ncol(data))
    diff <- list()
    for (t in 1:length(W)) {
      for (i in 1:ncol(data)) {
        quant[t, i] <- quantile(as.numeric(W_out[[(t)]][,
                                                        i]), probs = 0.05)
        diff[[t]] <- W_out[[(t)]] - colMeans(W_out[[(t)]])
        for (j in 1:length(W_out[[t]][, i])) {
          if (W_out[[(t)]][j, i] > quant[t, i]) {
            diff[[t]][j, i] = 0
          }
        }
      }
    }
    C_edc <- list()
    for (t in 1:length(W)) {
      aux <- matrix(0, nrow = ncol(data), ncol = ncol(data))
      for (i in 1:ncol(data)) {
        for (j in 1:ncol(data)) {
          aux[i, j] <- (mean(diff[[t]][, i] * diff[[t]][,
                                                        j]))
        }
      }
      C_edc[[t]] <- aux
    }
    r <- list()
    meanret <- list()
    stdev <- list()
    COV_edc <- list()
    for (t in 1:length(W_out)) {
      r[[t]] <- matrix(colMeans(ret[[t]]))
      meanret[[t]] <- sum(r[[t]])/n
      stdev[[t]] <- apply(W_out[[t]], 2, sd)
      stdev[[t]] <- matrix(stdev[[t]])
      COV_edc[[t]] <- nearPD(C_edc[[t]])$mat
    }
    tw <- length(W) - 1
    B <- list()
    f <- list()
    sol <- list()
    w <- list()
    Pvalue_edc <- NULL
    Pvalue_EDC <- list()
    oos.values <- list()
    for (t in 1:tw) {
      B[[t]] <- matrix(1, 1, n)
      B[[t]] <- rbind(B[[t]], t(r[[t]]), diag(n), -diag(n))
      f[[t]] <- c(1, meanret[[t]], rep(0, n), rep(-1,
                                                  n))
      sol[[t]] <- solve.QP(Dmat = COV_edc[[t]], dvec = 0 *
                             r[[t]], Amat = t(B[[t]]), bvec = f[[t]], meq = 2)
      w[[t]] <- round(sol[[t]]$solution, 6)
      Pvalue_edc[t] = sum(w[[t]] * W[[t]][(FT - 1), ])
      Pvalue_EDC[[t]] = (t(w[[t]]) %*% t(W[[t + 1]][((FT -
                                                        1) - ov):(FT - 1), ]) - Pvalue_edc[t])/Pvalue_edc[t]
      oos.values[[t]] <- Pvalue_EDC[[t]]
    }
  }
  if (ptf == "C_EDC") {
    for (t in 1:length(W)) {
      LR_cp_measure[[t]] <- matrix(unlist(LR_cp_measure[[t]]),
                                   nrow = 1, ncol = ncol(data))
    }
    ret <- list()
    y <- list()
    for (i in 1:ncol(data)) {
      for (t in 1:length(W)) {
        ret[[t]] <- matrix(0, nrow = nrow(W[[1]]) -
                             1, ncol = ncol(data))
        y[[t]] <- matrix(0, nrow = nrow(W[[1]]) - 1,
                         ncol = ncol(data))
      }
    }
    for (i in 1:ncol(data)) {
      for (t in 1:length(W)) {
        ret[[t]][, i] <- diff(log(W[[t]][, i]))
        y[[t]][, i] <- exp(ret[[t]][, i])
      }
    }
    W_out <- list()
    for (t in 1:length(W)) {
      W_out[[(t)]] = ret[[t]][((FT - 1) - ((ov + 1) *
                                             2)):(FT - 1), ]
    }
    quant <- matrix(0, nrow = length(W), ncol = ncol(data))
    diff <- list()
    for (t in 1:length(W)) {
      for (i in 1:ncol(data)) {
        quant[t, i] <- quantile(as.numeric(W_out[[(t)]][,
                                                        i]), probs = 0.05)
        diff[[t]] <- W_out[[(t)]] - colMeans(W_out[[(t)]])
        for (j in 1:length(W_out[[t]][, i])) {
          if (W_out[[(t)]][j, i] > quant[t, i]) {
            diff[[t]][j, i] = 0
          }
        }
      }
    }
    C_edc <- list()
    for (t in 1:length(W)) {
      aux <- matrix(0, nrow = ncol(data), ncol = ncol(data))
      for (i in 1:ncol(data)) {
        for (j in 1:ncol(data)) {
          aux[i, j] <- (mean(diff[[t]][, i] * diff[[t]][,
                                                        j]))
        }
      }
      C_edc[[t]] <- aux
    }
    r <- list()
    meanret <- list()
    stdev <- list()
    COV_edc <- list()
    for (t in 1:length(W_out)) {
      r[[t]] <- matrix(colMeans(ret[[t]]))
      meanret[[t]] <- sum(r[[t]])/n
      stdev[[t]] <- apply(W_out[[t]], 2, sd)
      stdev[[t]] <- matrix(stdev[[t]])
      COV_edc[[t]] <- nearPD(C_edc[[t]])$mat
    }
    B <- list()
    f <- list()
    sol <- list()
    w <- list()
    Pvalue_edc <- NULL
    Pvalue_mod_EDC <- list()
    lambda <- 1e-04
    oos.values <- list()
    tw <- length(W) - 1
    for (t in 1:tw) {
      B[[t]] <- matrix(1, 1, n)
      B[[t]] <- rbind(B[[t]], t(r[[t]]), diag(n), -diag(n))
      f[[t]] <- c(1, meanret[[t]], rep(0, n), rep(-1,
                                                  n))
      sol[[t]] <- solve.QP(Dmat = COV_edc[[t]], dvec = lambda *
                             LR_cp_measure[[t]], Amat = t(B[[t]]), bvec = f[[t]],
                           meq = 2)
      w[[t]] <- round(sol[[t]]$solution, 6)
      Pvalue_edc[t] = sum(w[[t]] * W[[t]][(FT - 1), ])
      Pvalue_mod_EDC[[t]] = (t(w[[t]]) %*% t(W[[t + 1]][((FT -
                                                            1) - ov):(FT - 1), ]) - Pvalue_edc[t])/Pvalue_edc[t]
      oos.values[[t]] <- Pvalue_mod_EDC[[t]]
    }
  }
  return(list(ptf.oos.value = oos.values, ptf.weights = w))
}
#'
#' Plot methods for a ptf_construction object
#'
#' @param ptf.oos.values object of class ptf_construction
#'
#' @return plot oos portfolio values
#'
#' @export
plot_ptf <- function(ptf.oos.values) {

  sample.values <- NULL
  icont <- 0
  count <- 1:length(ptf.oos.values)

  for (t in count) {
    for (j in 1:125) {
      icont <- icont + 1
      sample.values[icont] <- ptf.oos.values[[t]][1, j]
      }
  }

  # Calculate dynamic ylim range for each sample

    ylim <- range(sample.values, na.rm = TRUE)

  # Plotting

    boxplot(
    sample.values, ylim = ylim, outline = FALSE,
    main = ""
  )
  abline(h = mean(sample.values, na.rm = TRUE), col = 'red')  # Add mean line

}

#'
#' Plot methods for a ptf_construction object
#'
#' @param ptf.oos.values object of class ptf_construction
#'
#' @return summary of oos portfolio values
#'
#' @export
#'
summary_ptf <- function(ptf.oos.values) {

  sample.values <- NULL
  icont <- 0

  count <- 1:length(ptf.oos.values)
  for (t in count) {
    for (j in 1:125) {
      icont <- icont + 1
      sample.values[icont] <- ptf.oos.values[[t]][1, j]
    }
  }

  # Summary statistics
  oos_summary<-summary(sample.values)
  return(oos_summary)
}

