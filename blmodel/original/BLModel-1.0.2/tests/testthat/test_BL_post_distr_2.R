library(mvtnorm)

generate_data_ts <- function (means, covmat, num, PARAM)
{
  k = ncol(covmat)
  dfn = PARAM[1]
  D_disp = (dfn-2)/dfn*covmat

  sim_data =  rmvt (n=num, sigma=D_disp, df=dfn, delta=means, type="shifted")
  sim_data = matrix(num, k, data = sim_data)
  prob = matrix(1/num,num,1)
  mod_returns = cbind(sim_data, prob)
  return (mod_returns)
}


prepare_data_sim <- function ()
{


  w_m = c(0.146486521694715,0.183398991248777,0.387578736989452,0.0885946553732219,0.193941094693835)


  sample_cov = matrix(5,5, data = c(0.021204189795595,0.0141727898251296,0.0197603481450792,0.0206464173391891,0.0216768443446113,
                                    0.0141727898251296,0.0148456668040403,0.015499182880108,0.0138503998073033,0.0171510988328791,
                                    0.0197603481450792,0.015499182880108,0.032850023408191,0.0211583075709909,0.0234021692546761,
                                    0.0206464173391891,0.0138503998073033,0.0211583075709909,0.0351503966044722,0.0219059406230902,
                                    0.0216768443446113,0.0171510988328791,0.0234021692546761,0.0219059406230902,0.0312387894713628))

  sample_mean = c(0.09997125, 0.09276261, 0.19931271, 0.13299774, 0.11639217)
  q = c(-0.0177725,0.1073775,-0.4624975,0.2648525,0.0862525)


  return(list(w_m = w_m, lambda =0, sample_mean = sample_mean, sample_cov = sample_cov, q = q))

}

options(warn = -1)

data_nA = prepare_data_sim()
sample_cov = data_nA$sample_cov
sam_mean =  data_nA$sample_mean
k = ncol(sample_cov)

qe = (data_nA$q)
w_m = data_nA$w_m
Pe = .make_diag(rep(1,k))

alpha = 0.95  # quantile for CVaR inverse optimization 
tau = 0.05 # Black-Litterman tau
SR = 0.5 # market Sharpe ratio


num = 1000 # number of simulations for for market distribution

 

returns_coef = 1
returns = generate_data_ts(sam_mean/returns_coef, sample_cov/returns_coef, num, c(9,0)) # other data generators may be used

ret = BL_post_distr(returns, returns_coef, NULL, w_m, SR, Pe, qe, tau, "MAD", alpha, views_distr = observ_powerexp, "diag", sample_cov)

##########################
test_that("Output type",{
  expect_is(BL_post_distr(returns, returns_coef, NULL, w_m, SR, Pe, qe, tau, "MAD", alpha, views_distr = observ_powerexp, "diag", sample_cov),
            "list" )
})
##############################

##########################
test_that("Output type",{
  expect_error(BL_post_distr(returns, returns_coef, NULL, w_m, SR, Pe, qe, tau, "MAD", alpha, views_distr = observ_powerexp, "diag", sample_cov[-1,]) )
})
##############################
