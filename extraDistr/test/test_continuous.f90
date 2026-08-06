program test_continuous
  use extra_distr
  use test_support
  implicit none
  real(dp),parameter::p=0.37_dp,tol=3.0e-9_dp
  real(dp)::x

  x=qbetapr(p,2.2_dp,3.1_dp,1.4_dp);call assert_close(pbetapr(x,2.2_dp,3.1_dp,1.4_dp),p,tol,'beta prime')
  call assert_close(pbhatt(0.5_dp,0.5_dp,1.2_dp,0.7_dp),0.5_dp,1.0e-12_dp,'Bhattacharjee symmetry')
  x=qfatigue(p,0.8_dp,2.0_dp,-1.0_dp);call assert_close(pfatigue(x,0.8_dp,2.0_dp,-1.0_dp),p,tol,'fatigue')
  x=qfrechet(p,2.5_dp,-1.0_dp,3.0_dp);call assert_close(pfrechet(x,2.5_dp,-1.0_dp,3.0_dp),p,tol,'Frechet')
  x=qgev(p,0.2_dp,1.3_dp,0.15_dp);call assert_close(pgev(x,0.2_dp,1.3_dp,0.15_dp),p,tol,'GEV')
  x=qgev(p,0.2_dp,1.3_dp,0.0_dp);call assert_close(pgev(x,0.2_dp,1.3_dp,0.0_dp),p,tol,'Gumbel limit')
  x=qgompertz(p,1.2_dp,0.7_dp);call assert_close(pgompertz(x,1.2_dp,0.7_dp),p,tol,'Gompertz')
  x=qgpd(p,0.1_dp,1.4_dp,0.2_dp);call assert_close(pgpd(x,0.1_dp,1.4_dp,0.2_dp),p,tol,'GPD')
  x=qgumbel(p,-0.2_dp,2.0_dp);call assert_close(pgumbel(x,-0.2_dp,2.0_dp),p,tol,'Gumbel')
  x=qhcauchy(p,1.7_dp);call assert_close(phcauchy(x,1.7_dp),p,tol,'half Cauchy')
  x=qhnorm(p,1.7_dp);call assert_close(phnorm(x,1.7_dp),p,tol,'half normal')
  x=qht(p,6.0_dp,1.7_dp);call assert_close(pht(x,6.0_dp,1.7_dp),p,tol,'half t')
  x=qhuber(p,0.4_dp,1.2_dp,1.345_dp);call assert_close(phuber(x,0.4_dp,1.2_dp,1.345_dp),p,2.0e-8_dp,'Huber')
  x=qinvgamma(p,3.5_dp,2.0_dp);call assert_close(pinvgamma(x,3.5_dp,2.0_dp),p,tol,'inverse gamma')
  x=qinvchisq(p,8.0_dp,2.0_dp);call assert_close(pinvchisq(x,8.0_dp,2.0_dp),p,tol,'inverse chi-square')
  x=qkumar(p,2.0_dp,3.0_dp);call assert_close(pkumar(x,2.0_dp,3.0_dp),p,tol,'Kumaraswamy')
  x=qlaplace(p,0.4_dp,1.3_dp);call assert_close(plaplace(x,0.4_dp,1.3_dp),p,tol,'Laplace')
  x=qlst(p,7.0_dp,0.4_dp,1.3_dp);call assert_close(plst(x,7.0_dp,0.4_dp,1.3_dp),p,tol,'location-scale t')
  x=qlomax(p,1.1_dp,2.4_dp);call assert_close(plomax(x,1.1_dp,2.4_dp),p,tol,'Lomax')
  x=qnsbeta(p,2.0_dp,4.0_dp,-2.0_dp,3.0_dp);call assert_close(pnsbeta(x,2.0_dp,4.0_dp,-2.0_dp,3.0_dp),p,tol,'nonstandard beta')
  x=qpareto(p,2.3_dp,1.5_dp);call assert_close(ppareto(x,2.3_dp,1.5_dp),p,tol,'Pareto')
  x=qpower(p,2.3_dp,1.5_dp);call assert_close(ppower(x,2.3_dp,1.5_dp),p,tol,'power')
  x=qprop(p,10.0_dp,0.3_dp,0.5_dp);call assert_close(pprop(x,10.0_dp,0.3_dp,0.5_dp),p,tol,'proportion beta')
  x=qrayleigh(p,1.4_dp);call assert_close(prayleigh(x,1.4_dp),p,tol,'Rayleigh')
  call assert_close(psgomp(0.0_dp,1.2_dp,2.0_dp),0.0_dp,tol,'shifted Gompertz support')
  call assert_close(pslash(0.4_dp,0.4_dp,1.2_dp),0.5_dp,tol,'slash symmetry')
  x=qtriang(p,-2.0_dp,4.0_dp,1.0_dp);call assert_close(ptriang(x,-2.0_dp,4.0_dp,1.0_dp),p,tol,'triangular')
  x=qtnorm(p,0.3_dp,1.4_dp,-1.0_dp,2.0_dp);call assert_close(ptnorm(x,0.3_dp,1.4_dp,-1.0_dp,2.0_dp),p,tol,'truncated normal')
  call assert_close(qtlambda(0.8_dp,0.0_dp),log(4.0_dp),tol,'Tukey lambda logistic')
  call assert_true(pwald(1.0_dp,1.2_dp,2.0_dp)>0.0_dp.and.pwald(1.0_dp,1.2_dp,2.0_dp)<1.0_dp,'Wald CDF')

  call assert_true(dbetapr(1.0_dp,2.0_dp,3.0_dp)>0.0_dp,'positive beta-prime density')
  call assert_true(dfatigue(2.0_dp,0.8_dp)>0.0_dp,'positive fatigue density')
  call assert_true(dgev(0.0_dp)>0.0_dp,'positive GEV density')
  call assert_true(dhuber(0.0_dp)>0.0_dp,'positive Huber density')
  call assert_true(dwald(1.0_dp,1.2_dp,2.0_dp)>0.0_dp,'positive Wald density')

  call finish_tests('continuous distributions')
end program test_continuous
