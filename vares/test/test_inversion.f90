! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
program test_inversion
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use vares
  implicit none
  integer :: failures
  failures = 0
  call check_close('varexponential', pexponential(varexponential(0.37_dp, lambda=1.3_dp), lambda=1.3_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varkumexp', pkumexp(varkumexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varexpexp', pexpexp(varexpexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), lambda=1.3_dp, &
    & a=0.80000000000000004_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varinvexpexp', pinvexpexp(varinvexpexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & lambda=1.3_dp, a=0.80000000000000004_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetaexp', pbetaexp(varbetaexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varlogisexp', plogisexp(varlogisexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & lambda=1.3_dp, a=0.80000000000000004_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varexpext', pexpext(varexpext(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), lambda=1.3_dp, &
    & a=0.80000000000000004_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varmoexp', pmoexp(varmoexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), lambda=1.3_dp, &
    & a=0.80000000000000004_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varperks', pperks(varperks(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbeard', pbeard(varbeard(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & rho=0.59999999999999998_dp), a=0.80000000000000004_dp, b=1.7_dp, rho=0.59999999999999998_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('vargompertz', pgompertz(vargompertz(0.37_dp, b=1.7_dp, eta=1.1000000000000001_dp), b=1.7_dp, &
    & eta=1.1000000000000001_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetagompertz', pbetagompertz(varbetagompertz(0.37_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp, &
    & eta=1.1000000000000001_dp), b=1.7_dp, c=1.2_dp, d=1.5_dp, eta=1.1000000000000001_dp), 0.37_dp, 3.0e-7_dp, &
    & failures)
  call check_close('varlfr', plfr(varlfr(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), a=0.80000000000000004_dp, &
    & b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varpareto', ppareto(varpareto(0.37_dp, k=1.0_dp, c=1.2_dp), k=1.0_dp, c=1.2_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varkumpareto', pkumpareto(varkumpareto(0.37_dp, k=1.0_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & c=1.2_dp), k=1.0_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varF', pf(varf(0.37_dp, d1=5.0_dp, d2=8.0_dp), d1=5.0_dp, d2=8.0_dp), 0.37_dp, 3.0e-7_dp, &
    & failures)
  call check_close('vargenpareto', pgenpareto(vargenpareto(0.37_dp, k=0.40000000000000002_dp, c=1.2_dp), &
    & k=0.40000000000000002_dp, c=1.2_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetapareto', pbetapareto(varbetapareto(0.37_dp, k=1.0_dp, a=0.80000000000000004_dp, &
    & c=1.2_dp, d=1.5_dp), k=1.0_dp, a=0.80000000000000004_dp, c=1.2_dp, d=1.5_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varparetostable', pparetostable(varparetostable(0.37_dp, lambda=1.3_dp, nu=4.5_dp, &
    & sigma=1.3_dp), lambda=1.3_dp, nu=4.5_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargamma', pgamma(vargamma(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varkumgamma', pkumgamma(varkumgamma(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varnakagami', pnakagami(varnakagami(0.37_dp, m=2.0_dp, a=0.80000000000000004_dp), m=2.0_dp, &
    & a=0.80000000000000004_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varrgamma', prgamma(varrgamma(0.37_dp, a=0.80000000000000004_dp, theta=0.59999999999999998_dp, &
    & phi=1.2_dp), a=0.80000000000000004_dp, theta=0.59999999999999998_dp, phi=1.2_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varclg', pclg(varclg(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & theta=0.59999999999999998_dp), a=0.80000000000000004_dp, b=1.7_dp, theta=0.59999999999999998_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varloggamma', ploggamma(varloggamma(0.37_dp, a=0.80000000000000004_dp, &
    & r=1.3999999999999999_dp), a=0.80000000000000004_dp, r=1.3999999999999999_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varinvgamma', pinvgamma(varinvgamma(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varstacygamma', pstacygamma(varstacygamma(0.37_dp, gamma=0.69999999999999996_dp, c=1.2_dp, &
    & theta=0.59999999999999998_dp), gamma=0.69999999999999996_dp, c=1.2_dp, theta=0.59999999999999998_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetadist', pbetadist(varbetadist(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varuniform', puniform(varuniform(0.37_dp, a=0.0_dp, b=2.0_dp), a=0.0_dp, b=2.0_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('vargenunif', pgenunif(vargenunif(0.37_dp, a=0.0_dp, c=1.0_dp, h=0.80000000000000004_dp, &
    & k=1.2_dp), a=0.0_dp, c=1.0_dp, h=0.80000000000000004_dp, k=1.2_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varpower1', ppower1(varpower1(0.37_dp, a=0.80000000000000004_dp), a=0.80000000000000004_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varpower2', ppower2(varpower2(0.37_dp, b=1.7_dp), b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varlogbeta', plogbeta(varlogbeta(0.37_dp, a=0.5_dp, b=1.5_dp, c=1.2_dp, d=2.0_dp), a=0.5_dp, &
    & b=1.5_dp, c=1.2_dp, d=2.0_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varcompbeta', pcompbeta(varcompbeta(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varLNbeta', plnbeta(varlnbeta(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varMRbeta', pmrbeta(varmrbeta(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & r=1.3999999999999999_dp, q=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, r=1.3999999999999999_dp, q=1.3_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargenbeta', pgenbeta(vargenbeta(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vararcsine', parcsine(vararcsine(0.37_dp, a=0.0_dp, b=2.0_dp), a=0.0_dp, b=2.0_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('vartriangular', ptriangular(vartriangular(0.37_dp, a=0.0_dp, b=2.0_dp, &
    & c=0.80000000000000004_dp), a=0.0_dp, b=2.0_dp, c=0.80000000000000004_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargenbeta2', pgenbeta2(vargenbeta2(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varinvbeta', pinvbeta(varinvbeta(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargeninvbeta', pgeninvbeta(vargeninvbeta(0.37_dp, a=0.80000000000000004_dp, c=1.2_dp, &
    & d=1.5_dp), a=0.80000000000000004_dp, c=1.2_dp, d=1.5_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vartsp', ptsp(vartsp(0.37_dp, a=1.3_dp, theta=0.40000000000000002_dp), a=1.3_dp, &
    & theta=0.40000000000000002_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varkum', pkum(varkum(0.37_dp, a=1.2_dp, b=1.5_dp), a=1.2_dp, b=1.5_dp), 0.37_dp, 3.0e-7_dp, &
    & failures)
  call check_close('varnormal', pnormal(varnormal(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varkumnormal', pkumnormal(varkumnormal(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp, b=1.7_dp), mu=0.20000000000000001_dp, sigma=1.3_dp, a=0.80000000000000004_dp, &
    & b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varexppower', pexppower(varexppower(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp), mu=0.20000000000000001_dp, sigma=1.3_dp, a=0.80000000000000004_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varaep', paep(varaep(0.37_dp, q1=1.3999999999999999_dp, q2=1.8_dp, &
    & alpha=0.65000000000000002_dp), q1=1.3999999999999999_dp, q2=1.8_dp, alpha=0.65000000000000002_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varbetanorm', pbetanorm(varbetanorm(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp, b=1.7_dp), mu=0.20000000000000001_dp, sigma=1.3_dp, a=0.80000000000000004_dp, &
    & b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varhalfnorm', phalfnorm(varhalfnorm(0.37_dp, sigma=1.3_dp), sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, &
    & failures)
  call check_close('varkumhalfnorm', pkumhalfnorm(varkumhalfnorm(0.37_dp, sigma=1.3_dp, a=0.80000000000000004_dp, &
    & b=1.7_dp), sigma=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varT', pt(vart(0.37_dp, n=1.0_dp), n=1.0_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varast', past(varast(0.37_dp, nu1=4.5_dp, nu2=6.0_dp, alpha=0.65000000000000002_dp), &
    & nu1=4.5_dp, nu2=6.0_dp, alpha=0.65000000000000002_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varhalfT', phalft(varhalft(0.37_dp, n=1.0_dp), n=1.0_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varCauchy', pcauchy(varcauchy(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varlogcauchy', plogcauchy(varlogcauchy(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varhalfcauchy', phalfcauchy(varhalfcauchy(0.37_dp, sigma=1.3_dp), sigma=1.3_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varlaplace', plaplace(varlaplace(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varPCTAlaplace', ppctalaplace(varpctalaplace(0.37_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp), a=0.80000000000000004_dp, theta=0.59999999999999998_dp), 0.37_dp, 3.0e-7_dp, &
    & failures)
  call check_close('varHBlaplace', phblaplace(varhblaplace(0.37_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp, phi=1.2_dp), a=0.80000000000000004_dp, theta=0.59999999999999998_dp, &
    & phi=1.2_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varMlaplace', pmlaplace(varmlaplace(0.37_dp, theta=0.59999999999999998_dp, phi=1.2_dp, &
    & psi=0.90000000000000002_dp), theta=0.59999999999999998_dp, phi=1.2_dp, psi=0.90000000000000002_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varloglaplace', ploglaplace(varloglaplace(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & delta=1.0_dp), a=0.80000000000000004_dp, b=1.7_dp, delta=1.0_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varasylaplace', pasylaplace(varasylaplace(0.37_dp, tau=1.2_dp, kappa=1.1000000000000001_dp, &
    & theta=0.59999999999999998_dp), tau=1.2_dp, kappa=1.1000000000000001_dp, theta=0.59999999999999998_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varasypower', pasypower(varasypower(0.37_dp, a=0.80000000000000004_dp, lambda=1.3_dp, &
    & delta=1.0_dp), a=0.80000000000000004_dp, lambda=1.3_dp, delta=1.0_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varlogistic', plogistic(varlogistic(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varsecant', psecant(varsecant(0.37_dp)), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargenlogis', pgenlogis(vargenlogis(0.37_dp, a=0.80000000000000004_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), a=0.80000000000000004_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargenlogis3', pgenlogis3(vargenlogis3(0.37_dp, alpha=0.65000000000000002_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), alpha=0.65000000000000002_dp, mu=0.20000000000000001_dp, &
    & sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargenlogis4', pgenlogis4(vargenlogis4(0.37_dp, a=0.80000000000000004_dp, &
    & alpha=0.65000000000000002_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), a=0.80000000000000004_dp, &
    & alpha=0.65000000000000002_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varhalflogis', phalflogis(varhalflogis(0.37_dp, lambda=1.3_dp), lambda=1.3_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varloglogis', ploglogis(varloglogis(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varkumloglogis', pkumloglogis(varkumloglogis(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, beta=1.3999999999999999_dp), a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, beta=1.3999999999999999_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varexplogis', pexplogis(varexplogis(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varHlogis', phlogis(varhlogis(0.37_dp, k=0.40000000000000002_dp), k=0.40000000000000002_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varlognorm', plognorm(varlognorm(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetalognorm', pbetalognorm(varbetalognorm(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, mu=0.20000000000000001_dp, &
    & sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varburr', pburr(varburr(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), a=0.80000000000000004_dp, &
    & b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetaburr', pbetaburr(varbetaburr(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varburr7', pburr7(varburr7(0.37_dp, k=0.40000000000000002_dp, c=1.2_dp), &
    & k=0.40000000000000002_dp, c=1.2_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varkumburr7', pkumburr7(varkumburr7(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & k=0.40000000000000002_dp, c=1.2_dp), a=0.80000000000000004_dp, b=1.7_dp, k=0.40000000000000002_dp, c=1.2_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetaburr7', pbetaburr7(varbetaburr7(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & k=0.40000000000000002_dp), a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, k=0.40000000000000002_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('vardagum', pdagum(vardagum(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varlomax', plomax(varlomax(0.37_dp, a=0.80000000000000004_dp, lambda=1.3_dp), &
    & a=0.80000000000000004_dp, lambda=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetalomax', pbetalomax(varbetalomax(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, lambda=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, lambda=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargumbel', pgumbel(vargumbel(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varkumgumbel', pkumgumbel(varkumgumbel(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, mu=0.20000000000000001_dp, &
    & sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetagumbel', pbetagumbel(varbetagumbel(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, mu=0.20000000000000001_dp, &
    & sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargumbel2', pgumbel2(vargumbel2(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetagumbel2', pbetagumbel2(varbetagumbel2(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & c=1.2_dp, d=1.5_dp), a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varfrechet', pfrechet(varfrechet(0.37_dp, alpha=0.65000000000000002_dp, sigma=1.3_dp), &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetafrechet', pbetafrechet(varbetafrechet(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varWeibull', pweibull(varweibull(0.37_dp, alpha=0.65000000000000002_dp, sigma=1.3_dp), &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varkumweibull', pkumweibull(varkumweibull(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varlogisrayleigh', plogisrayleigh(varlogisrayleigh(0.37_dp, a=0.80000000000000004_dp, &
    & lambda=1.3_dp), a=0.80000000000000004_dp, lambda=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varmoweibull', pmoweibull(varmoweibull(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & lambda=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, lambda=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varbetaweibull', pbetaweibull(varbetaweibull(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vardweibull', pdweibull(vardweibull(0.37_dp, c=1.2_dp, mu=0.20000000000000001_dp, &
    & sigma=1.3_dp), c=1.2_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varexpweibull', pexpweibull(varexpweibull(0.37_dp, a=0.80000000000000004_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), a=0.80000000000000004_dp, alpha=0.65000000000000002_dp, &
    & sigma=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargenpowerweibull', pgenpowerweibull(vargenpowerweibull(0.37_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp), a=0.80000000000000004_dp, theta=0.59999999999999998_dp), 0.37_dp, 3.0e-7_dp, &
    & failures)
  call check_close('varchen', pchen(varchen(0.37_dp, b=1.7_dp, lambda=1.3_dp), b=1.7_dp, lambda=1.3_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varxie', pxie(varxie(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, lambda=1.3_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp, lambda=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varloglog', ploglog(varloglog(0.37_dp, a=0.80000000000000004_dp, lambda=1.3_dp), &
    & a=0.80000000000000004_dp, lambda=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varexplog', pexplog(varexplog(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varexpgeo', pexpgeo(varexpgeo(0.37_dp, theta=0.59999999999999998_dp, lambda=1.3_dp), &
    & theta=0.59999999999999998_dp, lambda=1.3_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varexppois', pexppois(varexppois(0.37_dp, b=1.7_dp, lambda=1.3_dp), b=1.7_dp, lambda=1.3_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varTL2', ptl2(vartl2(0.37_dp, b=1.7_dp), b=1.7_dp), 0.37_dp, 3.0e-7_dp, failures)
  call check_close('varquad', pquad(varquad(0.37_dp, a=0.0_dp, b=2.0_dp), a=0.0_dp, b=2.0_dp), 0.37_dp, 3.0e-7_dp, &
    & failures)
  call check_close('varschabe', pschabe(varschabe(0.37_dp, gamma=0.69999999999999996_dp, &
    & theta=0.59999999999999998_dp), gamma=0.69999999999999996_dp, theta=0.59999999999999998_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  call check_close('varBS', pbs(varbs(0.37_dp, gamma=0.69999999999999996_dp), gamma=0.69999999999999996_dp), &
    & 0.37_dp, 3.0e-7_dp, failures)
  call check_close('vargev', pgev(vargev(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & xi=0.20000000000000001_dp), mu=0.20000000000000001_dp, sigma=1.3_dp, xi=0.20000000000000001_dp), 0.37_dp, &
    & 3.0e-7_dp, failures)
  if (failures /= 0) error stop 'VaRES inversion tests failed'
  print '(a,i0,a)', 'VaRES CDF/quantile inversion checks passed: ', 110, '.'
contains
  subroutine check_close(label, actual, expected, tol, failures)
    character(*), intent(in) :: label
    real(dp), intent(in) :: actual, expected, tol
    integer, intent(inout) :: failures
    if (.not. ieee_is_finite(actual) .or. &
        abs(actual-expected) > tol*(1.0_dp+abs(expected))) then
      failures = failures + 1
      print '(a,2es24.14)', trim(label)//' failed: ', actual, expected
    end if
  end subroutine check_close
end program test_inversion
