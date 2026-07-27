! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
program demo_lgarch
  use lgarch_kinds, only : dp
  use lgarch_rng, only : seed_rng
  use lgarch_univariate, only : lgarch_simulate,fit_lgarch,lgarch_fit_result,LGARCH_ML
  use lgarch_multivariate, only : mlgarch_simulate,fit_mlgarch,mlgarch_fit_result
  implicit none
  integer,parameter::n=300,m=2
  real(dp)::y(n),ym(n,m),constant(m),arch(m,m),garch(m,m),cov(m,m)
  type(lgarch_fit_result)::ufit
  type(mlgarch_fit_result)::mfit
  integer::i
  call seed_rng(12345)
  call lgarch_simulate(n,-0.12_dp,[0.10_dp],[0.78_dp],y)
  call fit_lgarch(y,1,1,LGARCH_ML,ufit,compute_vcov=.false.,max_iter=2500)
  print '(a,l1)', 'Univariate optimizer converged: ',ufit%converged
  print '(a,*(1x,f10.5))','Univariate log-GARCH parameters:',ufit%lgarch_par
  print '(a,f12.4)','Univariate model log likelihood: ',ufit%loglik_model

  constant=[-0.15_dp,-0.10_dp]; arch=0.0_dp; garch=0.0_dp
  arch(1,1)=0.10_dp; arch(2,2)=0.08_dp; arch(1,2)=0.01_dp
  garch(1,1)=0.72_dp; garch(2,2)=0.75_dp
  cov=reshape([1.0_dp,0.25_dp,0.25_dp,1.0_dp],[m,m])
  call mlgarch_simulate(n,constant,arch,garch,ym,innovations_vcov=cov)
  call fit_mlgarch(ym,1,1,mfit,compute_vcov=.false.,max_iter=1800,tol=5.0e-6_dp)
  print '(a,l1)', 'Multivariate optimizer converged: ',mfit%converged
  print '(a,f12.4)','VARMA log likelihood: ',mfit%objective_varma
  print '(a)','First five fitted standard deviations:'
  do i=1,5
    print '(*(1x,f10.5))',mfit%fitted_sd(i,:)
  end do
end program demo_lgarch
