! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
program xreg_example
  use lgarch_kinds, only : dp
  use lgarch_rng, only : seed_rng
  use lgarch_univariate, only : lgarch_simulate,fit_lgarch,lgarch_fit_result,LGARCH_LS
  implicit none
  integer,parameter::n=250
  real(dp)::y(n),x(n,1),forcing(n)
  type(lgarch_fit_result)::fit
  integer::i
  do i=1,n
    x(i,1)=sin(2.0_dp*acos(-1.0_dp)*real(i,dp)/40.0_dp)
  end do
  forcing=0.08_dp*x(:,1)
  call seed_rng(90210)
  call lgarch_simulate(n,-0.2_dp,[0.12_dp],[0.74_dp],y,xreg=forcing)
  call fit_lgarch(y,1,1,LGARCH_LS,fit,xreg=x,compute_vcov=.true.)
  print '(a,*(1x,f10.5))','Estimated log-GARCH-X parameters:',fit%lgarch_par
  print '(a,f12.5)','ARMA residual sum of squares: ',fit%rss
end program xreg_example
