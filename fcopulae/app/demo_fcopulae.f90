! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program demo_fcopulae
  use fcopulae, only : dp, i8, copula_fit_result, archm_rng, archm_fit, archm_tau, &
    elliptical_rng, elliptical_fit, elliptical_tau, ev_rng, ev_fit, ev_tau
  implicit none
  real(dp),allocatable :: x(:,:),param(:)
  type(copula_fit_result) :: fit

  call archm_rng(180,2.0_dp,1,x,20260724_i8)
  fit=archm_fit(x(:,1),x(:,2),1,max_iter=350)
  write(*,'(a,f10.5)')'Clayton true alpha:      ',2.0_dp
  write(*,'(a,f10.5)')'Clayton fitted alpha:    ',fit%param(1)
  write(*,'(a,f10.5)')'Clayton Kendall tau:     ',archm_tau(fit%param(1),1)

  call elliptical_rng(180,0.55_dp,'norm',x,20260725_i8)
  fit=elliptical_fit(x(:,1),x(:,2),'norm',max_iter=350)
  write(*,'(a,f10.5)')'Gaussian true rho:       ',0.55_dp
  write(*,'(a,f10.5)')'Gaussian fitted rho:     ',fit%param(1)
  write(*,'(a,f10.5)')'Gaussian Kendall tau:    ',elliptical_tau(fit%param(1))

  allocate(param(1));param=2.0_dp
  call ev_rng(180,param,'gumbel',x,20260726_i8)
  fit=ev_fit(x(:,1),x(:,2),'gumbel',max_iter=350)
  write(*,'(a,f10.5)')'EV Gumbel true alpha:    ',2.0_dp
  write(*,'(a,f10.5)')'EV Gumbel fitted alpha:  ',fit%param(1)
  write(*,'(a,f10.5)')'EV Gumbel Kendall tau:   ',ev_tau(fit%param,'gumbel')
end program demo_fcopulae
