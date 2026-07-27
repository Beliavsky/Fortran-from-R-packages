! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program dependence_example
  use fcopulae, only : dp, archm_tau, archm_rho, archm_tail_coeff, elliptical_tau, &
    elliptical_rho, elliptical_tail_coeff, ev_tau, ev_rho, ev_tail_coeff
  implicit none
  real(dp)::lower,upper,param(1)
  param=2.0_dp
  call archm_tail_coeff(2.0_dp,1,lower,upper)
  write(*,'(a,3f12.6)')'Clayton tau/rho/lower-tail: ',archm_tau(2.0_dp,1),archm_rho(2.0_dp,1),lower
  call elliptical_tail_coeff(0.6_dp,'t',lower,upper,5.0_dp)
  write(*,'(a,3f12.6)')'Student-t tau/rho/tail:     ',elliptical_tau(0.6_dp), &
    elliptical_rho(0.6_dp,'t',5.0_dp),lower
  call ev_tail_coeff(param,'gumbel',lower,upper)
  write(*,'(a,3f12.6)')'EV Gumbel tau/rho/upper:    ',ev_tau(param,'gumbel'),ev_rho(param,'gumbel'),upper
end program dependence_example
