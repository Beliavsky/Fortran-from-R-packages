! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program demo_fmultivar
  use fmultivar, only : dp, i8, pnorm2d, pt2d, dnorm2d, dt2d, &
    mvsnorm_rng, fit_multivariate_normal, fit_skew_normal, skew_fit_result, &
    density2d, grid_data, square_binning, binning_result
  implicit none
  real(dp), allocatable :: x(:,:)
  real(dp) :: mu(1),omega(1,1),alpha(1)
  type(skew_fit_result) :: normal_fit, skew_fit
  type(grid_data) :: kde
  type(binning_result) :: bins
  logical :: ok

  write(*,'(a,es14.6)') 'Bivariate Normal density: ',dnorm2d(0.3_dp,-0.2_dp,0.5_dp)
  write(*,'(a,es14.6)') 'Bivariate Normal CDF:     ',pnorm2d(0.3_dp,-0.2_dp,0.5_dp)
  write(*,'(a,es14.6)') 'Bivariate t density:      ',dt2d(0.3_dp,-0.2_dp,0.5_dp,6.0_dp)
  write(*,'(a,es14.6)') 'Bivariate t CDF:          ',pt2d(0.3_dp,-0.2_dp,0.5_dp,6.0_dp)

  mu=[-0.25_dp];omega=reshape([1.2_dp],[1,1]);alpha=[2.5_dp]
  call mvsnorm_rng(500,mu,omega,alpha,x,20260723_i8,ok)
  if(.not.ok)error stop 'Skew-Normal simulation failed'
  normal_fit=fit_multivariate_normal(x)
  skew_fit=fit_skew_normal(x,1200,3.0e-6_dp)
  write(*,'(a,f10.5)') 'Normal log likelihood:      ',normal_fit%loglik
  write(*,'(a,f10.5)') 'Skew-Normal log likelihood: ',skew_fit%loglik
  write(*,'(a,f10.5)') 'Estimated shape:            ',skew_fit%alpha(1)

  kde=density2d(x(:,1),0.4_dp*x(:,1)+0.8_dp,25)
  bins=square_binning(x(:,1),0.4_dp*x(:,1)+0.8_dp,16,16)
  write(*,'(a,i0,a,i0)') 'KDE grid: ',size(kde%x),' x ',size(kde%y)
  write(*,'(a,i0,a,i0)') 'Nonempty square bins: ',size(bins%count), &
    '; retained observations: ',sum(bins%count)
end program demo_fmultivar
