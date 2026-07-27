! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program test_rolling_cleaning
  use kinds_mod, only: dp
  use rolling_mod
  use returns_mod, only: clean_boudt, loc_scale_robust, geltner_unsmooth
  use test_support_mod
  implicit none
  real(dp) :: r(8), out(8), out2(8), x(8), y(8), location, scale
  integer :: n
  r=[0.01_dp,0.02_dp,-0.01_dp,0.03_dp,0.00_dp,0.04_dp,-0.02_dp,0.01_dp]
  call rolling_statistic(r,3,'mean',out,n)
  call assert_true(n==6,'rolling count')
  call assert_close(out(1),sum(r(1:3))/3.0_dp,1.0e-13_dp,'rolling mean')
  call expanding_statistic(r,3,'sd',out,n)
  call assert_true(n==6,'expanding count')
  x=[(real(n,dp),n=1,8)];y=2.0_dp*x+1.0_dp
  call rolling_correlation(x,y,4,out,n)
  call assert_true(all(abs(out(:n)-1.0_dp)<1.0e-12_dp),'rolling correlation')

  x=[0.0_dp,0.01_dp,-0.01_dp,0.02_dp,-0.02_dp,0.01_dp,0.0_dp,1.0_dp]
  call loc_scale_robust(x,location,scale)
  call clean_boudt(x,0.01_dp,out)
  call assert_true(out(8)<x(8),'boudt cleaning')
  call geltner_unsmooth(r,out2)
  call assert_true(all(abs(out2)<1.0_dp),'geltner finite')
  write(*,'(a)')'Rolling and cleaning tests passed.'
end program test_rolling_cleaning
