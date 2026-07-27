! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program exact_moment_shrinkage
  use performanceanalytics_mod, only: dp, exact_shrinkage_result, exact_m2_shrinkage, &
    exact_m3_shrinkage, exact_m4_shrinkage
  implicit none
  integer,parameter::n=80,p=3
  real(dp)::r(n,p),factor(n,1)
  type(exact_shrinkage_result)::fit
  integer::t

  do t=1,n
    factor(t,1)=sin(0.17_dp*real(t,dp))+0.35_dp*cos(0.07_dp*real(t,dp))
    r(t,1)=0.75_dp*factor(t,1)+0.25_dp*sin(0.53_dp*real(t,dp))
    r(t,2)=-0.45_dp*factor(t,1)+0.30_dp*cos(0.41_dp*real(t,dp))
    r(t,3)=0.30_dp*factor(t,1)+0.20_dp*sin(0.29_dp*real(t,dp))**3
  end do

  call exact_m2_shrinkage(r,[1,2,3,4],fit,factor)
  write(*,'(a,*(1x,f8.5))')'M2 shrinkage weights:',fit%lambda

  call exact_m3_shrinkage(r,[1,2,3,4,5,6],fit,factor)
  write(*,'(a,*(1x,f8.5))')'M3 shrinkage weights:',fit%lambda

  call exact_m4_shrinkage(r,[1,2,3,4],fit,factor)
  write(*,'(a,*(1x,f8.5))')'M4 shrinkage weights:',fit%lambda
end program exact_moment_shrinkage
