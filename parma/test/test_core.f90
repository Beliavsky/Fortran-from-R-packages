! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
program test_core
   use parma
   implicit none
   real(dp), parameter :: tol = 1.0e-10_dp
   real(dp) :: w(2), data(4,2), r(4), expected, value, eps
   real(dp) :: m1(2), m2(2,2), grad(2), wp(2), wm(2)
   real(dp) :: x(4), lagged(4)
   character(len=10), allocatable :: dates(:)
   integer :: nout, info

   w = [0.6_dp,0.4_dp]
   data = reshape([0.01_dp,-0.01_dp,0.03_dp,-0.02_dp, &
                   0.02_dp,0.00_dp,-0.02_dp,0.01_dp],[4,2])
   r = matmul(data,w)
   expected = sum((r-sum(r)/4.0_dp)**2)/4.0_dp
   value = variance_risk(w,data)
   call assert_close(value,expected,tol,'variance risk')

   expected = sum(abs(r-sum(r)/4.0_dp))/4.0_dp
   call assert_close(mad_risk(w,data),expected,tol,'MAD risk')
   call assert_close(sentropy([0.5_dp,0.5_dp]),log(2.0_dp),tol,'Shannon entropy')
   call assert_close(empirical_quantile([1.0_dp,4.0_dp,2.0_dp,3.0_dp],0.5_dp), &
      2.5_dp,tol,'empirical quantile')
   call assert_close(lpm_risk(w,data,0.0_dp,1.0_dp), &
      sum(max(-r,0.0_dp))/4.0_dp,tol,'LPM')
   call assert_close(upm_risk(w,data,0.0_dp,1.0_dp), &
      sum(max(r,0.0_dp))/4.0_dp,tol,'UPM')

   x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   lagged = lag_vector(x,2,-1.0_dp)
   call assert_vector(lagged,[-1.0_dp,-1.0_dp,1.0_dp,2.0_dp],tol,'lag vector')

   call sequence_weekdays('2026-07-24','2026-07-28',dates,nout,info)
   if (info /= 0 .or. nout /= 3) error stop 'sequence_weekdays failed'
   if (dates(1) /= '2026-07-24' .or. dates(3) /= '2026-07-28') error stop 'weekday dates failed'

   m1 = [0.01_dp,0.02_dp]
   m2 = reshape([0.04_dp,0.01_dp,0.01_dp,0.09_dp],[2,2])
   grad = cara2_gradient(w,2.0_dp,m1,m2)
   eps = 1.0e-6_dp
   wp = w
   wm = w
   wp(1) = wp(1)+eps
   wm(1) = wm(1)-eps
   expected = (cara2_value(wp,2.0_dp,m1,m2)-cara2_value(wm,2.0_dp,m1,m2))/(2.0_dp*eps)
   call assert_close(grad(1),expected,1.0e-6_dp,'CARA2 gradient')

   print '(a)', 'test_core: PASS'

contains

   subroutine assert_close(actual,wanted,epsilonx,name)
      real(dp), intent(in) :: actual,wanted,epsilonx
      character(len=*), intent(in) :: name
      if (abs(actual-wanted) > epsilonx*max(1.0_dp,abs(wanted))) then
         write(*,'(a,2es24.14)') trim(name)//' failed: ',actual,wanted
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_vector(actual,wanted,epsilonx,name)
      real(dp), intent(in) :: actual(:),wanted(:),epsilonx
      character(len=*), intent(in) :: name
      if (size(actual) /= size(wanted) .or. any(abs(actual-wanted) > epsilonx)) then
         write(*,'(a)') trim(name)//' failed'
         error stop 1
      end if
   end subroutine assert_vector

end program test_core
