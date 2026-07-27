! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on MarkowitzR, copyright 2014-2020 Steven E. Pav.
program conditional_constraints
   use markowitzr, only: dp, markowitz_result, mp_vcov
   implicit none
   integer, parameter :: n = 200
   real(dp) :: x(n,3), feat(n,2), jmat(2,3), gmat(1,3)
   type(markowitz_result) :: unconstrained, constrained
   integer :: i

   do i = 1, n
      feat(i,1) = sin(0.07_dp*real(i,dp))
      feat(i,2) = cos(0.11_dp*real(i,dp))
      x(i,1) = 0.4_dp*feat(i,1)-0.2_dp*feat(i,2)+0.3_dp*sin(0.31_dp*real(i,dp))
      x(i,2) = -0.1_dp*feat(i,1)+0.5_dp*feat(i,2)+0.2_dp*cos(0.23_dp*real(i,dp))
      x(i,3) = 0.2_dp*feat(i,1)+0.1_dp*feat(i,2)+0.25_dp*sin(0.19_dp*real(i,dp))
   end do

   unconstrained = mp_vcov(x,feat=feat)
   if (unconstrained%status /= 0) error stop 1

   jmat = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.0_dp],[2,3])
   gmat = reshape([1.0_dp,0.0_dp,1.0_dp],[1,3])
   constrained = mp_vcov(x,feat=feat,jmat=jmat,gmat=gmat)
   if (constrained%status /= 0) error stop 1

   print '(a)', 'Unconstrained coefficient matrix:'
   do i = 1, 3
      print '(*(1x,es13.5))', unconstrained%w(i,:)
   end do
   print '(a)', 'Constrained and hedged coefficient matrix:'
   do i = 1, 3
      print '(*(1x,es13.5))', constrained%w(i,:)
   end do
end program conditional_constraints
