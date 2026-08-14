! SPDX-License-Identifier: GPL-2.0-only
program test_hypervolume
   use mco, only : dp, dominated_hypervolume
   implicit none
   real(dp) :: p1(3,1), p2(3,2), ref(3)
   ref=[4.0_dp,5.0_dp,6.0_dp]
   p1(:,1)=[1.0_dp,2.0_dp,3.0_dp]
   call close(dominated_hypervolume(p1,ref),27.0_dp,1.0e-12_dp,"single box")
   p2(:,1)=[1.0_dp,2.0_dp,3.0_dp]
   p2(:,2)=[2.0_dp,1.0_dp,2.0_dp]
   ! 27 + 32 - intersection [2,2,3] with volume 18 = 41.
   call close(dominated_hypervolume(p2,ref),41.0_dp,1.0e-12_dp,"two boxes")
   print '(a)', 'test_hypervolume: PASS'
contains
   subroutine close(x,y,tol,msg)
      real(dp),intent(in)::x,y,tol; character(*),intent(in)::msg
      if(abs(x-y)>tol) error stop msg
   end subroutine
end program
