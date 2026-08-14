! SPDX-License-Identifier: GPL-2.0-only
program test_functions
   use mco, only : dp, belegundu, belegundu_constraints, binh1, zdt1, zdt2, zdt3, vnt
   implicit none
   real(dp) :: x2(2), f3(3), f2(2), g2(2), x30(30)
   x2=[1.0_dp,2.0_dp]
   call belegundu(x2,f2); call close_vec(f2,[0.0_dp,4.0_dp],1.0e-12_dp,"belegundu")
   call belegundu_constraints(x2,g2); call close_vec(g2,[0.0_dp,4.0_dp],1.0e-12_dp,"constraints")
   call binh1(x2,f2); call close_vec(f2,[5.0_dp,25.0_dp],1.0e-12_dp,"binh1")
   x30=0.0_dp; x30(1)=0.25_dp
   call zdt1(x30,f2); call close_vec(f2,[0.25_dp,0.5_dp],1.0e-12_dp,"zdt1")
   call zdt2(x30,f2); call close_vec(f2,[0.25_dp,0.9375_dp],1.0e-12_dp,"zdt2")
   call zdt3(x30,f2); call close_vec(f2,[0.25_dp,0.5_dp-0.25_dp*sin(2.5_dp*acos(-1.0_dp))],1.0e-12_dp,"zdt3")
   call vnt([0.0_dp,0.0_dp],f3)
   call close_vec(f3,[0.0_dp,17.037037037037037_dp,-0.1_dp],1.0e-12_dp,"vnt")
   print '(a)', 'test_functions: PASS'
contains
   subroutine close_vec(x,y,tol,msg)
      real(dp),intent(in)::x(:),y(:),tol; character(*),intent(in)::msg
      if(maxval(abs(x-y))>tol) error stop msg
   end subroutine
end program
