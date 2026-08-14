! SPDX-License-Identifier: GPL-2.0-only
program test_pareto_quality
   use mco, only : dp, pareto_mask, pareto_filter, normalize_front, &
      generational_distance, generalized_spread, dominated_hypervolume, epsilon_indicator
   implicit none
   real(dp) :: a(2,4), ref(2), hv, b(2,3), n(2,3)
   logical :: mask(4)
   real(dp), allocatable :: f(:,:)
   a=reshape([1.0_dp,3.0_dp, 2.0_dp,2.0_dp, 3.0_dp,1.0_dp, 3.0_dp,3.0_dp],[2,4])
   mask=pareto_mask(a)
   call check(all(mask .eqv. [.true.,.true.,.true.,.false.]),"pareto mask")
   f=pareto_filter(a)
   call check(size(f,2)==3,"pareto filter size")
   ref=[4.0_dp,4.0_dp]
   hv=dominated_hypervolume(f,ref)
   call close(hv,6.0_dp,1.0e-12_dp,"2D hypervolume")
   b=f
   n=normalize_front(b,[1.0_dp,1.0_dp],[3.0_dp,3.0_dp])
   call close(n(1,1),0.0_dp,1.0e-12_dp,"normalize low")
   call close(n(2,1),1.0_dp,1.0e-12_dp,"normalize high")
   call close(generational_distance(f,f),0.0_dp,1.0e-12_dp,"GD identity")
   call close(generalized_spread(f,f),0.0_dp,1.0e-12_dp,"spread identity")
   call close(epsilon_indicator(f,f),0.0_dp,1.0e-12_dp,"epsilon identity")
   print '(a)', 'test_pareto_quality: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok; character(*),intent(in)::msg
      if(.not.ok) error stop msg
   end subroutine
   subroutine close(x,y,tol,msg)
      real(dp),intent(in)::x,y,tol; character(*),intent(in)::msg
      call check(abs(x-y)<=tol,msg)
   end subroutine
end program
