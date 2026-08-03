program test_ordering
   use portvine, only : dp, greedy_dvine_order
   implicit none
   integer, parameter :: n=120,d=4
   real(dp) :: u(d,n), t
   integer :: i,order(d),status
   do i=1,n
      t=real(i,dp)/real(n+1,dp)
      u(1,i)=t
      u(2,i)=max(1.0e-5_dp,min(1.0_dp-1.0e-5_dp,t+0.03_dp*sin(13.0_dp*t)))
      u(3,i)=real(mod(37*i,121)+1,dp)/122.0_dp
      u(4,i)=max(1.0e-5_dp,min(1.0_dp-1.0e-5_dp,u(3,i)+0.04_dp*cos(9.0_dp*t)))
   end do
   call greedy_dvine_order(u,order,[3],2,status)
   if(status/=0 .or. order(1)/=3 .or. .not.is_permutation(order))error stop 1
   call greedy_dvine_order(u,order,[3,4],2,status)
   if(status/=0 .or. any(order(1:2)/=[3,4]) .or. .not.is_permutation(order))error stop 2
   print '(a)', 'test_ordering: PASS'
contains
   logical function is_permutation(x)
      integer,intent(in)::x(:)
      integer::j
      is_permutation=.true.
      do j=1,size(x)
         if(count(x==j)/=1)is_permutation=.false.
      end do
   end function is_permutation
end program test_ordering
