program v03_features
   use rfast
   implicit none
   integer, parameter :: n=12
   real(dp) :: x(n,2), y(n)
   integer :: id(n),i
   type(random_intercept_result) :: ri
   type(selection_result) :: sel
   do i=1,n
      x(i,1)=real(i,dp)/real(n,dp)
      x(i,2)=sin(real(i,dp))
      id(i)=1+(i-1)/3
      y(i)=1.0_dp+2.0_dp*x(i,1)+0.25_dp*real(id(i)-2,dp)+0.05_dp*x(i,2)
   end do
   ri=rint_reg(y,x(:,1:1),id)
   sel=ompr(y,x,OMP_BIC)
   print '(a,2f10.5)','random-intercept beta: ',ri%beta
   print '(a,*(i0,1x))','OMP selected variables: ',sel%selected
end program v03_features
