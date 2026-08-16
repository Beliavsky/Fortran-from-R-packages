program basic_fit
   use powerlaw
   implicit none
   type(powerlaw_dist) :: m
   type(estimate_xmin_result) :: est
   real(dp) :: x(20)
   integer :: i
   do i=1,20
      x(i)=real(i,dp)
   end do
   m=displ(x)
   est=estimate_xmin(m)
   print '(a,f8.3)', "xmin  = ",est%xmin
   print '(a,f10.5)', "alpha = ",est%pars(1)
   print '(a,f10.6)', "KS    = ",est%gof
end program
