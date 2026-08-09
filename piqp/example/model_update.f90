program model_update
   use piqp
   implicit none
   real(dp) :: p(2,2), c(2)
   type(piqp_model_type) :: model
   p=0.0_dp; p(1,1)=2.0_dp; p(2,2)=2.0_dp; c=[-2.0_dp,-4.0_dp]
   call model%setup(p,c)
   call model%solve()
   print '(a,2f12.7)','first x  = ',model%result%x
   c=[-4.0_dp,-2.0_dp]
   call model%update(c=c)
   call model%solve()
   print '(a,2f12.7)','second x = ',model%result%x
end program model_update
