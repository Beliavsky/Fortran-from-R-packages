program basic
   use misc_tools
   implicit none
   real(dp), allocatable :: m(:,:),v(:),se(:)
   real(dp) :: cov(2,2)
   type(coef_table_result) :: ct

   call sym_matrix([2.0_dp,-1.0_dp,0.0_dp,2.0_dp,-1.0_dp,2.0_dp],m,3)
   print '(a,l1)', "positive semidefinite: ", semidefiniteness(m)
   call vecli(m,v)
   print '(a,6f8.3)', "independent elements: ", v

   cov = reshape([0.25_dp,0.0_dp,0.0_dp,0.04_dp],[2,2])
   call std_er(cov,se)
   ct = coef_table([1.0_dp,-0.5_dp],se,20.0_dp)
   print '(a,2f10.4)', "standard errors: ", se
   print '(a,2es12.4)', "two-sided p-values: ", ct%table(:,4)
end program basic
