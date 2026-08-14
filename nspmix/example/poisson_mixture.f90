program poisson_mixture
   use nspmix
   implicit none
   type(nsp_data) :: x
   type(nspmix_result) :: fit
   real(dp) :: v(9),w(9)
   integer :: i
   v=[(real(i-1,dp),i=1,9)]
   w=[34.0_dp,31.0_dp,25.0_dp,24.0_dp,21.0_dp,18.0_dp,14.0_dp,8.0_dp,5.0_dp]
   call make_nppois_data(v,x,w)
   call cnm(x,fit,maxit=30,ngrid=80,kmax=10)
   print *, "log-likelihood:", fit%ll
   print *, "support:", fit%mix%pt
   print *, "probabilities:", fit%mix%pr
end program
