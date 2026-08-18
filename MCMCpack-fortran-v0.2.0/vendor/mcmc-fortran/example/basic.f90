module basic_support
   use mcmc, only : dp
   implicit none
contains
   subroutine log_target(state,value,data)
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      value=-0.5_dp*sum(state*state)
      if (present(data)) value=value
   end subroutine log_target
end module basic_support

program basic
   use mcmc
   use basic_support
   implicit none
   type(metrop_result) :: r
   real(dp) :: mean(2)

   call set_mcmc_seed(1234)
   r=metrop(log_target,[0.0_dp,0.0_dp],500,blen=20,scale=scale_constant(1.0_dp))
   mean=sum(r%batch(101:,:),dim=1)/real(size(r%batch,1)-100,dp)
   print '(a,f8.4)', 'acceptance: ',r%accept
   print '(a,2f10.5)', 'estimated mean: ',mean
end program basic
