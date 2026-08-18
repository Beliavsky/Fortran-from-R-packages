module test_metrop_support
   use mcmc, only : dp
   implicit none
contains
   subroutine normal_lud(state,value,data)
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      value=-0.5_dp*sum(state*state)
      if (present(data)) value=value
   end subroutine normal_lud
   subroutine first_output(state,value,data)
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: value(:)
      class(*), intent(in), optional :: data
      value(1)=state(1)
      if (present(data)) value(1)=value(1)
   end subroutine first_output
end module test_metrop_support

program test_metrop
   use mcmc
   use test_metrop_support
   implicit none
   type(metrop_result) :: r
   type(mcmc_scale) :: sc
   real(dp) :: mean1
   integer :: fails
   fails=0
   call set_mcmc_seed(12345)
   sc=scale_constant(1.2_dp)
   r=metrop(normal_lud,[0.0_dp,0.0_dp],400,blen=20,nspac=1,scale=sc)
   if (r%status/=0) fails=fails+1
   if (r%accept<0.2_dp .or. r%accept>0.8_dp) fails=fails+1
   mean1=sum(r%batch(101:,1))/real(size(r%batch(101:,1)),dp)
   if (abs(mean1)>0.15_dp) fails=fails+1
   if (abs(sum(r%accept_batch)/real(size(r%accept_batch),dp)-r%accept)>1.0e-13_dp) fails=fails+1

   call set_mcmc_seed(6789)
   r=metrop(normal_lud,[0.0_dp,0.0_dp],100,blen=1,scale=scale_diagonal([0.5_dp,1.0_dp]),debug=.true.)
   if (r%status/=0 .or. .not. allocated(r%current)) fails=fails+1
   if (size(r%current,1)/=100) fails=fails+1

   call set_mcmc_seed(42)
   r=metrop(normal_lud,[0.0_dp,0.0_dp],100,blen=10,scale=scale_constant(0.8_dp), &
            out_dim=1,outfun=first_output)
   if (r%status/=0 .or. size(r%batch,2)/=1) fails=fails+1

   call set_mcmc_seed(43)
   r=metrop(normal_lud,[0.0_dp,0.0_dp],50,blen=5, &
            scale=scale_full(reshape([0.7_dp,0.2_dp,0.0_dp,0.6_dp],[2,2])))
   if (r%status/=0 .or. r%accept<=0.0_dp .or. r%accept>=1.0_dp) fails=fails+1

   if (fails/=0) then
      print *,"test_metrop: FAIL",fails,r%accept
      error stop 1
   end if
   print *,"test_metrop: PASS"
end program test_metrop
