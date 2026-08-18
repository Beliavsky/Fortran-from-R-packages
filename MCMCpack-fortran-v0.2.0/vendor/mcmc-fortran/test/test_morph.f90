module test_morph_support
   use mcmc, only : dp
   implicit none
contains
   subroutine t_lud(state,value,data)
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      value=-0.5_dp*sum(state*state)
      if (present(data)) value=value
   end subroutine t_lud
end module test_morph_support

program test_morph
   use mcmc
   use test_morph_support
   implicit none
   type(morph_transform) :: mo,id
   type(metrop_result) :: r
   real(dp) :: x(3),y(3),z(3),lj
   integer :: fails
   fails=0
   x=[1.2_dp,-0.7_dp,0.4_dp]
   id=morph_identity(3)
   call id%transform(x,y); call id%inverse(y,z)
   if (maxval(abs(z-x))>1.0e-14_dp) fails=fails+1
   if (abs(id%log_jacobian(y))>1.0e-14_dp) fails=fails+1

   mo=morph_create(3,b=0.7_dp,r=0.5_dp,power=3.5_dp,center=[0.2_dp])
   call mo%transform(x,y); call mo%inverse(y,z)
   if (maxval(abs(z-x))>2.0e-10_dp) fails=fails+1
   lj=mo%log_jacobian(y)
   if (abs(lj)>=huge(1.0_dp)/2.0_dp) fails=fails+1

   call set_mcmc_seed(31415)
   r=morph_metrop(t_lud,mo,[0.0_dp,0.0_dp,0.0_dp],200,blen=10,scale=scale_constant(0.8_dp))
   if (r%status/=0) fails=fails+1
   if (r%accept<=0.0_dp .or. r%accept>=1.0_dp) fails=fails+1
   if (maxval(abs(sum(r%batch(51:,:),dim=1)/real(size(r%batch(51:,1)),dp)))>0.5_dp) fails=fails+1

   if (fails/=0) then
      print *,"test_morph: FAIL",fails,r%status,r%accept
      error stop 1
   end if
   print *,"test_morph: PASS"
end program test_morph
