module test_temper_support
   use mcmc, only : dp
   implicit none
contains
   subroutine temper_lud(state,value,data)
      real(dp), intent(in) :: state(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      integer :: i
      real(dp) :: sigma
      i=nint(state(1))
      if (i==1) then
         sigma=1.0_dp
      else if (i==2) then
         sigma=2.0_dp
      else
         value=-huge(1.0_dp); return
      end if
      value=-0.5_dp*sum((state(2:)/sigma)**2)-real(size(state)-1,dp)*log(sigma)
      if (present(data)) value=value
   end subroutine temper_lud
end module test_temper_support

program test_temper
   use mcmc
   use test_temper_support
   implicit none
   logical :: neigh(2,2)
   type(mcmc_scale) :: scales(2)
   type(temper_serial_result) :: sr
   type(temper_parallel_result) :: pr
   real(dp) :: initp(2,1),occ(2),mean1,mean2
   integer :: fails
   fails=0
   neigh=.false.; neigh(1,2)=.true.; neigh(2,1)=.true.
   scales(1)=scale_constant(1.0_dp); scales(2)=scale_constant(2.0_dp)

   call set_mcmc_seed(111)
   sr=temper_serial(temper_lud,[1.0_dp,0.0_dp],neigh,500,blen=20,scales=scales)
   if (sr%status/=0) fails=fails+1
   occ=sum(sr%ibatch(101:,:),dim=1)/real(size(sr%ibatch(101:,1)),dp)
   if (maxval(abs(occ-0.5_dp))>0.12_dp) fails=fails+1
   if (any(sr%acceptx<0.0_dp) .or. any(sr%acceptx>1.0_dp)) fails=fails+1

   initp=0.0_dp
   call set_mcmc_seed(222)
   pr=temper_parallel(temper_lud,initp,neigh,400,blen=15,scales=scales)
   if (pr%status/=0) fails=fails+1
   mean1=sum(pr%batch(101:,1,1))/real(size(pr%batch(101:,1,1)),dp)
   mean2=sum(pr%batch(101:,2,1))/real(size(pr%batch(101:,2,1)),dp)
   if (abs(mean1)>0.25_dp .or. abs(mean2)>0.5_dp) fails=fails+1
   if (any(pr%acceptx<0.0_dp) .or. any(pr%acceptx>1.0_dp)) fails=fails+1

   if (fails/=0) then
      print *,"test_temper: FAIL",fails,sr%status,pr%status,occ,mean1,mean2
      error stop 1
   end if
   print *,"test_temper: PASS"
end program test_temper
