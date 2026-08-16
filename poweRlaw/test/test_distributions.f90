program test_distributions
   use powerlaw
   implicit none
   type(powerlaw_dist) :: m
   real(dp) :: ref
   integer :: fails
   fails=0
   m=displ([1.0_dp,1.0_dp,3.0_dp])
   call m%set_pars([2.0_dp])
   if(abs(m%cdf(1.0_dp)-0.6079271018544_dp)>2.0e-10_dp) fails=fails+1
   if(abs(m%cdf(3.0_dp)-0.8274563330796_dp)>2.0e-10_dp) fails=fails+1
   call m%set_xmin(2.0_dp)
   if(abs(m%cdf(3.0_dp)-0.5599194238202_dp)>2.0e-10_dp) fails=fails+1

   m=conpl([2.0_dp,2.0_dp,4.0_dp])
   call m%set_xmin(1.0_dp); call m%set_pars([2.0_dp])
   if(abs(m%cdf(2.0_dp)-0.5_dp)>1.0e-12_dp) fails=fails+1
   if(abs(m%pdf(2.0_dp)-0.25_dp)>1.0e-12_dp) fails=fails+1

   m=disexp([1.0_dp,1.0_dp,3.0_dp]); call m%set_pars([1.0_dp])
   if(abs(m%cdf(1.0_dp)-(1.0_dp-exp(-1.0_dp)))>1.0e-12_dp) fails=fails+1

   m=conexp([2.0_dp,2.0_dp,4.0_dp])
   call m%set_xmin(1.0_dp); call m%set_pars([1.0_dp])
   if(abs(m%cdf(2.0_dp)-(1.0_dp-exp(-1.0_dp)))>1.0e-12_dp) fails=fails+1

   m=dispois([1.0_dp,1.0_dp,3.0_dp]); call m%set_pars([1.0_dp])
   ref=exp(-1.0_dp)/(1.0_dp-exp(-1.0_dp))
   if(abs(m%cdf(1.0_dp)-ref)>2.0e-12_dp) fails=fails+1

   m=conweibull([1.0_dp,2.0_dp,3.0_dp])
   call m%set_xmin(1.0_dp); call m%set_pars([2.0_dp,2.0_dp])
   ref=(exp(-0.25_dp)-exp(-1.0_dp))/exp(-0.25_dp)
   if(abs(m%cdf(2.0_dp)-ref)>2.0e-12_dp) fails=fails+1

   if(fails/=0) then
      print *,"test_distributions: FAIL",fails
      error stop 1
   end if
   print *,"test_distributions: PASS"
end program
