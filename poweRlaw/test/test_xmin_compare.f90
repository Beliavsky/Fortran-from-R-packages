program test_xmin_compare
   use powerlaw
   implicit none
   type(powerlaw_dist) :: m1,m2
   type(estimate_xmin_result) :: e
   type(compare_result) :: c
   real(dp) :: x(10), cand(101,1), xmins(8)
   integer :: i,fails
   fails=0
   do i=1,10
      x(i)=real(i,dp)
   end do
   do i=1,101
      cand(i,1)=2.0_dp+0.01_dp*real(i-1,dp)
   end do
   do i=1,8
      xmins(i)=real(i,dp)
   end do
   m1=displ(x)
   e=estimate_xmin(m1,xmins=xmins,candidates=cand,distance="ks")
   if(abs(e%xmin-5.0_dp)>1.0e-12_dp) fails=fails+1
   if(abs(e%pars(1)-3.0_dp)>0.11_dp) fails=fails+1

   m1=conpl([1.0_dp,1.2_dp,1.5_dp,2.0_dp,3.0_dp,5.0_dp])
   call m1%set_xmin(1.0_dp); call m1%set_pars([2.0_dp])
   m2=conexp(m1%data)
   call m2%set_xmin(1.0_dp); call m2%set_pars([0.5_dp])
   c=compare_distributions(m1,m2)
   if(size(c%ratio)/=6) fails=fails+1
   if(c%p_two_sided<0.0_dp .or. c%p_two_sided>1.0_dp) fails=fails+1

   if(fails/=0) then
      print *,"test_xmin_compare: FAIL",fails,e%xmin,e%pars
      error stop 1
   end if
   print *,"test_xmin_compare: PASS"
end program
