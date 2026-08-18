program test_rank
   use suppdists
   implicit none
   integer :: fails
   real(dp), allocatable :: s(:)
   fails=0
   call check(pkendall(-0.4_dp,5),0.24166666666666667_dp,1e-14_dp,'pkendall')
   call check(dkendall(-0.4_dp,5),0.125_dp,1e-14_dp,'dkendall')
   call check(pkruskalwallis(5.0_dp,4,30,0.6_dp),0.8273930563325306_dp,2e-12_dp,'pkw')
   call check(dkruskalwallis(5.0_dp,4,30,0.6_dp),0.07981613495549193_dp,2e-12_dp,'dkw')
   if(abs(pfriedman(qfriedman(0.7_dp,5,4),5,4)-0.7_dp)>0.2_dp)fails=fails+1
   if(abs(pspearman(0.0_dp,10)-0.5_dp)>0.15_dp)fails=fails+1
   call norm_order(9,s)
   if(size(s)/=9 .or. abs(s(5))>1e-15_dp .or. maxval(abs(s+s(9:1:-1)))>1e-12_dp)fails=fails+1
   if(fails==0)then;print '(a)','test_rank: PASS';else;error stop 1;end if
contains
   subroutine check(got,want,tol,name)
      real(dp),intent(in)::got,want,tol;character(*),intent(in)::name
      if(abs(got-want)>tol)then
         print '(a,2es24.14)',trim(name)//' FAIL ',got,want;fails=fails+1
      end if
   end subroutine
end program test_rank
