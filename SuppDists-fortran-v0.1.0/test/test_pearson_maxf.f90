program test_pearson_maxf
   use suppdists
   implicit none
   integer :: fails
   real(dp) :: p,x
   fails=0
   call check(dpearson(0.3_dp,10,0.0_dp),0.8242182812499996_dp,2e-10_dp,'dpearson')
   call check(ppearson(0.3_dp,10,0.0_dp),0.800154265625_dp,3e-8_dp,'ppearson')
   x=qpearson(0.8_dp,10,0.0_dp);call check(x,0.29981286864575285_dp,2e-6_dp,'qpearson')
   call check(pmaxfratio(3.0_dp,5,4),0.348035997791578_dp,2e-6_dp,'pmaxfratio')
   call check(dmaxfratio(3.0_dp,5,4),0.2146393307605734_dp,3e-6_dp,'dmaxfratio')
   p=rpearson(10,0.2_dp);if(abs(p)>1.0_dp)fails=fails+1
   p=rmaxfratio(5,4);if(p<1.0_dp)fails=fails+1
   if(fails==0)then;print '(a)','test_pearson_maxf: PASS';else;error stop 1;end if
contains
   subroutine check(got,want,tol,name)
      real(dp),intent(in)::got,want,tol;character(*),intent(in)::name
      if(abs(got-want)>tol)then
         print '(a,2es24.14)',trim(name)//' FAIL ',got,want;fails=fails+1
      end if
   end subroutine
end program test_pearson_maxf
