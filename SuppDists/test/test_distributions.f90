program test_distributions
   use suppdists
   implicit none
   integer :: fails
   real(dp) :: x,p
   fails=0
   call check(dinvgauss(1.5_dp,2.0_dp,3.0_dp),0.3533380431253714_dp,2e-12_dp,'dinvgauss')
   call check(pinvgauss(1.5_dp,2.0_dp,3.0_dp),0.49569012484162966_dp,2e-12_dp,'pinvgauss')
   call check(qinvgauss(0.7_dp,2.0_dp,3.0_dp),2.253125822882857_dp,2e-10_dp,'qinvgauss')
   call check(dghyper(3,5.0_dp,7.0_dp,20.0_dp),0.17608359133126936_dp,2e-13_dp,'dghyper classic')
   call check(pghyper(3,5.0_dp,7.0_dp,20.0_dp),0.9692982456140351_dp,2e-13_dp,'pghyper classic')
   x=rinvgauss(2.0_dp,3.0_dp)
   if(x<=0.0_dp)fails=fails+1
   p=real(rghyper(5.0_dp,7.0_dp,20.0_dp),dp)
   if(p<0.0_dp .or. p>5.0_dp)fails=fails+1
   if(fails==0)then;print '(a)','test_distributions: PASS';else;error stop 1;end if
contains
   subroutine check(got,want,tol,name)
      real(dp),intent(in)::got,want,tol;character(*),intent(in)::name
      if(abs(got-want)>tol)then
         print '(a,2es24.14)',trim(name)//' FAIL ',got,want;fails=fails+1
      end if
   end subroutine
end program test_distributions
