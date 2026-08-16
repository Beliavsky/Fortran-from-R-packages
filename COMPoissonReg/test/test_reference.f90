program test_reference
   use compoissonreg
   implicit none
   integer :: fails
   fails=0
   call check_close('Z 3,0.8',ncmp(3.0_dp,0.8_dp),36.054407304327465_dp,3.0e-6_dp,fails)
   call check_close('d(3) 3,0.8',dcmp(3,3.0_dp,0.8_dp),0.17860121262347748_dp,3.0e-6_dp,fails)
   call check_close('P(4) 3,0.8',pcmp(4,3.0_dp,0.8_dp),0.6096645544097626_dp,3.0e-6_dp,fails)
   call check_close('mean 3,0.8',ecmp(3.0_dp,0.8_dp),4.083321469318807_dp,3.0e-6_dp,fails)
   call check_close('var 3,0.8',vcmp(3.0_dp,0.8_dp),4.919377995256668_dp,2.0e-5_dp,fails)
   call check_close('Z .5,2',ncmp(0.5_dp,2.0_dp),1.5660829297563505_dp,3.0e-6_dp,fails)
   call check_close('d0 .5,2',dcmp(0,0.5_dp,2.0_dp),0.6385357895163182_dp,3.0e-6_dp,fails)
   call check_close('geometric mean',ecmp(0.25_dp,0.0_dp),1.0_dp/3.0_dp,1.0e-12_dp,fails)
   call check_close('geometric var',vcmp(0.25_dp,0.0_dp),4.0_dp/9.0_dp,1.0e-12_dp,fails)
   if(fails==0)then;print *,'test_reference: PASS';else;error stop 1;end if
contains
   subroutine check_close(name,x,y,tol,fails)
      character(*),intent(in)::name
      real(dp),intent(in)::x,y,tol
      integer,intent(inout)::fails
      if(abs(x-y)>tol*max(1.0_dp,abs(y)))then
         print *,'FAIL ',trim(name),x,y
         fails=fails+1
      end if
   end subroutine check_close
end program test_reference
