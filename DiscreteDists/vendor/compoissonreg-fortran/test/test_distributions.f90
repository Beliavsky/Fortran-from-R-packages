program test_distributions
   use compoissonreg
   implicit none
   integer :: fails, q
   real(dp) :: got, ref, p
   fails=0

   ref=exp(-2.0_dp)*2.0_dp**3/6.0_dp
   got=dcmp(3,2.0_dp,1.0_dp)
   call check_close('CMP Poisson density',got,ref,5.0e-7_dp,fails)
   call check_close('CMP Poisson mean',ecmp(10.0_dp,1.0_dp),10.0_dp,1.0e-10_dp,fails)
   call check_close('CMP Poisson variance',vcmp(10.0_dp,1.0_dp),10.0_dp,1.0e-10_dp,fails)
   call check_close('CMP Bernoulli limit',ecmp(0.5_dp,100.0_dp),1.0_dp/3.0_dp,1.0e-2_dp,fails)
   call check_close('CMP geometric limit',ecmp(0.25_dp,0.0001_dp),1.0_dp/3.0_dp,1.0e-2_dp,fails)

   p=pcmp(4,3.0_dp,0.8_dp)
   q=qcmp(p-1.0e-12_dp,3.0_dp,0.8_dp)
   if(q/=4)then;print *, 'FAIL quantile inversion',q;fails=fails+1;end if

   call check_close('ZICMP zero',dzicmp(0,2.0_dp,1.0_dp,0.25_dp), &
      0.25_dp+0.75_dp*exp(-2.0_dp),5.0e-7_dp,fails)
   call check_close('ZICMP mean',ezicmp(2.0_dp,1.0_dp,0.25_dp),1.5_dp,1.0e-10_dp,fails)
   call check_close('ZIP zero',dzip(0,5.0_dp,0.25_dp),0.25_dp+0.75_dp*exp(-5.0_dp),1.0e-12_dp,fails)
   call check_close('ZIP mean',ezip(5.0_dp,0.25_dp),3.75_dp,1.0e-12_dp,fails)

   if(fails==0)then;print *, 'test_distributions: PASS';else;error stop 1;end if
contains
   subroutine check_close(name,x,y,tol,fails)
      character(*),intent(in)::name
      real(dp),intent(in)::x,y,tol
      integer,intent(inout)::fails
      if(abs(x-y)>tol*max(1.0_dp,abs(y)))then
         print *, 'FAIL ',trim(name),x,y
         fails=fails+1
      end if
   end subroutine check_close
end program test_distributions
