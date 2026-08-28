module banded_problem
   use nleqslv_fortran, only : dp
   implicit none
contains
   subroutine brdban(x,y)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      integer::k,j,k1,k2,n
      real(dp)::temp
      n=size(x)
      do k=1,n
         k1=max(1,k-2); k2=min(n,k+2); temp=0.0_dp
         do j=k1,k2
            if(j/=k) temp=temp+x(j)*(1.0_dp+x(j))
         end do
         y(k)=x(k)*(2.0_dp+5.0_dp*x(k)**2)+1.0_dp-temp
      end do
   end subroutine
end module
program test_banded
   use nleqslv_fortran
   use banded_problem
   implicit none
   type(nleq_options)::o
   type(nleq_result)::a,b
   real(dp)::x0(10)
   integer::fails
   fails=0; x0=-1.0_dp
   o=nleq_options(); o%method=NLEQ_NEWTON
   call solve_nleqslv(x0,brdban,a,o)
   o%dsub=2; o%dsuper=2
   call solve_nleqslv(x0,brdban,b,o)
   if(a%termcd/=1 .or. b%termcd/=1) fails=fails+1
   if(maxval(abs(a%x-b%x))>1.0e-9_dp) fails=fails+1
   o=nleq_options(); o%method=NLEQ_BROYDEN; o%dsub=2; o%dsuper=2
   call solve_nleqslv(x0,brdban,b,o)
   if(b%termcd/=1) fails=fails+1
   if(maxval(abs(a%x-b%x))>1.0e-7_dp) fails=fails+1
   if(fails==0) then
      print *, 'test_banded: PASS'
   else
      print *, 'test_banded: FAIL',fails; error stop 1
   end if
end program
