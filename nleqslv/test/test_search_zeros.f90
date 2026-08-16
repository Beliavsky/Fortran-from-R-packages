module search_problem
   use nleqslv_fortran, only : dp
   implicit none
contains
   subroutine dsln(x,y)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      y(1)=x(1)**2+x(2)**2-2.0_dp
      y(2)=exp(x(1)-1.0_dp)+x(2)**3-2.0_dp
   end subroutine
end module
program test_search_zeros
   use nleqslv_fortran
   use search_problem
   implicit none
   type(search_zeros_result)::r
   type(nleq_options)::o
   real(dp)::starts(25,2), f(2)
   integer::i,j,k,fails
   fails=0; k=0
   do i=0,4
      do j=0,4
         k=k+1
         starts(k,1)=-2.0_dp+real(i,dp)
         starts(k,2)=-2.0_dp+real(j,dp)
      end do
   end do
   o=nleq_options(); o%method=NLEQ_BROYDEN; o%global=NLEQ_DBLDOG
   call search_zeros(starts,dsln,r,6,o)
   if(size(r%x,1)<2) fails=fails+1
   do i=1,size(r%x,1)
      call dsln(r%x(i,:),f)
      if(maxval(abs(f))>1.0e-6_dp) fails=fails+1
   end do
   if(fails==0) then
      print *, 'test_search_zeros: PASS, roots=',size(r%x,1)
   else
      print *, 'test_search_zeros: FAIL',fails,' roots=',size(r%x,1); error stop 1
   end if
end program
