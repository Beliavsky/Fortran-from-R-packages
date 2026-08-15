program test_core
   use benfordtests
   implicit none
   real(dp), allocatable :: p(:),q(:),x(:)
   integer, allocatable :: d(:),seq(:)
   integer :: fails
   fails=0
   p=pbenf(1);q=qbenf(1);seq=signifd_seq(2)
   call check(abs(sum(p)-1.0_dp)<1.0d-14,'pbenf sums to one')
   call check(abs(p(1)-log10(2.0_dp))<1.0d-14,'first Benford probability')
   call check(abs(q(9)-1.0_dp)<1.0d-14,'qbenf endpoint')
   call check(size(seq)==90 .and. seq(1)==10 .and. seq(90)==99,'signifd sequence')
   allocate(x(5));x=[123.4_dp,-0.00456_dp,99.0_dp,1.0_dp,0.0789_dp]
   d=significant_digits(x,2)
   call check(all(d==[12,45,99,10,78]),'significant digits')
   deallocate(x);allocate(x(20000));call rbenf(x)
   call check(all(x>=1.0_dp .and. x<10.0_dp),'rbenf support')
   call check(abs(real(count(significant_digits(x,1)==1),dp)/size(x)-p(1))<0.02_dp,'rbenf digit frequency')
   if(fails==0)then
      print '(a)','test_core: PASS'
   else
      print '(a,i0)','test_core: FAIL ',fails
      error stop 1
   end if
contains
   subroutine check(ok,name)
      logical,intent(in)::ok;character(len=*),intent(in)::name
      if(.not.ok)then;fails=fails+1;print '(a,a)','FAIL: ',trim(name);end if
   end subroutine
end program
