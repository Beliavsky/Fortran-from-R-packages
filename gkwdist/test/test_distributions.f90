program test_distributions
   use gkwdist
   implicit none
   integer :: fails
   real(dp) :: x,p,q,dens
   fails=0

   x=0.4_dp
   call check_close(dbeta_(x,2.0_dp,2.0_dp),1.728_dp,2.0e-12_dp,'beta density',fails)
   call check_close(pbeta_(x,2.0_dp,2.0_dp),0.5248_dp,2.0e-12_dp,'beta cdf',fails)
   call check_close(qbeta_(0.5248_dp,2.0_dp,2.0_dp),x,2.0e-11_dp,'beta quantile',fails)

   dens=6.0_dp*x*(1.0_dp-x*x)**2
   p=1.0_dp-(1.0_dp-x*x)**3
   call check_close(dkw(x,2.0_dp,3.0_dp),dens,2.0e-12_dp,'kw density',fails)
   call check_close(pkw(x,2.0_dp,3.0_dp),p,2.0e-12_dp,'kw cdf',fails)
   call check_close(qkw(p,2.0_dp,3.0_dp),x,2.0e-11_dp,'kw quantile',fails)

   x=0.37_dp
   call check_close(dgkw(x,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp), &
      2.0820234137809606_dp,2.0e-11_dp,'gkw density reference',fails)
   p=pgkw(x,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp)
   call check_close(p,0.39173287381033073_dp,2.0e-11_dp,'gkw cdf reference',fails)
   q=qgkw(p,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp)
   call check_close(q,x,5.0e-11_dp,'gkw inversion',fails)

   call check_close(dgkw(x,1.7_dp,2.4_dp,1.0_dp,0.0_dp,1.0_dp), &
      dkw(x,1.7_dp,2.4_dp),1.0e-12_dp,'kw nesting',fails)
   call check_close(dgkw(x,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.0_dp), &
      dbkw(x,1.7_dp,2.4_dp,1.3_dp,0.8_dp),1.0e-12_dp,'bkw nesting',fails)
   call check_close(dgkw(x,1.7_dp,2.4_dp,1.0_dp,0.8_dp,1.2_dp), &
      dkkw(x,1.7_dp,2.4_dp,0.8_dp,1.2_dp),1.0e-12_dp,'kkw nesting',fails)
   call check_close(dgkw(x,1.7_dp,2.4_dp,1.0_dp,0.0_dp,1.2_dp), &
      dekw(x,1.7_dp,2.4_dp,1.2_dp),1.0e-12_dp,'ekw nesting',fails)
   call check_close(dgkw(x,1.0_dp,1.0_dp,1.3_dp,0.8_dp,1.2_dp), &
      dmc(x,1.3_dp,0.8_dp,1.2_dp),1.0e-12_dp,'mc nesting',fails)
   call check_close(dgkw(x,1.0_dp,1.0_dp,1.3_dp,0.8_dp,1.0_dp), &
      dbeta_(x,1.3_dp,0.8_dp),1.0e-12_dp,'beta nesting',fails)

   call check_close(exp(pgkw(x,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp,log_p=.true.)), &
      p,2.0e-14_dp,'log cdf',fails)
   call check_close(pgkw(x,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp)+ &
      pgkw(x,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp,lower_tail=.false.), &
      1.0_dp,3.0e-14_dp,'tail complement',fails)

   if(fails==0) then
      print '(a)','test_distributions: PASS'
   else
      print '(a,i0)','test_distributions: FAIL ',fails
      error stop 1
   end if
contains
   subroutine check_close(a,b,tol,label,fails)
      real(dp),intent(in)::a,b,tol
      character(len=*),intent(in)::label
      integer,intent(inout)::fails
      if(abs(a-b)>tol .or. a/=a) then
         print '(a,2es24.15)','FAIL '//trim(label)//': ',a,b
         fails=fails+1
      end if
   end subroutine
end program test_distributions
