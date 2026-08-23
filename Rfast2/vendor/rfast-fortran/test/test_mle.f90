program test_mle
   use rfast
   implicit none
   real(dp) :: x(5), g(6), b(6), p(5)
   integer :: k(6)
   type(mle_result) :: fit

   x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   fit = normal_mle(x)
   call assert_close(fit%param(1),3.0_dp,1e-12_dp,'normal mean')
   call assert_close(fit%param(2),2.0_dp,1e-12_dp,'normal var')
   fit = exponential_mle(x)
   call assert_close(fit%param(1),3.0_dp,1e-12_dp,'exp scale')
   fit = laplace_mle(x)
   call assert_close(fit%param(1),3.0_dp,1e-12_dp,'laplace location')
   call assert_close(fit%param(2),1.2_dp,1e-12_dp,'laplace scale')

   g = [0.8_dp,1.1_dp,1.5_dp,2.0_dp,2.4_dp,3.2_dp]
   fit = gamma_mle(g)
   call assert_true(fit%status==0 .and. all(fit%param>0), 'gamma mle')
   b = [0.1_dp,0.2_dp,0.35_dp,0.55_dp,0.7_dp,0.85_dp]
   fit = beta_mle(b)
   call assert_true(fit%status==0 .and. all(fit%param>0), 'beta mle')
   p = [1.0_dp,1.5_dp,2.0_dp,4.0_dp,8.0_dp]
   fit = pareto_mle(p)
   call assert_close(fit%param(1),1.0_dp,1e-12_dp,'pareto scale')

   k = [0,1,2,2,3,4]
   fit = poisson_mle(k)
   call assert_close(fit%param(1),2.0_dp,1e-12_dp,'poisson mean')
   fit = geometric_mle(k)
   call assert_close(fit%param(1),1.0_dp/3.0_dp,1e-12_dp,'geometric')

   print *, 'test_mle: PASS'
contains
   subroutine assert_close(got,want,tol,msg)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::msg
      if(abs(got-want)>tol)then
         print *, 'FAIL ',trim(msg),got,want
         error stop 1
      end if
   end subroutine
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         print *, 'FAIL ',trim(msg)
         error stop 1
      end if
   end subroutine
end program test_mle
