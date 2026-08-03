! SPDX-License-Identifier: GPL-3.0-or-later
program test_pin_factorizations
   use pinstimation
   implicit none
   type(trade_counts) :: data
   type(pin_parameters) :: p
   real(dp), allocatable :: posterior(:,:)
   real(dp) :: l0, le, llk, leho
   integer :: i

   data%buys = [12_i8, 20_i8, 7_i8, 30_i8, 15_i8, 8_i8]
   data%sells = [10_i8, 8_i8, 18_i8, 12_i8, 14_i8, 22_i8]
   p = pin_parameters(0.35_dp, 0.45_dp, 9.0_dp, 11.0_dp, 10.0_dp)

   l0 = pin_loglik(data, p)
   le = pin_loglik_e(data, p)
   llk = pin_loglik_lk(data, p)
   leho = pin_loglik_eho(data, p)
   call assert_close(l0, le, 1.0e-10_dp, 'direct and E likelihoods')
   call assert_close(l0, llk, 1.0e-10_dp, 'direct and LK likelihoods')
   call assert_close(l0, leho, 1.0e-9_dp, 'direct and EHO likelihoods')
   call pin_posteriors(data, p, posterior)
   do i = 1, size(posterior,1)
      call assert_close(sum(posterior(i,:)), 1.0_dp, 1.0e-12_dp, 'posterior normalization')
   end do
   call assert_close(pin_value(p), p%alpha*p%mu/(p%alpha*p%mu+p%eps_b+p%eps_s), 1.0e-14_dp, 'PIN formula')
   print '(a)', 'test_pin_factorizations: PASS'
contains
   subroutine assert_close(x,y,tol,message)
      real(dp),intent(in)::x,y,tol
      character(len=*),intent(in)::message
      if (abs(x-y)>tol*(1.0_dp+abs(y))) then
         write(*,*) 'FAIL: ',trim(message),x,y
         error stop 1
      end if
   end subroutine assert_close
end program test_pin_factorizations
