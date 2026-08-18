program test_multivariate
   use biasedurn
   implicit none
   integer :: m2(2), x2(2), comb, failures
   integer :: m3(3), x3(3), a, b, c
   real(dp) :: w2(2), w3(3), mean2(2), var2(2), s

   failures = 0
   m2 = [10, 15]
   x2 = [4, 4]
   w2 = [2.5_dp, 1.0_dp]

   call check_close(dmfnchypergeo(x2, m2, 8, w2), &
      dfnchypergeo(4, 10, 15, 8, 2.5_dp), 5.0e-12_dp, &
      'multivariate Fisher two-color reduction', failures)
   call check_close(dmwnchypergeo(x2, m2, 8, w2), &
      dwnchypergeo(4, 10, 15, 8, 2.5_dp), 5.0e-12_dp, &
      'multivariate Wallenius two-color reduction', failures)

   call momentsmfnchypergeo(m2, 8, w2, mean2, var2, combinations=comb)
   call check_close(mean2(1), meanfnchypergeo(10, 15, 8, 2.5_dp), &
      5.0e-11_dp, 'multivariate Fisher mean', failures)
   call check_close(var2(1), varfnchypergeo(10, 15, 8, 2.5_dp), &
      5.0e-11_dp, 'multivariate Fisher variance', failures)
   if (comb /= 9) failures = failures + 1

   call momentsmwnchypergeo(m2, 8, w2, mean2, var2, combinations=comb)
   call check_close(mean2(1), meanwnchypergeo(10, 15, 8, 2.5_dp), &
      5.0e-10_dp, 'multivariate Wallenius mean', failures)
   call check_close(var2(1), varwnchypergeo(10, 15, 8, 2.5_dp), &
      5.0e-10_dp, 'multivariate Wallenius variance', failures)

   ! Three-color probability normalization on a small support.
   m3 = [4, 5, 3]
   w3 = [0.7_dp, 2.0_dp, 1.3_dp]
   s = 0.0_dp
   do a = 0, min(m3(1), 5)
      do b = 0, min(m3(2), 5 - a)
         c = 5 - a - b
         if (c < 0 .or. c > m3(3)) cycle
         x3 = [a, b, c]
         s = s + dmfnchypergeo(x3, m3, 5, w3)
      end do
   end do
   call check_close(s, 1.0_dp, 1.0e-11_dp, &
      'multivariate Fisher normalization', failures)

   s = 0.0_dp
   do a = 0, min(m3(1), 5)
      do b = 0, min(m3(2), 5 - a)
         c = 5 - a - b
         if (c < 0 .or. c > m3(3)) cycle
         x3 = [a, b, c]
         s = s + dmwnchypergeo(x3, m3, 5, w3)
      end do
   end do
   call check_close(s, 1.0_dp, 2.0e-9_dp, &
      'multivariate Wallenius normalization', failures)

   if (failures == 0) then
      print *, 'test_multivariate: PASS'
   else
      print *, 'test_multivariate: FAIL', failures
      error stop 1
   end if

contains
   subroutine check_close(got, expected, atol, label, failures)
      real(dp), intent(in) :: got, expected, atol
      character(*), intent(in) :: label
      integer, intent(inout) :: failures
      if (abs(got - expected) > atol) then
         print *, trim(label), ' got=', got, ' expected=', expected
         failures = failures + 1
      end if
   end subroutine check_close
end program test_multivariate
