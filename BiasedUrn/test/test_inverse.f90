program test_inverse
   use biasedurn
   implicit none
   integer :: failures
   integer :: m(3)
   real(dp) :: mu(3), odds(3), got(3), urn(3), pair(2)

   failures = 0
   call check_close(oddsfnchypergeo(4.0_dp, 10, 15, 8), &
      4.0_dp * 11.0_dp / (6.0_dp * 4.0_dp), 1.0e-14_dp, &
      'odds Fisher', failures)
   call check_close(oddswnchypergeo(4.0_dp, 10, 15, 8), &
      log(0.6_dp) / log(11.0_dp / 15.0_dp), 1.0e-14_dp, &
      'odds Wallenius', failures)

   pair = numfnchypergeo(4.0_dp, 8, 25, 2.5_dp)
   call check_close(sum(pair), 25.0_dp, 1.0e-12_dp, 'num Fisher sum', failures)
   pair = numwnchypergeo(4.0_dp, 8, 25, 2.5_dp)
   call check_close(sum(pair), 25.0_dp, 1.0e-9_dp, 'num Wallenius sum', failures)

   m = [20, 30, 50]
   mu = [3.0_dp, 7.0_dp, 10.0_dp]
   got = oddsmfnchypergeo(mu, m, 20)
   if (.not. any(abs(got - 1.0_dp) < 1.0e-14_dp)) failures = failures + 1
   got = oddsmwnchypergeo(mu, m, 20)
   if (.not. any(abs(got - 1.0_dp) < 1.0e-14_dp)) failures = failures + 1

   odds = [0.5_dp, 1.0_dp, 2.0_dp]
   urn = nummfnchypergeo(mu, 20, 100, odds)
   call check_close(sum(urn), 100.0_dp, 1.0e-8_dp, 'num multi Fisher sum', failures)
   urn = nummwnchypergeo(mu, 20, 100, odds)
   call check_close(sum(urn), 100.0_dp, 1.0e-7_dp, 'num multi Wallenius sum', failures)

   if (failures == 0) then
      print *, 'test_inverse: PASS'
   else
      print *, 'test_inverse: FAIL', failures
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
end program test_inverse
