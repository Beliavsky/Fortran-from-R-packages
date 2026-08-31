! SPDX-License-Identifier: GPL-3.0-only
program test_rng_linalg
   use pan, only : dp
   use pan_linalg, only : invwishart_draw, is_spd, spd_inverse
   use pan_rng, only : rng_seed, rng_state, rng_uniform
   implicit none

   type(rng_state) :: rng1
   type(rng_state) :: rng2
   integer :: info
   real(dp) :: a(2, 2)
   real(dp) :: ainv(2, 2)
   real(dp) :: iw(2, 2)
   real(dp) :: u1
   real(dp) :: u2

   call rng_seed(rng1, 13579)
   call rng_seed(rng2, 13579)
   u1 = rng_uniform(rng1)
   u2 = rng_uniform(rng2)
   call assert_close(u1, u2, 0.0_dp, "seeded uniforms are reproducible")
   call assert_close(u1, 0.10627426817374037_dp, 1.0e-15_dp, "Park-Miller first uniform")

   a = reshape([2.0_dp, 0.5_dp, 0.5_dp, 1.5_dp], [2, 2])
   call spd_inverse(a, ainv, info)
   call assert_true(info == 0, "SPD inverse succeeds")
   call assert_close(ainv(1, 1), 0.5454545454545454_dp, 1.0e-12_dp, "inverse (1,1)")
   call assert_close(ainv(1, 2), -0.1818181818181818_dp, 1.0e-12_dp, "inverse (1,2)")

   call invwishart_draw(rng1, a, 6.0_dp, iw, info)
   call assert_true(info == 0, "inverse-Wishart draw succeeds")
   call assert_true(is_spd(iw), "inverse-Wishart draw is positive definite")
   call assert_close(iw(1, 2), iw(2, 1), 1.0e-14_dp, "inverse-Wishart draw is symmetric")

   call finish_tests()

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Test predicate expected to be true.
      character(len=*), intent(in) :: message !! Diagnostic printed if the predicate is false.

      if (.not. condition) then
         write (*, '(a)') "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tol, message)
      real(dp), intent(in) :: actual !! Computed scalar value.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tol !! Maximum allowed absolute error.
      character(len=*), intent(in) :: message !! Diagnostic printed if the comparison fails.

      if (abs(actual - expected) > tol) then
         write (*, '(a,2(es16.8,1x))') "FAIL: " // trim(message) // " actual/ref: ", actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine finish_tests()
   end subroutine finish_tests

end program test_rng_linalg
