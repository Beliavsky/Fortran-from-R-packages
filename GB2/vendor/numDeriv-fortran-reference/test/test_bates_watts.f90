! SPDX-License-Identifier: GPL-2.0-or-later
module test_bates_watts_functions
   use numderiv, only : dp
   implicit none
   real(dp), parameter :: concentration(12) = [ &
      0.02_dp, 0.02_dp, 0.06_dp, 0.06_dp, 0.11_dp, 0.11_dp, &
      0.22_dp, 0.22_dp, 0.56_dp, 0.56_dp, 1.10_dp, 1.10_dp]
contains
   function puromycin_residuals(theta_local) result(value)
      real(dp), intent(in) :: theta_local(:)
      real(dp), allocatable :: value(:)
      real(dp), parameter :: response(12) = [ &
         76.0_dp, 47.0_dp, 97.0_dp, 107.0_dp, 123.0_dp, 139.0_dp, &
         159.0_dp, 152.0_dp, 191.0_dp, 201.0_dp, 207.0_dp, 200.0_dp]
      value = theta_local(1) * concentration / (theta_local(2) + concentration) - response
   end function puromycin_residuals
end module test_bates_watts_functions

program test_bates_watts
   use test_bates_watts_functions
   use numderiv, only : dp, deriv_options, gend_result, nd_success, gend
   implicit none

   real(dp) :: theta(2), analytic(12, 5), relative_error
   type(deriv_options) :: options
   type(gend_result) :: result
   integer :: i

   theta = [212.7_dp, 0.0641_dp]
   options = deriv_options(d=0.01_dp)
   call gend(puromycin_residuals, theta, result, options)
   call check(result%status == nd_success, 'Bates-Watts genD status')

   do i = 1, 12
      analytic(i, 1) = concentration(i) / (theta(2) + concentration(i))
      analytic(i, 2) = -theta(1) * concentration(i) / &
         (theta(2) + concentration(i)) ** 2
      analytic(i, 3) = 0.0_dp
      analytic(i, 4) = -concentration(i) / (theta(2) + concentration(i)) ** 2
      analytic(i, 5) = 2.0_dp * theta(1) * concentration(i) / &
         (theta(2) + concentration(i)) ** 3
   end do

   relative_error = maxval(abs((result%dmat - analytic) / (analytic + 1.0e-4_dp)))
   call check(relative_error < 1.0e-6_dp, 'Bates-Watts relative accuracy')

   print '(a)', 'test_bates_watts: PASS'

contains
   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check
end program test_bates_watts
