! SPDX-License-Identifier: GPL-2.0-or-later
program test_apca
   use fints
   implicit none

   integer, parameter :: observations = 30, series = 12
   real(dp) :: x(observations, series), factor(observations), loading
   type(apca_result) :: fit
   integer :: i, j

   do i = 1, observations
      factor(i) = sin(0.21_dp * real(i, dp)) + 0.4_dp * cos(0.47_dp * real(i, dp))
   end do
   do j = 1, series
      loading = 0.3_dp + 0.12_dp * real(j, dp)
      do i = 1, observations
         x(i, j) = 0.02_dp * real(j, dp) + loading * factor(i) + &
            0.03_dp * sin(0.31_dp * real(i * j, dp))
      end do
   end do

   call apca(x, 1, fit)
   call check(fit%status == fints_ok, 'apca status')
   call check(size(fit%factors, 1) == observations .and. size(fit%factors, 2) == 1, &
      'apca factor dimensions')
   call check(size(fit%loadings, 1) == series, 'apca loading dimensions')
   call check(minval(fit%r_squared) > 0.90_dp, 'apca high r squared')
   call check(all(fit%eigenvalues(1:size(fit%eigenvalues) - 1) >= &
      fit%eigenvalues(2:size(fit%eigenvalues))), 'apca eigenvalue ordering')
   call check(abs(dot_product(fit%factors(:, 1), fit%factors(:, 1)) - 1.0_dp) < 1.0e-10_dp, &
      'apca normalized factor')

   print '(a)', 'test_apca: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_apca
