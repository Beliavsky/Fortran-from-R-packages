! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
program test_factor_simulation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ftsa, only : dp, hurst_result, ftsa_er_gr, ftsa_h_estimate_rs
   use ftsa, only : ftsa_sim_farma, ftsa_stop_time_sim_data
   implicit none
   type(hurst_result) :: hres
   real(dp) :: x(30, 20), fun(5, 20)
   real(dp), allocatable :: eig(:), sim(:, :), grid(:), stopped(:, :)
   integer :: i, j, k_er, k_gr, status

   do i = 1, 30
      do j = 1, 20
         x(i, j) = sin(6.283185307179586_dp * real(i * j, dp) / 31.0_dp) + &
            0.2_dp * cos(0.13_dp * real(i + 2 * j, dp))
      end do
   end do
   call ftsa_er_gr(x, k_er, k_gr, eig, status)
   call check(status == 0 .and. k_er >= 1 .and. k_gr >= 1, 'ER_GR')
   do j = 1, 20
      do i = 1, 5
         fun(i, j) = sin(0.25_dp * real(j, dp)) * real(i, dp) + 0.1_dp * cos(real(i * j, dp))
      end do
   end do
   call ftsa_h_estimate_rs(fun, hres, band_const=0.5_dp, status=status)
   call check(status == 0 .and. hres%ncomp >= 1 .and. all(ieee_is_finite(hres%score)), 'H estimate')
   call ftsa_sim_farma(6, 0.1_dp, 8, 7, .false., sim, status)
   call check(status == 0 .and. all(shape(sim) == [8, 6]) .and. all(ieee_is_finite(sim)), 'FAR simulation')
   call ftsa_stop_time_sim_data(8, 0.2_dp, 9, grid, stopped, status)
   call check(status == 0 .and. all(shape(stopped) == [101, 7]), 'stop-time simulation')
   call check(abs(grid(1)) < 1.0e-12_dp .and. abs(grid(101) - 1.0_dp) < 1.0e-12_dp, 'stop-time grid')

contains
   subroutine check(condition, message)
      logical, intent(in) :: condition !! Assertion condition expected to be true.
      character(len=*), intent(in) :: message !! Diagnostic emitted if the assertion fails.
      if (.not. condition) error stop message
   end subroutine check
end program test_factor_simulation
