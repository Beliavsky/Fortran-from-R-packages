! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_simulation_mod
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fracdiff_kinds, only : dp
   use fracdiff_status, only : fd_ok, fd_invalid_input, fd_gamma_error, fd_unstable_ar, &
                               fracdiff_status_message
   use fracdiff_types, only : fracdiff_simulation
   use fracdiff_rng, only : fracdiff_rng_state, seed_rng, fill_normal
   use fracdiff_polynomial, only : minimum_ar_root_modulus
   implicit none
   private

   public :: fracdiff_sim, fractional_arma_filter_simulation

contains

   subroutine fractional_arma_filter_simulation(innovations, ar, ma, d, mean_value, series, status)
      real(dp), intent(in) :: innovations(:)
      real(dp), intent(in) :: ar(:), ma(:)
      real(dp), intent(in) :: d, mean_value
      real(dp), intent(out) :: series(:)
      integer, intent(out) :: status

      real(dp), allocatable :: y(:), phi(:)
      real(dp) :: vk, temp, amk, dk1, dk1d, dj, value
      integer :: n, p, q, total, i, j, k, jmax

      n = size(series)
      p = size(ar)
      q = size(ma)
      total = n + q
      status = fd_ok
      series = 0.0_dp
      if (size(innovations) < total .or. n < 1 .or. d < -0.5_dp .or. d > 0.5_dp) then
         status = fd_invalid_input
         return
      end if
      if (d >= 0.5_dp) then
         status = fd_gamma_error
         return
      end if

      allocate(y(total), phi(max(1,total)))
      y = innovations(1:total)
      phi = 0.0_dp

      vk = gamma(1.0_dp - 2.0_dp*d)/(gamma(1.0_dp - d)**2)
      if (.not. ieee_is_finite(vk) .or. vk <= 0.0_dp) then
         status = fd_gamma_error
         return
      end if
      y(1) = y(1)*sqrt(vk)

      if (total >= 2) then
         temp = d/(1.0_dp - d)
         vk = vk*(1.0_dp - temp*temp)
         amk = temp*y(1)
         phi(1) = temp
         y(2) = amk + y(2)*sqrt(max(vk,0.0_dp))
      end if

      do k = 3, total
         dk1 = real(k - 1, dp)
         dk1d = dk1 - d
         do j = 1, k - 2
            dj = dk1 - real(j,dp)
            phi(j) = phi(j)*dk1*(dj - d)/(dk1d*dj)
         end do
         temp = d/dk1d
         phi(k - 1) = temp
         vk = vk*(1.0_dp - temp*temp)
         amk = 0.0_dp
         do j = 1, k - 1
            amk = amk + phi(j)*y(k - j)
         end do
         y(k) = amk + y(k)*sqrt(max(vk,0.0_dp))
      end do

      do k = 1, n
         value = 0.0_dp
         jmax = min(p, k - 1)
         do i = 1, jmax
            value = value + ar(i)*series(k - i)
         end do
         do j = 1, q
            value = value - ma(j)*y(k + q - j)
         end do
         series(k) = value + y(k + q)
      end do
      if (abs(mean_value) > 0.0_dp) series = series + mean_value
   end subroutine fractional_arma_filter_simulation

   function fracdiff_sim(n, d, ar, ma, mean_value, n_start, back_compatible, &
                         allow_zero_start, innovations, start_innovations, seed, &
                         standard_deviation) result(simulation)
      integer, intent(in) :: n
      real(dp), intent(in) :: d
      real(dp), intent(in), optional :: ar(:), ma(:)
      real(dp), intent(in), optional :: mean_value
      integer, intent(in), optional :: n_start
      logical, intent(in), optional :: back_compatible, allow_zero_start
      real(dp), intent(in), optional :: innovations(:), start_innovations(:)
      integer(int64), intent(in), optional :: seed
      real(dp), intent(in), optional :: standard_deviation
      type(fracdiff_simulation) :: simulation

      real(dp), allocatable :: ar_local(:), ma_local(:), all_innovations(:), generated(:), full_series(:)
      real(dp) :: mu, minimum_root, sd
      integer :: p, q, burn, total_n, status, root_status, start_index
      integer(int64) :: seed_value
      logical :: back_comp, allow_zero
      type(fracdiff_rng_state) :: rng

      simulation%status = fd_ok
      if (present(ar)) then
         allocate(ar_local(size(ar)))
         ar_local = ar
      else
         allocate(ar_local(0))
      end if
      if (present(ma)) then
         allocate(ma_local(size(ma)))
         ma_local = ma
      else
         allocate(ma_local(0))
      end if
      p = size(ar_local)
      q = size(ma_local)
      mu = 0.0_dp
      if (present(mean_value)) mu = mean_value
      sd = 1.0_dp
      if (present(standard_deviation)) sd = standard_deviation
      back_comp = .true.
      if (present(back_compatible)) back_comp = back_compatible
      allow_zero = .false.
      if (present(allow_zero_start)) allow_zero = allow_zero_start

      if (n < 1 .or. d < -0.5_dp .or. d > 0.5_dp .or. sd < 0.0_dp) then
         simulation%status = fd_invalid_input
         simulation%message = fracdiff_status_message(simulation%status)
         allocate(simulation%series(0), simulation%ar(p), simulation%ma(q))
         simulation%ar = ar_local
         simulation%ma = ma_local
         return
      end if

      if (present(n_start)) then
         burn = n_start
      else
         if (p > 0) then
            minimum_root = minimum_ar_root_modulus(ar_local, root_status)
            if (root_status /= 0 .or. minimum_root <= 1.0_dp) then
               minimum_root = 1.01_dp
               simulation%status = fd_unstable_ar
            end if
            burn = p + q + ceiling(6.0_dp/log(minimum_root))
         else
            burn = p + q
         end if
      end if

      if (burn < 0 .or. (burn < p + q .and. .not. allow_zero)) then
         simulation%status = fd_invalid_input
         simulation%message = fracdiff_status_message(simulation%status)
         allocate(simulation%series(0), simulation%ar(p), simulation%ma(q))
         simulation%ar = ar_local
         simulation%ma = ma_local
         return
      end if

      total_n = n + burn
      allocate(all_innovations(total_n + q), generated(total_n + q), full_series(total_n))
      seed_value = 123456789_int64
      if (present(seed)) seed_value = seed
      call seed_rng(rng, seed_value)

      if (burn > 0) then
         if (present(start_innovations)) then
            if (size(start_innovations) < burn) then
               simulation%status = fd_invalid_input
               simulation%message = "start_innovations is shorter than n_start"
               allocate(simulation%series(0), simulation%ar(p), simulation%ma(q))
               simulation%ar = ar_local
               simulation%ma = ma_local
               return
            end if
            all_innovations(1:burn) = start_innovations(1:burn)
         else
            call fill_normal(rng, all_innovations(1:burn), sd)
         end if
      end if

      if (present(innovations)) then
         if (size(innovations) < n + q) then
            simulation%status = fd_invalid_input
            simulation%message = "innovations is shorter than n + q"
            allocate(simulation%series(0), simulation%ar(p), simulation%ma(q))
            simulation%ar = ar_local
            simulation%ma = ma_local
            return
         end if
         all_innovations(burn+1:) = innovations(1:n+q)
      else
         call fill_normal(rng, all_innovations(burn+1:), sd)
      end if
      generated = all_innovations

      call fractional_arma_filter_simulation(generated, ar_local, ma_local, d, mu, full_series, status)
      if (status /= fd_ok) then
         simulation%status = status
         simulation%message = fracdiff_status_message(status)
         allocate(simulation%series(0), simulation%ar(p), simulation%ma(q))
         simulation%ar = ar_local
         simulation%ma = ma_local
         return
      end if

      if (back_comp) then
         start_index = burn + 1
      else
         start_index = burn - q + 1
      end if
      if (start_index < 1 .or. start_index + n - 1 > total_n) then
         simulation%status = fd_invalid_input
         simulation%message = "requested compatibility slice is outside generated series"
         allocate(simulation%series(0), simulation%ar(p), simulation%ma(q))
         simulation%ar = ar_local
         simulation%ma = ma_local
         return
      end if

      allocate(simulation%series(n), simulation%ar(p), simulation%ma(q))
      simulation%series = full_series(start_index:start_index+n-1)
      simulation%ar = ar_local
      simulation%ma = ma_local
      simulation%d = d
      simulation%mean = mu
      simulation%n_start = burn
      if (simulation%status == fd_ok) then
         simulation%message = "ok"
      else
         simulation%message = fracdiff_status_message(simulation%status)
      end if
   end function fracdiff_sim

end module fracdiff_simulation_mod
