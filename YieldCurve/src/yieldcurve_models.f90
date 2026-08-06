! SPDX-License-Identifier: GPL-2.0-or-later
module yieldcurve_models
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use yieldcurve_kinds, only : dp
   use yieldcurve_status, only : yc_success, yc_invalid_argument, yc_dimension_error, &
      yc_rank_deficient, yc_no_solution, set_yc_status
   use yieldcurve_factors, only : factor_beta1, factor_beta2, beta1_spot, beta2_spot, &
      beta1_forward, beta2_forward
   use yieldcurve_linalg, only : least_squares
   use yieldcurve_optimization, only : maximize_factor_beta2, maximize_beta2_spot
   implicit none
   private

   public :: nelson_siegel_fit, svensson_fit
   public :: ns_rates, svensson_rates

   interface nelson_siegel_fit
      module procedure nelson_siegel_fit_vector
      module procedure nelson_siegel_fit_matrix
   end interface nelson_siegel_fit

   interface svensson_fit
      module procedure svensson_fit_vector
      module procedure svensson_fit_matrix
   end interface svensson_fit

   interface ns_rates
      module procedure ns_rates_vector
      module procedure ns_rates_matrix
   end interface ns_rates

   interface svensson_rates
      module procedure svensson_rates_vector
      module procedure svensson_rates_matrix
   end interface svensson_rates

contains

   subroutine ns_rates_vector(coeff, maturity, curve, stat, message)
      real(dp), intent(in) :: coeff(:), maturity(:)
      real(dp), intent(out) :: curve(:)
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      if (size(coeff) /= 4 .or. size(curve) /= size(maturity)) then
         call set_yc_status(yc_dimension_error, 'NS coefficients require four values and matching output.', stat, message)
         return
      end if
      if (.not. valid_maturities(maturity) .or. .not. all(ieee_is_finite(coeff))) then
         call set_yc_status(yc_invalid_argument, 'Maturities must be positive/increasing and coefficients finite.', stat, message)
         return
      end if

      curve = coeff(1) + coeff(2) * factor_beta1(coeff(4), maturity) + &
         coeff(3) * factor_beta2(coeff(4), maturity)
      call set_yc_status(yc_success, '', stat, message)
   end subroutine ns_rates_vector

   subroutine ns_rates_matrix(coeff, maturity, curve, stat, message)
      real(dp), intent(in) :: coeff(:, :), maturity(:)
      real(dp), intent(out) :: curve(:, :)
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      integer :: i, local_stat
      character(len=160) :: local_message

      if (size(coeff, 2) /= 4 .or. size(curve, 1) /= size(coeff, 1) .or. &
          size(curve, 2) /= size(maturity)) then
         call set_yc_status(yc_dimension_error, 'NS coefficient/output matrix dimensions are inconsistent.', stat, message)
         return
      end if

      do i = 1, size(coeff, 1)
         call ns_rates_vector(coeff(i, :), maturity, curve(i, :), local_stat, local_message)
         if (local_stat /= yc_success) then
            call set_yc_status(local_stat, trim(local_message), stat, message)
            return
         end if
      end do
      call set_yc_status(yc_success, '', stat, message)
   end subroutine ns_rates_matrix

   subroutine svensson_rates_vector(coeff, maturity, curve, rate_type, stat, message)
      real(dp), intent(in) :: coeff(:), maturity(:)
      real(dp), intent(out) :: curve(:)
      character(len=*), intent(in), optional :: rate_type
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      character(len=16) :: kind

      if (size(coeff) /= 6 .or. size(curve) /= size(maturity)) then
         call set_yc_status(yc_dimension_error, 'Svensson coefficients require six values and matching output.', stat, message)
         return
      end if
      if (.not. valid_maturities(maturity) .or. .not. all(ieee_is_finite(coeff)) .or. &
          coeff(5) <= 0.0_dp .or. coeff(6) <= 0.0_dp) then
         call set_yc_status(yc_invalid_argument, 'Maturities and tau values must be positive and finite.', stat, message)
         return
      end if

      kind = 'forward'
      if (present(rate_type)) kind = lower_ascii(adjustl(rate_type))

      select case (trim(kind))
      case ('forward')
         curve = coeff(1) + coeff(2) * beta1_forward(maturity, coeff(5)) + &
            coeff(3) * beta2_forward(maturity, coeff(5)) + &
            coeff(4) * beta2_forward(maturity, coeff(6))
      case ('spot')
         curve = coeff(1) + coeff(2) * beta1_spot(maturity, coeff(5)) + &
            coeff(3) * beta2_spot(maturity, coeff(5)) + &
            coeff(4) * beta2_spot(maturity, coeff(6))
      case default
         call set_yc_status(yc_invalid_argument, 'rate_type must be "spot" or "forward".', stat, message)
         return
      end select
      call set_yc_status(yc_success, '', stat, message)
   end subroutine svensson_rates_vector

   subroutine svensson_rates_matrix(coeff, maturity, curve, rate_type, stat, message)
      real(dp), intent(in) :: coeff(:, :), maturity(:)
      real(dp), intent(out) :: curve(:, :)
      character(len=*), intent(in), optional :: rate_type
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      integer :: i, local_stat
      character(len=160) :: local_message

      if (size(coeff, 2) /= 6 .or. size(curve, 1) /= size(coeff, 1) .or. &
          size(curve, 2) /= size(maturity)) then
         call set_yc_status(yc_dimension_error, 'Svensson coefficient/output matrix dimensions are inconsistent.', stat, message)
         return
      end if

      do i = 1, size(coeff, 1)
         call svensson_rates_vector(coeff(i, :), maturity, curve(i, :), rate_type, local_stat, local_message)
         if (local_stat /= yc_success) then
            call set_yc_status(local_stat, trim(local_message), stat, message)
            return
         end if
      end do
      call set_yc_status(yc_success, '', stat, message)
   end subroutine svensson_rates_matrix

   subroutine nelson_siegel_fit_vector(rate, maturity, coeff, stat, message)
      real(dp), intent(in) :: rate(:), maturity(:)
      real(dp), intent(out) :: coeff(:)
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: lambda_values(:), residual(:)
      real(dp) :: beta(3), candidate(4), lambda_temp, ssr, best_ssr
      integer :: i, local_stat, valid_count

      coeff = 0.0_dp
      if (size(coeff) /= 4 .or. size(rate) /= size(maturity)) then
         call set_yc_status(yc_dimension_error, 'Rate, maturity, and four-coefficient output dimensions must match.', stat, message)
         return
      end if
      if (size(rate) < 3 .or. .not. valid_maturities(maturity)) then
         call set_yc_status(yc_invalid_argument, &
            'Nelson-Siegel fitting needs at least three positive increasing maturities.', &
            stat, message)
         return
      end if

      call make_sequence(maturity(1), maturity(size(maturity)), 0.5_dp, lambda_values)
      allocate(residual(size(rate)))
      best_ssr = huge(1.0_dp)
      valid_count = 0

      do i = 1, size(lambda_values)
         lambda_temp = maximize_factor_beta2(lambda_values(i), 0.001_dp, 1.0_dp)
         call ns_estimator(rate, maturity, lambda_temp, beta, residual, local_stat)
         if (local_stat == yc_success .and. beta(1) > 0.0_dp .and. beta(1) < 20.0_dp) then
            ssr = sum(residual**2, mask=ieee_is_finite(residual))
            candidate = [beta, lambda_temp]
            valid_count = valid_count + 1
         else
            ssr = 1.0e5_dp
            candidate = [beta, lambda_values(i)]
         end if
         if (ssr < best_ssr) then
            best_ssr = ssr
            coeff = candidate
         end if
      end do

      if (valid_count == 0) then
         call set_yc_status(yc_no_solution, 'No full-rank Nelson-Siegel candidate satisfied 0 < beta_0 < 20.', stat, message)
      else
         call set_yc_status(yc_success, '', stat, message)
      end if
   end subroutine nelson_siegel_fit_vector

   subroutine nelson_siegel_fit_matrix(rate, maturity, coeff, stat, message)
      real(dp), intent(in) :: rate(:, :), maturity(:)
      real(dp), intent(out) :: coeff(:, :)
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      integer :: i, local_stat
      character(len=160) :: local_message

      if (size(rate, 2) /= size(maturity) .or. size(coeff, 1) /= size(rate, 1) .or. &
          size(coeff, 2) /= 4) then
         call set_yc_status(yc_dimension_error, 'Nelson-Siegel matrix dimensions are inconsistent.', stat, message)
         return
      end if

      do i = 1, size(rate, 1)
         call nelson_siegel_fit_vector(rate(i, :), maturity, coeff(i, :), local_stat, local_message)
         if (local_stat /= yc_success) then
            call set_yc_status(local_stat, trim(local_message), stat, message)
            return
         end if
      end do
      call set_yc_status(yc_success, '', stat, message)
   end subroutine nelson_siegel_fit_matrix

   subroutine svensson_fit_vector(rate, maturity, coeff, stat, message)
      real(dp), intent(in) :: rate(:), maturity(:)
      real(dp), intent(out) :: coeff(:)
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: tau1_values(:), tau2_values(:), residual(:)
      real(dp) :: beta(4), inner_best(7), candidate(7)
      real(dp) :: tau1_temp, tau2_temp, ssr, inner_best_ssr, best_ssr, middle
      integer :: i, j, local_stat, valid_count

      coeff = 0.0_dp
      if (size(coeff) /= 6 .or. size(rate) /= size(maturity)) then
         call set_yc_status(yc_dimension_error, 'Rate, maturity, and six-coefficient output dimensions must match.', stat, message)
         return
      end if
      if (size(rate) < 4 .or. .not. valid_maturities(maturity)) then
         call set_yc_status(yc_invalid_argument, &
            'Svensson fitting needs at least four positive increasing maturities.', &
            stat, message)
         return
      end if

      middle = median_sorted(maturity)
      call make_sequence(maturity(1), middle, 1.0_dp, tau1_values)
      call make_sequence(middle, maturity(size(maturity)), 1.5_dp, tau2_values)
      allocate(residual(size(rate)))
      best_ssr = huge(1.0_dp)
      valid_count = 0

      do i = 1, size(tau1_values)
         tau1_temp = maximize_beta2_spot(tau1_values(i), 0.001_dp, maxval(tau1_values))
         inner_best_ssr = huge(1.0_dp)
         inner_best = 0.0_dp

         do j = 1, size(tau2_values)
            tau2_temp = maximize_beta2_spot(tau2_values(j), 0.001_dp, maturity(size(maturity)))
            call nss_estimator(rate, maturity, tau1_temp, tau2_temp, beta, residual, local_stat)
            if (local_stat == yc_success) then
               ssr = sum(residual**2, mask=ieee_is_finite(residual))
               candidate = [beta, tau1_temp, tau2_temp, ssr]
               valid_count = valid_count + 1
            else
               ssr = huge(1.0_dp)
               candidate = [beta, tau1_temp, tau2_temp, ssr]
            end if
            if (ssr < inner_best_ssr) then
               inner_best_ssr = ssr
               inner_best = candidate
            end if
         end do

         if (inner_best_ssr < best_ssr) then
            best_ssr = inner_best_ssr
            coeff = inner_best(1:6)
         end if
      end do

      if (valid_count == 0) then
         call set_yc_status(yc_no_solution, 'No full-rank Svensson candidate was available.', stat, message)
      else
         call set_yc_status(yc_success, '', stat, message)
      end if
   end subroutine svensson_fit_vector

   subroutine svensson_fit_matrix(rate, maturity, coeff, stat, message)
      real(dp), intent(in) :: rate(:, :), maturity(:)
      real(dp), intent(out) :: coeff(:, :)
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      integer :: i, local_stat
      character(len=160) :: local_message

      if (size(rate, 2) /= size(maturity) .or. size(coeff, 1) /= size(rate, 1) .or. &
          size(coeff, 2) /= 6) then
         call set_yc_status(yc_dimension_error, 'Svensson matrix dimensions are inconsistent.', stat, message)
         return
      end if

      do i = 1, size(rate, 1)
         call svensson_fit_vector(rate(i, :), maturity, coeff(i, :), local_stat, local_message)
         if (local_stat /= yc_success) then
            call set_yc_status(local_stat, trim(local_message), stat, message)
            return
         end if
      end do
      call set_yc_status(yc_success, '', stat, message)
   end subroutine svensson_fit_matrix

   subroutine ns_estimator(rate, maturity, lambda, beta, residual, stat)
      real(dp), intent(in) :: rate(:), maturity(:), lambda
      real(dp), intent(out) :: beta(3), residual(:)
      integer, intent(out) :: stat

      real(dp), allocatable :: a(:, :), b(:), fit_residual(:), x(:)
      integer, allocatable :: index(:)
      integer :: i, nvalid

      beta = 0.0_dp
      residual = 0.0_dp
      nvalid = count(ieee_is_finite(rate))
      if (nvalid < 3) then
         stat = yc_rank_deficient
         return
      end if

      allocate(index(nvalid), a(nvalid, 3), b(nvalid), fit_residual(nvalid), x(3))
      index = pack([(i, i = 1, size(rate))], ieee_is_finite(rate))
      do i = 1, nvalid
         a(i, :) = [1.0_dp, factor_beta1(lambda, maturity(index(i))), &
            factor_beta2(lambda, maturity(index(i)))]
         b(i) = rate(index(i))
      end do
      call least_squares(a, b, x, fit_residual, stat)
      if (stat /= yc_success) return
      beta = x
      residual = 0.0_dp
      do i = 1, nvalid
         residual(index(i)) = fit_residual(i)
      end do
   end subroutine ns_estimator

   subroutine nss_estimator(rate, maturity, tau1, tau2, beta, residual, stat)
      real(dp), intent(in) :: rate(:), maturity(:), tau1, tau2
      real(dp), intent(out) :: beta(4), residual(:)
      integer, intent(out) :: stat

      real(dp), allocatable :: a(:, :), b(:), fit_residual(:), x(:)
      integer, allocatable :: index(:)
      integer :: i, nvalid

      beta = 0.0_dp
      residual = 0.0_dp
      nvalid = count(ieee_is_finite(rate))
      if (nvalid < 4) then
         stat = yc_rank_deficient
         return
      end if

      allocate(index(nvalid), a(nvalid, 4), b(nvalid), fit_residual(nvalid), x(4))
      index = pack([(i, i = 1, size(rate))], ieee_is_finite(rate))
      do i = 1, nvalid
         a(i, :) = [1.0_dp, beta1_spot(maturity(index(i)), tau1), &
            beta2_spot(maturity(index(i)), tau1), beta2_spot(maturity(index(i)), tau2)]
         b(i) = rate(index(i))
      end do
      call least_squares(a, b, x, fit_residual, stat)
      if (stat /= yc_success) return
      beta = x
      residual = 0.0_dp
      do i = 1, nvalid
         residual(index(i)) = fit_residual(i)
      end do
   end subroutine nss_estimator

   pure logical function valid_maturities(maturity) result(valid)
      real(dp), intent(in) :: maturity(:)
      integer :: i

      valid = size(maturity) > 0 .and. all(ieee_is_finite(maturity)) .and. all(maturity > 0.0_dp)
      if (.not. valid) return
      do i = 2, size(maturity)
         if (maturity(i) <= maturity(i - 1)) then
            valid = .false.
            return
         end if
      end do
   end function valid_maturities

   pure function median_sorted(values) result(median)
      real(dp), intent(in) :: values(:)
      real(dp) :: median
      integer :: n

      n = size(values)
      if (mod(n, 2) == 1) then
         median = values((n + 1) / 2)
      else
         median = 0.5_dp * (values(n / 2) + values(n / 2 + 1))
      end if
   end function median_sorted

   subroutine make_sequence(first, last, step, values)
      real(dp), intent(in) :: first, last, step
      real(dp), allocatable, intent(out) :: values(:)
      integer :: i, n

      n = int(floor((last - first) / step + 1.0e-12_dp)) + 1
      n = max(1, n)
      allocate(values(n))
      do i = 1, n
         values(i) = first + real(i - 1, dp) * step
      end do
   end subroutine make_sequence

   pure function lower_ascii(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code

      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lower_ascii

end module yieldcurve_models
