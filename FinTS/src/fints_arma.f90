! SPDX-License-Identifier: GPL-2.0-or-later
module fints_arma
   use fints_kinds, only : dp, pi_dp
   use fints_status, only : fints_ok, fints_invalid_input, fints_nonstationary, &
      fints_iteration_limit, fints_numerical_failure
   use fints_types, only : arma_acf_result
   use fints_time_series, only : pacf_from_acf
   implicit none
   private
   public :: arma_true_acf, find_conjugates, polynomial_roots

contains

   subroutine arma_true_acf(ar, ma, lag_max, result, partial, complex_tolerance)
      real(dp), intent(in) :: ar(:), ma(:)
      integer, intent(in) :: lag_max
      type(arma_acf_result), intent(out) :: result
      logical, intent(in), optional :: partial
      real(dp), intent(in), optional :: complex_tolerance
      real(dp), allocatable :: polynomial(:), psi(:), correlations(:), partial_values(:)
      complex(dp), allocatable :: raw_roots(:), conjugates(:)
      real(dp) :: tolerance, energy, tail_scale, root_modulus
      integer :: p, q, degree, stat, i, j, k, max_terms, last_term, quiet_count
      logical :: want_partial

      result = arma_acf_result()
      p = size(ar)
      q = size(ma)
      if (lag_max < 0) then
         result%status = fints_invalid_input
         return
      end if
      want_partial = .false.
      if (present(partial)) want_partial = partial
      tolerance = 1000.0_dp * epsilon(1.0_dp)
      if (present(complex_tolerance)) tolerance = complex_tolerance

      if (p > 0) then
         degree = p
         allocate(polynomial(0:degree))
         polynomial(0) = 1.0_dp
         polynomial(1:degree) = -ar
         call polynomial_roots(polynomial, raw_roots, stat)
         if (stat /= fints_ok) then
            result%status = stat
            return
         end if
         allocate(result%roots(size(raw_roots)))
         do i = 1, size(raw_roots)
            if (abs(raw_roots(i)) <= tiny(1.0_dp)) then
               result%roots(i) = cmplx(huge(1.0_dp), 0.0_dp, kind=dp)
            else
               result%roots(i) = 1.0_dp / raw_roots(i)
            end if
         end do
         call sort_complex_roots(result%roots)
         result%stationary = all(abs(result%roots) < 1.0_dp - 100.0_dp * epsilon(1.0_dp))
      else
         allocate(result%roots(0))
         result%stationary = .true.
      end if

      if (.not. result%stationary) then
         result%status = fints_nonstationary
         return
      end if

      max_terms = max(20000, lag_max + 5000)
      allocate(psi(0:max_terms))
      psi = 0.0_dp
      psi(0) = 1.0_dp
      quiet_count = 0
      energy = 1.0_dp
      last_term = max(lag_max, q)
      do j = 1, max_terms
         if (j <= q) psi(j) = ma(j)
         do i = 1, min(p, j)
            psi(j) = psi(j) + ar(i) * psi(j - i)
         end do
         energy = energy + psi(j) ** 2
         tail_scale = sqrt(energy) * 1.0e-13_dp
         if (j > lag_max + 100 .and. abs(psi(j)) <= tail_scale) then
            quiet_count = quiet_count + 1
         else
            quiet_count = 0
         end if
         last_term = j
         if (quiet_count >= 100) exit
      end do
      if (last_term == max_terms .and. abs(psi(last_term)) > 1.0e-8_dp * sqrt(energy)) then
         result%status = fints_iteration_limit
         return
      end if

      allocate(correlations(0:lag_max))
      do k = 0, lag_max
         correlations(k) = dot_product(psi(0:last_term - k), psi(k:last_term))
      end do
      if (correlations(0) <= 0.0_dp) then
         result%status = fints_numerical_failure
         return
      end if
      correlations = correlations / correlations(0)

      result%partial = want_partial
      if (want_partial) then
         if (lag_max == 0) then
            allocate(result%value(0))
         else
            call pacf_from_acf(correlations, partial_values, stat)
            if (stat /= fints_ok) then
               result%status = stat
               return
            end if
            allocate(result%value(lag_max))
            result%value = partial_values
         end if
      else
         allocate(result%value(lag_max + 1))
         result%value = correlations
      end if

      call find_conjugates(result%roots, conjugates, tolerance)
      allocate(result%damping(size(conjugates)), result%period(size(conjugates)))
      do i = 1, size(conjugates)
         root_modulus = abs(conjugates(i))
         result%damping(i) = root_modulus
         if (root_modulus > 0.0_dp) then
            result%period(i) = 2.0_dp * pi_dp / &
               acos(max(-1.0_dp, min(1.0_dp, real(conjugates(i), dp) / root_modulus)))
         else
            result%period(i) = huge(1.0_dp)
         end if
      end do
      result%status = fints_ok
   end subroutine arma_true_acf

   subroutine find_conjugates(x, representatives, complex_tolerance)
      complex(dp), intent(in) :: x(:)
      complex(dp), allocatable, intent(out) :: representatives(:)
      real(dp), intent(in), optional :: complex_tolerance
      complex(dp), allocatable :: work(:)
      real(dp) :: tolerance, scale
      integer :: i, j, count
      logical :: conjugate_pair, distinct

      tolerance = epsilon(1.0_dp)
      if (present(complex_tolerance)) tolerance = complex_tolerance
      allocate(work(size(x)))
      count = 0
      do i = 2, size(x)
         do j = 1, i - 1
            scale = max(abs(x(i)), abs(x(j)))
            if (scale <= 0.0_dp) cycle
            conjugate_pair = abs(x(i) - conjg(x(j))) / scale < tolerance
            distinct = abs(x(i) - x(j)) / scale > tolerance
            if (conjugate_pair .and. distinct) then
               count = count + 1
               work(count) = x(i)
               exit
            end if
         end do
      end do
      allocate(representatives(count))
      if (count > 0) representatives = work(1:count)
   end subroutine find_conjugates

   subroutine polynomial_roots(coefficients, roots, status, tolerance, max_iterations)
      real(dp), intent(in) :: coefficients(0:)
      complex(dp), allocatable, intent(out) :: roots(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      complex(dp) :: numerator, denominator, update
      real(dp) :: tol, radius, max_update
      integer :: degree, iterations, iter, i, j

      degree = ubound(coefficients, 1)
      allocate(roots(max(0, degree)))
      if (degree < 1 .or. abs(coefficients(degree)) <= tiny(1.0_dp)) then
         status = fints_invalid_input
         return
      end if
      if (degree == 1) then
         roots(1) = cmplx(-coefficients(0) / coefficients(1), 0.0_dp, kind=dp)
         status = fints_ok
         return
      end if

      tol = 1.0e-13_dp
      if (present(tolerance)) tol = tolerance
      iterations = 5000
      if (present(max_iterations)) iterations = max_iterations
      radius = 1.0_dp + maxval(abs(coefficients(0:degree - 1) / coefficients(degree)))
      do i = 1, degree
         roots(i) = radius * exp(cmplx(0.0_dp, &
            2.0_dp * pi_dp * (real(i - 1, dp) + 0.25_dp) / real(degree, dp), kind=dp))
      end do

      status = fints_iteration_limit
      do iter = 1, iterations
         max_update = 0.0_dp
         do i = 1, degree
            numerator = evaluate_polynomial(coefficients, roots(i))
            denominator = cmplx(1.0_dp, 0.0_dp, kind=dp)
            do j = 1, degree
               if (j /= i) denominator = denominator * (roots(i) - roots(j))
            end do
            if (abs(denominator) <= tiny(1.0_dp)) then
               roots(i) = roots(i) + cmplx(tol * real(i, dp), tol, kind=dp)
               cycle
            end if
            update = numerator / denominator
            roots(i) = roots(i) - update
            max_update = max(max_update, abs(update))
         end do
         if (max_update <= tol * max(1.0_dp, maxval(abs(roots)))) then
            status = fints_ok
            exit
         end if
      end do
   end subroutine polynomial_roots

   pure complex(dp) function evaluate_polynomial(coefficients, z) result(value)
      real(dp), intent(in) :: coefficients(0:)
      complex(dp), intent(in) :: z
      integer :: i

      value = cmplx(coefficients(ubound(coefficients, 1)), 0.0_dp, kind=dp)
      do i = ubound(coefficients, 1) - 1, 0, -1
         value = value * z + coefficients(i)
      end do
   end function evaluate_polynomial

   subroutine sort_complex_roots(roots)
      complex(dp), intent(inout) :: roots(:)
      complex(dp) :: temp
      integer :: i, j, k

      do i = 1, size(roots) - 1
         k = i
         do j = i + 1, size(roots)
            if (root_less(roots(j), roots(k))) k = j
         end do
         if (k /= i) then
            temp = roots(i)
            roots(i) = roots(k)
            roots(k) = temp
         end if
      end do
   end subroutine sort_complex_roots

   pure logical function root_less(a, b) result(is_less)
      complex(dp), intent(in) :: a, b
      real(dp) :: aa, ab

      aa = abs(a)
      ab = abs(b)
      if (abs(aa - ab) > 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, aa, ab)) then
         is_less = aa < ab
      else
         is_less = aimag(a) < aimag(b)
      end if
   end function root_less

end module fints_arma
