module bivkld_exact
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_positive_inf, ieee_quiet_nan, ieee_value
   use bivkld_kinds, only : dp
   implicit none
   private

   public :: biv_kld_discrete
   public :: biv_kld_discrete_matrix
   public :: biv_kld_discrete_vector
   public :: biv_kld_independent_weibull
   public :: biv_kld_normal
   public :: biv_kld_pareto2

   interface biv_kld_discrete
      module procedure biv_kld_discrete_vector
      module procedure biv_kld_discrete_matrix
   end interface biv_kld_discrete

contains

   pure function biv_kld_discrete_vector(p, q, normalize, tolerance) result(value)
      real(dp), intent(in) :: p(:) !! First probability vector; entries must be finite and non-negative.
      real(dp), intent(in) :: q(:) !! Second probability vector; it must have the same length as p.
      logical, intent(in), optional :: normalize !! If true, rescale each input to sum to one before computing KL.
      real(dp), intent(in), optional :: tolerance !! Absolute sum tolerance when normalization is disabled.
      real(dp) :: value

      logical :: do_normalize
      real(dp) :: ptot, qtot, tol, pi, qi
      integer :: i

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(p) == 0 .or. size(q) /= size(p)) return
      if (any(.not. ieee_is_finite(p)) .or. any(.not. ieee_is_finite(q))) return
      if (any(p < 0.0_dp) .or. any(q < 0.0_dp)) return

      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize
      tol = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) tol = tolerance
      if (.not. ieee_is_finite(tol) .or. tol < 0.0_dp) return

      ptot = sum(p)
      qtot = sum(q)
      if (ptot <= 0.0_dp .or. qtot <= 0.0_dp) return
      if (.not. do_normalize) then
         if (abs(ptot - 1.0_dp) > tol .or. abs(qtot - 1.0_dp) > tol) return
      end if

      value = 0.0_dp
      do i = 1, size(p)
         if (do_normalize) then
            pi = p(i)/ptot
            qi = q(i)/qtot
         else
            pi = p(i)
            qi = q(i)
         end if
         if (pi > 0.0_dp) then
            if (qi <= 0.0_dp) then
               value = ieee_value(0.0_dp, ieee_positive_inf)
               return
            end if
            value = value + pi*(log(pi) - log(qi))
         end if
      end do
      if (value < 0.0_dp .and. abs(value) < 100.0_dp*epsilon(1.0_dp)) value = 0.0_dp
   end function biv_kld_discrete_vector

   pure function biv_kld_discrete_matrix(p, q, normalize, tolerance) result(value)
      real(dp), intent(in) :: p(:, :) !! First probability table; entries must be finite and non-negative.
      real(dp), intent(in) :: q(:, :) !! Second probability table; its shape must match p exactly.
      logical, intent(in), optional :: normalize !! If true, rescale each table to sum to one before computing KL.
      real(dp), intent(in), optional :: tolerance !! Absolute sum tolerance when normalization is disabled.
      real(dp) :: value

      logical :: do_normalize
      real(dp) :: ptot, qtot, tol, pi, qi
      integer :: i, j

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(p) == 0 .or. any(shape(q) /= shape(p))) return
      if (any(.not. ieee_is_finite(p)) .or. any(.not. ieee_is_finite(q))) return
      if (any(p < 0.0_dp) .or. any(q < 0.0_dp)) return

      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize
      tol = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) tol = tolerance
      if (.not. ieee_is_finite(tol) .or. tol < 0.0_dp) return

      ptot = sum(p)
      qtot = sum(q)
      if (ptot <= 0.0_dp .or. qtot <= 0.0_dp) return
      if (.not. do_normalize) then
         if (abs(ptot - 1.0_dp) > tol .or. abs(qtot - 1.0_dp) > tol) return
      end if

      value = 0.0_dp
      do j = 1, size(p, 2)
         do i = 1, size(p, 1)
            if (do_normalize) then
               pi = p(i, j)/ptot
               qi = q(i, j)/qtot
            else
               pi = p(i, j)
               qi = q(i, j)
            end if
            if (pi > 0.0_dp) then
               if (qi <= 0.0_dp) then
                  value = ieee_value(0.0_dp, ieee_positive_inf)
                  return
               end if
               value = value + pi*(log(pi) - log(qi))
            end if
         end do
      end do
      if (value < 0.0_dp .and. abs(value) < 100.0_dp*epsilon(1.0_dp)) value = 0.0_dp
   end function biv_kld_discrete_matrix

   pure function biv_kld_normal(mean1, sigma1, mean2, sigma2) result(value)
      real(dp), intent(in) :: mean1(2) !! Mean vector of the first bivariate normal distribution.
      real(dp), intent(in) :: sigma1(2, 2) !! Symmetric positive-definite covariance matrix of the first distribution.
      real(dp), intent(in) :: mean2(2) !! Mean vector of the second bivariate normal distribution.
      real(dp), intent(in) :: sigma2(2, 2) !! Symmetric positive-definite covariance matrix of the second distribution.
      real(dp) :: value

      real(dp) :: det1, det2, inverse2(2, 2), delta(2), trace_term, quad_term

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (any(.not. ieee_is_finite(mean1)) .or. any(.not. ieee_is_finite(mean2))) return
      if (.not. is_spd_2x2(sigma1) .or. .not. is_spd_2x2(sigma2)) return

      det1 = determinant_2x2(sigma1)
      det2 = determinant_2x2(sigma2)
      call inverse_2x2(sigma2, inverse2)
      delta = mean2 - mean1
      trace_term = inverse2(1, 1)*sigma1(1, 1) + inverse2(1, 2)*sigma1(2, 1) &
                   + inverse2(2, 1)*sigma1(1, 2) + inverse2(2, 2)*sigma1(2, 2)
      quad_term = delta(1)*(inverse2(1, 1)*delta(1) + inverse2(1, 2)*delta(2)) &
                  + delta(2)*(inverse2(2, 1)*delta(1) + inverse2(2, 2)*delta(2))
      value = 0.5_dp*(log(det2) - log(det1) - 2.0_dp + trace_term + quad_term)
      if (value < 0.0_dp .and. abs(value) < 100.0_dp*epsilon(1.0_dp)) value = 0.0_dp
   end function biv_kld_normal

   pure elemental function biv_kld_pareto2(alpha, beta) result(value)
      real(dp), intent(in) :: alpha !! Positive shape parameter of the first bivariate Pareto type-II model.
      real(dp), intent(in) :: beta !! Positive shape parameter of the second bivariate Pareto type-II model.
      real(dp) :: value

      if (.not. ieee_is_finite(alpha) .or. .not. ieee_is_finite(beta) &
          .or. alpha <= 0.0_dp .or. beta <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      value = (beta - alpha)*(1.0_dp + 2.0_dp*alpha)/(alpha*(1.0_dp + alpha)) &
              + log(alpha*(1.0_dp + alpha)/(beta*(1.0_dp + beta)))
      if (value < 0.0_dp .and. abs(value) < 100.0_dp*epsilon(1.0_dp)) value = 0.0_dp
   end function biv_kld_pareto2

   pure function biv_kld_independent_weibull(alpha, beta) result(value)
      real(dp), intent(in) :: alpha(2) !! Positive rate parameters of the first independent bivariate Weibull model.
      real(dp), intent(in) :: beta(2) !! Positive rate parameters of the second independent bivariate Weibull model.
      real(dp) :: value

      if (any(.not. ieee_is_finite(alpha)) .or. any(.not. ieee_is_finite(beta)) &
          .or. any(alpha <= 0.0_dp) .or. any(beta <= 0.0_dp)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      value = sum(beta/alpha) - 2.0_dp + log(product(alpha)/product(beta))
      if (value < 0.0_dp .and. abs(value) < 100.0_dp*epsilon(1.0_dp)) value = 0.0_dp
   end function biv_kld_independent_weibull

   pure function determinant_2x2(a) result(det_a)
      real(dp), intent(in) :: a(2, 2) !! Two-by-two matrix whose determinant is required.
      real(dp) :: det_a

      det_a = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
   end function determinant_2x2

   pure subroutine inverse_2x2(a, inverse_a)
      real(dp), intent(in) :: a(2, 2) !! Nonsingular two-by-two matrix to invert.
      real(dp), intent(out) :: inverse_a(2, 2) !! Matrix inverse of a.

      real(dp) :: det_a

      det_a = determinant_2x2(a)
      inverse_a(1, 1) = a(2, 2)/det_a
      inverse_a(1, 2) = -a(1, 2)/det_a
      inverse_a(2, 1) = -a(2, 1)/det_a
      inverse_a(2, 2) = a(1, 1)/det_a
   end subroutine inverse_2x2

   pure function is_spd_2x2(a) result(ok)
      real(dp), intent(in) :: a(2, 2) !! Candidate symmetric positive-definite two-by-two matrix.
      logical :: ok

      real(dp) :: symmetry_tolerance

      symmetry_tolerance = sqrt(epsilon(1.0_dp))
      ok = all(ieee_is_finite(a))
      if (.not. ok) return
      ok = abs(a(1, 2) - a(2, 1)) <= symmetry_tolerance &
           .and. a(1, 1) > 0.0_dp .and. a(2, 2) > 0.0_dp &
           .and. determinant_2x2(a) > 0.0_dp
   end function is_spd_2x2

end module bivkld_exact
