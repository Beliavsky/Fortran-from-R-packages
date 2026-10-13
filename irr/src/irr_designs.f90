! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation of irr 0.85: N.cohen.kappa, N2.cohen.kappa,
! and relInterIntra. Upstream authors: Ian Fellows; Puspendra Singh and
! Jim Lemon; Tore Wentzel-Larsen, with output formatting by Jim Lemon.
module irr_designs
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use r_kinds, only: dp
   use r_distributions, only: r_qnorm, r_pf, r_qf
   use irr_common, only: irr_ok, irr_invalid, irr_degenerate, irr_resource, max_raters, max_cells
   use lpsolve, only: solve_lp, lp_result, lp_control, LP_MAX, LP_EQ, LP_OPTIMAL, LP_INFEASIBLE, LP_TIMEOUT
   implicit none
   private
   public :: n_cohen_kappa, n2_cohen_kappa, rel_inter_intra, rel_inter_intra_result

   type :: rel_inter_intra_result
      integer :: status = irr_invalid
      integer :: subjects = 0
      integer :: raters = 0
      integer :: repetitions = 0
      real(dp), allocatable :: rohat(:) !! Inter/intra random, inter/intra fixed, then per-rater random and fixed.
      real(dp), allocatable :: fs(:) !! Interrater F, pooled intrarater F, then per-rater intrarater F.
      real(dp), allocatable :: pvalue(:) !! Upper-tail probabilities corresponding to fs.
      real(dp), allocatable :: lowvalue(:) !! Random/fixed interrater lower bounds, pooled and per-rater intra bounds.
      real(dp), allocatable :: sem(:) !! Pooled and per-rater intra SEM, then random/fixed interrater SEM.
      real(dp), allocatable :: mean_squares(:) !! MSS, MSR, MSSR, MSE, followed by each rater's MSE.
      real(dp) :: variance_random(4) = 0.0_dp !! Subject, rater, interaction, and residual variance estimates.
      real(dp) :: variance_fixed(4) = 0.0_dp !! Subject, rater, interaction, and residual variance estimates.
      real(dp) :: nu(2) = 0.0_dp !! Approximate denominator degrees of freedom for interrater bounds.
   end type rel_inter_intra_result

contains

   pure subroutine n_cohen_kappa(rate1, rate2, k1, k0, n, status, alpha, power, twosided, raw_n)
      real(dp), intent(in) :: rate1 !! First rater's positive probability, in [0,1].
      real(dp), intent(in) :: rate2 !! Second rater's positive probability, in [0,1].
      real(dp), intent(in) :: k1 !! Alternative kappa, in [-1,1] and feasible for the two marginals.
      real(dp), intent(in) :: k0 !! Null kappa, distinct from k1 and feasible for the two marginals.
      real(dp), intent(out) :: n !! Required sample size rounded upward, represented as real to avoid integer overflow.
      integer, intent(out) :: status !! irr_ok on success; invalid, degenerate, or resource otherwise.
      real(dp), intent(in), optional :: alpha !! Type-I error probability strictly between zero and one; default 0.05.
      real(dp), intent(in), optional :: power !! Desired power strictly between zero and one; default 0.8.
      logical, intent(in), optional :: twosided !! Whether to use a two-sided test; default false.
      real(dp), intent(out), optional :: raw_n !! Unrounded sample size, or NaN if calculation failed.
      real(dp) :: a, b, side, q0, q1, raw
      integer :: s0, s1

      n = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(raw_n)) raw_n = n
      call design_options(k1, k0, alpha, power, twosided, a, b, side, status)
      if (status /= irr_ok) return
      if (.not. all(ieee_is_finite([rate1, rate2]))) then
         status = irr_invalid
         return
      end if
      if (min(rate1, rate2) < 0.0_dp .or. max(rate1, rate2) > 1.0_dp) then
         status = irr_invalid
         return
      end if
      call binary_variance(rate1, rate2, k0, q0, s0)
      call binary_variance(rate1, rate2, k1, q1, s1)
      status = max(s0, s1)
      if (status /= irr_ok) return
      call finish_sample_size(q0, q1, k1 - k0, a, b, side, n, raw, status)
      if (present(raw_n)) raw_n = raw
   end subroutine n_cohen_kappa

   subroutine n2_cohen_kappa(mrg, k1, k0, n, status, alpha, power, twosided, raw_n)
      real(dp), intent(in) :: mrg(:) !! Shared marginal probabilities for 2 to 10 categories; nonnegative and sum to one.
      real(dp), intent(in) :: k1 !! Alternative kappa in [-1,1], feasible for the specified marginals.
      real(dp), intent(in) :: k0 !! Null kappa in [-1,1], distinct from k1 and feasible for the marginals.
      real(dp), intent(out) :: n !! Upward-rounded sample size as real(dp), or NaN on failure.
      integer, intent(out) :: status !! irr_ok, irr_invalid, irr_degenerate, or irr_resource.
      real(dp), intent(in), optional :: alpha !! Type-I error probability in (0,1); default 0.05.
      real(dp), intent(in), optional :: power !! Desired power in (0,1); default 0.8.
      logical, intent(in), optional :: twosided !! Use two tails when true; default false.
      real(dp), intent(out), optional :: raw_n !! Unrounded required sample size, or NaN on failure.
      real(dp) :: a, b, side, q0, q1, raw
      integer :: s0, s1

      n = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(raw_n)) raw_n = n
      call design_options(k1, k0, alpha, power, twosided, a, b, side, status)
      if (status /= irr_ok) return
      status = irr_invalid
      if (size(mrg) < 2) return
      if (size(mrg) > 10) then
         status = irr_resource
         return
      end if
      if (.not. all(ieee_is_finite(mrg))) return
      if (any(mrg < 0.0_dp) .or. any(mrg > 1.0_dp)) return
      if (abs(sum(mrg) - 1.0_dp) > 64.0_dp * epsilon(1.0_dp)) return
      call maximum_kappa_variance(mrg, k0, q0, s0)
      call maximum_kappa_variance(mrg, k1, q1, s1)
      status = max(s0, s1)
      if (status /= irr_ok) return
      call finish_sample_size(q0, q1, k1 - k0, a, b, side, n, raw, status)
      if (present(raw_n)) raw_n = raw
   end subroutine n2_cohen_kappa

   pure subroutine design_options(k1, k0, alpha, power, twosided, a, b, side, status)
      real(dp), intent(in) :: k1 !! Alternative kappa, required finite in [-1,1].
      real(dp), intent(in) :: k0 !! Null kappa, required finite in [-1,1] and unequal to k1.
      real(dp), intent(in), optional :: alpha !! Optional type-I error probability, default 0.05.
      real(dp), intent(in), optional :: power !! Optional desired power, default 0.8.
      logical, intent(in), optional :: twosided !! Optional two-tailed test flag, default false.
      real(dp), intent(out) :: a !! Effective type-I error probability.
      real(dp), intent(out) :: b !! Effective power probability.
      real(dp), intent(out) :: side !! Number of test tails, one or two.
      integer, intent(out) :: status !! Validation status using irr status constants.

      a = 0.05_dp
      if (present(alpha)) a = alpha
      b = 0.8_dp
      if (present(power)) b = power
      side = 1.0_dp
      if (present(twosided)) then
         if (twosided) side = 2.0_dp
      end if
      status = irr_invalid
      if (.not. all(ieee_is_finite([k0, k1, a, b]))) return
      if (max(abs(k0), abs(k1)) > 1.0_dp .or. k0 == k1) return
      if (a <= 0.0_dp .or. a >= 1.0_dp .or. b <= 0.0_dp .or. b >= 1.0_dp) return
      status = irr_ok
   end subroutine design_options

   pure subroutine binary_variance(rate1, rate2, kappa, variance, status)
      real(dp), intent(in) :: rate1 !! First rater's positive probability in [0,1].
      real(dp), intent(in) :: rate2 !! Second rater's positive probability in [0,1].
      real(dp), intent(in) :: kappa !! Kappa for which the multinomial asymptotic variance is required.
      real(dp), intent(out) :: variance !! Sample-size-independent variance coefficient of kappa.
      integer, intent(out) :: status !! Feasibility and numerical status.
      real(dp) :: pe, po, cells(4), a, b, numerator

      variance = ieee_value(0.0_dp, ieee_quiet_nan)
      a = 1.0_dp - rate1
      b = 1.0_dp - rate2
      pe = rate1 * rate2 + a * b
      status = irr_degenerate
      if (pe >= 1.0_dp) return
      po = kappa * (1.0_dp - pe) + pe
      cells(4) = (po - rate1 + b) / 2.0_dp
      cells(1) = po - cells(4)
      cells(2) = rate1 - cells(1)
      cells(3) = rate2 - cells(1)
      status = irr_invalid
      if (any(cells < -64.0_dp * epsilon(1.0_dp))) return
      cells = max(cells, 0.0_dp)
      numerator = cells(1) * (1.0_dp - pe - (rate2 + rate1) * (1.0_dp - po)) ** 2 &
         + cells(4) * (1.0_dp - pe - (b + a) * (1.0_dp - po)) ** 2 &
         + (1.0_dp - po) ** 2 * (cells(2) * (rate2 + a) ** 2 + cells(3) * (b + rate1) ** 2) &
         - (po * pe - 2.0_dp * pe + po) ** 2
      status = irr_degenerate
      if (numerator < -128.0_dp * epsilon(1.0_dp)) return
      variance = max(0.0_dp, numerator) / (1.0_dp - pe) ** 4
      if (.not. ieee_is_finite(variance)) return
      status = irr_ok
   end subroutine binary_variance

   subroutine maximum_kappa_variance(mrg, kappa, variance, status)
      real(dp), intent(in) :: mrg(:) !! Shared category marginals, already validated, for at most ten categories.
      real(dp), intent(in) :: kappa !! Feasible null or alternative kappa.
      real(dp), intent(out) :: variance !! Largest multinomial kappa variance coefficient under these marginals.
      integer, intent(out) :: status !! irr_ok or a validation, degeneracy, or LP resource failure.
      real(dp), allocatable :: objective(:), constraints(:, :), rhs(:)
      integer, allocatable :: sense(:)
      real(dp) :: pe, po, numerator
      integer :: categories, i, j, k
      type(lp_result) :: optimum
      type(lp_control) :: control

      variance = ieee_value(0.0_dp, ieee_quiet_nan)
      categories = size(mrg)
      pe = sum(mrg ** 2)
      po = kappa + pe * (1.0_dp - kappa)
      status = irr_degenerate
      if (pe >= 1.0_dp) return
      status = irr_invalid
      if (po <= 0.0_dp .or. po > 1.0_dp + 64.0_dp * epsilon(1.0_dp)) return
      if (po < max(0.0_dp, 2.0_dp * maxval(mrg) - 1.0_dp) - 64.0_dp * epsilon(1.0_dp)) return
      allocate(objective(categories ** 2), constraints(2 * categories + 2, categories ** 2))
      allocate(rhs(2 * categories + 2), sense(2 * categories + 2))
      constraints = 0.0_dp
      do j = 1, categories
         do i = 1, categories
            k = i + (j - 1) * categories
            if (i == j) then
               objective(k) = -(1.0_dp - po) * 2.0_dp * mrg(i) &
                  * (2.0_dp * (1.0_dp - pe) - (1.0_dp - po) * 2.0_dp * mrg(i))
               constraints(2 * categories + 2, k) = 1.0_dp
            else
               objective(k) = (1.0_dp - po) ** 2 * (mrg(i) + mrg(j)) ** 2
            end if
            constraints(i, k) = 1.0_dp
            constraints(categories + j, k) = 1.0_dp
            constraints(2 * categories + 1, k) = 1.0_dp
         end do
      end do
      rhs(1:categories) = mrg
      rhs(categories + 1:2 * categories) = mrg
      rhs(2 * categories + 1:) = [1.0_dp, po]
      sense = LP_EQ
      control%feasibility_tol = 1.0e-11_dp
      control%optimality_tol = 1.0e-12_dp
      control%max_simplex_iter = 200000
      call solve_lp(LP_MAX, objective, constraints, sense, rhs, optimum, control)
      if (optimum%status /= LP_OPTIMAL) then
         status = irr_degenerate
         if (optimum%status == LP_INFEASIBLE) status = irr_invalid
         if (optimum%status == LP_TIMEOUT) status = irr_resource
         return
      end if
      numerator = po * (1.0_dp - pe) ** 2 + optimum%objective &
         - (po * pe - 2.0_dp * pe + po) ** 2
      status = irr_degenerate
      if (numerator <= 0.0_dp) return
      variance = numerator / (1.0_dp - pe) ** 4
      if (.not. ieee_is_finite(variance)) return
      status = irr_ok
   end subroutine maximum_kappa_variance

   pure subroutine finish_sample_size(q0, q1, delta, alpha, power, side, n, raw_n, status)
      real(dp), intent(in) :: q0 !! Null variance coefficient, finite and nonnegative.
      real(dp), intent(in) :: q1 !! Alternative variance coefficient, finite and nonnegative.
      real(dp), intent(in) :: delta !! Nonzero alternative-minus-null kappa difference.
      real(dp), intent(in) :: alpha !! Type-I error probability in (0,1).
      real(dp), intent(in) :: power !! Desired power in (0,1).
      real(dp), intent(in) :: side !! One or two test tails.
      real(dp), intent(out) :: n !! Upward-rounded sample size as a real number.
      real(dp), intent(out) :: raw_n !! Unrounded sample size, NaN if numerical range is exceeded.
      integer, intent(out) :: status !! irr_ok or irr_resource for a nonrepresentable sample size.
      real(dp) :: ratio

      n = ieee_value(0.0_dp, ieee_quiet_nan)
      raw_n = n
      status = irr_resource
      ratio = (r_qnorm(alpha / side, lower_tail = .false.) * sqrt(q0) + r_qnorm(power) * sqrt(q1)) / delta
      if (.not. ieee_is_finite(ratio)) return
      if (abs(ratio) > sqrt(huge(1.0_dp))) return
      raw_n = ratio ** 2
      n = aint(raw_n)
      if (n < raw_n) n = n + 1.0_dp
      status = irr_ok
   end subroutine finish_sample_size

   pure subroutine rel_inter_intra(x, nrater, result, rho0inter, rho0intra, conf_level)
      real(dp), intent(in) :: x(:, :) !! Finite scores: subjects by contiguous rater blocks of repeated measurements.
      integer, intent(in) :: nrater !! Number of raters, at least two; must divide size(x,2), with at least two repeats.
      type(rel_inter_intra_result), intent(out) :: result !! Reliability estimates, tests, lower bounds, SEMs, and status.
      real(dp), intent(in), optional :: rho0inter !! Null interrater reliability in [0,1); default 0.6.
      real(dp), intent(in), optional :: rho0intra !! Null intrarater reliability in [0,1); default 0.8.
      real(dp), intent(in), optional :: conf_level !! One-sided confidence level strictly between zero and one; default 0.95.
      real(dp), allocatable :: cell_means(:, :), subject_means(:), rater_means(:), part(:)
      real(dp) :: n, t, m, rho_inter, rho_intra, confidence, grand, mss, msr, mssr, mse
      real(dp) :: sr, rr, interaction, sf, total_random, total_fixed, intra_random, intra_fixed
      real(dp) :: rh, numerator, denominator, f1, f2, f3, df1, df2, nan
      integer :: subjects, repetitions, i, j, first, last

      result%status = irr_invalid
      subjects = size(x, 1)
      result%subjects = subjects
      result%raters = nrater
      if (subjects < 2 .or. nrater < 2) return
      if (nrater > max_raters .or. size(x, kind = 8) > int(max_cells, 8)) then
         result%status = irr_resource
         return
      end if
      if (mod(size(x, 2), nrater) /= 0) return
      repetitions = size(x, 2) / nrater
      result%repetitions = repetitions
      if (repetitions < 2) return
      if (.not. all(ieee_is_finite(x))) return
      rho_inter = 0.6_dp
      if (present(rho0inter)) rho_inter = rho0inter
      rho_intra = 0.8_dp
      if (present(rho0intra)) rho_intra = rho0intra
      confidence = 0.95_dp
      if (present(conf_level)) confidence = conf_level
      if (.not. all(ieee_is_finite([rho_inter, rho_intra, confidence]))) return
      if (min(rho_inter, rho_intra) < 0.0_dp .or. max(rho_inter, rho_intra) >= 1.0_dp) return
      if (confidence <= 0.0_dp .or. confidence >= 1.0_dp) return
      n = real(subjects, dp)
      t = real(nrater, dp)
      m = real(repetitions, dp)
      allocate(cell_means(subjects, nrater), subject_means(subjects), rater_means(nrater), part(nrater))
      allocate(result%rohat(4 + 2 * nrater), result%fs(2 + nrater), result%pvalue(2 + nrater))
      allocate(result%lowvalue(3 + nrater), result%sem(3 + nrater), result%mean_squares(4 + nrater))
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      result%rohat = nan
      result%fs = nan
      result%pvalue = nan
      result%lowvalue = nan
      result%sem = nan
      result%mean_squares = nan
      result%variance_random = nan
      result%variance_fixed = nan
      result%nu = nan
      do j = 1, nrater
         first = (j - 1) * repetitions + 1
         last = j * repetitions
         cell_means(:, j) = sum(x(:, first:last) / m, dim = 2)
         part(j) = 0.0_dp
         do i = 1, subjects
            part(j) = part(j) + sum((x(i, first:last) - cell_means(i, j)) ** 2)
         end do
         part(j) = part(j) / (n * (m - 1.0_dp))
      end do
      subject_means = sum(cell_means / t, dim = 2)
      rater_means = sum(cell_means / n, dim = 1)
      grand = sum(subject_means / n)
      mss = m * t * sum((subject_means - grand) ** 2) / (n - 1.0_dp)
      msr = m * n * sum((rater_means - grand) ** 2) / (t - 1.0_dp)
      mssr = 0.0_dp
      do j = 1, nrater
         mssr = mssr + sum((cell_means(:, j) - subject_means - rater_means(j) + grand) ** 2)
      end do
      mssr = m * mssr / ((n - 1.0_dp) * (t - 1.0_dp))
      mse = sum(part / t)
      result%mean_squares = [mss, msr, mssr, mse, part]
      result%status = irr_resource
      if (.not. all(ieee_is_finite(result%mean_squares))) return
      sr = (mss - mssr) / (m * t)
      rr = (msr - mssr) / (m * n)
      interaction = (mssr - mse) / m
      sf = (mss - mse) / (m * t)
      result%variance_random = [sr, rr, interaction, mse]
      result%variance_fixed = [sf, rr, interaction, mse]
      total_random = sr + rr + interaction + mse
      total_fixed = sf + (t - 1.0_dp) * interaction / t + mse
      intra_random = sr + rr + interaction
      intra_fixed = sf + (t - 1.0_dp) * interaction / t
      result%rohat(1:4) = [quotient(sr, total_random), quotient(intra_random, total_random), &
         quotient(sf - interaction / t, total_fixed), quotient(intra_fixed, total_fixed)]
      result%rohat(5:4 + nrater) = quotient(intra_random, intra_random + part)
      result%rohat(5 + nrater:) = quotient(intra_fixed, intra_fixed + part)
      result%fs(1) = quotient((1.0_dp - rho_inter) * mss, (1.0_dp + (t - 1.0_dp) * rho_inter) * mssr)
      result%fs(2) = quotient((1.0_dp - rho_intra) * mss, (1.0_dp + (m - 1.0_dp) * rho_intra) * mse * t)
      result%fs(3:) = quotient((1.0_dp - rho_intra) * mss, (1.0_dp + (m - 1.0_dp) * rho_intra) * part * t)
      df1 = n - 1.0_dp
      df2 = (n - 1.0_dp) * (t - 1.0_dp)
      result%pvalue(1) = r_pf(result%fs(1), df1, df2, lower_tail = .false.)
      result%pvalue(2:) = r_pf(result%fs(2:), df1, n * (m - 1.0_dp), lower_tail = .false.)
      rh = result%rohat(1)
      numerator = df2 * (t * rh * (msr - mssr) + n * (1.0_dp + (t - 1.0_dp) * rh) * mssr &
         + n * t * (m - 1.0_dp) * rh * mse) ** 2
      denominator = (n - 1.0_dp) * (t * rh) ** 2 * msr ** 2 &
         + (n * (1.0_dp + (t - 1.0_dp) * rh) - t * rh) ** 2 * mssr ** 2 &
         + df2 * n * t * (m - 1.0_dp) * rh ** 2 * mse ** 2
      result%nu(1) = quotient(numerator, denominator)
      rh = result%rohat(3)
      numerator = df2 * (n * (1.0_dp + (t - 1.0_dp) * rh) * mssr &
         + n * t * (m - 1.0_dp) * rh * mse) ** 2
      denominator = (n * (1.0_dp + (t - 1.0_dp) * rh)) ** 2 * mssr ** 2 &
         + df2 * n * t * (m - 1.0_dp) * rh ** 2 * mse ** 2
      result%nu(2) = quotient(numerator, denominator)
      f1 = r_qf(confidence, df1, result%nu(1))
      f2 = r_qf(confidence, df1, result%nu(2))
      f3 = r_qf(confidence, df1, n * (m - 1.0_dp))
      result%lowvalue(1) = quotient(n * (mss - f1 * mssr), n * mss &
         + f1 * (t * (msr - mssr) + n * (t - 1.0_dp) * mssr + n * t * (m - 1.0_dp) * mse))
      result%lowvalue(2) = quotient(n * (mss - f2 * mssr), n * mss &
         + f2 * (n * (t - 1.0_dp) * mssr + n * t * (m - 1.0_dp) * mse))
      result%lowvalue(3) = quotient(mss / t - f3 * mse, mss / t + f3 * (m - 1.0_dp) * mse)
      result%lowvalue(4:) = quotient(mss / t - f3 * part, mss / t + f3 * (m - 1.0_dp) * part)
      do i = 1, size(result%lowvalue)
         if (result%lowvalue(i) > 1.0_dp) result%lowvalue(i) = 1.0_dp
      end do
      result%sem(1) = sqrt(mse)
      result%sem(2:1 + nrater) = sqrt(part)
      result%sem(2 + nrater) = nonnegative_sqrt(rr + interaction + mse)
      result%sem(3 + nrater) = nonnegative_sqrt(interaction + mse)
      result%status = irr_degenerate
      if (.not. all(ieee_is_finite(result%rohat))) return
      if (.not. all(ieee_is_finite(result%fs))) return
      if (.not. all(ieee_is_finite(result%pvalue))) return
      if (.not. all(ieee_is_finite(result%lowvalue))) return
      if (.not. all(ieee_is_finite(result%sem))) return
      if (.not. all(ieee_is_finite(result%nu))) return
      result%status = irr_ok
   end subroutine rel_inter_intra

   pure elemental real(dp) function quotient(numerator, denominator) result(value)
      real(dp), intent(in) :: numerator !! Numerator of a reliability formula, possibly nonfinite after a degenerate model.
      real(dp), intent(in) :: denominator !! Denominator; zero yields a signed infinity or NaN for zero divided by zero.

      if (denominator == 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         if (numerator /= 0.0_dp .and. ieee_is_finite(numerator)) then
            value = sign(ieee_value(0.0_dp, ieee_positive_inf), numerator)
         end if
      else
         value = numerator / denominator
      end if
   end function quotient

   pure elemental real(dp) function nonnegative_sqrt(variance) result(value)
      real(dp), intent(in) :: variance !! Estimated variance; a negative estimate yields NaN as in R's sqrt.

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (variance >= 0.0_dp) value = sqrt(variance)
   end function nonnegative_sqrt
end module irr_designs
