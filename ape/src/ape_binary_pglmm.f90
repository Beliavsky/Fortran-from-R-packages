! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Binary phylogenetic GLMM translated from ape R/binaryPGLMM.R.
! Formula/model-frame plumbing and the stochastic simulation wrapper are omitted.
module ape_binary_pglmm
   use r_kinds, only : dp
   use r_linalg, only : solve_spd, spd_inverse_logdet
   use ape_types, only : phylo_tree
   use ape_topology, only : phylogenetic_vcv
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   real(dp), parameter :: pi = 3.1415926535897932384626433832795029_dp
   real(dp), parameter :: two_pi = 2.0_dp * pi

   type, public :: binary_pglmm_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: coefficient_se(:)
      real(dp), allocatable :: coefficient_covariance(:, :)
      real(dp), allocatable :: z_score(:)
      real(dp), allocatable :: p_value(:)
      real(dp), allocatable :: fitted_probability(:)
      real(dp), allocatable :: random_effect(:)
      real(dp), allocatable :: working_residual(:)
      real(dp), allocatable :: phylogenetic_covariance(:, :)
      real(dp), allocatable :: working_covariance(:, :)
      real(dp) :: s2 = 0.0_dp
      real(dp) :: s2_p_value = 1.0_dp
      real(dp) :: conditional_reml_log_likelihood = -huge(1.0_dp)
      real(dp) :: conditional_reml_log_likelihood_null = -huge(1.0_dp)
      real(dp) :: convergence_s2 = huge(1.0_dp)
      real(dp) :: convergence_coefficients = huge(1.0_dp)
      integer :: iterations = 0
      integer :: covariance_resets = 0
      logical :: converged = .false.
   end type binary_pglmm_result

   public :: binary_pglmm_fit
   public :: binary_pglmm_reml_objective

contains

   subroutine binary_pglmm_fit(y, design, tree, result, info, s2_init, coefficients_init, tol_pql, &
      maxit_pql, maxit_reml, s2_upper)
      !! Fits ape's binary phylogenetic generalized linear mixed model by penalized quasi-likelihood.
      integer, intent(in) :: y(:) !! Binary response vector in numeric tip order; every value must be zero or one.
      real(dp), intent(in) :: design(:, :) !! Fixed-effect design matrix with species in rows and coefficients in columns.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogeny with branch lengths and one tip for every response row.
      type(binary_pglmm_result), intent(out) :: result !! Fixed/random effects, phylogenetic variance, tests, and convergence state.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid input, linear algebra, or optimizer failure.
      real(dp), intent(in), optional :: s2_init !! Initial phylogenetic variance; default `0.1` as in upstream ape.
      real(dp), intent(in), optional :: coefficients_init(:) !! Optional initial fixed-effect coefficients; otherwise logistic IRLS.
      real(dp), intent(in), optional :: tol_pql !! PQL coefficient/variance convergence tolerance; default `1e-6`.
      integer, intent(in), optional :: maxit_pql !! Maximum outer and inner PQL iterations; default 200.
      integer, intent(in), optional :: maxit_reml !! Maximum bounded-BFGS iterations for each variance update; default 100.
      real(dp), intent(in), optional :: s2_upper !! Numerical upper bound for the nonnegative variance search; default 1000.
      real(dp), allocatable :: b(:)
      real(dp), allocatable :: b_old(:)
      real(dp), allocatable :: c(:, :)
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: coefficients_old(:)
      real(dp), allocatable :: coefficients_inner_old(:)
      real(dp), allocatable :: denom(:, :)
      real(dp), allocatable :: denom_inverse(:, :)
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: inv_v(:, :)
      real(dp), allocatable :: inv_weight(:)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: normal_inverse(:, :)
      real(dp), allocatable :: normal_matrix(:, :)
      real(dp), allocatable :: num(:)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: vphy(:, :)
      real(dp), allocatable :: working_response(:)
      real(dp), allocatable :: xb(:)
      real(dp) :: logdet_design
      real(dp) :: logdet_normal
      real(dp) :: logdet_v
      real(dp) :: objective
      real(dp) :: objective_null
      real(dp) :: s2
      real(dp) :: s2_old
      real(dp) :: search_upper
      real(dp) :: tolerance
      real(dp) :: w
      integer :: i
      integer :: inner
      integer :: inner_limit
      integer :: n
      integer :: outer
      integer :: outer_limit
      integer :: p
      integer :: reml_limit
      integer :: status
      logical :: inner_converged

      result = binary_pglmm_result()
      info = 0
      n = size(y)
      p = size(design, 2)
      if (n < 2 .or. p < 1 .or. size(design, 1) /= n .or. tree%n_tip /= n .or. p >= n) then
         info = 1
         return
      end if
      if (any((y /= 0) .and. (y /= 1)) .or. all(y == 0) .or. all(y == 1)) then
         info = 2
         return
      end if
      if (any(.not. ieee_is_finite(design))) then
         info = 3
         return
      end if
      tolerance = 1.0e-6_dp
      if (present(tol_pql)) tolerance = tol_pql
      outer_limit = 200
      if (present(maxit_pql)) outer_limit = maxit_pql
      inner_limit = outer_limit
      reml_limit = 100
      if (present(maxit_reml)) reml_limit = maxit_reml
      search_upper = 1000.0_dp
      if (present(s2_upper)) search_upper = s2_upper
      if (tolerance <= 0.0_dp .or. outer_limit < 1 .or. reml_limit < 1 .or. search_upper <= 0.0_dp) then
         info = 4
         return
      end if

      call phylogenetic_vcv(tree, vphy, status)
      if (status /= 0 .or. maxval(vphy) <= 0.0_dp) then
         info = 10 + max(1, status)
         return
      end if
      ! Upstream evaluates, but does not assign, a determinant-normalization expression here.
      vphy = vphy / maxval(vphy)

      allocate(coefficients(p), coefficients_old(p), coefficients_inner_old(p))
      if (present(coefficients_init)) then
         if (size(coefficients_init) /= p .or. any(.not. ieee_is_finite(coefficients_init))) then
            info = 5
            return
         end if
         coefficients = coefficients_init
      else
         call logistic_initial_fit(y, design, coefficients, status)
         if (status /= 0) coefficients = 0.001_dp
      end if
      s2 = 0.1_dp
      if (present(s2_init)) s2 = abs(s2_init)
      s2 = min(search_upper, s2)

      allocate(b(n), b_old(n), mu(n), xb(n), inv_weight(n), working_response(n), h(n))
      b = 0.0_dp
      xb = matmul(design, coefficients)
      call logistic_vector(xb, mu)
      allocate(c(n, n))
      c = s2 * vphy
      coefficients_old = 1.0e6_dp
      s2_old = 1.0e6_dp
      result%converged = .false.

      do outer = 1, outer_limit
         result%iterations = outer
         s2_old = s2
         coefficients_old = coefficients
         coefficients_inner_old = 1.0e6_dp
         inner_converged = .false.

         do inner = 1, inner_limit
            if (sqrt(sum((coefficients - coefficients_inner_old)**2)) / sqrt(real(p, dp)) <= tolerance) then
               inner_converged = .true.
               exit
            end if
            coefficients_inner_old = coefficients
            do i = 1, n
               w = max(mu(i) * (1.0_dp - mu(i)), 1.0e-12_dp)
               inv_weight(i) = 1.0_dp / w
            end do
            allocate(v(n, n))
            v = c
            do i = 1, n
               v(i, i) = v(i, i) + inv_weight(i)
            end do
            call spd_inverse_logdet(v, inv_v, logdet_v, status)
            deallocate(v)
            if (status /= 0) then
               result%covariance_resets = result%covariance_resets + 1
               coefficients = 0.001_dp
               b = 0.0_dp
               xb = matmul(design, coefficients)
               call logistic_vector(xb, mu)
               if (result%covariance_resets >= 3) then
                  info = 20 + status
                  return
               end if
               cycle
            end if
            xb = matmul(design, coefficients)
            working_response = xb + b + (real(y, dp) - mu) * inv_weight
            allocate(denom(p, p), num(p))
            denom = matmul(transpose(design), matmul(inv_v, design))
            num = matmul(transpose(design), matmul(inv_v, working_response))
            call solve_spd(denom, num, coefficients, status)
            deallocate(denom, num)
            if (status /= 0) then
               info = 30 + status
               return
            end if
            xb = matmul(design, coefficients)
            b = matmul(c, matmul(inv_v, working_response - xb))
            call logistic_vector(xb + b, mu)
         end do
         if (.not. inner_converged .and. inner == inner_limit + 1) then
            ! Match upstream by continuing the outer iteration with the last mean-component estimates.
         end if

         h = working_response - matmul(design, coefficients)
         call minimize_reml_variance(inv_weight, h, vphy, design, search_upper, reml_limit, tolerance, &
            s2, objective, status)
         if (status /= 0) then
            info = 40 + status
            return
         end if
         c = s2 * vphy

         result%convergence_s2 = abs(s2 - s2_old)
         result%convergence_coefficients = sqrt(sum((coefficients - coefficients_old)**2)) / real(p, dp)
         if (result%convergence_s2 <= tolerance .and. &
            sqrt(sum((coefficients - coefficients_old)**2)) / sqrt(real(p, dp)) <= tolerance) then
            result%converged = .true.
            exit
         end if
      end do

      ! Recompute the final PQL mean and covariance quantities as upstream does after its outer loop.
      do i = 1, n
         w = max(mu(i) * (1.0_dp - mu(i)), 1.0e-12_dp)
         inv_weight(i) = 1.0_dp / w
      end do
      allocate(v(n, n))
      v = c
      do i = 1, n
         v(i, i) = v(i, i) + inv_weight(i)
      end do
      call spd_inverse_logdet(v, inv_v, logdet_v, status)
      if (status /= 0) then
         info = 50 + status
         return
      end if
      xb = matmul(design, coefficients)
      working_response = xb + b + (real(y, dp) - mu) * inv_weight
      allocate(denom(p, p), num(p))
      denom = matmul(transpose(design), matmul(inv_v, design))
      num = matmul(transpose(design), matmul(inv_v, working_response))
      call solve_spd(denom, num, coefficients, status)
      if (status /= 0) then
         info = 60 + status
         return
      end if
      xb = matmul(design, coefficients)
      b = matmul(c, matmul(inv_v, working_response - xb))
      call logistic_vector(xb + b, mu)
      h = working_response - xb

      call spd_inverse_logdet(denom, denom_inverse, logdet_normal, status)
      if (status /= 0) then
         info = 70 + status
         return
      end if
      allocate(result%coefficients(p), result%coefficient_se(p), result%coefficient_covariance(p, p))
      allocate(result%z_score(p), result%p_value(p))
      result%coefficients = coefficients
      result%coefficient_covariance = denom_inverse
      do i = 1, p
         result%coefficient_se(i) = sqrt(max(0.0_dp, denom_inverse(i, i)))
         if (result%coefficient_se(i) > 0.0_dp) then
            result%z_score(i) = coefficients(i) / result%coefficient_se(i)
            result%p_value(i) = erfc(abs(result%z_score(i)) / sqrt(2.0_dp))
         else
            result%z_score(i) = 0.0_dp
            result%p_value(i) = 1.0_dp
         end if
      end do

      ! Upstream reports the last variance-optimization objective even after recomputing the final mean component.
      objective_null = binary_pglmm_reml_objective(0.0_dp, inv_weight, h, vphy, design, status)
      if (status /= 0) then
         info = 90 + status
         return
      end if
      normal_matrix = matmul(transpose(design), design)
      call spd_inverse_logdet(normal_matrix, normal_inverse, logdet_design, status)
      if (status /= 0) then
         info = 100 + status
         return
      end if
      result%conditional_reml_log_likelihood = -0.5_dp * real(n - p, dp) * log(two_pi) &
         + 0.5_dp * logdet_design - 0.5_dp * objective
      result%conditional_reml_log_likelihood_null = -0.5_dp * real(n - p, dp) * log(two_pi) &
         + 0.5_dp * logdet_design - 0.5_dp * objective_null
      result%s2 = s2
      result%s2_p_value = 0.5_dp * erfc(sqrt(max(0.0_dp, result%conditional_reml_log_likelihood &
         - result%conditional_reml_log_likelihood_null)))
      result%fitted_probability = mu
      result%random_effect = b
      result%working_residual = h
      result%phylogenetic_covariance = vphy
      result%working_covariance = v
      if (.not. result%converged) info = 110
   end subroutine binary_pglmm_fit

   function binary_pglmm_reml_objective(s2, inv_weight, h, vphy, design, info) result(value)
      !! Evaluates the upstream `pglmm.reml` objective for a fixed nonnegative phylogenetic variance.
      real(dp), intent(in) :: s2 !! Nonnegative phylogenetic random-effect variance.
      real(dp), intent(in) :: inv_weight(:) !! Diagonal of inverse Bernoulli working weights.
      real(dp), intent(in) :: h(:) !! PQL working residual vector `Z - X*B`.
      real(dp), intent(in) :: vphy(:, :) !! Normalized phylogenetic covariance matrix used by the PQL iteration.
      real(dp), intent(in) :: design(:, :) !! Fixed-effect design matrix.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid dimensions/non-SPD covariance matrices.
      real(dp) :: value
      real(dp), allocatable :: inv_v(:, :)
      real(dp), allocatable :: normal_inverse(:, :)
      real(dp), allocatable :: normal_matrix(:, :)
      real(dp), allocatable :: v(:, :)
      real(dp) :: logdet_normal
      real(dp) :: logdet_v
      integer :: i
      integer :: n
      integer :: status

      value = huge(1.0_dp) / 100.0_dp
      info = 0
      n = size(h)
      if (n < 1 .or. size(inv_weight) /= n .or. size(vphy, 1) /= n .or. size(vphy, 2) /= n &
         .or. size(design, 1) /= n .or. size(design, 2) < 1 .or. s2 < 0.0_dp) then
         info = 1
         return
      end if
      if (any(inv_weight <= 0.0_dp) .or. any(.not. ieee_is_finite(inv_weight)) &
         .or. any(.not. ieee_is_finite(h)) .or. any(.not. ieee_is_finite(vphy)) &
         .or. any(.not. ieee_is_finite(design))) then
         info = 2
         return
      end if
      allocate(v(n, n))
      v = s2 * vphy
      do i = 1, n
         v(i, i) = v(i, i) + inv_weight(i)
      end do
      call spd_inverse_logdet(v, inv_v, logdet_v, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      normal_matrix = matmul(transpose(design), matmul(inv_v, design))
      call spd_inverse_logdet(normal_matrix, normal_inverse, logdet_normal, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      value = logdet_v + dot_product(h, matmul(inv_v, h)) + logdet_normal
   end function binary_pglmm_reml_objective

   subroutine minimize_reml_variance(inv_weight, h, vphy, design, upper, max_iter, tolerance, s2, objective, info)
      !! Minimizes the one-dimensional PQL REML variance objective by deterministic golden-section search.
      real(dp), intent(in) :: inv_weight(:) !! Diagonal of inverse Bernoulli working weights.
      real(dp), intent(in) :: h(:) !! Current PQL working residual vector.
      real(dp), intent(in) :: vphy(:, :) !! Normalized phylogenetic covariance matrix.
      real(dp), intent(in) :: design(:, :) !! Fixed-effect design matrix.
      real(dp), intent(in) :: upper !! Inclusive numerical upper bound for the phylogenetic variance.
      integer, intent(in) :: max_iter !! Maximum golden-section iterations.
      real(dp), intent(in) :: tolerance !! Relative interval-width convergence tolerance.
      real(dp), intent(out) :: s2 !! Estimated nonnegative phylogenetic variance.
      real(dp), intent(out) :: objective !! REML objective at the returned variance.
      integer, intent(out) :: info !! Zero on success or nonzero if all evaluated objectives are invalid.
      real(dp), parameter :: golden = 0.6180339887498948482045868343656381_dp
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: fc
      real(dp) :: fd
      real(dp) :: f0
      real(dp) :: fu
      integer :: iteration
      integer :: status

      info = 0
      a = 0.0_dp
      b = upper
      f0 = binary_pglmm_reml_objective(a, inv_weight, h, vphy, design, status)
      if (status /= 0) f0 = huge(1.0_dp) / 100.0_dp
      fu = binary_pglmm_reml_objective(b, inv_weight, h, vphy, design, status)
      if (status /= 0) fu = huge(1.0_dp) / 100.0_dp
      c = b - golden * (b - a)
      d = a + golden * (b - a)
      fc = binary_pglmm_reml_objective(c, inv_weight, h, vphy, design, status)
      if (status /= 0) fc = huge(1.0_dp) / 100.0_dp
      fd = binary_pglmm_reml_objective(d, inv_weight, h, vphy, design, status)
      if (status /= 0) fd = huge(1.0_dp) / 100.0_dp
      do iteration = 1, max_iter
         if (b - a <= tolerance * (1.0_dp + 0.5_dp * (a + b))) exit
         if (fc <= fd) then
            b = d
            fu = fd
            d = c
            fd = fc
            c = b - golden * (b - a)
            fc = binary_pglmm_reml_objective(c, inv_weight, h, vphy, design, status)
            if (status /= 0) fc = huge(1.0_dp) / 100.0_dp
         else
            a = c
            f0 = fc
            c = d
            fc = fd
            d = a + golden * (b - a)
            fd = binary_pglmm_reml_objective(d, inv_weight, h, vphy, design, status)
            if (status /= 0) fd = huge(1.0_dp) / 100.0_dp
         end if
      end do
      if (f0 <= fc .and. f0 <= fd .and. f0 <= fu) then
         s2 = 0.0_dp
         objective = f0
      else if (fu <= fc .and. fu <= fd) then
         s2 = upper
         objective = fu
      else if (fc <= fd) then
         s2 = c
         objective = fc
      else
         s2 = d
         objective = fd
      end if
      if (.not. ieee_is_finite(objective) .or. objective >= huge(1.0_dp) / 1000.0_dp) info = 1
   end subroutine minimize_reml_variance

   subroutine logistic_initial_fit(y, design, coefficients, info)
      !! Computes deterministic no-phylogeny logistic-regression starting values by Fisher scoring/IRLS.
      integer, intent(in) :: y(:) !! Binary response vector.
      real(dp), intent(in) :: design(:, :) !! Fixed-effect design matrix.
      real(dp), intent(out) :: coefficients(:) !! Logistic-regression coefficient estimate used to initialize PQL.
      integer, intent(out) :: info !! Zero on convergence or nonzero if a weighted normal equation is singular.
      real(dp), allocatable :: eta(:)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: normal(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: trial(:)
      real(dp), allocatable :: weighted_design(:, :)
      real(dp), allocatable :: working(:)
      real(dp), allocatable :: weight(:)
      integer :: i
      integer :: iteration
      integer :: n
      integer :: p
      integer :: status

      n = size(y)
      p = size(design, 2)
      info = 0
      coefficients = 0.0_dp
      allocate(eta(n), mu(n), working(n), weight(n), weighted_design(n, p), normal(p, p), rhs(p), trial(p))
      do iteration = 1, 100
         eta = matmul(design, coefficients)
         call logistic_vector(eta, mu)
         do i = 1, n
            weight(i) = max(mu(i) * (1.0_dp - mu(i)), 1.0e-10_dp)
            working(i) = eta(i) + (real(y(i), dp) - mu(i)) / weight(i)
            weighted_design(i, :) = weight(i) * design(i, :)
         end do
         normal = matmul(transpose(design), weighted_design)
         rhs = matmul(transpose(weighted_design), working)
         call solve_spd(normal, rhs, trial, status)
         if (status /= 0) then
            info = 10 + status
            return
         end if
         if (maxval(abs(trial - coefficients)) <= 1.0e-10_dp * (1.0_dp + maxval(abs(coefficients)))) then
            coefficients = trial
            return
         end if
         coefficients = trial
      end do
   end subroutine logistic_initial_fit

   pure subroutine logistic_vector(eta, probability)
      !! Evaluates the logistic inverse link without overflowing for large positive or negative predictors.
      real(dp), intent(in) :: eta(:) !! Linear predictors.
      real(dp), intent(out) :: probability(:) !! Logistic probabilities in the open unit interval.
      real(dp) :: t
      integer :: i

      do i = 1, size(eta)
         if (eta(i) >= 0.0_dp) then
            probability(i) = 1.0_dp / (1.0_dp + exp(-min(eta(i), 700.0_dp)))
         else
            t = exp(max(eta(i), -700.0_dp))
            probability(i) = t / (1.0_dp + t)
         end if
         probability(i) = max(1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, probability(i)))
      end do
   end subroutine logistic_vector

end module ape_binary_pglmm
