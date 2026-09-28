module psych_omega_sem
   use r_kinds, only : dp
   use r_linalg, only : spd_inverse_logdet, symmetric_eigen
   use psych_computational_extended, only : omega_from_cfa_loadings
   use psych_omega_extended, only : schmid_leiman_decomposition
   implicit none
   private
   public :: omega_sem_bifactor

contains

   subroutine omega_sem_bifactor(r, n_factors, cfa_loadings, uniqueness, omega_h, omega_total, score_r2, &
      omega_group, fit_discrepancy, iterations, info, n_obs, chi_square, max_iter, tolerance)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix fitted by the orthogonal bifactor CFA.
      integer, intent(in) :: n_factors !! Number of first-order/group factors requested; one gives a one-factor CFA.
      real(dp), allocatable, intent(out) :: cfa_loadings(:, :) !! Standardized CFA loadings; general factor is column one.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! Positive residual variances from the fitted CFA model.
      real(dp), intent(out) :: omega_h !! Omega hierarchical calculated from the fitted general-factor loadings.
      real(dp), intent(out) :: omega_total !! Omega total calculated from all fitted common-factor loadings.
      real(dp), allocatable, intent(out) :: score_r2(:) !! Factor-score determinacy-like diagonal used by omegaFromSem.
      real(dp), allocatable, intent(out) :: omega_group(:, :) !! Total/general/group omega table for the fitted bifactor model.
      real(dp), intent(out) :: fit_discrepancy !! ML discrepancy log|Sigma| + tr(R Sigma^-1) - log|R| - p.
      integer, intent(out) :: iterations !! Number of BFGS iterations used by the CFA fit.
      integer, intent(out) :: info !! Zero on convergence; nonzero for dimensions, initialization, SPD, or convergence failure.
      integer, intent(in), optional :: n_obs !! Sample size used only to scale the optional chi-square statistic.
      real(dp), intent(out), optional :: chi_square !! Approximate normal-theory chi-square `(n_obs-1)*fit_discrepancy`.
      integer, intent(in), optional :: max_iter !! Maximum BFGS iterations; default 500.
      real(dp), intent(in), optional :: tolerance !! Maximum absolute gradient convergence tolerance; default 1e-7.
      real(dp), allocatable :: sl(:, :), general(:), groups(:, :), phi(:, :), higher(:), h2(:), u2(:), p2(:)
      real(dp), allocatable :: complexity(:), start_loadings(:, :), params(:), keyed_r(:, :)
      integer, allocatable :: assignment(:)
      real(dp) :: logdet_r, raw_fit, tol
      integer :: p, k, n_groups, i, f, local_info, fit_info, maximum_iterations

      omega_h = 0.0_dp
      omega_total = 0.0_dp
      fit_discrepancy = huge(1.0_dp)
      iterations = 0
      info = 0
      if (present(chi_square)) chi_square = huge(1.0_dp)
      allocate(cfa_loadings(0, 0), uniqueness(0), score_r2(0), omega_group(0, 0))

      if (size(r, 1) /= size(r, 2) .or. size(r, 1) < 2 .or. n_factors < 1) then
         info = 1
         return
      end if
      p = size(r, 1)
      k = min(n_factors, p)
      maximum_iterations = 500
      if (present(max_iter)) maximum_iterations = max(1, max_iter)
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))

      call observed_logdet(r, logdet_r, local_info)
      if (local_info /= 0) then
         info = 2
         return
      end if

      if (k == 1) then
         n_groups = 0
         call one_factor_start(r, start_loadings, uniqueness, local_info)
         if (local_info /= 0) then
            info = 10 + local_info
            return
         end if
         allocate(assignment(p))
         assignment = 0
      else
         n_groups = k
         call schmid_leiman_decomposition(r, k, sl, general, groups, phi, higher, h2, u2, p2, complexity, local_info)
         if (local_info /= 0) then
            info = 20 + local_info
            return
         end if
         if (allocated(uniqueness)) deallocate(uniqueness)
         allocate(start_loadings(p, k + 1), assignment(p), uniqueness(p))
         start_loadings = 0.0_dp
         start_loadings(:, 1) = sl(:, 1)
         do i = 1, p
            assignment(i) = maxloc(abs(sl(i, 2:k + 1)), dim=1)
            start_loadings(i, assignment(i) + 1) = sl(i, assignment(i) + 1)
            uniqueness(i) = max(0.05_dp, r(i, i) - start_loadings(i, 1) ** 2 - &
               start_loadings(i, assignment(i) + 1) ** 2)
         end do
      end if

      call pack_parameters(start_loadings, uniqueness, assignment, n_groups, params)
      call fit_bifactor_ml(r, assignment, n_groups, params, raw_fit, iterations, local_info, maximum_iterations, tol)
      fit_info = local_info
      if (fit_info > 1) then
         info = 30 + fit_info
         return
      end if
      call unpack_parameters(params, assignment, n_groups, cfa_loadings, uniqueness)
      call align_bifactor_signs(r, cfa_loadings, keyed_r)
      call omega_from_cfa_loadings(keyed_r, cfa_loadings, omega_h, omega_total, score_r2, omega_group, local_info)
      if (local_info /= 0) then
         info = 40 + local_info
         return
      end if

      fit_discrepancy = max(0.0_dp, raw_fit - logdet_r - real(p, dp))
      if (present(chi_square)) then
         if (present(n_obs)) then
            if (n_obs > 1) then
               chi_square = real(n_obs - 1, dp) * fit_discrepancy
            else
               chi_square = fit_discrepancy
            end if
         else
            chi_square = 499.0_dp * fit_discrepancy
         end if
      end if
      if (fit_info == 0) then
         info = 0
      else
         info = 3
      end if

      do f = 2, size(cfa_loadings, 2)
         if (sum(cfa_loadings(:, f)) < 0.0_dp) cfa_loadings(:, f) = -cfa_loadings(:, f)
      end do
   end subroutine omega_sem_bifactor

   subroutine one_factor_start(r, loadings, uniqueness, info)
      real(dp), intent(in) :: r(:, :) !! Positive-definite correlation matrix for the one-factor initialization.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Item-by-one starting loading matrix.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! Positive starting residual variances.
      integer, intent(out) :: info !! Zero on success; nonzero when the leading eigenpair cannot be obtained.
      real(dp), allocatable :: values(:), vectors(:, :)
      integer :: i

      call symmetric_eigen(r, values, vectors, info, descending=.true.)
      if (info /= 0) then
         allocate(loadings(0, 0), uniqueness(0))
         return
      end if
      allocate(loadings(size(r, 1), 1), uniqueness(size(r, 1)))
      loadings(:, 1) = vectors(:, 1) * sqrt(max(values(1), 0.0_dp))
      if (sum(loadings(:, 1)) < 0.0_dp) loadings(:, 1) = -loadings(:, 1)
      do i = 1, size(r, 1)
         uniqueness(i) = max(0.05_dp, r(i, i) - loadings(i, 1) ** 2)
      end do
   end subroutine one_factor_start


   pure subroutine pack_parameters(loadings, uniqueness, assignment, n_groups, params)
      real(dp), intent(in) :: loadings(:, :) !! Starting CFA loadings with general factor in column one.
      real(dp), intent(in) :: uniqueness(:) !! Positive starting residual variances.
      integer, intent(in) :: assignment(:) !! Per-item group-factor index; zero denotes a one-factor model.
      integer, intent(in) :: n_groups !! Number of group-factor columns retained in the CFA model.
      real(dp), allocatable, intent(out) :: params(:) !! Packed free loadings followed by log residual variances.
      integer :: p, i, offset

      p = size(loadings, 1)
      if (n_groups == 0) then
         allocate(params(2 * p))
         params(1:p) = loadings(:, 1)
         params(p + 1:2 * p) = log(max(uniqueness, 1.0e-10_dp))
      else
         allocate(params(3 * p))
         params(1:p) = loadings(:, 1)
         do i = 1, p
            params(p + i) = loadings(i, assignment(i) + 1)
         end do
         offset = 2 * p
         params(offset + 1:offset + p) = log(max(uniqueness, 1.0e-10_dp))
      end if
   end subroutine pack_parameters

   pure subroutine unpack_parameters(params, assignment, n_groups, loadings, uniqueness)
      real(dp), intent(in) :: params(:) !! Packed fitted loadings and log residual variances.
      integer, intent(in) :: assignment(:) !! Per-item group-factor index; zero denotes a one-factor model.
      integer, intent(in) :: n_groups !! Number of group-factor columns retained in the CFA model.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Expanded sparse CFA loading matrix.
      real(dp), allocatable, intent(out) :: uniqueness(:) !! Positive fitted residual variances.
      integer :: p, i, offset

      p = size(assignment)
      if (n_groups == 0) then
         allocate(loadings(p, 1), uniqueness(p))
         loadings(:, 1) = params(1:p)
         uniqueness = exp(clamp_vector(params(p + 1:2 * p), -30.0_dp, 30.0_dp))
      else
         allocate(loadings(p, n_groups + 1), uniqueness(p))
         loadings = 0.0_dp
         loadings(:, 1) = params(1:p)
         do i = 1, p
            loadings(i, assignment(i) + 1) = params(p + i)
         end do
         offset = 2 * p
         uniqueness = exp(clamp_vector(params(offset + 1:offset + p), -30.0_dp, 30.0_dp))
      end if
   end subroutine unpack_parameters

   pure subroutine fit_bifactor_ml(r, assignment, n_groups, params, objective, iterations, info, max_iter, tolerance)
      real(dp), intent(in) :: r(:, :) !! Positive-definite target correlation matrix.
      integer, intent(in) :: assignment(:) !! Per-item group-factor index; zero denotes a one-factor model.
      integer, intent(in) :: n_groups !! Number of group factors in the fitted CFA model.
      real(dp), intent(inout) :: params(:) !! Starting parameters on input and fitted parameters on output.
      real(dp), intent(out) :: objective !! Final ML objective without constants from the observed matrix.
      integer, intent(out) :: iterations !! Number of BFGS iterations performed.
      integer, intent(out) :: info !! Zero on convergence; one for iteration limit; two for numerical line-search failure.
      integer, intent(in) :: max_iter !! Maximum number of BFGS iterations.
      real(dp), intent(in) :: tolerance !! Infinity-norm gradient convergence tolerance.
      real(dp), allocatable :: gradient(:), new_gradient(:), direction(:), trial(:), step_vector(:), y(:)
      real(dp), allocatable :: hessian_inverse(:, :), hy(:)
      real(dp) :: f, f_new, directional_derivative, step, ys, yhy, coefficient
      real(dp) :: relative_change
      integer :: m, i, iter, line_iter, eval_info

      m = size(params)
      allocate(gradient(m), new_gradient(m), direction(m), trial(m), step_vector(m), y(m), hy(m))
      allocate(hessian_inverse(m, m))
      hessian_inverse = 0.0_dp
      do i = 1, m
         hessian_inverse(i, i) = 1.0_dp
      end do
      call bifactor_objective_gradient(r, assignment, n_groups, params, f, gradient, eval_info)
      if (eval_info /= 0) then
         objective = huge(1.0_dp)
         iterations = 0
         info = 2
         return
      end if
      info = 1
      iterations = 0

      do iter = 1, max_iter
         iterations = iter
         if (maxval(abs(gradient)) <= tolerance) then
            info = 0
            exit
         end if
         direction = -matmul(hessian_inverse, gradient)
         directional_derivative = dot_product(gradient, direction)
         if (directional_derivative >= -epsilon(1.0_dp) * max(1.0_dp, dot_product(gradient, gradient))) then
            hessian_inverse = 0.0_dp
            do i = 1, m
               hessian_inverse(i, i) = 1.0_dp
            end do
            direction = -gradient
            directional_derivative = -dot_product(gradient, gradient)
         end if

         step = 1.0_dp
         eval_info = 1
         do line_iter = 1, 40
            trial = params + step * direction
            call bifactor_objective_gradient(r, assignment, n_groups, trial, f_new, new_gradient, eval_info)
            if (eval_info == 0) then
               if (f_new <= f + 1.0e-4_dp * step * directional_derivative) exit
            end if
            step = 0.5_dp * step
         end do
         if (eval_info /= 0 .or. line_iter > 40) then
            info = 2
            exit
         end if

         step_vector = trial - params
         y = new_gradient - gradient
         ys = dot_product(y, step_vector)
         relative_change = abs(f - f_new) / max(1.0_dp, abs(f))
         params = trial
         f = f_new
         gradient = new_gradient

         if (maxval(abs(gradient)) <= tolerance .or. &
            (relative_change <= tolerance * 0.1_dp .and. maxval(abs(step_vector)) <= sqrt(tolerance))) then
            info = 0
            exit
         end if

         if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, sqrt(dot_product(y, y) * dot_product(step_vector, step_vector)))) then
            hy = matmul(hessian_inverse, y)
            yhy = dot_product(y, hy)
            coefficient = (ys + yhy) / (ys * ys)
            hessian_inverse = hessian_inverse + coefficient * outer_product(step_vector, step_vector) - &
               (outer_product(hy, step_vector) + outer_product(step_vector, hy)) / ys
         else
            hessian_inverse = 0.0_dp
            do i = 1, m
               hessian_inverse(i, i) = 1.0_dp
            end do
         end if
      end do
      objective = f
   end subroutine fit_bifactor_ml

   pure subroutine bifactor_objective_gradient(r, assignment, n_groups, params, objective, gradient, info)
      real(dp), intent(in) :: r(:, :) !! Positive-definite target correlation matrix.
      integer, intent(in) :: assignment(:) !! Per-item group-factor index; zero denotes a one-factor model.
      integer, intent(in) :: n_groups !! Number of group factors represented by the sparse loading structure.
      real(dp), intent(in) :: params(:) !! Packed free loadings and log residual variances.
      real(dp), intent(out) :: objective !! ML objective `log|Sigma| + tr(R Sigma^-1)`.
      real(dp), intent(out) :: gradient(:) !! Analytic gradient with respect to packed parameters.
      integer, intent(out) :: info !! Zero on success; nonzero when the implied covariance is numerically invalid.
      real(dp), allocatable :: loadings(:, :), uniqueness(:), sigma(:, :), sigma_inverse(:, :), gmat(:, :), gl(:, :)
      real(dp) :: logdet, trace_term
      integer :: p, k, i, offset

      p = size(r, 1)
      call unpack_parameters(params, assignment, n_groups, loadings, uniqueness)
      k = size(loadings, 2)
      allocate(sigma(p, p))
      sigma = matmul(loadings, transpose(loadings))
      do i = 1, p
         sigma(i, i) = sigma(i, i) + uniqueness(i)
      end do
      call spd_inverse_logdet(sigma, sigma_inverse, logdet, info)
      if (info /= 0) then
         objective = huge(1.0_dp)
         gradient = 0.0_dp
         return
      end if
      trace_term = sum(transpose(r) * sigma_inverse)
      objective = logdet + trace_term
      if (objective >= huge(1.0_dp) * 0.1_dp) then
         info = 2
         gradient = 0.0_dp
         return
      end if

      allocate(gmat(p, p), gl(p, k))
      gmat = sigma_inverse - matmul(matmul(sigma_inverse, r), sigma_inverse)
      gl = 2.0_dp * matmul(gmat, loadings)
      gradient = 0.0_dp
      gradient(1:p) = gl(:, 1)
      if (k == 1) then
         offset = p
      else
         do i = 1, p
            gradient(p + i) = gl(i, assignment(i) + 1)
         end do
         offset = 2 * p
      end if
      do i = 1, p
         gradient(offset + i) = gmat(i, i) * uniqueness(i)
      end do
      info = 0
   end subroutine bifactor_objective_gradient

   pure subroutine align_bifactor_signs(r, loadings, keyed_r)
      real(dp), intent(in) :: r(:, :) !! Original correlation matrix before item-key sign alignment.
      real(dp), intent(inout) :: loadings(:, :) !! Fitted loadings, row-keyed and factor-sign aligned on output.
      real(dp), allocatable, intent(out) :: keyed_r(:, :) !! Correlation matrix after omegaFromSem-style item flips.
      integer :: i, f

      allocate(keyed_r(size(r, 1), size(r, 2)))
      keyed_r = r
      do i = 1, size(loadings, 1)
         if (loadings(i, 1) < 0.0_dp) then
            loadings(i, :) = -loadings(i, :)
            keyed_r(i, :) = -keyed_r(i, :)
            keyed_r(:, i) = -keyed_r(:, i)
         end if
      end do
      do f = 2, size(loadings, 2)
         if (sum(loadings(:, f)) < 0.0_dp) loadings(:, f) = -loadings(:, f)
      end do
   end subroutine align_bifactor_signs

   pure subroutine observed_logdet(a, logdet, info)
      real(dp), intent(in) :: a(:, :) !! Positive-definite observed correlation matrix.
      real(dp), intent(out) :: logdet !! Natural logarithm of the observed matrix determinant.
      integer, intent(out) :: info !! Zero when the shared SPD factorization succeeds.
      real(dp), allocatable :: inverse(:, :)

      call spd_inverse_logdet(a, inverse, logdet, info)
   end subroutine observed_logdet


   pure function outer_product(a, b) result(matrix)
      real(dp), intent(in) :: a(:) !! Left vector in the rank-one outer product.
      real(dp), intent(in) :: b(:) !! Right vector in the rank-one outer product.
      real(dp) :: matrix(size(a), size(b))
      integer :: i

      do i = 1, size(a)
         matrix(i, :) = a(i) * b
      end do
   end function outer_product

   pure function clamp_vector(x, lower, upper) result(y)
      real(dp), intent(in) :: x(:) !! Values to clamp elementwise.
      real(dp), intent(in) :: lower !! Inclusive lower bound applied to every element.
      real(dp), intent(in) :: upper !! Inclusive upper bound applied to every element.
      real(dp) :: y(size(x))

      y = min(upper, max(lower, x))
   end function clamp_vector

end module psych_omega_sem
