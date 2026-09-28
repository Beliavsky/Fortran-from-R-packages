! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_extended
   use r_kinds, only : dp
   use r_linalg, only : least_squares_svd, solve_system, symmetric_eigen, thin_svd
   use ftsa_factor, only : ftsa_er_gr
   use ftsa_models, only : ftsa_forecast_ftsm, ftsa_ftsm
   use ftsa_statistics, only : ftsa_invalid_input, ftsa_ok, ftsa_unsupported
   use ftsa_types, only : ftsm_fit
   implicit none
   private

   public :: ftsa_dynamic_flr, ftsa_dynupdate, ftsa_forecast_fplsr, ftsa_fplsr
   public :: ftsa_ftsm_weight_select, ftsa_isfe
   public :: ftsa_maf_multivariate, ftsa_pcscore_bootstrap_data, ftsa_stop_time_detect
   public :: ftsa_simpls

contains

   subroutine ftsa_simpls(xtrain, ytrain, ncomp, coefficients, prediction, status, xtest, weight, beta)
      real(dp), intent(in) :: xtrain(:, :) !! Predictor training matrix with observations in rows and variables in columns.
      real(dp), intent(in) :: ytrain(:, :) !! Response training matrix with the same observation count as `xtrain`.
      integer, intent(in) :: ncomp !! Number of SIMPLS components; must be between one and `min(n-1,p)`.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Predictor-by-response regression coefficient matrix.
      real(dp), allocatable, intent(out) :: prediction(:, :) !! Optional-test predictions; zero rows when `xtest` is absent.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or a failed SVD.
      real(dp), intent(in), optional :: xtest(:, :) !! Test predictors with the same variable count as `xtrain`.
      logical, intent(in), optional :: weight !! Apply upstream exponentially decreasing observation weights when true.
      real(dp), intent(in), optional :: beta !! Exponential-weight parameter in `(0,1]`; default is `0.1`.
      real(dp), allocatable :: p_load(:, :), q_load(:, :), r_load(:, :), s_cross(:, :)
      real(dp), allocatable :: t_scores(:, :), u(:, :), singular(:), vt(:, :), v_load(:, :)
      real(dp), allocatable :: x(:, :), y(:, :), xcenter(:, :), ycenter(:, :)
      real(dp) :: beta_work, norm_t, norm_v
      real(dp) :: mean_x(size(xtrain, 2)), mean_y(size(ytrain, 2))
      real(dp) :: pp(size(xtrain, 2)), qq(size(ytrain, 2)), rr(size(xtrain, 2))
      real(dp) :: tt(size(xtrain, 1)), uu(size(xtrain, 1)), vv(size(xtrain, 2))
      integer :: a, info, m, n, p
      logical :: use_weight

      n = size(xtrain, 1)
      p = size(xtrain, 2)
      m = size(ytrain, 2)
      status = ftsa_ok
      allocate(coefficients(0, 0), prediction(0, 0))
      if (n < 2 .or. p < 1 .or. m < 1 .or. size(ytrain, 1) /= n) then
         status = ftsa_invalid_input
         return
      end if
      if (ncomp < 1 .or. ncomp > min(n - 1, p)) then
         status = ftsa_invalid_input
         return
      end if
      if (present(xtest)) then
         if (size(xtest, 2) /= p) then
            status = ftsa_invalid_input
            return
         end if
      end if

      use_weight = .false.
      if (present(weight)) use_weight = weight
      beta_work = 0.1_dp
      if (present(beta)) beta_work = beta
      if (use_weight .and. (beta_work <= 0.0_dp .or. beta_work > 1.0_dp)) then
         status = ftsa_invalid_input
         return
      end if

      mean_x = sum(xtrain, dim=1) / real(n, dp)
      mean_y = sum(ytrain, dim=1) / real(n, dp)
      allocate(xcenter(n, p), ycenter(n, m), x(n, p), y(n, m))
      xcenter = xtrain - spread(mean_x, 1, n)
      ycenter = ytrain - spread(mean_y, 1, n)
      x = xcenter
      y = ycenter
      if (use_weight) then
         do a = 1, n
            x(a, :) = x(a, :) * beta_work * (1.0_dp - beta_work)**real(n - a, dp)
            y(a, :) = y(a, :) * beta_work * (1.0_dp - beta_work)**real(n - a, dp)
         end do
      end if

      allocate(s_cross(p, m), r_load(p, ncomp), p_load(p, ncomp), q_load(m, ncomp))
      allocate(t_scores(n, ncomp), v_load(p, ncomp))
      r_load = 0.0_dp
      p_load = 0.0_dp
      q_load = 0.0_dp
      t_scores = 0.0_dp
      v_load = 0.0_dp
      s_cross = matmul(transpose(x), y)

      do a = 1, ncomp
         call thin_svd(s_cross, u, singular, vt, info)
         if (info /= 0 .or. size(vt, 1) < 1) then
            status = merge(info, ftsa_invalid_input, info /= 0)
            return
         end if
         qq = vt(1, :)
         rr = matmul(s_cross, qq)
         tt = matmul(x, rr)
         tt = tt - sum(tt) / real(n, dp)
         norm_t = sqrt(sum(tt * tt))
         if (norm_t <= tiny(1.0_dp)) then
            status = ftsa_invalid_input
            return
         end if
         tt = tt / norm_t
         rr = rr / norm_t
         pp = matmul(transpose(x), tt)
         qq = matmul(transpose(y), tt)
         uu = matmul(y, qq)
         vv = pp
         if (a > 1) then
            vv = vv - matmul(v_load(:, 1:a - 1), matmul(transpose(v_load(:, 1:a - 1)), pp))
            uu = uu - matmul(t_scores(:, 1:a - 1), matmul(transpose(t_scores(:, 1:a - 1)), uu))
         end if
         norm_v = sqrt(sum(vv * vv))
         if (norm_v <= tiny(1.0_dp)) then
            status = ftsa_invalid_input
            return
         end if
         vv = vv / norm_v
         s_cross = s_cross - spread(vv, 2, m) * spread(matmul(vv, s_cross), 1, p)
         r_load(:, a) = rr
         p_load(:, a) = pp
         q_load(:, a) = qq
         t_scores(:, a) = tt
         v_load(:, a) = vv
      end do

      deallocate(coefficients)
      allocate(coefficients(p, m))
      coefficients = matmul(r_load, transpose(q_load))
      if (present(xtest)) then
         deallocate(prediction)
         allocate(prediction(size(xtest, 1), m))
         prediction = matmul(xtest - spread(mean_x, 1, size(xtest, 1)), coefficients) + &
            spread(mean_y, 1, size(xtest, 1))
      end if
   end subroutine ftsa_simpls

   subroutine ftsa_fplsr(y, order, prediction, coefficients, fitted, residuals, status, weight, beta)
      real(dp), intent(in) :: y(:, :) !! Functional observations with grid points in rows and successive times in columns.
      integer, intent(in) :: order !! Number of SIMPLS components used for the lag-one functional regression.
      real(dp), allocatable, intent(out) :: prediction(:) !! One-step-ahead functional prediction on the input grid.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Grid-by-grid SIMPLS coefficient matrix.
      real(dp), allocatable, intent(out) :: fitted(:, :) !! Grid-by-training-time fitted response curves.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Grid-by-training-time response residual curves.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or failed SIMPLS fitting.
      logical, intent(in), optional :: weight !! Apply upstream exponentially decreasing SIMPLS weights when true.
      real(dp), intent(in), optional :: beta !! Exponential-weight parameter in `(0,1]`; default is `0.1`.
      real(dp), allocatable :: pred_matrix(:, :), xtrain(:, :), ytrain(:, :)
      real(dp) :: mean_x(size(y, 1)), mean_y(size(y, 1))
      integer :: n, p

      p = size(y, 1)
      n = size(y, 2)
      allocate(prediction(0), coefficients(0, 0), fitted(0, 0), residuals(0, 0))
      if (p < 1 .or. n < 3 .or. order < 1 .or. order > min(n - 2, p)) then
         status = ftsa_invalid_input
         return
      end if
      allocate(xtrain(n - 1, p), ytrain(n - 1, p))
      xtrain = transpose(y(:, 1:n - 1))
      ytrain = transpose(y(:, 2:n))
      if (present(weight) .and. present(beta)) then
         call ftsa_simpls(xtrain, ytrain, order, coefficients, pred_matrix, status, &
            xtest=transpose(y(:, n:n)), weight=weight, beta=beta)
      else if (present(weight)) then
         call ftsa_simpls(xtrain, ytrain, order, coefficients, pred_matrix, status, &
            xtest=transpose(y(:, n:n)), weight=weight)
      else if (present(beta)) then
         call ftsa_simpls(xtrain, ytrain, order, coefficients, pred_matrix, status, &
            xtest=transpose(y(:, n:n)), beta=beta)
      else
         call ftsa_simpls(xtrain, ytrain, order, coefficients, pred_matrix, status, xtest=transpose(y(:, n:n)))
      end if
      if (status /= ftsa_ok) return

      mean_x = sum(xtrain, dim=1) / real(n - 1, dp)
      mean_y = sum(ytrain, dim=1) / real(n - 1, dp)
      deallocate(prediction, fitted, residuals)
      allocate(prediction(p), fitted(p, n - 1), residuals(p, n - 1))
      prediction = pred_matrix(1, :)
      fitted = transpose(matmul(xtrain - spread(mean_x, 1, n - 1), coefficients) + &
         spread(mean_y, 1, n - 1))
      residuals = y(:, 2:n) - fitted
   end subroutine ftsa_fplsr

   subroutine ftsa_forecast_fplsr(y, components, horizon, forecasts, status, weight, beta)
      real(dp), intent(in) :: y(:, :) !! Initial functional observations with grid points in rows and times in columns.
      integer, intent(in) :: components !! SIMPLS component count refitted after each generated forecast.
      integer, intent(in) :: horizon !! Number of recursive one-step functional forecasts to generate.
      real(dp), allocatable, intent(out) :: forecasts(:, :) !! Grid-by-horizon recursively generated forecasts.
      integer, intent(out) :: status !! Zero when every recursive fit succeeds.
      logical, intent(in), optional :: weight !! Apply exponentially decreasing SIMPLS weights when true.
      real(dp), intent(in), optional :: beta !! Exponential-weight parameter in `(0,1]`; default is `0.1`.
      real(dp), allocatable :: coefficients(:, :), current(:, :), fitted(:, :), next(:, :), prediction(:), residuals(:, :)
      integer :: i

      allocate(forecasts(0, 0))
      if (size(y, 1) < 1 .or. size(y, 2) < 3 .or. horizon < 1) then
         status = ftsa_invalid_input
         return
      end if
      allocate(current(size(y, 1), size(y, 2)))
      current = y
      deallocate(forecasts)
      allocate(forecasts(size(y, 1), horizon))
      do i = 1, horizon
         if (present(weight) .and. present(beta)) then
            call ftsa_fplsr(current, components, prediction, coefficients, fitted, residuals, status, weight, beta)
         else if (present(weight)) then
            call ftsa_fplsr(current, components, prediction, coefficients, fitted, residuals, status, weight=weight)
         else if (present(beta)) then
            call ftsa_fplsr(current, components, prediction, coefficients, fitted, residuals, status, beta=beta)
         else
            call ftsa_fplsr(current, components, prediction, coefficients, fitted, residuals, status)
         end if
         if (status /= ftsa_ok) return
         forecasts(:, i) = prediction
         allocate(next(size(current, 1), size(current, 2) + 1))
         next(:, 1:size(current, 2)) = current
         next(:, size(current, 2) + 1) = prediction
         call move_alloc(next, current)
      end do
   end subroutine ftsa_forecast_fplsr

   subroutine ftsa_maf_multivariate(data, threshold, maf, loadings, reconstruction, residuals, &
      ncomp_threshold, ncomp_eigen_ratio, status)
      real(dp), intent(in) :: data(:, :) !! Variable-by-time matrix used by upstream `MAF_multivariate`.
      real(dp), intent(in) :: threshold !! Cumulative generalized-eigenvalue fraction in `(0,1]` for threshold selection.
      real(dp), allocatable, intent(out) :: maf(:, :) !! Time-by-component normalized maximum-autocorrelation factors.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Component-by-variable normalized reconstruction loadings.
      real(dp), allocatable, intent(out) :: reconstruction(:, :) !! Time-by-variable centered reconstruction.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Centered-data reconstruction residuals.
      integer, intent(out) :: ncomp_threshold !! Components required to reach the requested cumulative eigenvalue fraction.
      integer, intent(out) :: ncomp_eigen_ratio !! ER/GR-style component choice derived from generalized eigenvalues.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid input or failed linear algebra.
      real(dp), allocatable :: diff_x(:, :), eigenvalues(:), eigenvectors(:, :), inv_sqrt(:, :)
      real(dp), allocatable :: root(:, :), s_diff(:, :), s_z(:, :), temp(:, :), v(:, :), x(:, :), y(:, :), z(:, :)
      real(dp), allocatable :: basis(:, :), generalized(:, :), gen_values(:), gen_vectors(:, :)
      real(dp), allocatable :: er(:), gr(:)
      real(dp) :: cumulative, norm_const, total, u1, u2, v1, v2
      real(dp) :: means(size(data, 1))
      integer :: i, info, j, k, kmax, kstar, m, n, p, rank

      p = size(data, 1)
      n = size(data, 2)
      m = min(n, p)
      ncomp_threshold = 0
      ncomp_eigen_ratio = 0
      allocate(maf(0, 0), loadings(0, 0), reconstruction(0, 0), residuals(0, 0))
      if (p < 2 .or. n < 3 .or. threshold <= 0.0_dp .or. threshold > 1.0_dp) then
         status = ftsa_invalid_input
         return
      end if

      means = sum(data, dim=2) / real(n, dp)
      allocate(z(n, p))
      z = transpose(data - spread(means, 2, n))
      allocate(s_z(p, p))
      s_z = matmul(transpose(z), z) / real(n - 1, dp)
      call symmetric_eigen(s_z, eigenvalues, eigenvectors, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      allocate(root(p, p), inv_sqrt(p, p))
      root = 0.0_dp
      inv_sqrt = 0.0_dp
      do j = 1, p
         if (eigenvalues(j) > epsilon(1.0_dp) * max(1.0_dp, eigenvalues(1)) * real(p, dp)) then
            root = root + sqrt(eigenvalues(j)) * spread(eigenvectors(:, j), 2, p) * &
               spread(eigenvectors(:, j), 1, p)
            inv_sqrt = inv_sqrt + 1.0_dp / sqrt(eigenvalues(j)) * spread(eigenvectors(:, j), 2, p) * &
               spread(eigenvectors(:, j), 1, p)
         end if
      end do
      allocate(x(n, p), diff_x(n - 1, p), s_diff(p, p))
      x = matmul(z, root)
      diff_x = x(2:n, :) - x(1:n - 1, :)
      temp = diff_x - spread(sum(diff_x, dim=1) / real(n - 1, dp), 1, n - 1)
      s_diff = matmul(transpose(temp), temp) / real(max(1, n - 2), dp)
      call symmetric_eigen(s_diff, eigenvalues, v, info, descending=.false.)
      if (info /= 0) then
         status = info
         return
      end if
      allocate(y(n, p))
      y = matmul(z, matmul(root, v))
      allocate(basis(p, p))
      call least_squares_svd(y, z, basis, rank, info)
      if (info /= 0) then
         status = info
         return
      end if

      deallocate(maf, loadings, reconstruction, residuals)
      allocate(maf(n, p), loadings(p, p), reconstruction(n, p), residuals(n, p))
      loadings = basis
      maf = y
      do i = 1, p
         norm_const = sum(loadings(i, :))
         if (abs(norm_const) > sqrt(epsilon(1.0_dp))) then
            loadings(i, :) = loadings(i, :) / norm_const
            maf(:, i) = maf(:, i) * norm_const
         end if
      end do
      reconstruction = matmul(maf, loadings)
      residuals = z - reconstruction

      allocate(generalized(p, p))
      generalized = matmul(inv_sqrt, matmul(s_diff, inv_sqrt))
      generalized = 0.5_dp * (generalized + transpose(generalized))
      call symmetric_eigen(generalized, gen_values, gen_vectors, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      gen_values = max(gen_values, 0.0_dp)
      total = sum(gen_values(1:m))
      if (total <= tiny(1.0_dp)) then
         ncomp_threshold = 1
         ncomp_eigen_ratio = 1
         status = ftsa_ok
         return
      end if
      cumulative = 0.0_dp
      ncomp_threshold = m
      do k = 1, m
         cumulative = cumulative + gen_values(k)
         if (cumulative / total >= threshold) then
            ncomp_threshold = k
            exit
         end if
      end do

      kstar = count(gen_values(1:m) > total / real(m, dp))
      kmax = min(kstar, max(2, int(0.1_dp * real(m, dp))))
      if (kmax < 2) then
         ncomp_eigen_ratio = 1
      else
         allocate(er(kmax - 1), gr(kmax - 1))
         do k = 1, kmax - 1
            er(k) = gen_values(k) / max(tiny(1.0_dp), gen_values(k + 1))
            v1 = sum(gen_values(k + 1:m))
            v2 = sum(gen_values(k + 2:m))
            if (v1 <= tiny(1.0_dp) .or. v2 <= tiny(1.0_dp)) then
               gr(k) = 0.0_dp
            else
               u1 = gen_values(k) / v1
               u2 = gen_values(k + 1) / v2
               gr(k) = log(1.0_dp + u1) / max(tiny(1.0_dp), log(1.0_dp + u2))
            end if
         end do
         ncomp_eigen_ratio = max(maxloc(er, dim=1), maxloc(gr, dim=1))
      end if
      status = ftsa_ok
   end subroutine ftsa_maf_multivariate

   subroutine ftsa_pcscore_bootstrap_data(data, bootrep, statistic, bootdata, centers, status)
      real(dp), intent(in) :: data(:, :) !! Grid-by-time functional matrix whose principal-score rows are resampled.
      integer, intent(in) :: bootrep !! Number of bootstrap replicates; must be positive.
      character(len=*), intent(in) :: statistic !! Center statistic: `mean` or `variance`.
      real(dp), allocatable, intent(out) :: bootdata(:, :, :) !! Grid-by-time-by-replicate bootstrap curves.
      real(dp), allocatable, intent(out) :: centers(:, :) !! Grid-by-replicate bootstrap center functions.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions/statistic or failed SVD.
      real(dp), allocatable :: centered(:, :), reconstructed(:, :), singular(:), u(:, :), vt(:, :)
      real(dp), allocatable :: draw_u(:, :), diag_vt(:, :)
      real(dp) :: mean_curve(size(data, 1)), random_value
      character(len=:), allocatable :: key
      integer :: b, i, info, j, n, p, r, sample_index

      p = size(data, 1)
      n = size(data, 2)
      allocate(bootdata(0, 0, 0), centers(0, 0))
      key = lower_ascii(trim(statistic))
      if (p < 1 .or. n < 2 .or. bootrep < 1 .or. (key /= 'mean' .and. key /= 'variance')) then
         status = ftsa_invalid_input
         return
      end if
      mean_curve = sum(data, dim=2) / real(n, dp)
      allocate(centered(n, p))
      centered = transpose(data - spread(mean_curve, 2, n))
      call thin_svd(centered, u, singular, vt, info)
      if (info /= 0) then
         status = info
         return
      end if
      r = size(singular)
      allocate(draw_u(n, r), diag_vt(r, p), reconstructed(n, p))
      do i = 1, r
         diag_vt(i, :) = singular(i) * vt(i, :)
      end do
      deallocate(bootdata, centers)
      allocate(bootdata(p, n, bootrep), centers(p, bootrep))
      do b = 1, bootrep
         do i = 1, n
            call random_number(random_value)
            sample_index = min(n, 1 + int(random_value * real(n, dp)))
            draw_u(i, :) = u(sample_index, :)
         end do
         reconstructed = matmul(draw_u, diag_vt) + spread(mean_curve, 1, n)
         bootdata(:, :, b) = transpose(reconstructed)
         select case (key)
         case ('mean')
            centers(:, b) = sum(bootdata(:, :, b), dim=2) / real(n, dp)
         case ('variance')
            do j = 1, p
               centers(j, b) = sum((bootdata(j, :, b) - sum(bootdata(j, :, b)) / real(n, dp))**2) / real(n, dp)
            end do
         end select
      end do
      status = ftsa_ok
   end subroutine ftsa_pcscore_bootstrap_data

   subroutine ftsa_ftsm_weight_select(y, ncomp, ntestyear, criterion, beta_opt, status)
      real(dp), intent(in) :: y(:, :) !! Grid-by-time functional observations used for rolling validation.
      integer, intent(in) :: ncomp !! Principal-component count fitted at every validation origin.
      integer, intent(in) :: ntestyear !! Number of final observations reserved for sequential validation.
      character(len=*), intent(in) :: criterion !! Validation loss: `mae`, `mse`, or `mape`.
      real(dp), intent(out) :: beta_opt !! Approximately optimal exponential-weight parameter in `(0,1)`.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid settings or failed model fits.
      real(dp) :: a, b, c, d, fc, fd, golden
      character(len=:), allocatable :: key
      integer :: iter

      key = lower_ascii(trim(criterion))
      beta_opt = 0.0_dp
      if (size(y, 2) < 4 .or. ntestyear < 1 .or. ntestyear >= size(y, 2) - 1 .or. ncomp < 0 .or. &
         (key /= 'mae' .and. key /= 'mse' .and. key /= 'mape')) then
         status = ftsa_invalid_input
         return
      end if
      a = 1.0e-5_dp
      b = 1.0_dp - 1.0e-5_dp
      golden = (sqrt(5.0_dp) - 1.0_dp) / 2.0_dp
      c = b - golden * (b - a)
      d = a + golden * (b - a)
      call weighted_validation_loss(y, ncomp, ntestyear, key, c, fc, status)
      if (status /= ftsa_ok) return
      call weighted_validation_loss(y, ncomp, ntestyear, key, d, fd, status)
      if (status /= ftsa_ok) return
      do iter = 1, 48
         if (fc <= fd) then
            b = d
            d = c
            fd = fc
            c = b - golden * (b - a)
            call weighted_validation_loss(y, ncomp, ntestyear, key, c, fc, status)
         else
            a = c
            c = d
            fc = fd
            d = a + golden * (b - a)
            call weighted_validation_loss(y, ncomp, ntestyear, key, d, fd, status)
         end if
         if (status /= ftsa_ok) return
         if (b - a < 1.0e-5_dp) exit
      end do
      beta_opt = 0.5_dp * (a + b)
      status = ftsa_ok
   end subroutine ftsa_ftsm_weight_select

   subroutine ftsa_isfe(y, max_order, initial_n, horizons, values, status, drift)
      real(dp), intent(in) :: y(:, :) !! Grid-by-time functional observations for rolling forecast evaluation.
      integer, intent(in) :: max_order !! Largest principal-component order evaluated, including order zero.
      integer, intent(in) :: initial_n !! First rolling forecast origin; must leave room for `max(horizons)` observations.
      integer, intent(in) :: horizons(:) !! Positive forecast horizons to evaluate.
      real(dp), allocatable, intent(out) :: values(:, :) !! `(max_order+1)`-by-horizon average integrated squared errors.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid inputs or a failed FTSM fit/forecast.
      logical, intent(in), optional :: drift !! Use random-walk-with-drift component forecasts when true.
      type(ftsm_fit) :: fit
      real(dp), allocatable :: forecasts(:, :)
      real(dp) :: err
      integer :: hmax, i, j, n, n_origins, order, origin

      n = size(y, 2)
      allocate(values(0, 0))
      if (max_order < 0 .or. initial_n < 2 .or. size(horizons) < 1 .or. any(horizons < 1)) then
         status = ftsa_invalid_input
         return
      end if
      hmax = maxval(horizons)
      if (initial_n + hmax > n .or. max_order > min(size(y, 1), initial_n)) then
         status = ftsa_invalid_input
         return
      end if
      n_origins = n - hmax - initial_n + 1
      if (n_origins < 1) then
         status = ftsa_invalid_input
         return
      end if
      deallocate(values)
      allocate(values(max_order + 1, size(horizons)))
      values = 0.0_dp
      do order = 0, max_order
         do origin = initial_n, n - hmax
            call ftsa_ftsm(y(:, 1:origin), order, fit, status)
            if (status /= ftsa_ok) return
            if (present(drift)) then
               call ftsa_forecast_ftsm(fit, hmax, forecasts, status, drift)
            else
               call ftsa_forecast_ftsm(fit, hmax, forecasts, status)
            end if
            if (status /= ftsa_ok) return
            do j = 1, size(horizons)
               i = horizons(j)
               err = sum((forecasts(:, i) - y(:, origin + i))**2) / real(size(y, 1), dp)
               values(order + 1, j) = values(order + 1, j) + err
            end do
         end do
      end do
      values = values / real(n_origins, dp)
      status = ftsa_ok
   end subroutine ftsa_isfe

   subroutine ftsa_dynupdate(y, newdata, order, method, lambda, forecast, status, drift)
      real(dp), intent(in) :: y(:, :) !! Grid-by-time historical functional observations used to fit the FTSM basis.
      real(dp), intent(in) :: newdata(:) !! Observed prefix of the next curve, ordered from the first grid point.
      integer, intent(in) :: order !! Principal-component count in the historical FTSM.
      character(len=*), intent(in) :: method !! Update method: `ts`, `ols`, `pls`, or `ridge`.
      real(dp), intent(in) :: lambda !! Nonnegative ridge/PLS penalty; ignored for `ts` and `ols`.
      real(dp), allocatable, intent(out) :: forecast(:) !! Updated forecast for the unobserved suffix of the next curve.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid inputs or failed linear algebra.
      logical, intent(in), optional :: drift !! Use random-walk-with-drift score prediction when true.
      type(ftsm_fit) :: fit
      real(dp), allocatable :: base(:, :), base1(:, :), base2(:, :), coefficients(:), forecast_curve(:, :)
      real(dp), allocatable :: lhs(:, :), rhs(:)
      real(dp) :: future_coef(order + 1)
      character(len=:), allocatable :: key
      integer :: i, info, n2, p, rank

      p = size(y, 1)
      n2 = size(newdata)
      key = lower_ascii(trim(method))
      allocate(forecast(0))
      if (p < 2 .or. size(y, 2) < 2 .or. n2 < 0 .or. n2 >= p .or. order < 0 .or. lambda < 0.0_dp .or. &
         (key /= 'ts' .and. key /= 'ols' .and. key /= 'pls' .and. key /= 'ridge')) then
         status = ftsa_invalid_input
         return
      end if
      call ftsa_ftsm(y, order, fit, status)
      if (status /= ftsa_ok) return
      if (present(drift)) then
         call ftsa_forecast_ftsm(fit, 1, forecast_curve, status, drift)
      else
         call ftsa_forecast_ftsm(fit, 1, forecast_curve, status)
      end if
      if (status /= ftsa_ok) return
      allocate(base(p, order + 1))
      base(:, 1) = fit%mean_curve
      if (order > 0) base(:, 2:order + 1) = fit%basis(:, 1:order)
      future_coef(1) = 1.0_dp
      if (order > 0) then
         future_coef(2:order + 1) = matmul(transpose(fit%basis(:, 1:order)), forecast_curve(:, 1) - fit%mean_curve)
      end if
      if (key == 'ts' .or. n2 == 0) then
         deallocate(forecast)
         allocate(forecast(p - n2))
         forecast = forecast_curve(n2 + 1:p, 1)
         status = ftsa_ok
         return
      end if
      allocate(base1(n2, order + 1), base2(p - n2, order + 1), coefficients(order + 1))
      base1 = base(1:n2, :)
      base2 = base(n2 + 1:p, :)
      select case (key)
      case ('ols')
         call least_squares_svd(base1, newdata, coefficients, rank, info)
      case ('ridge', 'pls')
         allocate(lhs(order + 1, order + 1), rhs(order + 1))
         lhs = matmul(transpose(base1), base1)
         do i = 1, order + 1
            lhs(i, i) = lhs(i, i) + lambda
         end do
         rhs = matmul(transpose(base1), newdata)
         if (key == 'pls') rhs = rhs + lambda * future_coef
         call solve_system(lhs, rhs, coefficients, info)
      case default
         info = ftsa_unsupported
      end select
      if (info /= 0) then
         status = info
         return
      end if
      deallocate(forecast)
      allocate(forecast(p - n2))
      forecast = matmul(base2, coefficients)
      status = ftsa_ok
   end subroutine ftsa_dynupdate


   subroutine ftsa_dynamic_flr(dat, newdata, order_k_fraction, order_m_fraction, forecast, order_k, order_m, status)
      real(dp), intent(in) :: dat(:, :) !! Grid-by-time historical curves used to estimate prefix-to-suffix covariance regression.
      real(dp), intent(in) :: newdata(:) !! Observed prefix of the next curve; its length partitions the functional grid.
      real(dp), intent(in) :: order_k_fraction !! Cumulative prefix-eigenvalue fraction in `(0,1]` used to choose `order_k`.
      real(dp), intent(in) :: order_m_fraction !! Cumulative suffix-eigenvalue fraction in `(0,1]` used to choose `order_m`.
      real(dp), allocatable, intent(out) :: forecast(:) !! Point forecast for the unobserved suffix of the next curve.
      integer, intent(out) :: order_k !! Selected prefix principal-component dimension.
      integer, intent(out) :: order_m !! Selected suffix principal-component dimension.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid inputs or failed eigensolutions.
      real(dp), allocatable :: cross_cov(:, :), data_first(:, :), data_last(:, :), eig_k(:), eig_m(:)
      real(dp), allocatable :: phi(:, :), psi(:, :), vec_k(:, :), vec_m(:, :), final_beta(:, :)
      real(dp) :: centered_new(size(newdata)), lambda_km, prefix_mean(size(newdata))
      real(dp), allocatable :: suffix_mean(:)
      integer :: i, info, j, k, m, n, p

      p = size(dat, 1)
      n = size(dat, 2)
      k = size(newdata)
      m = p - k
      order_k = 0
      order_m = 0
      allocate(forecast(0))
      if (p < 2 .or. n < 2 .or. k < 1 .or. m < 1 .or. &
         order_k_fraction <= 0.0_dp .or. order_k_fraction > 1.0_dp .or. &
         order_m_fraction <= 0.0_dp .or. order_m_fraction > 1.0_dp) then
         status = ftsa_invalid_input
         return
      end if

      allocate(data_first(k, n), data_last(m, n), suffix_mean(m))
      prefix_mean = sum(dat(1:k, :), dim=2) / real(n, dp)
      suffix_mean = sum(dat(k + 1:p, :), dim=2) / real(n, dp)
      data_first = dat(1:k, :) - spread(prefix_mean, 2, n)
      data_last = dat(k + 1:p, :) - spread(suffix_mean, 2, n)
      centered_new = newdata - prefix_mean
      allocate(cross_cov(m, k))
      cross_cov = matmul(data_last, transpose(data_first))

      call symmetric_eigen(matmul(data_first, transpose(data_first)), eig_k, vec_k, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      call symmetric_eigen(matmul(data_last, transpose(data_last)), eig_m, vec_m, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      order_k = cumulative_dimension(eig_k, order_k_fraction)
      order_m = cumulative_dimension(eig_m, order_m_fraction)
      if (order_k < 1 .or. order_m < 1) then
         status = ftsa_invalid_input
         return
      end if
      allocate(phi(k, order_k), psi(m, order_m), final_beta(k, m))
      phi = vec_k(:, 1:order_k)
      psi = vec_m(:, 1:order_m)
      final_beta = 0.0_dp
      do j = 1, order_k
         if (eig_k(j) <= tiny(1.0_dp)) cycle
         do i = 1, order_m
            lambda_km = dot_product(psi(:, i), matmul(cross_cov, phi(:, j)))
            final_beta = final_beta + lambda_km / eig_k(j) * &
               spread(phi(:, j), 2, m) * spread(psi(:, i), 1, k)
         end do
      end do
      deallocate(forecast)
      allocate(forecast(m))
      forecast = suffix_mean + matmul(centered_new, final_beta)
      status = ftsa_ok
   end subroutine ftsa_dynamic_flr

   subroutine ftsa_stop_time_detect(y, forward_error, backward_error, forward_components, backward_components, &
      break_forward, break_backward, status, drift)
      real(dp), intent(in) :: y(:, :) !! Grid-by-time functional observations used for forward and reverse one-step monitoring.
      real(dp), allocatable, intent(out) :: forward_error(:) !! Integrated squared one-step errors for times `3:n`.
      real(dp), allocatable, intent(out) :: backward_error(:) !! Reverse-time integrated squared errors for times `1:n-2`.
      integer, allocatable, intent(out) :: forward_components(:) !! ER/GR-selected component counts at forward origins.
      integer, allocatable, intent(out) :: backward_components(:) !! ER/GR-selected component counts at reverse origins.
      integer, intent(out) :: break_forward !! One-based time index of the minimum-SSE forward mean-shift breakpoint.
      integer, intent(out) :: break_backward !! One-based time index of the minimum-SSE backward mean-shift breakpoint.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or unrecoverable model failures.
      logical, intent(in), optional :: drift !! Use random-walk-with-drift score forecasts when true.
      type(ftsm_fit) :: fit
      real(dp), allocatable :: eigenvalues(:), forecast(:, :), reversed(:, :), training(:, :)
      integer :: i, k_er, k_gr, local_status, n, order, p

      p = size(y, 1)
      n = size(y, 2)
      break_forward = 0
      break_backward = 0
      allocate(forward_error(0), backward_error(0), forward_components(0), backward_components(0))
      if (p < 2 .or. n < 5) then
         status = ftsa_invalid_input
         return
      end if
      deallocate(forward_error, backward_error, forward_components, backward_components)
      allocate(forward_error(n - 2), backward_error(n - 2), forward_components(n - 2), backward_components(n - 2))

      do i = 2, n - 1
         call ftsa_er_gr(transpose(y(:, 1:i)), k_er, k_gr, eigenvalues, local_status)
         if (local_status == ftsa_ok) then
            order = max(k_er, k_gr)
         else
            order = min(1, min(p, i))
         end if
         forward_components(i - 1) = order
         call ftsa_ftsm(y(:, 1:i), order, fit, local_status)
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
         if (present(drift)) then
            call ftsa_forecast_ftsm(fit, 1, forecast, local_status, drift)
         else
            call ftsa_forecast_ftsm(fit, 1, forecast, local_status)
         end if
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
         forward_error(i - 1) = sum((y(:, i + 1) - forecast(:, 1))**2)
      end do

      do i = 2, n - 1
         allocate(reversed(p, n - i + 1))
         reversed = y(:, n:i:-1)
         call ftsa_er_gr(transpose(reversed), k_er, k_gr, eigenvalues, local_status)
         if (local_status == ftsa_ok) then
            order = max(k_er, k_gr)
         else
            order = min(1, min(p, size(reversed, 2)))
         end if
         backward_components(i - 1) = order
         allocate(training(p, size(reversed, 2)))
         training = reversed
         call ftsa_ftsm(training, order, fit, local_status)
         deallocate(reversed, training)
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
         if (present(drift)) then
            call ftsa_forecast_ftsm(fit, 1, forecast, local_status, drift)
         else
            call ftsa_forecast_ftsm(fit, 1, forecast, local_status)
         end if
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
         backward_error(i - 1) = sum((y(:, i - 1) - forecast(:, 1))**2)
      end do

      break_forward = 2 + single_break_index(forward_error)
      break_backward = single_break_index(backward_error)
      status = ftsa_ok
   end subroutine ftsa_stop_time_detect

   subroutine weighted_validation_loss(y, ncomp, ntestyear, criterion, beta, loss, status)
      real(dp), intent(in) :: y(:, :) !! Grid-by-time observations for sequential holdout validation.
      integer, intent(in) :: ncomp !! Principal-component count.
      integer, intent(in) :: ntestyear !! Number of final holdout times.
      character(len=*), intent(in) :: criterion !! Scalar loss name: `mae`, `mse`, or `mape`.
      real(dp), intent(in) :: beta !! Exponential-weight parameter in `(0,1)`.
      real(dp), intent(out) :: loss !! Mean validation loss over all grid points and holdout times.
      integer, intent(out) :: status !! Zero on success; nonzero when a fit or forecast fails.
      type(ftsm_fit) :: fit
      real(dp), allocatable :: forecasts(:, :), weights(:)
      real(dp) :: error_value, denom
      integer :: origin, start_origin

      start_origin = size(y, 2) - ntestyear
      loss = 0.0_dp
      do origin = start_origin, size(y, 2) - 1
         allocate(weights(origin))
         call exponential_weights(origin, beta, weights)
         call ftsa_ftsm(y(:, 1:origin), ncomp, fit, status, weights)
         deallocate(weights)
         if (status /= ftsa_ok) return
         call ftsa_forecast_ftsm(fit, 1, forecasts, status)
         if (status /= ftsa_ok) return
         select case (criterion)
         case ('mae')
            error_value = sum(abs(forecasts(:, 1) - y(:, origin + 1))) / real(size(y, 1), dp)
         case ('mse')
            error_value = sum((forecasts(:, 1) - y(:, origin + 1))**2) / real(size(y, 1), dp)
         case ('mape')
            denom = sum(abs(y(:, origin + 1)))
            if (denom <= tiny(1.0_dp)) then
               error_value = huge(1.0_dp)
            else
               error_value = sum(abs((forecasts(:, 1) - y(:, origin + 1)) / &
                  max(abs(y(:, origin + 1)), tiny(1.0_dp)))) / real(size(y, 1), dp)
            end if
         case default
            status = ftsa_invalid_input
            return
         end select
         loss = loss + error_value
      end do
      loss = loss / real(ntestyear, dp)
      status = ftsa_ok
   end subroutine weighted_validation_loss

   pure subroutine exponential_weights(n, beta, weights)
      integer, intent(in) :: n !! Number of time observations receiving exponential weights.
      real(dp), intent(in) :: beta !! Exponential-weight parameter in `(0,1)`.
      real(dp), intent(out) :: weights(:) !! Chronological weights with largest weight on the most recent observation.
      integer :: i

      do i = 1, n
         weights(i) = beta * (1.0_dp - beta)**real(n - i, dp)
      end do
   end subroutine exponential_weights


   pure integer function cumulative_dimension(values, fraction) result(dimension)
      real(dp), intent(in) :: values(:) !! Descending nonnegative eigenvalues.
      real(dp), intent(in) :: fraction !! Target cumulative fraction in `(0,1]`.
      real(dp) :: cumulative, total
      integer :: i

      total = sum(max(values, 0.0_dp))
      if (total <= tiny(1.0_dp)) then
         dimension = 0
         return
      end if
      cumulative = 0.0_dp
      dimension = size(values)
      do i = 1, size(values)
         cumulative = cumulative + max(values(i), 0.0_dp)
         if (cumulative / total >= fraction) then
            dimension = i
            return
         end if
      end do
   end function cumulative_dimension

   pure integer function single_break_index(values) result(index)
      real(dp), intent(in) :: values(:) !! Scalar sequence in which one mean-shift breakpoint is selected.
      real(dp) :: left_mean, right_mean, sse, best_sse
      integer :: k, n

      n = size(values)
      if (n < 4) then
         index = max(1, n / 2)
         return
      end if
      index = 2
      best_sse = huge(1.0_dp)
      do k = 2, n - 2
         left_mean = sum(values(1:k)) / real(k, dp)
         right_mean = sum(values(k + 1:n)) / real(n - k, dp)
         sse = sum((values(1:k) - left_mean)**2) + sum((values(k + 1:n) - right_mean)**2)
         if (sse < best_sse) then
            best_sse = sse
            index = k
         end if
      end do
   end function single_break_index

   pure function lower_ascii(text) result(lower)
      character(len=*), intent(in) :: text !! ASCII string to convert to lower case.
      character(len=len(text)) :: lower
      integer :: code, i

      lower = text
      do i = 1, len(text)
         code = iachar(lower(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lower_ascii

end module ftsa_extended
