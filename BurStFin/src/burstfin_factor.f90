! Public-domain Burns Statistics Ltd. algorithms, translated to modern Fortran.
module burstfin_factor
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_quantiles, only: r_quantile_type7
   use r_linalg, only: thin_svd, symmetric_eigen
   implicit none
   private
   public :: statistical_factor_model, factor_model_stat, factor_covariance

   type :: statistical_factor_model
      real(dp), allocatable :: loadings(:, :), uniquenesses(:), sdev(:)
      real(dp), allocatable :: cumulative_variance_fraction(:)
      logical, allocatable :: constant_assets(:)
      integer :: iterations = 0
      logical :: converged = .false.
   end type statistical_factor_model
contains
   subroutine factor_model_stat(x, model, status, weights, center, center_values, frac_var, iter_max, &
         nfac_miss, full_min, reg_min, sd_min, quan_sd, tol, zero_load, range_factors, &
         constant_returns_okay, specific_floor, floor_type)
      real(dp), intent(in) :: x(:, :) !! Returns with observations in rows; all nonfinite entries are treated as missing.
      type(statistical_factor_model), intent(out) :: model !! Factor loadings, uniquenesses, volatilities, and diagnostics.
      integer, intent(out) :: status !! 0 success, 1 invalid option/shape, 2 insufficient data, 3 degenerate data, 4 solver failure.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative row weights; scalar means equal weights, default linear 0.5 to 1.5.
      logical, intent(in), optional :: center !! Subtract upstream weighted centers when true; default true.
      real(dp), intent(in), optional :: center_values(:) !! Explicit finite asset centers, overriding the logical center option.
      real(dp), intent(in), optional :: frac_var !! SVD cumulative-variance target in [0,1]; default 0.5.
      integer, intent(in), optional :: iter_max !! Maximum principal-factor refinements; default one, zero disables refinement.
      integer, intent(in), optional :: nfac_miss(:) !! Factors used for each missing-count value; last entry extends, default one.
      integer, intent(in), optional :: full_min !! Minimum complete assets when missing values remain; default 20.
      integer, intent(in), optional :: reg_min !! Regress only assets with strictly more observed rows; default 40.
      integer, intent(in), optional :: sd_min !! Minimum observations to estimate incomplete-asset volatility; default 20.
      real(dp), intent(in), optional :: quan_sd !! Type-7 replacement-volatility quantile in [0,1]; default 0.9.
      real(dp), intent(in), optional :: tol !! Absolute uniqueness-change convergence tolerance; default 0.001.
      logical, intent(in), optional :: zero_load !! Use zero rather than mean loadings for sparse assets; default false.
      integer, intent(in), optional :: range_factors(2) !! Nonnegative factor-count lower/upper bounds; default [0,min(n,p)].
      logical, intent(in), optional :: constant_returns_okay !! Permit zero-volatility assets with zero covariances; default false.
      real(dp), intent(in), optional :: specific_floor !! Positive uniqueness floor or quantile probability; default 0.1.
      character(*), intent(in), optional :: floor_type !! Exact choice 'quantile' (default) or 'fraction' for the specific floor.
      real(dp), allocatable :: z(:, :), w(:), scales(:), centers(:), u(:, :), d(:), vt(:, :)
      real(dp), allocatable :: full_load(:, :), cm(:, :), reduced(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: unique_old(:), unique_new(:), scores(:, :), singular_squared(:)
      logical, allocatable :: keep(:), observed(:, :), complete(:), valid_sd(:)
      integer, allocatable :: rows(:), columns(:), missing(:)
      real(dp) :: fv, qs, tolerance, sf, maximum, denom, comm, uf, total
      integer :: n, p, nr, nc, i, j, k, m, nfac, limit, fmin, rmin, smin, nf, info, bounds(2)
      logical :: do_center, allow_constants, zeros
      character(:), allocatable :: ft
      status = 1
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1) return
      fv = 0.5_dp
      qs = 0.9_dp
      tolerance = 0.001_dp
      sf = 0.1_dp
      limit = 1
      fmin = 20
      rmin = 40
      smin = 20
      do_center = .true.
      allow_constants = .false.
      zeros = .false.
      ft = 'quantile'
      if (present(frac_var)) fv = frac_var
      if (present(quan_sd)) qs = quan_sd
      if (present(tol)) tolerance = tol
      if (present(specific_floor)) sf = specific_floor
      if (present(iter_max)) limit = iter_max
      if (present(full_min)) fmin = full_min
      if (present(reg_min)) rmin = reg_min
      if (present(sd_min)) smin = sd_min
      if (present(center)) do_center = center
      if (present(constant_returns_okay)) allow_constants = constant_returns_okay
      if (present(zero_load)) zeros = zero_load
      if (present(floor_type)) ft = floor_type
      if (.not. all(ieee_is_finite([fv, qs, tolerance, sf]))) return
      if (min(fv, qs, tolerance) < 0.0_dp .or. max(fv, qs) > 1.0_dp) return
      if (min(limit, fmin, rmin) < 0 .or. smin < 1) return
      if (ft /= 'quantile' .and. ft /= 'fraction') return
      if (ft == 'quantile' .and. sf > 1.0_dp) return
      if (present(nfac_miss)) then
         if (size(nfac_miss) == 0) return
         if (any(nfac_miss < 0)) return
      end if
      keep = any(ieee_is_finite(x), dim = 2)
      nr = count(keep)
      status = 2
      if (nr < 2) return
      rows = pack([(i, i = 1, n)], keep)
      z = x(rows, :)
      observed = ieee_is_finite(z)
      where (.not. observed) z = 0.0_dp
      missing = nr - count(observed, dim = 1)
      complete = missing == 0
      nc = count(complete)
      if (any(missing > 0) .and. nc < fmin) return
      if (nc == 0) return
      columns = pack([(j, j = 1, p)], complete)
      status = 1
      bounds = [0, min(nr, nc)]
      if (present(range_factors)) bounds = [minval(range_factors), maxval(range_factors)]
      if (bounds(1) < 0) return
      allocate(w(nr), scales(p), centers(p), model%sdev(p), model%constant_assets(p), valid_sd(p))
      if (present(weights)) then
         if (.not. all(ieee_is_finite(weights))) return
         if (any(weights < 0.0_dp)) return
         if (size(weights) == 1) then
            if (weights(1) <= 0.0_dp) return
            w = 1.0_dp
         else if (size(weights) == n) then
            w = weights(rows)
         else if (size(weights) == nr) then
            w = weights
         else
            return
         end if
      else
         do i = 1, nr
            w(i) = 0.5_dp + real(i - 1, dp) / real(nr - 1, dp)
         end do
      end if
      maximum = maxval(w)
      if (maximum <= 0.0_dp) return
      w = w / maximum
      w = w / sum(w)
      if (present(center_values)) then
         if (size(center_values) /= p) return
         if (.not. all(ieee_is_finite(center_values))) return
      end if
      status = 3
      model%sdev = ieee_value(0.0_dp, ieee_quiet_nan)
      model%constant_assets = .false.
      valid_sd = .false.
      do j = 1, p
         scales(j) = maxval(abs(z(:, j)), mask = observed(:, j))
         if (count(observed(:, j)) == 0) scales(j) = 1.0_dp
         if (present(center_values)) scales(j) = max(scales(j), abs(center_values(j)))
         if (scales(j) == 0.0_dp) scales(j) = 1.0_dp
         where (observed(:, j)) z(:, j) = z(:, j) / scales(j)
         centers(j) = 0.0_dp
         if (do_center) centers(j) = sum(w * z(:, j), mask = observed(:, j))
         if (present(center_values)) centers(j) = center_values(j) / scales(j)
         where (observed(:, j)) z(:, j) = z(:, j) - centers(j)
         if (missing(j) > 0 .and. count(observed(:, j)) < smin) cycle
         denom = 1.0_dp
         if (missing(j) > 0) denom = sum(w, mask = observed(:, j))
         if (denom <= 0.0_dp) cycle
         maximum = sqrt(sum(w * z(:, j) ** 2, mask = observed(:, j)) / denom)
         model%sdev(j) = maximum * scales(j)
         valid_sd(j) = .true.
         if (maximum <= 0.0_dp) then
            if (.not. allow_constants) return
            model%constant_assets(j) = .true.
            model%sdev(j) = 1.0e-16_dp
         end if
      end do
      if (.not. any(valid_sd)) return
      if (any(.not. valid_sd)) then
         uf = r_quantile_type7(pack(model%sdev, valid_sd), qs)
         where (.not. valid_sd) model%sdev = uf
      end if
      if (any(model%sdev <= 0.0_dp) .or. .not. all(ieee_is_finite(model%sdev))) return
      do j = 1, p
         ! Divide by the volatility in normalized units before applying the row weight.
         where (observed(:, j)) z(:, j) = z(:, j) / (model%sdev(j) / scales(j)) * sqrt(w)
      end do
      call thin_svd(z(:, columns), u, d, vt, info)
      status = 4
      if (info /= 0) return
      status = 3
      total = sum(d ** 2)
      if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) return
      allocate(model%cumulative_variance_fraction(size(d)))
      total = 0.0_dp
      do i = 1, size(d)
         total = total + d(i) ** 2
         model%cumulative_variance_fraction(i) = total
      end do
      model%cumulative_variance_fraction = model%cumulative_variance_fraction / total
      nfac = count(model%cumulative_variance_fraction < fv) + 1
      nfac = min(size(d), max(bounds(1), min(bounds(2), nfac)))
      allocate(full_load(nc, nfac))
      do k = 1, nfac
         full_load(:, k) = vt(k, :) * d(k)
      end do
      singular_squared = d(:nfac) ** 2
      if (limit > 0 .and. nfac > 0) then
         cm = matmul(transpose(z(:, columns)), z(:, columns))
         unique_old = min(1.0_dp, max(0.0_dp, 1.0_dp - sum(full_load ** 2, dim = 2)))
         do m = 1, limit
            reduced = cm
            do j = 1, nc
               reduced(j, j) = reduced(j, j) - unique_old(j)
            end do
            call symmetric_eigen(reduced, eigenvalues, eigenvectors, info, descending = .true.)
            status = 4
            if (info /= 0) return
            do k = 1, nfac
               full_load(:, k) = eigenvectors(:, k) * sqrt(max(0.0_dp, eigenvalues(k)))
            end do
            unique_new = min(1.0_dp, max(0.0_dp, 1.0_dp - sum(full_load ** 2, dim = 2)))
            model%iterations = m
            model%converged = all(abs(unique_new - unique_old) < tolerance)
            if (model%converged) exit
            unique_old = unique_new
         end do
      end if
      allocate(model%loadings(p, nfac))
      model%loadings = 0.0_dp
      if (.not. zeros) then
         do k = 1, nfac
            model%loadings(:, k) = sum(full_load(:, k)) / real(nc, dp)
         end do
      end if
      model%loadings(columns, :) = full_load
      if (any(missing > 0) .and. nfac > 0) then
         scores = matmul(z(:, columns), full_load)
         do j = 1, p
            if (missing(j) == 0 .or. nr - missing(j) <= rmin) cycle
            nf = 1
            if (present(nfac_miss)) nf = nfac_miss(min(missing(j), size(nfac_miss)))
            nf = min(nf, nfac)
            if (nf == 0) cycle
            denom = sum(z(:, j) ** 2, mask = observed(:, j))
            status = 3
            if (denom <= 0.0_dp .or. any(singular_squared(:nf) <= tiny(1.0_dp))) return
            do k = 1, nf
               ! Deliberately preserve lsfit(asset, scores), not the reverse regression.
               model%loadings(j, k) = sum(z(:, j) * scores(:, k), mask = observed(:, j)) / &
                  denom / singular_squared(k)
            end do
         end do
      end if
      allocate(model%uniquenesses(p))
      do j = 1, p
         comm = sum(model%loadings(j, :) ** 2)
         if (comm > 1.0_dp) model%loadings(j, :) = model%loadings(j, :) / sqrt(comm)
         model%uniquenesses(j) = 1.0_dp - min(1.0_dp, comm)
      end do
      if (sf > 0.0_dp) then
         uf = sf
         if (ft == 'quantile') uf = r_quantile_type7(model%uniquenesses, sf)
         model%uniquenesses = max(uf, model%uniquenesses)
      end if
      where (model%constant_assets) model%sdev = 0.0_dp
      status = 0
   end subroutine factor_model_stat

   pure subroutine factor_covariance(model, covariance, status, output)
      type(statistical_factor_model), intent(in) :: model !! Successfully fitted model from factor_model_stat.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Square covariance in original asset order.
      integer, intent(out) :: status !! Zero on success; one for invalid model or output choice; three for overflow.
      character(*), intent(in), optional :: output !! Exact 'full' (default), 'systematic', or 'specific' covariance.
      character(:), allocatable :: mode
      integer :: p, i, j
      status = 1
      mode = 'full'
      if (present(output)) mode = output
      if (mode /= 'full' .and. mode /= 'systematic' .and. mode /= 'specific') return
      if (.not. allocated(model%sdev) .or. .not. allocated(model%loadings)) return
      if (.not. allocated(model%uniquenesses)) return
      p = size(model%sdev)
      if (size(model%loadings, 1) /= p .or. size(model%uniquenesses) /= p) return
      allocate(covariance(p, p))
      covariance = 0.0_dp
      if (mode /= 'specific') then
         covariance = matmul(model%loadings, transpose(model%loadings))
         do j = 1, p
            do i = 1, p
               covariance(i, j) = (covariance(i, j) * model%sdev(i)) * model%sdev(j)
            end do
         end do
      end if
      if (mode /= 'systematic') then
         do i = 1, p
            covariance(i, i) = covariance(i, i) + (model%uniquenesses(i) * model%sdev(i)) * model%sdev(i)
         end do
      end if
      status = 3
      if (.not. all(ieee_is_finite(covariance))) return
      status = 0
   end subroutine factor_covariance
end module burstfin_factor
