module psych_simulation_algorithms
   use r_kinds, only : dp, r_pi
   use r_linalg, only : symmetric_eigen
   use psych_stats, only : correlation_matrix, scale_columns
   use psych_distributions, only : normal_cdf
   implicit none
   private
   public :: simulate_structural_model, simulate_hierarchical_data, simulate_minor_factors
   public :: simulate_spherical_items, simulate_multilevel_data, simulate_anova_data
   public :: simulate_poly_ideal, simulate_poly_ideal_logistic, simulate_poly_ideal_normal
contains

   subroutine simulate_structural_model(loadings, model, reliability, phi, n_cases, factor_means, observed, theta, &
      sample_correlation, items, low, high, categories, threshold, info)
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loading matrix defining the common-score structure.
      real(dp), allocatable, intent(out) :: model(:, :) !! Implied item correlation matrix with unit diagonal.
      real(dp), allocatable, intent(out) :: reliability(:) !! Item communalities from the supplied loading matrix.
      real(dp), intent(in), optional :: phi(:, :) !! Factor correlation matrix; identity is used when absent.
      integer, intent(in), optional :: n_cases !! Number of random observations to generate; default zero.
      real(dp), intent(in), optional :: factor_means(:) !! Factor means mapped to item means through `loadings`.
      real(dp), allocatable, intent(out), optional :: observed(:, :) !! Simulated observation-by-item matrix.
      real(dp), allocatable, intent(out), optional :: theta(:, :) !! Simulated correlated factor scores before mean shifts.
      real(dp), allocatable, intent(out), optional :: sample_correlation(:, :) !! Correlation matrix of simulated observations.
      logical, intent(in), optional :: items !! Rescale continuous observations to integer-like item categories when true.
      real(dp), intent(in), optional :: low !! Lower continuous clipping bound used before item rescaling; default -2.
      real(dp), intent(in), optional :: high !! Upper continuous clipping bound used before item rescaling; default 2.
      integer, intent(in), optional :: categories !! Number of item categories used when `items` is true; default 5.
      real(dp), intent(in), optional :: threshold(:) !! Optional item thresholds converting final observations to zero/one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or eigensolver failure.
      real(dp), allocatable :: common(:, :), local_phi(:, :), local_theta(:, :), error(:, :), means(:), work(:, :)
      real(dp), allocatable :: roots(:, :), values(:), vectors(:, :), uniqueness(:)
      real(dp) :: lo, hi, xmin, xmax
      integer :: i, j, nf, nobs, p, ncat, local_info
      logical :: make_items

      p = size(loadings, 1)
      nf = size(loadings, 2)
      nobs = 0
      if (present(n_cases)) nobs = max(0, n_cases)
      local_info = 0
      if (p < 1 .or. nf < 1) local_info = 1
      if (present(phi)) then
         if (size(phi, 1) /= nf .or. size(phi, 2) /= nf) local_info = 1
      end if
      if (present(factor_means)) then
         if (size(factor_means) /= nf) local_info = 1
      end if
      if (present(threshold)) then
         if (size(threshold) /= p) local_info = 1
      end if
      if (local_info /= 0) then
         allocate(model(0, 0), reliability(0))
         if (present(observed)) allocate(observed(0, 0))
         if (present(theta)) allocate(theta(0, 0))
         if (present(sample_correlation)) allocate(sample_correlation(0, 0))
         if (present(info)) info = local_info
         return
      end if

      allocate(local_phi(nf, nf))
      local_phi = 0.0_dp
      do i = 1, nf
         local_phi(i, i) = 1.0_dp
      end do
      if (present(phi)) local_phi = phi

      allocate(common(p, p), model(p, p), reliability(p), uniqueness(p), means(p))
      common = matmul(matmul(loadings, local_phi), transpose(loadings))
      model = common
      do i = 1, p
         model(i, i) = 1.0_dp
         reliability(i) = sum(loadings(i, :) ** 2)
         uniqueness(i) = sqrt(max(0.0_dp, 1.0_dp - common(i, i)))
      end do
      means = 0.0_dp
      if (present(factor_means)) means = matmul(loadings, factor_means)

      if (nobs > 0 .and. (present(observed) .or. present(theta) .or. present(sample_correlation))) then
         call symmetric_eigen(local_phi, values, vectors, local_info, descending=.true.)
         if (local_info /= 0) then
            if (present(observed)) allocate(observed(0, 0))
            if (present(theta)) allocate(theta(0, 0))
            if (present(sample_correlation)) allocate(sample_correlation(0, 0))
            if (present(info)) info = local_info
            return
         end if
         allocate(roots(nf, nf), local_theta(nobs, nf), error(nobs, p), work(nobs, p))
         roots = vectors
         do j = 1, nf
            roots(:, j) = roots(:, j) * sqrt(max(values(j), 0.0_dp))
         end do
         call fill_normal_matrix(local_theta)
         local_theta = matmul(local_theta, transpose(roots))
         call fill_normal_matrix(error)
         work = matmul(local_theta, transpose(loadings))
         do j = 1, p
            work(:, j) = work(:, j) + uniqueness(j) * error(:, j) + means(j)
         end do

         make_items = .false.
         if (present(items)) make_items = items
         lo = -2.0_dp
         hi = 2.0_dp
         ncat = 5
         if (present(low)) lo = low
         if (present(high)) hi = high
         if (present(categories)) ncat = max(2, categories)
         if (make_items) then
            work = max(lo, min(hi, work))
            xmin = minval(work)
            xmax = maxval(work)
            if (xmax > xmin) then
               work = 1.0_dp + real(ncat - 1, dp) * (work - xmin) / (xmax - xmin)
               work = real(nint(work), dp)
            else
               work = 1.0_dp
            end if
         end if
         if (present(threshold)) then
            do j = 1, p
               work(:, j) = merge(1.0_dp, 0.0_dp, work(:, j) > threshold(j))
            end do
         end if
         if (present(observed)) observed = work
         if (present(theta)) theta = local_theta
         if (present(sample_correlation)) then
            allocate(sample_correlation(p, p))
            call correlation_matrix(work, sample_correlation)
         end if
      else
         if (present(observed)) allocate(observed(0, p))
         if (present(theta)) allocate(theta(0, nf))
         if (present(sample_correlation)) allocate(sample_correlation(0, 0))
      end if
      if (present(info)) info = local_info
   end subroutine simulate_structural_model

   subroutine simulate_hierarchical_data(gload, fload, model, sl_loadings, n_cases, factor_means, observed, theta, &
      sample_correlation, categorical, low, high, threshold, info)
      real(dp), intent(in) :: gload(:) !! Higher-order general-factor loadings for the group factors.
      real(dp), intent(in) :: fload(:, :) !! Item-by-group-factor loading matrix.
      real(dp), allocatable, intent(out) :: model(:, :) !! Population item correlation matrix.
      real(dp), allocatable, intent(out) :: sl_loadings(:, :) !! Item-by-(general plus group) Schmid-Leiman loadings.
      integer, intent(in), optional :: n_cases !! Number of observations to simulate; default zero.
      real(dp), intent(in), optional :: factor_means(:) !! Means of orthogonal general/group scores; default zero.
      real(dp), allocatable, intent(out), optional :: observed(:, :) !! Simulated observation-by-item matrix.
      real(dp), allocatable, intent(out), optional :: theta(:, :) !! Orthogonal general and group factor scores.
      real(dp), allocatable, intent(out), optional :: sample_correlation(:, :) !! Sample item correlation matrix.
      logical, intent(in), optional :: categorical !! Round and clamp observations when true; default false.
      real(dp), intent(in), optional :: low !! Lower categorical bound; default -3.
      real(dp), intent(in), optional :: high !! Upper categorical bound; default 3.
      real(dp), intent(in), optional :: threshold(:) !! Optional per-item thresholds converting observations to zero/one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for dimension errors.
      real(dp), allocatable :: fcor(:, :), primeload(:, :), fstar(:, :), gstar(:), common(:, :), uniqueness(:)
      real(dp), allocatable :: factors(:, :), error(:, :), work(:, :), means(:)
      real(dp) :: lo, hi
      integer :: i, j, nf, nobs, p, local_info
      logical :: make_cat

      p = size(fload, 1)
      nf = size(fload, 2)
      nobs = 0
      if (present(n_cases)) nobs = max(0, n_cases)
      local_info = 0
      if (p < 1 .or. nf < 1 .or. size(gload) /= nf) local_info = 1
      if (present(factor_means)) then
         if (size(factor_means) /= nf + 1) local_info = 1
      end if
      if (present(threshold)) then
         if (size(threshold) /= p) local_info = 1
      end if
      if (local_info /= 0) then
         allocate(model(0, 0), sl_loadings(0, 0))
         if (present(observed)) allocate(observed(0, 0))
         if (present(theta)) allocate(theta(0, 0))
         if (present(sample_correlation)) allocate(sample_correlation(0, 0))
         if (present(info)) info = local_info
         return
      end if

      allocate(fcor(nf, nf), primeload(p, nf), fstar(p, nf), gstar(p), sl_loadings(p, nf + 1))
      fcor = spread(gload, 2, nf) * spread(gload, 1, nf)
      do i = 1, nf
         fcor(i, i) = 1.0_dp
      end do
      model = matmul(matmul(fload, fcor), transpose(fload))
      common = model
      do i = 1, p
         model(i, i) = 1.0_dp
      end do
      do j = 1, nf
         primeload(:, j) = fload(:, j) * gload(j)
      end do
      fstar = sqrt(max(0.0_dp, fload * fload - primeload * primeload))
      gstar = matmul(fload, gload)
      sl_loadings(:, 1) = gstar
      sl_loadings(:, 2:nf + 1) = fstar

      if (nobs > 0 .and. (present(observed) .or. present(theta) .or. present(sample_correlation))) then
         allocate(factors(nobs, nf + 1), error(nobs, p), work(nobs, p), uniqueness(p), means(nf + 1))
         call fill_normal_matrix(factors)
         call fill_normal_matrix(error)
         means = 0.0_dp
         if (present(factor_means)) means = factor_means
         do j = 1, nf + 1
            factors(:, j) = factors(:, j) + means(j)
         end do
         work = matmul(factors, transpose(sl_loadings))
         do j = 1, p
            uniqueness(j) = sqrt(max(0.0_dp, 1.0_dp - common(j, j)))
            work(:, j) = work(:, j) + uniqueness(j) * error(:, j)
         end do
         make_cat = .false.
         if (present(categorical)) make_cat = categorical
         lo = -3.0_dp
         hi = 3.0_dp
         if (present(low)) lo = low
         if (present(high)) hi = high
         if (make_cat) work = max(lo, min(hi, real(nint(work), dp)))
         if (present(threshold)) then
            do j = 1, p
               work(:, j) = merge(1.0_dp, 0.0_dp, work(:, j) > threshold(j))
            end do
         end if
         if (present(observed)) observed = work
         if (present(theta)) theta = factors
         if (present(sample_correlation)) then
            allocate(sample_correlation(p, p))
            call correlation_matrix(work, sample_correlation)
         end if
      else
         if (present(observed)) allocate(observed(0, p))
         if (present(theta)) allocate(theta(0, nf + 1))
         if (present(sample_correlation)) allocate(sample_correlation(0, 0))
      end if
      if (present(info)) info = local_info
   end subroutine simulate_hierarchical_data

   subroutine simulate_minor_factors(n_variables, n_factors, loadings, model, n_cases, observed, sample_correlation, &
      general_loadings, big_loadings, small_values, n_small, phi, bipolar, info)
      integer, intent(in) :: n_variables !! Number of observed variables when `big_loadings` is absent.
      integer, intent(in) :: n_factors !! Number of major group factors when `big_loadings` is absent.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Generated major/general/minor loading matrix.
      real(dp), allocatable, intent(out) :: model(:, :) !! Population item correlation matrix.
      integer, intent(in), optional :: n_cases !! Number of observations to simulate; default zero.
      real(dp), allocatable, intent(out), optional :: observed(:, :) !! Optional simulated observation-by-item matrix.
      real(dp), allocatable, intent(out), optional :: sample_correlation(:, :) !! Optional sample correlation matrix.
      real(dp), intent(in), optional :: general_loadings(:) !! General-factor loadings; values may be recycled randomly.
      real(dp), intent(in), optional :: big_loadings(:, :) !! Explicit item-by-major-factor loadings replacing defaults.
      real(dp), intent(in), optional :: small_values(:) !! Candidate minor-factor loadings; default `[-0.2, 0.2, 0, ...]`.
      integer, intent(in), optional :: n_small !! Number of minor factors; default floor(number of items/4).
      real(dp), intent(in), optional :: phi(:, :) !! Correlation matrix for all generated factors.
      logical, intent(in), optional :: bipolar !! Randomly flip default major-loading signs when true; default true.
      integer, intent(out), optional :: info !! Zero on success; nonzero for incompatible dimensions or h2 above one.
      real(dp), allocatable :: major(:, :), small(:, :), candidates(:), local_phi(:, :), reliability(:), temp_obs(:, :)
      real(dp), allocatable :: temp_cor(:, :), g(:)
      real(dp) :: u
      integer :: group_size, i, j, k, ns, p, nf, nobs, local_info, pick
      logical :: use_bipolar

      p = n_variables
      nf = n_factors
      nobs = 0
      if (present(n_cases)) nobs = max(0, n_cases)
      local_info = 0
      if (present(big_loadings)) then
         p = size(big_loadings, 1)
         nf = size(big_loadings, 2)
      end if
      if (p < 1 .or. nf < 1) local_info = 1
      if (.not. present(big_loadings)) then
         if (mod(p, nf) /= 0) local_info = 1
      end if
      if (local_info /= 0) then
         allocate(loadings(0, 0), model(0, 0))
         if (present(observed)) allocate(observed(0, 0))
         if (present(sample_correlation)) allocate(sample_correlation(0, 0))
         if (present(info)) info = local_info
         return
      end if

      allocate(major(p, nf))
      if (present(big_loadings)) then
         major = big_loadings
      else
         major = 0.0_dp
         group_size = p / nf
         use_bipolar = .true.
         if (present(bipolar)) use_bipolar = bipolar
         do j = 1, nf
            do i = (j - 1) * group_size + 1, j * group_size
               call random_number(u)
               if (u < 0.5_dp) then
                  major(i, j) = 0.8_dp
               else
                  major(i, j) = 0.6_dp
               end if
               if (use_bipolar) then
                  call random_number(u)
                  if (u < 0.5_dp) major(i, j) = -major(i, j)
               end if
            end do
         end do
      end if

      ns = p / 4
      if (present(n_small)) ns = max(0, n_small)
      if (present(general_loadings)) then
         allocate(g(p))
         do i = 1, p
            if (size(general_loadings) >= p) then
               g(i) = general_loadings(i)
            else
               call random_number(u)
               pick = min(size(general_loadings), int(u * real(size(general_loadings), dp)) + 1)
               g(i) = general_loadings(pick)
            end if
         end do
      else
         allocate(g(0))
      end if

      if (ns > 0) then
         if (present(small_values)) then
            allocate(candidates(size(small_values) + ns))
            candidates(1:size(small_values)) = small_values
            candidates(size(small_values) + 1:) = 0.0_dp
         else
            allocate(candidates(ns + 2))
            candidates(1:2) = [-0.2_dp, 0.2_dp]
            candidates(3:) = 0.0_dp
         end if
         allocate(small(p, ns))
         do j = 1, ns
            do i = 1, p
               call random_number(u)
               pick = min(size(candidates), int(u * real(size(candidates), dp)) + 1)
               small(i, j) = candidates(pick)
            end do
         end do
      else
         allocate(small(p, 0))
      end if

      allocate(loadings(p, nf + ns + merge(1, 0, size(g) > 0)))
      k = 0
      if (size(g) > 0) then
         loadings(:, 1) = g
         k = 1
      end if
      loadings(:, k + 1:k + nf) = major
      if (ns > 0) loadings(:, k + nf + 1:) = small
      if (maxval(sum(loadings * loadings, dim=2)) > 1.0_dp + 1.0e-12_dp) then
         allocate(model(0, 0))
         if (present(observed)) allocate(observed(0, 0))
         if (present(sample_correlation)) allocate(sample_correlation(0, 0))
         local_info = 2
         if (present(info)) info = local_info
         return
      end if

      if (present(phi)) then
         if (size(phi, 1) /= size(loadings, 2) .or. size(phi, 2) /= size(loadings, 2)) then
            allocate(model(0, 0))
            if (present(observed)) allocate(observed(0, 0))
            if (present(sample_correlation)) allocate(sample_correlation(0, 0))
            local_info = 3
            if (present(info)) info = local_info
            return
         end if
         local_phi = phi
      else
         allocate(local_phi(size(loadings, 2), size(loadings, 2)))
         local_phi = 0.0_dp
         do i = 1, size(local_phi, 1)
            local_phi(i, i) = 1.0_dp
         end do
      end if
      call simulate_structural_model(loadings, model, reliability, phi=local_phi, n_cases=nobs, observed=temp_obs, &
         sample_correlation=temp_cor, info=local_info)
      if (present(observed)) observed = temp_obs
      if (present(sample_correlation)) sample_correlation = temp_cor
      if (present(info)) info = local_info
   end subroutine simulate_minor_factors

   subroutine simulate_spherical_items(items, simple, nx, ny, n_subjects, x_loading, y_loading, z_loading, g_loading, &
      x_bias, y_bias, z_bias, categorical, low, high, threshold, loadings, latent, info)
      real(dp), allocatable, intent(out) :: items(:, :) !! Simulated subject-by-item response matrix.
      logical, intent(in), optional :: simple !! Use the four-direction simple structure instead of spherical geometry.
      integer, intent(in), optional :: nx !! Number of polar x-coordinate levels for the spherical design; default 7.
      integer, intent(in), optional :: ny !! Number of azimuth divisions; default 12.
      integer, intent(in), optional :: n_subjects !! Number of simulated subjects; default 500.
      real(dp), intent(in), optional :: x_loading !! Spherical x-coordinate loading scale; default 0.55.
      real(dp), intent(in), optional :: y_loading !! Spherical y-coordinate loading scale; default 0.55.
      real(dp), intent(in), optional :: z_loading !! Spherical z-coordinate loading scale; default 0.55.
      real(dp), intent(in), optional :: g_loading !! General-factor loading shared by items; default zero.
      real(dp), intent(in), optional :: x_bias !! Mean shift of the x latent coordinate; default zero.
      real(dp), intent(in), optional :: y_bias !! Mean shift of the y latent coordinate; default zero.
      real(dp), intent(in), optional :: z_bias !! Mean shift of the z latent coordinate; default zero.
      logical, intent(in), optional :: categorical !! Round and clamp responses when true; default false.
      real(dp), intent(in), optional :: low !! Lower categorical bound; default -3.
      real(dp), intent(in), optional :: high !! Upper categorical bound; default 3.
      real(dp), intent(in), optional :: threshold(:) !! Optional per-item thresholds returning dichotomous responses.
      real(dp), allocatable, intent(out), optional :: loadings(:, :) !! Item-by-(g,x,y,z) generating loading matrix.
      real(dp), allocatable, intent(out), optional :: latent(:, :) !! Subject-by-(g,x,y,z) latent score matrix.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid sizes or impossible loadings.
      real(dp), allocatable :: f(:, :), true_scores(:, :), error(:, :)
      real(dp) :: xl, yl, zl, gl, xb, yb, zb, lo, hi, error_weight, polar, azimuth
      integer :: ix, iy, j, nsub, nxv, nyv, nvar, local_info
      logical :: use_simple, make_cat

      use_simple = .false.
      if (present(simple)) use_simple = simple
      nxv = 7
      nyv = 12
      nsub = 500
      if (present(nx)) nxv = nx
      if (present(ny)) nyv = ny
      if (present(n_subjects)) nsub = n_subjects
      xl = 0.55_dp
      yl = 0.55_dp
      zl = 0.55_dp
      gl = 0.0_dp
      xb = 0.0_dp
      yb = 0.0_dp
      zb = 0.0_dp
      lo = -3.0_dp
      hi = 3.0_dp
      if (present(x_loading)) xl = x_loading
      if (present(y_loading)) yl = y_loading
      if (present(z_loading)) zl = z_loading
      if (present(g_loading)) gl = g_loading
      if (present(x_bias)) xb = x_bias
      if (present(y_bias)) yb = y_bias
      if (present(z_bias)) zb = z_bias
      if (present(low)) lo = low
      if (present(high)) hi = high
      local_info = 0
      if (nxv < 2 .or. nyv < 2 .or. nsub < 1) local_info = 1
      if (xl * xl + yl * yl + zl * zl + gl * gl > 1.0_dp + 1.0e-12_dp) local_info = 2
      if (local_info /= 0) then
         allocate(items(0, 0))
         if (present(loadings)) allocate(loadings(0, 0))
         if (present(latent)) allocate(latent(0, 0))
         if (present(info)) info = local_info
         return
      end if

      if (use_simple) then
         nvar = 4 * nyv
      else
         nvar = nxv * (nyv - 1)
      end if
      if (present(threshold)) then
         if (size(threshold) /= nvar) then
            allocate(items(0, 0))
            if (present(loadings)) allocate(loadings(0, 0))
            if (present(latent)) allocate(latent(0, 0))
            if (present(info)) info = 3
            return
         end if
      end if
      allocate(f(nvar, 4))
      f = 0.0_dp
      f(:, 1) = gl
      if (use_simple) then
         do j = 1, nvar
            select case (mod(j - 1, 4))
            case (0)
               f(j, 2) = xl
               f(j, 3) = yl
            case (1)
               f(j, 4) = zl
            case (2)
               f(j, 2) = -xl
               f(j, 3) = -yl
            case default
               f(j, 4) = -zl
            end select
         end do
      else
         j = 0
         do ix = 1, nxv
            polar = r_pi * real(ix - 1, dp) / real(nxv - 1, dp)
            do iy = 1, nyv - 1
               j = j + 1
               azimuth = 2.0_dp * r_pi * real(iy - 1, dp) / real(nyv - 1, dp)
               f(j, 2) = cos(polar) * xl
               f(j, 3) = cos(azimuth) * yl * (1.0_dp - f(j, 2) ** 2)
               f(j, 4) = sin(azimuth) * zl * (1.0_dp - f(j, 2) ** 2)
            end do
         end do
      end if

      allocate(true_scores(nsub, 4), error(nsub, nvar), items(nsub, nvar))
      call fill_normal_matrix(true_scores)
      call fill_normal_matrix(error)
      true_scores(:, 2) = true_scores(:, 2) + xb
      true_scores(:, 3) = true_scores(:, 3) + yb
      true_scores(:, 4) = true_scores(:, 4) + zb
      error_weight = sqrt(max(0.0_dp, 1.0_dp - xl * xl - yl * yl - zl * zl - gl * gl))
      items = matmul(true_scores, transpose(f)) + error_weight * error
      make_cat = .false.
      if (present(categorical)) make_cat = categorical
      if (make_cat) items = max(lo, min(hi, real(nint(items), dp)))
      if (present(threshold)) then
         do j = 1, nvar
            items(:, j) = merge(1.0_dp, 0.0_dp, items(:, j) > threshold(j))
         end do
      end if
      if (present(loadings)) loadings = f
      if (present(latent)) latent = true_scores
      if (present(info)) info = 0
   end subroutine simulate_spherical_items

   subroutine simulate_multilevel_data(r_within, r_between, eta, n_groups, n_cases, within_scores, between_scores, &
      group, observations, info)
      real(dp), intent(in) :: r_within(:, :) !! Target within-group correlation matrix.
      real(dp), intent(in) :: r_between(:, :) !! Target between-group correlation matrix.
      real(dp), intent(in) :: eta(:) !! Per-variable within-level weights; squared values weight within scores.
      integer, intent(in) :: n_groups !! Number of groups.
      integer, intent(in) :: n_cases !! Total case count, divisible by `n_groups`.
      real(dp), allocatable, intent(out) :: within_scores(:, :) !! Correlated within-level case scores.
      real(dp), allocatable, intent(out) :: between_scores(:, :) !! Group-level scores replicated to cases.
      integer, allocatable, intent(out) :: group(:) !! One-based group identifier for each case.
      real(dp), allocatable, intent(out) :: observations(:, :) !! Weighted within-plus-between simulated scores.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or eigensolver failure.
      real(dp), allocatable :: z(:, :), z_group(:, :), means(:), sds(:), roots(:, :), values(:), vectors(:, :)
      real(dp), allocatable :: roots_bg(:, :), values_bg(:), vectors_bg(:, :), group_scores(:, :)
      integer :: cases_per_group, g, i, j, p, local_info

      p = size(r_within, 1)
      local_info = 0
      if (p < 1 .or. size(r_within, 2) /= p .or. size(r_between, 1) /= p .or. size(r_between, 2) /= p) local_info = 1
      if (size(eta) /= p .or. n_groups < 1 .or. n_cases < 1 .or. mod(n_cases, n_groups) /= 0) local_info = 1
      if (local_info /= 0) then
         allocate(within_scores(0, 0), between_scores(0, 0), group(0), observations(0, 0))
         if (present(info)) info = local_info
         return
      end if

      call symmetric_eigen(r_within, values, vectors, local_info, descending=.true.)
      if (local_info /= 0) then
         allocate(within_scores(0, 0), between_scores(0, 0), group(0), observations(0, 0))
         if (present(info)) info = local_info
         return
      end if
      call symmetric_eigen(r_between, values_bg, vectors_bg, local_info, descending=.true.)
      if (local_info /= 0) then
         allocate(within_scores(0, 0), between_scores(0, 0), group(0), observations(0, 0))
         if (present(info)) info = local_info
         return
      end if
      allocate(roots(p, p), roots_bg(p, p))
      roots = vectors
      roots_bg = vectors_bg
      do j = 1, p
         roots(:, j) = roots(:, j) * sqrt(max(values(j), 0.0_dp))
         roots_bg(:, j) = roots_bg(:, j) * sqrt(max(values_bg(j), 0.0_dp))
      end do

      allocate(z(n_cases, p), means(p), sds(p), within_scores(n_cases, p))
      call fill_normal_matrix(z)
      call scale_columns(z, within_scores, means, sds)
      within_scores = matmul(within_scores, transpose(roots))

      allocate(z_group(n_groups, p), group_scores(n_groups, p))
      call fill_normal_matrix(z_group)
      call scale_columns(z_group, group_scores, means, sds)
      group_scores = matmul(group_scores, transpose(roots_bg))

      allocate(between_scores(n_cases, p), group(n_cases), observations(n_cases, p))
      cases_per_group = n_cases / n_groups
      do g = 1, n_groups
         do i = 1, cases_per_group
            j = (g - 1) * cases_per_group + i
            group(j) = g
            between_scores(j, :) = group_scores(g, :)
         end do
      end do
      do j = 1, p
         observations(:, j) = within_scores(:, j) * eta(j) ** 2 + &
            between_scores(:, j) * (1.0_dp - eta(j) ** 2)
      end do
      if (present(info)) info = 0
   end subroutine simulate_multilevel_data

   subroutine simulate_anova_data(n_per_cell, n1, n2, n3, observations, iv1, iv2, iv3, es1, es2, es3, es12, es13, &
      es23, es123, es11, es22, es33, standardize, center, info)
      integer, intent(in) :: n_per_cell !! Replicates per factorial cell when all factors have positive level counts.
      integer, intent(in) :: n1 !! Number of levels for factor one; zero uses a random continuous predictor.
      integer, intent(in) :: n2 !! Number of levels for factor two; zero uses a random continuous predictor.
      integer, intent(in) :: n3 !! Number of levels for factor three; zero uses a random continuous predictor.
      real(dp), allocatable, intent(out) :: observations(:) !! Simulated dependent variable.
      real(dp), allocatable, intent(out) :: iv1(:) !! First predictor after optional standardization/centering.
      real(dp), allocatable, intent(out) :: iv2(:) !! Second predictor after optional standardization/centering.
      real(dp), allocatable, intent(out) :: iv3(:) !! Third predictor after optional standardization/centering.
      real(dp), intent(in), optional :: es1 !! Main-effect coefficient for factor one; default zero.
      real(dp), intent(in), optional :: es2 !! Main-effect coefficient for factor two; default zero.
      real(dp), intent(in), optional :: es3 !! Main-effect coefficient for factor three; default zero.
      real(dp), intent(in), optional :: es12 !! Two-way interaction coefficient for factors one and two; default zero.
      real(dp), intent(in), optional :: es13 !! Two-way interaction coefficient for factors one and three; default zero.
      real(dp), intent(in), optional :: es23 !! Two-way interaction coefficient for factors two and three; default zero.
      real(dp), intent(in), optional :: es123 !! Three-way interaction coefficient; default zero.
      real(dp), intent(in), optional :: es11 !! Quadratic coefficient for factor one; default zero.
      real(dp), intent(in), optional :: es22 !! Quadratic coefficient for factor two; default zero.
      real(dp), intent(in), optional :: es33 !! Quadratic coefficient for factor three; default zero.
      logical, intent(in), optional :: standardize !! Standardize predictors to sample SD one when true; default true.
      logical, intent(in), optional :: center !! Keep zero-centered contrast coding when true; default true.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid sample sizes.
      real(dp), allocatable :: c1(:), c2(:), c3(:), noise(:)
      real(dp) :: b1, b2, b3, b12, b13, b23, b123, b11, b22, b33
      real(dp) :: m, s
      integer :: i1, i2, i3, rep, idx, l1, l2, l3, n
      logical :: do_std, do_center

      if (n_per_cell < 1 .or. n1 < 0 .or. n2 < 0 .or. n3 < 0) then
         allocate(observations(0), iv1(0), iv2(0), iv3(0))
         if (present(info)) info = 1
         return
      end if
      l1 = max(1, n1)
      l2 = max(1, n2)
      l3 = max(1, n3)
      n = n_per_cell * l1 * l2 * l3
      allocate(iv1(n), iv2(n), iv3(n), observations(n), noise(n))
      if (n1 > 0) call contrast_values(n1, c1)
      if (n2 > 0) call contrast_values(n2, c2)
      if (n3 > 0) call contrast_values(n3, c3)
      idx = 0
      do i3 = 1, l3
         do i2 = 1, l2
            do i1 = 1, l1
               do rep = 1, n_per_cell
                  idx = idx + 1
                  if (n1 > 0) then
                     iv1(idx) = c1(i1)
                  else
                     call fill_normal_scalar(iv1(idx))
                  end if
                  if (n2 > 0) then
                     iv2(idx) = c2(i2)
                  else
                     call fill_normal_scalar(iv2(idx))
                  end if
                  if (n3 > 0) then
                     iv3(idx) = c3(i3)
                  else
                     call fill_normal_scalar(iv3(idx))
                  end if
               end do
            end do
         end do
      end do
      do_std = .true.
      if (present(standardize)) do_std = standardize
      if (do_std) then
         call standardize_vector(iv1, m, s)
         call standardize_vector(iv2, m, s)
         call standardize_vector(iv3, m, s)
      end if
      call fill_normal_vector(noise)
      b1 = 0.0_dp
      b2 = 0.0_dp
      b3 = 0.0_dp
      b12 = 0.0_dp
      b13 = 0.0_dp
      b23 = 0.0_dp
      b123 = 0.0_dp
      b11 = 0.0_dp
      b22 = 0.0_dp
      b33 = 0.0_dp
      if (present(es1)) b1 = es1
      if (present(es2)) b2 = es2
      if (present(es3)) b3 = es3
      if (present(es12)) b12 = es12
      if (present(es13)) b13 = es13
      if (present(es23)) b23 = es23
      if (present(es123)) b123 = es123
      if (present(es11)) b11 = es11
      if (present(es22)) b22 = es22
      if (present(es33)) b33 = es33
      observations = b1 * iv1 + b2 * iv2 + b3 * iv3 + b12 * iv1 * iv2 + b13 * iv1 * iv3 + b23 * iv2 * iv3 + &
         b123 * iv1 * iv2 * iv3 + b11 * iv1 * iv1 + b22 * iv2 * iv2 + b33 * iv3 * iv3 + noise
      do_center = .true.
      if (present(center)) do_center = center
      if (.not. do_center) then
         iv1 = iv1 - minval(iv1) + 1.0_dp
         iv2 = iv2 - minval(iv2) + 1.0_dp
         iv3 = iv3 - minval(iv3) + 1.0_dp
      end if
      if (present(info)) info = 0
   end subroutine simulate_anova_data

   subroutine simulate_poly_ideal_logistic(n_cases, n_variables, items, probability, theta, discrimination, difficulty, &
      gamma, zeta, low, high, categories, mu, sigma, info)
      integer, intent(in) :: n_cases !! Number of simulated subjects.
      integer, intent(in) :: n_variables !! Number of ideal-point items.
      integer, allocatable, intent(out) :: items(:, :) !! Simulated binomial category counts from zero to `categories-1`.
      real(dp), allocatable, intent(out) :: probability(:, :) !! Ideal-point logistic response probabilities.
      real(dp), intent(in), optional :: theta(:) !! Optional supplied latent scores; otherwise normal scores are simulated.
      real(dp), intent(in), optional :: discrimination(:) !! Per-item discriminations; defaults to one.
      real(dp), intent(in), optional :: difficulty(:) !! Per-item ideal locations; defaults to an equally spaced sequence.
      real(dp), intent(in), optional :: gamma !! Lower asymptote parameter; default zero.
      real(dp), intent(in), optional :: zeta !! Upper asymptote parameter; default one.
      real(dp), intent(in), optional :: low !! Lowest default item location; default -2.
      real(dp), intent(in), optional :: high !! Highest default item location; default 2.
      integer, intent(in), optional :: categories !! Number of response categories; default five.
      real(dp), intent(in), optional :: mu !! Mean of simulated latent scores; default zero.
      real(dp), intent(in), optional :: sigma !! SD of simulated latent scores; default one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions.
      real(dp), allocatable :: latent(:), a(:), d(:)
      real(dp) :: c, z, lo, hi, mean_theta, sd_theta, q
      integer :: i, j, ncat, trials

      if (n_cases < 1 .or. n_variables < 1) then
         allocate(items(0, 0), probability(0, 0))
         if (present(info)) info = 1
         return
      end if
      if (present(theta)) then
         if (size(theta) /= n_cases) then
            allocate(items(0, 0), probability(0, 0))
            if (present(info)) info = 1
            return
         end if
      end if
      if (present(discrimination)) then
         if (size(discrimination) /= n_variables) then
            allocate(items(0, 0), probability(0, 0))
            if (present(info)) info = 1
            return
         end if
      end if
      if (present(difficulty)) then
         if (size(difficulty) /= n_variables) then
            allocate(items(0, 0), probability(0, 0))
            if (present(info)) info = 1
            return
         end if
      end if
      allocate(latent(n_cases), a(n_variables), d(n_variables), items(n_cases, n_variables), &
         probability(n_cases, n_variables))
      mean_theta = 0.0_dp
      sd_theta = 1.0_dp
      c = 0.0_dp
      z = 1.0_dp
      lo = -2.0_dp
      hi = 2.0_dp
      ncat = 5
      if (present(mu)) mean_theta = mu
      if (present(sigma)) sd_theta = sigma
      if (present(gamma)) c = gamma
      if (present(zeta)) z = zeta
      if (present(low)) lo = low
      if (present(high)) hi = high
      if (present(categories)) ncat = max(2, categories)
      if (present(theta)) then
         latent = theta
      else
         call fill_normal_vector(latent)
         latent = mean_theta + sd_theta * latent
      end if
      a = 1.0_dp
      if (present(discrimination)) a = discrimination
      call default_difficulty(lo, hi, d)
      if (present(difficulty)) d = difficulty
      do j = 1, n_variables
         do i = 1, n_cases
            q = exp(max(-700.0_dp, min(700.0_dp, a(j) * (d(j) - latent(i)))))
            probability(i, j) = 2.0_dp * (c + (z - c) * q / ((1.0_dp + q) ** 2))
            probability(i, j) = max(0.0_dp, min(1.0_dp, probability(i, j)))
         end do
      end do
      trials = ncat - 1
      call draw_binomial_matrix(probability, trials, items)
      if (present(info)) info = 0
   end subroutine simulate_poly_ideal_logistic

   subroutine simulate_poly_ideal_normal(n_cases, n_variables, items, probability, theta, discrimination, difficulty, &
      gamma, zeta, low, high, categories, mu, sigma, info)
      integer, intent(in) :: n_cases !! Number of simulated subjects.
      integer, intent(in) :: n_variables !! Number of ideal-point items.
      integer, allocatable, intent(out) :: items(:, :) !! Simulated binomial category counts from zero to `categories-1`.
      real(dp), allocatable, intent(out) :: probability(:, :) !! Upstream normal-ogive ideal-point probabilities.
      real(dp), intent(in), optional :: theta(:) !! Optional supplied latent scores; otherwise normal scores are simulated.
      real(dp), intent(in), optional :: discrimination(:) !! Per-item discriminations; defaults to one.
      real(dp), intent(in), optional :: difficulty(:) !! Per-item ideal locations; defaults to an equally spaced sequence.
      real(dp), intent(in), optional :: gamma !! Lower asymptote parameter; default zero.
      real(dp), intent(in), optional :: zeta !! Upper asymptote parameter; default one.
      real(dp), intent(in), optional :: low !! Lowest default item location; default -2.
      real(dp), intent(in), optional :: high !! Highest default item location; default 2.
      integer, intent(in), optional :: categories !! Number of response categories; default five.
      real(dp), intent(in), optional :: mu !! Mean of simulated latent scores; default zero.
      real(dp), intent(in), optional :: sigma !! SD of simulated latent scores; default one.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions.
      real(dp), allocatable :: latent(:), a(:), d(:)
      real(dp) :: c, z, lo, hi, mean_theta, sd_theta
      integer :: i, j, ncat, trials

      if (n_cases < 1 .or. n_variables < 1) then
         allocate(items(0, 0), probability(0, 0))
         if (present(info)) info = 1
         return
      end if
      if (present(theta)) then
         if (size(theta) /= n_cases) then
            allocate(items(0, 0), probability(0, 0))
            if (present(info)) info = 1
            return
         end if
      end if
      if (present(discrimination)) then
         if (size(discrimination) /= n_variables) then
            allocate(items(0, 0), probability(0, 0))
            if (present(info)) info = 1
            return
         end if
      end if
      if (present(difficulty)) then
         if (size(difficulty) /= n_variables) then
            allocate(items(0, 0), probability(0, 0))
            if (present(info)) info = 1
            return
         end if
      end if
      allocate(latent(n_cases), a(n_variables), d(n_variables), items(n_cases, n_variables), &
         probability(n_cases, n_variables))
      mean_theta = 0.0_dp
      sd_theta = 1.0_dp
      c = 0.0_dp
      z = 1.0_dp
      lo = -2.0_dp
      hi = 2.0_dp
      ncat = 5
      if (present(mu)) mean_theta = mu
      if (present(sigma)) sd_theta = sigma
      if (present(gamma)) c = gamma
      if (present(zeta)) z = zeta
      if (present(low)) lo = low
      if (present(high)) hi = high
      if (present(categories)) ncat = max(2, categories)
      if (present(theta)) then
         latent = theta
      else
         call fill_normal_vector(latent)
         latent = mean_theta + sd_theta * latent
      end if
      a = 1.0_dp
      if (present(discrimination)) a = discrimination
      call default_difficulty(lo, hi, d)
      if (present(difficulty)) d = difficulty
      do j = 1, n_variables
         do i = 1, n_cases
            probability(i, j) = c + (z - c) * normal_cdf(abs(a(j) * (latent(i) - d(j))))
            probability(i, j) = max(0.0_dp, min(1.0_dp, probability(i, j)))
         end do
      end do
      trials = ncat - 1
      call draw_binomial_matrix(probability, trials, items)
      if (present(info)) info = 0
   end subroutine simulate_poly_ideal_normal

   subroutine simulate_poly_ideal(n_cases, n_variables, items, probability, normal_model, theta, discrimination, &
      difficulty, gamma, zeta, low, high, categories, mu, sigma, info)
      integer, intent(in) :: n_cases !! Number of simulated subjects.
      integer, intent(in) :: n_variables !! Number of ideal-point items.
      integer, allocatable, intent(out) :: items(:, :) !! Simulated response-category counts.
      real(dp), allocatable, intent(out) :: probability(:, :) !! Ideal-point response probabilities.
      logical, intent(in), optional :: normal_model !! Use the normal-ogive variant when true; logistic otherwise.
      real(dp), intent(in), optional :: theta(:) !! Optional supplied latent scores.
      real(dp), intent(in), optional :: discrimination(:) !! Optional per-item discriminations.
      real(dp), intent(in), optional :: difficulty(:) !! Optional per-item ideal locations.
      real(dp), intent(in), optional :: gamma !! Lower asymptote; default zero.
      real(dp), intent(in), optional :: zeta !! Upper asymptote; default one.
      real(dp), intent(in), optional :: low !! Lowest default item location; default -2.
      real(dp), intent(in), optional :: high !! Highest default item location; default 2.
      integer, intent(in), optional :: categories !! Number of response categories; default five.
      real(dp), intent(in), optional :: mu !! Mean of simulated latent scores; default zero.
      real(dp), intent(in), optional :: sigma !! SD of simulated latent scores; default one.
      integer, intent(out), optional :: info !! Zero on success; propagated status from the selected simulator.
      logical :: use_normal
      integer :: local_info

      use_normal = .false.
      if (present(normal_model)) use_normal = normal_model
      if (use_normal) then
         call simulate_poly_ideal_normal(n_cases, n_variables, items, probability, theta, discrimination, difficulty, &
            gamma, zeta, low, high, categories, mu, sigma, local_info)
      else
         call simulate_poly_ideal_logistic(n_cases, n_variables, items, probability, theta, discrimination, difficulty, &
            gamma, zeta, low, high, categories, mu, sigma, local_info)
      end if
      if (present(info)) info = local_info
   end subroutine simulate_poly_ideal

   subroutine contrast_values(n_levels, values)
      integer, intent(in) :: n_levels !! Number of equally spaced contrast levels.
      real(dp), allocatable, intent(out) :: values(:) !! Symmetric integer contrast sequence matching upstream coding.
      integer :: i

      allocate(values(n_levels))
      if (mod(n_levels, 2) == 1) then
         do i = 1, n_levels
            values(i) = real(i - 1 - n_levels / 2, dp)
         end do
      else
         do i = 1, n_levels
            values(i) = real(-(n_levels - 1) + 2 * (i - 1), dp)
         end do
      end if
   end subroutine contrast_values

   subroutine standardize_vector(x, mean_value, sd_value)
      real(dp), intent(inout) :: x(:) !! Vector standardized in place.
      real(dp), intent(out) :: mean_value !! Sample mean before standardization.
      real(dp), intent(out) :: sd_value !! Sample standard deviation before standardization.

      mean_value = sum(x) / real(size(x), dp)
      if (size(x) > 1) then
         sd_value = sqrt(sum((x - mean_value) ** 2) / real(size(x) - 1, dp))
      else
         sd_value = 0.0_dp
      end if
      if (sd_value > tiny(1.0_dp)) then
         x = x / sd_value
      else
         x = 0.0_dp
      end if
   end subroutine standardize_vector

   subroutine default_difficulty(low, high, difficulty)
      real(dp), intent(in) :: low !! Lower endpoint of the equally spaced sequence.
      real(dp), intent(in) :: high !! Upper endpoint of the equally spaced sequence.
      real(dp), intent(out) :: difficulty(:) !! Equally spaced difficulty/ideal-point locations.
      integer :: j

      if (size(difficulty) == 1) then
         difficulty(1) = low
      else
         do j = 1, size(difficulty)
            difficulty(j) = low + (high - low) * real(j - 1, dp) / real(size(difficulty) - 1, dp)
         end do
      end if
   end subroutine default_difficulty

   subroutine draw_binomial_matrix(probability, trials, counts)
      real(dp), intent(in) :: probability(:, :) !! Success probabilities for each simulated subject and item.
      integer, intent(in) :: trials !! Number of independent Bernoulli trials per matrix element.
      integer, intent(out) :: counts(:, :) !! Binomial success counts.
      real(dp) :: u
      integer :: i, j, k

      counts = 0
      do j = 1, size(counts, 2)
         do i = 1, size(counts, 1)
            do k = 1, trials
               call random_number(u)
               if (u < probability(i, j)) counts(i, j) = counts(i, j) + 1
            end do
         end do
      end do
   end subroutine draw_binomial_matrix

   subroutine fill_normal_matrix(x)
      real(dp), intent(out) :: x(:, :) !! Matrix filled with independent standard-normal pseudorandom draws.
      integer :: j

      do j = 1, size(x, 2)
         call fill_normal_vector(x(:, j))
      end do
   end subroutine fill_normal_matrix

   subroutine fill_normal_vector(x)
      real(dp), intent(out) :: x(:) !! Vector filled with independent standard-normal pseudorandom draws.
      real(dp) :: u1, u2, radius, angle
      integer :: i

      i = 1
      do while (i <= size(x))
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         radius = sqrt(-2.0_dp * log(u1))
         angle = 2.0_dp * r_pi * u2
         x(i) = radius * cos(angle)
         if (i + 1 <= size(x)) x(i + 1) = radius * sin(angle)
         i = i + 2
      end do
   end subroutine fill_normal_vector

   subroutine fill_normal_scalar(x)
      real(dp), intent(out) :: x !! Scalar standard-normal pseudorandom draw.
      real(dp) :: one(1)

      call fill_normal_vector(one)
      x = one(1)
   end subroutine fill_normal_scalar

end module psych_simulation_algorithms
