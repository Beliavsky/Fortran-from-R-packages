! Modern Fortran translation of the GEE engine in src/gee2.cc and src/geesubs.cc.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_gee
   use r_kinds, only : dp
   use geepack_status, only : GEE_OK, GEE_ERR_SHAPE, GEE_ERR_ARGUMENT, GEE_ERR_SINGULAR, &
      GEE_ERR_MAXITER, GEE_ERR_INVALID_MEAN, GEE_ERR_CORRELATION
   use geepack_links, only : LINK_IDENTITY, link_function, link_inverse, link_derivative, &
      variance_function, variance_derivative, valid_mean
   use geepack_correlations, only : COR_INDEPENDENCE, COR_EXCHANGEABLE, COR_AR1, COR_UNSTRUCTURED, &
      COR_USERDEFINED, COR_FIXED, working_correlation, upper_triangle, pair_products
   use geepack_design, only : gen_zcor
   use geepack_matrix, only : inverse_checked, solve_checked, outer_product, diagonal_matrix, max_abs
   implicit none
   private

   type, public :: gee_spec
      integer :: corstr = COR_INDEPENDENCE
      integer :: corr_link = LINK_IDENTITY
      logical :: scale_fixed = .false.
      real(dp) :: scale_value = 1.0_dp
      real(dp) :: tolerance = 1.0e-4_dp
      integer :: max_iterations = 25
      logical :: approximate_jackknife = .false.
      logical :: one_step_jackknife = .false.
      logical :: fully_iterated_jackknife = .false.
      integer, allocatable :: mean_links(:)
      integer, allocatable :: variance_codes(:)
      integer, allocatable :: scale_links(:)
   end type gee_spec

   type, public :: gee_result
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: vbeta(:, :)
      real(dp), allocatable :: valpha(:, :)
      real(dp), allocatable :: vgamma(:, :)
      real(dp), allocatable :: vbeta_naive(:, :)
      real(dp), allocatable :: valpha_naive(:, :)
      real(dp), allocatable :: valpha_stable(:, :)
      real(dp), allocatable :: vbeta_ajs(:, :)
      real(dp), allocatable :: valpha_ajs(:, :)
      real(dp), allocatable :: vgamma_ajs(:, :)
      real(dp), allocatable :: vbeta_j1s(:, :)
      real(dp), allocatable :: valpha_j1s(:, :)
      real(dp), allocatable :: vgamma_j1s(:, :)
      real(dp), allocatable :: vbeta_fij(:, :)
      real(dp), allocatable :: valpha_fij(:, :)
      real(dp), allocatable :: vgamma_fij(:, :)
      real(dp), allocatable :: influence(:, :)
      integer :: iterations = 0
      integer :: error = GEE_OK
   end type gee_result

   public :: fit_geese

contains

   subroutine fit_geese(y, x, cluster_sizes, beta_initial, spec, result, alpha_initial, gamma_initial, &
      offset, scale_offset, weights, waves, zsca, zcor, cor_param)
      real(dp), intent(in) :: y(:) !! Response vector ordered by contiguous clusters.
      real(dp), intent(in) :: x(:, :) !! Mean-model design matrix, one row per response.
      integer, intent(in) :: cluster_sizes(:) !! Number of observations in each contiguous cluster.
      real(dp), intent(in) :: beta_initial(:) !! Starting mean-regression coefficients.
      type(gee_spec), intent(in) :: spec !! Links, variance family, working correlation, and controls.
      type(gee_result), intent(out) :: result !! Fitted parameters, covariance estimates, influences, and status.
      real(dp), optional, intent(in) :: alpha_initial(:) !! Starting association parameters; defaults to zero, or one for fixed.
      real(dp), optional, intent(in) :: gamma_initial(:) !! Starting scale-model coefficients.
      real(dp), optional, intent(in) :: offset(:) !! Mean-model offset; defaults to zero.
      real(dp), optional, intent(in) :: scale_offset(:) !! Scale-model offset; defaults to zero.
      real(dp), optional, intent(in) :: weights(:) !! Observation weights; defaults to one.
      integer, optional, intent(in) :: waves(:) !! One-based wave indices; defaults to within-cluster order.
      real(dp), optional, intent(in) :: zsca(:, :) !! Scale-model design matrix; defaults to an intercept column.
      real(dp), optional, intent(in) :: zcor(:, :) !! Association design matrix; generated from waves when absent.
      real(dp), optional, intent(in) :: cor_param(:) !! Known correlation metadata; defaults numerically to waves.
      real(dp), allocatable :: off(:)
      real(dp), allocatable :: soff(:)
      real(dp), allocatable :: wt(:)
      real(dp), allocatable :: zsca_work(:, :)
      real(dp), allocatable :: zcor_work(:, :)
      real(dp), allocatable :: corp(:)
      integer, allocatable :: wave_work(:)
      integer, allocatable :: mean_links(:)
      integer, allocatable :: variance_codes(:)
      integer, allocatable :: scale_links(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: alpha(:)
      integer :: status
      integer :: n
      integer :: r
      integer :: q
      integer :: iterations_work

      call initialize_empty_result(result, size(beta_initial), 0, 0, size(cluster_sizes))
      n = size(y)
      if (size(x, 1) /= n .or. size(x, 2) /= size(beta_initial) .or. sum(cluster_sizes) /= n .or. &
          any(cluster_sizes < 1)) then
         result%error = GEE_ERR_SHAPE
         return
      end if

      call setup_vectors(n, cluster_sizes, offset, scale_offset, weights, waves, off, soff, wt, wave_work, status)
      if (status /= GEE_OK) then
         result%error = status
         return
      end if
      call setup_links(spec, maxval(wave_work), mean_links, variance_codes, scale_links, status)
      if (status /= GEE_OK) then
         result%error = status
         return
      end if

      if (present(zsca)) then
         if (size(zsca, 1) /= n) then
            result%error = GEE_ERR_SHAPE
            return
         end if
         allocate(zsca_work(size(zsca, 1), size(zsca, 2)))
         zsca_work = zsca
      else
         allocate(zsca_work(n, 1))
         zsca_work = 1.0_dp
      end if
      r = size(zsca_work, 2)

      if (present(zcor)) then
         allocate(zcor_work(size(zcor, 1), size(zcor, 2)))
         zcor_work = zcor
      else
         call gen_zcor(cluster_sizes, wave_work, spec%corstr, zcor_work, status)
         if (status /= GEE_OK) then
            result%error = status
            return
         end if
      end if
      q = size(zcor_work, 2)

      allocate(corp(n))
      if (present(cor_param)) then
         if (size(cor_param) /= n) then
            result%error = GEE_ERR_SHAPE
            return
         end if
         corp = cor_param
      else
         corp = real(wave_work, dp)
      end if

      allocate(beta(size(beta_initial)))
      beta = beta_initial
      allocate(gamma(r))
      if (present(gamma_initial)) then
         if (size(gamma_initial) /= r) then
            result%error = GEE_ERR_SHAPE
            return
         end if
         gamma = gamma_initial
      else
         gamma = 0.0_dp
         if (r >= 1) gamma(1) = link_function(spec%scale_value, scale_links(1))
      end if
      allocate(alpha(q))
      if (q > 0) then
         if (present(alpha_initial)) then
            if (size(alpha_initial) /= q) then
               result%error = GEE_ERR_SHAPE
               return
            end if
            alpha = alpha_initial
         else if (spec%corstr == COR_FIXED) then
            alpha = 1.0_dp
         else
            alpha = 0.0_dp
         end if
      end if

      call validate_zcor_shape(cluster_sizes, spec%corstr, zcor_work, status)
      if (status /= GEE_OK) then
         result%error = status
         return
      end if
      call estimate_parameters(y, x, cluster_sizes, off, soff, wt, wave_work, zsca_work, zcor_work, corp, &
         mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, 0, spec%max_iterations, &
         iterations_work, status)

      call initialize_empty_result(result, size(beta), size(alpha), size(gamma), size(cluster_sizes))
      result%beta = beta
      result%alpha = alpha
      result%gamma = gamma
      result%iterations = max(iterations_work, 0)
      result%error = status
      if (status /= GEE_OK .and. status /= GEE_ERR_MAXITER) return

      call compute_covariance(y, x, cluster_sizes, off, soff, wt, wave_work, zsca_work, zcor_work, corp, &
         mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, result, status)
      if (status /= GEE_OK) then
         result%error = status
         return
      end if
      if (spec%approximate_jackknife) then
         call approximate_jackknife(y, x, cluster_sizes, off, soff, wt, wave_work, zsca_work, zcor_work, corp, &
            mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, result, status)
         if (status /= GEE_OK) result%error = status
      end if
      if (spec%one_step_jackknife) then
         call refit_jackknife(y, x, cluster_sizes, off, soff, wt, wave_work, zsca_work, zcor_work, corp, &
            mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, 1, &
            result%vbeta_j1s, result%vgamma_j1s, result%valpha_j1s, status)
         if (status /= GEE_OK) result%error = status
      end if
      if (spec%fully_iterated_jackknife) then
         call refit_jackknife(y, x, cluster_sizes, off, soff, wt, wave_work, zsca_work, zcor_work, corp, &
            mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, spec%max_iterations, &
            result%vbeta_fij, result%vgamma_fij, result%valpha_fij, status)
         if (status /= GEE_OK) result%error = status
      end if
   end subroutine fit_geese

   subroutine initialize_empty_result(result, p, q, r, ncluster)
      type(gee_result), intent(out) :: result !! Result object to allocate and zero.
      integer, intent(in) :: p !! Number of mean parameters.
      integer, intent(in) :: q !! Number of association parameters.
      integer, intent(in) :: r !! Number of scale parameters.
      integer, intent(in) :: ncluster !! Number of clusters, used for the influence matrix.

      allocate(result%beta(p), result%alpha(q), result%gamma(r))
      allocate(result%vbeta(p, p), result%valpha(q, q), result%vgamma(r, r))
      allocate(result%vbeta_naive(p, p), result%valpha_naive(q, q), result%valpha_stable(q, q))
      allocate(result%vbeta_ajs(p, p), result%valpha_ajs(q, q), result%vgamma_ajs(r, r))
      allocate(result%vbeta_j1s(p, p), result%valpha_j1s(q, q), result%vgamma_j1s(r, r))
      allocate(result%vbeta_fij(p, p), result%valpha_fij(q, q), result%vgamma_fij(r, r))
      allocate(result%influence(p + r + q, ncluster))
      result%beta = 0.0_dp
      result%alpha = 0.0_dp
      result%gamma = 0.0_dp
      result%vbeta = 0.0_dp
      result%valpha = 0.0_dp
      result%vgamma = 0.0_dp
      result%vbeta_naive = 0.0_dp
      result%valpha_naive = 0.0_dp
      result%valpha_stable = 0.0_dp
      result%vbeta_ajs = 0.0_dp
      result%valpha_ajs = 0.0_dp
      result%vgamma_ajs = 0.0_dp
      result%vbeta_j1s = 0.0_dp
      result%valpha_j1s = 0.0_dp
      result%vgamma_j1s = 0.0_dp
      result%vbeta_fij = 0.0_dp
      result%valpha_fij = 0.0_dp
      result%vgamma_fij = 0.0_dp
      result%influence = 0.0_dp
      result%iterations = 0
      result%error = GEE_OK
   end subroutine initialize_empty_result

   subroutine setup_vectors(n, cluster_sizes, offset, scale_offset, weights, waves, off, soff, wt, wave_work, status)
      integer, intent(in) :: n !! Number of observations.
      integer, intent(in) :: cluster_sizes(:) !! Number of observations in each cluster.
      real(dp), optional, intent(in) :: offset(:) !! Optional mean offsets.
      real(dp), optional, intent(in) :: scale_offset(:) !! Optional scale offsets.
      real(dp), optional, intent(in) :: weights(:) !! Optional observation weights.
      integer, optional, intent(in) :: waves(:) !! Optional one-based wave indices.
      real(dp), allocatable, intent(out) :: off(:) !! Materialized mean offsets.
      real(dp), allocatable, intent(out) :: soff(:) !! Materialized scale offsets.
      real(dp), allocatable, intent(out) :: wt(:) !! Materialized weights.
      integer, allocatable, intent(out) :: wave_work(:) !! Materialized wave indices.
      integer, intent(out) :: status !! GEE_OK on success or an error code.
      integer :: g
      integer :: i
      integer :: pos

      allocate(off(n), soff(n), wt(n), wave_work(n))
      off = 0.0_dp
      soff = 0.0_dp
      wt = 1.0_dp
      if (present(offset)) then
         if (size(offset) /= n) then
            status = GEE_ERR_SHAPE
            return
         end if
         off = offset
      end if
      if (present(scale_offset)) then
         if (size(scale_offset) /= n) then
            status = GEE_ERR_SHAPE
            return
         end if
         soff = scale_offset
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         wt = weights
      end if
      if (present(waves)) then
         if (size(waves) /= n .or. any(waves < 1)) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         wave_work = waves
      else
         pos = 0
         do g = 1, size(cluster_sizes)
            do i = 1, cluster_sizes(g)
               wave_work(pos + i) = i
            end do
            pos = pos + cluster_sizes(g)
         end do
      end if
      status = GEE_OK
   end subroutine setup_vectors

   subroutine setup_links(spec, maxwave, mean_links, variance_codes, scale_links, status)
      type(gee_spec), intent(in) :: spec !! Model specification containing optional per-wave codes.
      integer, intent(in) :: maxwave !! Largest one-based wave index in the data.
      integer, allocatable, intent(out) :: mean_links(:) !! Mean-link code for each wave.
      integer, allocatable, intent(out) :: variance_codes(:) !! Variance-function code for each wave.
      integer, allocatable, intent(out) :: scale_links(:) !! Scale-link code for each wave.
      integer, intent(out) :: status !! GEE_OK on success or an error code.

      call expand_codes(spec%mean_links, maxwave, LINK_IDENTITY, mean_links, status)
      if (status /= GEE_OK) return
      call expand_codes(spec%variance_codes, maxwave, 1, variance_codes, status)
      if (status /= GEE_OK) return
      call expand_codes(spec%scale_links, maxwave, LINK_IDENTITY, scale_links, status)
   end subroutine setup_links

   subroutine expand_codes(source, n, default_code, target, status)
      integer, allocatable, intent(in) :: source(:) !! Optional scalar or per-wave code vector.
      integer, intent(in) :: n !! Number of wave positions required.
      integer, intent(in) :: default_code !! Code used when source is unallocated.
      integer, allocatable, intent(out) :: target(:) !! Expanded code vector of length n.
      integer, intent(out) :: status !! GEE_OK on success or an error code.

      allocate(target(n))
      if (.not. allocated(source)) then
         target = default_code
      else if (size(source) == 1) then
         target = source(1)
      else if (size(source) == n) then
         target = source
      else
         target = default_code
         status = GEE_ERR_ARGUMENT
         return
      end if
      status = GEE_OK
   end subroutine expand_codes

   subroutine validate_zcor_shape(cluster_sizes, corstr, zcor, status)
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes used to determine required association rows.
      integer, intent(in) :: corstr !! Correlation structure identifier.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix to validate.
      integer, intent(out) :: status !! GEE_OK on success or an error code.
      integer :: expected

      select case (corstr)
      case (COR_INDEPENDENCE)
         expected = 0
      case (COR_EXCHANGEABLE, COR_AR1)
         expected = size(cluster_sizes)
      case default
         expected = sum(cluster_sizes * (cluster_sizes - 1) / 2)
      end select
      if (size(zcor, 1) /= expected) then
         status = GEE_ERR_SHAPE
      else
         status = GEE_OK
      end if
   end subroutine validate_zcor_shape

   subroutine estimate_parameters(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
      mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, omit_cluster, max_iterations, iterations, status)
      real(dp), intent(in) :: y(:) !! Response vector.
      real(dp), intent(in) :: x(:, :) !! Mean design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      real(dp), intent(in) :: offset(:) !! Mean offsets.
      real(dp), intent(in) :: scale_offset(:) !! Scale offsets.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: zsca(:, :) !! Scale design matrix.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix.
      real(dp), intent(in) :: cor_param(:) !! Correlation metadata, retained for upstream API parity.
      integer, intent(in) :: mean_links(:) !! Mean-link codes by wave.
      integer, intent(in) :: variance_codes(:) !! Variance codes by wave.
      integer, intent(in) :: scale_links(:) !! Scale-link codes by wave.
      type(gee_spec), intent(in) :: spec !! Model and iteration controls.
      real(dp), intent(inout) :: beta(:) !! Mean coefficients, updated in place.
      real(dp), intent(inout) :: gamma(:) !! Scale coefficients, updated in place.
      real(dp), intent(inout) :: alpha(:) !! Association coefficients, updated in place.
      integer, intent(in) :: omit_cluster !! Cluster index to omit, or zero for all clusters.
      integer, intent(in) :: max_iterations !! Maximum number of alternating updates.
      integer, intent(out) :: iterations !! Number of update sweeps attempted.
      integer, intent(out) :: status !! GEE_OK, GEE_ERR_MAXITER, or a numerical error code.
      real(dp), allocatable :: pr(:)
      real(dp), allocatable :: phi(:)
      real(dp) :: db
      real(dp) :: dg
      real(dp) :: da
      real(dp) :: delta
      integer :: iter

      allocate(pr(size(y)), phi(size(y)))
      status = GEE_OK
      iterations = 0
      do iter = 1, max_iterations
         call compute_phi(scale_offset, zsca, waves, scale_links, gamma, spec, phi)
         call update_beta(y, x, cluster_sizes, offset, weights, waves, zcor, cor_param, mean_links, &
            variance_codes, spec, phi, alpha, beta, omit_cluster, db, status)
         if (status /= GEE_OK) return
         call compute_pr(y, x, offset, waves, mean_links, variance_codes, beta, pr, status)
         if (status /= GEE_OK) return
         call update_gamma(pr, weights, cluster_sizes, waves, scale_offset, zsca, scale_links, spec, gamma, &
            omit_cluster, dg, status)
         if (status /= GEE_OK) return
         call compute_phi(scale_offset, zsca, waves, scale_links, gamma, spec, phi)
         call update_alpha(pr, phi, weights, cluster_sizes, waves, zcor, cor_param, spec, alpha, omit_cluster, da, status)
         if (status /= GEE_OK) return
         delta = max(db, max(dg, da))
         iterations = iter
         if (delta <= spec%tolerance) return
      end do
      status = GEE_ERR_MAXITER
   end subroutine estimate_parameters

   subroutine update_beta(y, x, cluster_sizes, offset, weights, waves, zcor, cor_param, mean_links, variance_codes, &
      spec, phi, alpha, beta, omit_cluster, delta, status)
      real(dp), intent(in) :: y(:) !! Response vector.
      real(dp), intent(in) :: x(:, :) !! Mean design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      real(dp), intent(in) :: offset(:) !! Mean offsets.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix.
      real(dp), intent(in) :: cor_param(:) !! Correlation metadata retained for parity.
      integer, intent(in) :: mean_links(:) !! Mean-link codes by wave.
      integer, intent(in) :: variance_codes(:) !! Variance codes by wave.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(in) :: phi(:) !! Current scale values.
      real(dp), intent(in) :: alpha(:) !! Current association coefficients.
      real(dp), intent(inout) :: beta(:) !! Mean coefficients, updated in place.
      integer, intent(in) :: omit_cluster !! Cluster to omit, or zero.
      real(dp), intent(out) :: delta !! Maximum absolute coefficient change.
      integer, intent(out) :: status !! GEE_OK or a numerical error code.
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: g(:)
      real(dp), allocatable :: step(:)
      real(dp), allocatable :: rmat(:, :)
      real(dp), allocatable :: rinv(:, :)
      real(dp), allocatable :: pr(:)
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: v(:)
      real(dp), allocatable :: vmu(:)
      real(dp), allocatable :: dw(:, :)
      real(dp), allocatable :: prw(:)
      real(dp), allocatable :: beta_new(:)
      integer :: gidx
      integer :: obs0
      integer :: z0
      integer :: s
      integer :: npair
      integer :: info
      integer :: tries

      allocate(h(size(beta), size(beta)), g(size(beta)), step(size(beta)), beta_new(size(beta)))
      h = 0.0_dp
      g = 0.0_dp
      obs0 = 0
      z0 = 0
      do gidx = 1, size(cluster_sizes)
         s = cluster_sizes(gidx)
         npair = s * (s - 1) / 2
         if (gidx == omit_cluster) then
            obs0 = obs0 + s
            z0 = z0 + association_rows(spec%corstr, s)
            cycle
         end if
         allocate(pr(s), d(s, size(beta)), mu(s), v(s), vmu(s), rmat(s, s), rinv(s, s), dw(s, size(beta)), prw(s))
         call mean_prep(y(obs0 + 1:obs0 + s), x(obs0 + 1:obs0 + s, :), offset(obs0 + 1:obs0 + s), &
            waves(obs0 + 1:obs0 + s), mean_links, variance_codes, beta, pr, d, mu, v, vmu, status)
         if (status /= GEE_OK) return
         call cluster_correlation(spec, alpha, zcor, z0, waves(obs0 + 1:obs0 + s), &
            cor_param(obs0 + 1:obs0 + s), npair, rmat, status)
         if (status /= GEE_OK) return
         call inverse_checked(rmat, rinv, info)
         if (info /= 0) then
            status = GEE_ERR_SINGULAR
            return
         end if
         dw = d
         prw = pr
         call row_scale_matrix(dw, sqrt(weights(obs0 + 1:obs0 + s) / phi(obs0 + 1:obs0 + s)))
         prw = prw * sqrt(weights(obs0 + 1:obs0 + s) / phi(obs0 + 1:obs0 + s))
         h = h + matmul(transpose(dw), matmul(rinv, dw))
         g = g + matmul(transpose(dw), matmul(rinv, prw))
         deallocate(pr, d, mu, v, vmu, rmat, rinv, dw, prw)
         obs0 = obs0 + s
         z0 = z0 + association_rows(spec%corstr, s)
      end do
      call solve_checked(h, g, step, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      beta_new = beta + step
      do tries = 1, 60
         if (all_valid_fitted(x, offset, waves, mean_links, variance_codes, beta_new)) exit
         step = 0.5_dp * step
         beta_new = beta + step
      end do
      if (.not. all_valid_fitted(x, offset, waves, mean_links, variance_codes, beta_new)) then
         status = GEE_ERR_INVALID_MEAN
         return
      end if
      beta = beta_new
      delta = max_abs(step)
      status = GEE_OK
   end subroutine update_beta

   subroutine update_gamma(pr, weights, cluster_sizes, waves, scale_offset, zsca, scale_links, spec, gamma, &
      omit_cluster, delta, status)
      real(dp), intent(in) :: pr(:) !! Pearson residuals standardized by the mean variance function.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: scale_offset(:) !! Scale offsets.
      real(dp), intent(in) :: zsca(:, :) !! Scale design matrix.
      integer, intent(in) :: scale_links(:) !! Scale-link codes by wave.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(inout) :: gamma(:) !! Scale coefficients, updated in place.
      integer, intent(in) :: omit_cluster !! Cluster to omit, or zero.
      real(dp), intent(out) :: delta !! Maximum absolute coefficient change.
      integer, intent(out) :: status !! GEE_OK or a numerical error code.
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: g(:)
      real(dp), allocatable :: step(:)
      real(dp), allocatable :: phi(:)
      real(dp), allocatable :: d2(:, :)
      real(dp), allocatable :: wv(:)
      integer :: obs0
      integer :: s
      integer :: gidx
      integer :: info

      delta = 0.0_dp
      status = GEE_OK
      if (spec%scale_fixed .or. size(gamma) == 0) return
      allocate(h(size(gamma), size(gamma)), g(size(gamma)), step(size(gamma)))
      h = 0.0_dp
      g = 0.0_dp
      obs0 = 0
      do gidx = 1, size(cluster_sizes)
         s = cluster_sizes(gidx)
         if (gidx == omit_cluster) then
            obs0 = obs0 + s
            cycle
         end if
         allocate(phi(s), d2(s, size(gamma)), wv(s))
         call scale_prep(scale_offset(obs0 + 1:obs0 + s), zsca(obs0 + 1:obs0 + s, :), &
            waves(obs0 + 1:obs0 + s), scale_links, gamma, phi, d2)
         wv = weights(obs0 + 1:obs0 + s) / (2.0_dp * phi)
         h = h + matmul(transpose(d2), row_scaled_copy(d2, wv))
         g = g + matmul(transpose(d2), wv * (pr(obs0 + 1:obs0 + s) ** 2 - phi))
         deallocate(phi, d2, wv)
         obs0 = obs0 + s
      end do
      call solve_checked(h, g, step, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      gamma = gamma + step
      delta = max_abs(step)
   end subroutine update_gamma

   subroutine update_alpha(pr, phi, weights, cluster_sizes, waves, zcor, cor_param, spec, alpha, omit_cluster, &
      delta, status)
      real(dp), intent(in) :: pr(:) !! Pearson residuals standardized by the mean variance function.
      real(dp), intent(in) :: phi(:) !! Current scale values.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix.
      real(dp), intent(in) :: cor_param(:) !! Correlation metadata retained for upstream parity.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(inout) :: alpha(:) !! Association coefficients, updated in place.
      integer, intent(in) :: omit_cluster !! Cluster to omit, or zero.
      real(dp), intent(out) :: delta !! Maximum absolute association-parameter change.
      integer, intent(out) :: status !! GEE_OK or a numerical error code.
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: g(:)
      real(dp), allocatable :: step(:)
      real(dp), allocatable :: rmat(:, :)
      real(dp), allocatable :: e(:, :)
      real(dp), allocatable :: zi(:)
      real(dp), allocatable :: rhoi(:)
      real(dp), allocatable :: pairw(:)
      real(dp), allocatable :: spr(:)
      integer :: gidx
      integer :: obs0
      integer :: z0
      integer :: s
      integer :: npair
      integer :: info

      delta = 0.0_dp
      status = GEE_OK
      if (.not. association_active(spec%corstr, size(alpha))) return
      allocate(h(size(alpha), size(alpha)), g(size(alpha)), step(size(alpha)))
      h = 0.0_dp
      g = 0.0_dp
      obs0 = 0
      z0 = 0
      do gidx = 1, size(cluster_sizes)
         s = cluster_sizes(gidx)
         npair = s * (s - 1) / 2
         if (gidx == omit_cluster .or. s == 1) then
            obs0 = obs0 + s
            z0 = z0 + association_rows(spec%corstr, s)
            cycle
         end if
         allocate(rmat(s, s), e(npair, size(alpha)), zi(npair), rhoi(npair), pairw(npair), spr(s))
         call cluster_correlation_and_e(spec, alpha, zcor, z0, waves(obs0 + 1:obs0 + s), &
            cor_param(obs0 + 1:obs0 + s), npair, rmat, e, status)
         if (status /= GEE_OK) return
         spr = pr(obs0 + 1:obs0 + s) / sqrt(phi(obs0 + 1:obs0 + s))
         call pair_products(spr, zi)
         call upper_triangle(rmat, rhoi)
         call pair_products(weights(obs0 + 1:obs0 + s), pairw)
         h = h + matmul(transpose(e), row_scaled_copy(e, pairw))
         g = g + matmul(transpose(e), pairw * (zi - rhoi))
         deallocate(rmat, e, zi, rhoi, pairw, spr)
         obs0 = obs0 + s
         z0 = z0 + association_rows(spec%corstr, s)
      end do
      call solve_checked(h, g, step, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      alpha = alpha + step
      delta = max_abs(step)
   end subroutine update_alpha

   subroutine compute_covariance(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
      mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, result, status)
      real(dp), intent(in) :: y(:) !! Response vector.
      real(dp), intent(in) :: x(:, :) !! Mean design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      real(dp), intent(in) :: offset(:) !! Mean offsets.
      real(dp), intent(in) :: scale_offset(:) !! Scale offsets.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: zsca(:, :) !! Scale design matrix.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix.
      real(dp), intent(in) :: cor_param(:) !! Correlation metadata retained for parity.
      integer, intent(in) :: mean_links(:) !! Mean-link codes.
      integer, intent(in) :: variance_codes(:) !! Variance codes.
      integer, intent(in) :: scale_links(:) !! Scale-link codes.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(in) :: beta(:) !! Fitted mean coefficients.
      real(dp), intent(in) :: gamma(:) !! Fitted scale coefficients.
      real(dp), intent(in) :: alpha(:) !! Fitted association coefficients.
      type(gee_result), intent(inout) :: result !! Result receiving covariance and influence estimates.
      integer, intent(out) :: status !! GEE_OK or a numerical error code.
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: hinv(:, :)
      real(dp), allocatable :: meat(:, :)
      real(dp), allocatable :: hi(:, :)
      real(dp), allocatable :: score(:)
      real(dp), allocatable :: stable(:, :)
      integer :: p
      integer :: r
      integer :: q
      integer :: l
      integer :: gidx
      integer :: info

      p = size(beta)
      r = size(gamma)
      q = size(alpha)
      l = p + r + q
      allocate(h(l, l), hinv(l, l), meat(l, l), hi(l, l), score(l))
      h = 0.0_dp
      meat = 0.0_dp
      do gidx = 1, size(cluster_sizes)
         call cluster_hg(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
            mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, gidx, hi, score, status)
         if (status /= GEE_OK) return
         h = h + hi
         meat = meat + outer_product(score, score)
         result%influence(:, gidx) = score
      end do
      call regularize_inactive_blocks(h, p, r, q, spec)
      call inverse_checked(h, hinv, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      result%influence = matmul(hinv, result%influence)
      stable = matmul(hinv, matmul(meat, transpose(hinv)))
      result%vbeta = stable(1:p, 1:p)
      result%vbeta_naive = hinv(1:p, 1:p)
      if (r > 0) result%vgamma = stable(p + 1:p + r, p + 1:p + r)
      if (q > 0) then
         result%valpha = stable(p + r + 1:l, p + r + 1:l)
         if (association_active(spec%corstr, q)) then
            result%valpha_naive = hinv(p + r + 1:l, p + r + 1:l)
         end if
      end if
      if (q > 0 .and. association_active(spec%corstr, q)) then
         call alpha_stable_covariance(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
            mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, result%valpha_stable, status)
      else
         status = GEE_OK
      end if
   end subroutine compute_covariance

   subroutine alpha_stable_covariance(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
      mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, covariance, status)
      real(dp), intent(in) :: y(:) !! Response vector.
      real(dp), intent(in) :: x(:, :) !! Mean design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      real(dp), intent(in) :: offset(:) !! Mean offsets.
      real(dp), intent(in) :: scale_offset(:) !! Scale offsets.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: zsca(:, :) !! Scale design matrix.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix.
      real(dp), intent(in) :: cor_param(:) !! Correlation metadata retained for parity.
      integer, intent(in) :: mean_links(:) !! Mean-link codes.
      integer, intent(in) :: variance_codes(:) !! Variance codes.
      integer, intent(in) :: scale_links(:) !! Scale-link codes.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(in) :: beta(:) !! Fitted mean coefficients.
      real(dp), intent(in) :: gamma(:) !! Fitted scale coefficients.
      real(dp), intent(in) :: alpha(:) !! Fitted association coefficients.
      real(dp), intent(out) :: covariance(:, :) !! Stable alpha sandwich using only the association equations.
      integer, intent(out) :: status !! GEE_OK or a numerical error code.
      real(dp), allocatable :: f(:, :)
      real(dp), allocatable :: finv(:, :)
      real(dp), allocatable :: l33(:, :)
      real(dp), allocatable :: hi(:, :)
      real(dp), allocatable :: score(:)
      integer :: p
      integer :: r
      integer :: q
      integer :: l
      integer :: gidx
      integer :: info

      p = size(beta)
      r = size(gamma)
      q = size(alpha)
      l = p + r + q
      allocate(f(q, q), finv(q, q), l33(q, q), hi(l, l), score(l))
      f = 0.0_dp
      l33 = 0.0_dp
      do gidx = 1, size(cluster_sizes)
         call cluster_hg(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
            mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, gidx, hi, score, status)
         if (status /= GEE_OK) return
         f = f + hi(p + r + 1:l, p + r + 1:l)
         l33 = l33 + outer_product(score(p + r + 1:l), score(p + r + 1:l))
      end do
      call inverse_checked(f, finv, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      covariance = matmul(finv, matmul(l33, transpose(finv)))
      status = GEE_OK
   end subroutine alpha_stable_covariance

   subroutine cluster_hg(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
      mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, cluster_index, h, score, status)
      real(dp), intent(in) :: y(:) !! Response vector.
      real(dp), intent(in) :: x(:, :) !! Mean design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      real(dp), intent(in) :: offset(:) !! Mean offsets.
      real(dp), intent(in) :: scale_offset(:) !! Scale offsets.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: zsca(:, :) !! Scale design matrix.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix.
      real(dp), intent(in) :: cor_param(:) !! Correlation metadata retained for parity.
      integer, intent(in) :: mean_links(:) !! Mean-link codes.
      integer, intent(in) :: variance_codes(:) !! Variance codes.
      integer, intent(in) :: scale_links(:) !! Scale-link codes.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(in) :: beta(:) !! Fitted mean coefficients.
      real(dp), intent(in) :: gamma(:) !! Fitted scale coefficients.
      real(dp), intent(in) :: alpha(:) !! Fitted association coefficients.
      integer, intent(in) :: cluster_index !! One-based cluster whose sensitivity and score are returned.
      real(dp), intent(out) :: h(:, :) !! Lower-block sensitivity matrix in beta,gamma,alpha order.
      real(dp), intent(out) :: score(:) !! Cluster estimating-equation score vector.
      integer, intent(out) :: status !! GEE_OK or a numerical error code.
      real(dp), allocatable :: pr(:)
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: v(:)
      real(dp), allocatable :: vmu(:)
      real(dp), allocatable :: phi(:)
      real(dp), allocatable :: d2(:, :)
      real(dp), allocatable :: rmat(:, :)
      real(dp), allocatable :: rinv(:, :)
      real(dp), allocatable :: e(:, :)
      real(dp), allocatable :: dw(:, :)
      real(dp), allocatable :: prw(:)
      real(dp), allocatable :: wv(:)
      real(dp), allocatable :: sbeta(:, :)
      real(dp), allocatable :: zi(:)
      real(dp), allocatable :: rhoi(:)
      real(dp), allocatable :: spr(:)
      real(dp), allocatable :: pairw(:)
      real(dp), allocatable :: zbeta(:, :)
      real(dp), allocatable :: zgamma(:, :)
      integer :: p
      integer :: r
      integer :: q
      integer :: l
      integer :: obs0
      integer :: z0
      integer :: s
      integer :: npair
      integer :: info
      integer :: gidx

      p = size(beta)
      r = size(gamma)
      q = size(alpha)
      l = p + r + q
      h = 0.0_dp
      score = 0.0_dp
      if (cluster_index > 1) then
         obs0 = sum(cluster_sizes(1:cluster_index - 1))
      else
         obs0 = 0
      end if
      z0 = 0
      do gidx = 1, cluster_index - 1
         z0 = z0 + association_rows(spec%corstr, cluster_sizes(gidx))
      end do
      s = cluster_sizes(cluster_index)
      npair = s * (s - 1) / 2
      allocate(pr(s), d(s, p), mu(s), v(s), vmu(s), phi(s), d2(s, r), rmat(s, s), rinv(s, s))
      call mean_prep(y(obs0 + 1:obs0 + s), x(obs0 + 1:obs0 + s, :), offset(obs0 + 1:obs0 + s), &
         waves(obs0 + 1:obs0 + s), mean_links, variance_codes, beta, pr, d, mu, v, vmu, status)
      if (status /= GEE_OK) return
      call scale_prep(scale_offset(obs0 + 1:obs0 + s), zsca(obs0 + 1:obs0 + s, :), &
         waves(obs0 + 1:obs0 + s), scale_links, gamma, phi, d2)
      if (spec%scale_fixed) then
         phi = spec%scale_value
         d2 = 0.0_dp
      end if
      call cluster_correlation(spec, alpha, zcor, z0, waves(obs0 + 1:obs0 + s), &
            cor_param(obs0 + 1:obs0 + s), npair, rmat, status)
      if (status /= GEE_OK) return
      call inverse_checked(matmul(diagonal_matrix(sqrt(phi)), &
         matmul(rmat, diagonal_matrix(sqrt(phi)))), rinv, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      allocate(dw(s, p), prw(s))
      dw = d
      prw = pr
      call row_scale_matrix(dw, sqrt(weights(obs0 + 1:obs0 + s)))
      prw = prw * sqrt(weights(obs0 + 1:obs0 + s))
      h(1:p, 1:p) = matmul(transpose(dw), matmul(rinv, dw))
      score(1:p) = matmul(transpose(dw), matmul(rinv, prw))

      if (.not. spec%scale_fixed .and. r > 0) then
         allocate(wv(s), sbeta(s, p))
         wv = weights(obs0 + 1:obs0 + s) / (2.0_dp * phi)
         h(p + 1:p + r, p + 1:p + r) = matmul(transpose(d2), row_scaled_copy(d2, wv))
         call get_s_beta(d, pr, v, vmu, sbeta)
         h(p + 1:p + r, 1:p) = -matmul(transpose(d2), row_scaled_copy(sbeta, wv))
         score(p + 1:p + r) = matmul(transpose(d2), wv * (pr ** 2 - phi))
      end if

      if (association_active(spec%corstr, q) .and. s > 1) then
         allocate(e(npair, q), zi(npair), rhoi(npair), spr(s), pairw(npair), zbeta(npair, p))
         call cluster_correlation_and_e(spec, alpha, zcor, z0, waves(obs0 + 1:obs0 + s), &
            cor_param(obs0 + 1:obs0 + s), npair, rmat, e, status)
         if (status /= GEE_OK) return
         spr = pr / sqrt(phi)
         call pair_products(spr, zi)
         call upper_triangle(rmat, rhoi)
         call pair_products(weights(obs0 + 1:obs0 + s), pairw)
         h(p + r + 1:l, p + r + 1:l) = matmul(transpose(e), row_scaled_copy(e, pairw))
         call get_z_beta(d, pr, v, vmu, phi, zbeta)
         h(p + r + 1:l, 1:p) = -matmul(transpose(e), row_scaled_copy(zbeta, pairw))
         score(p + r + 1:l) = matmul(transpose(e), pairw * (zi - rhoi))
         if (.not. spec%scale_fixed .and. r > 0) then
            allocate(zgamma(npair, r))
            call get_z_gamma(d2, pr, phi, zi, zgamma)
            h(p + r + 1:l, p + 1:p + r) = -matmul(transpose(e), row_scaled_copy(zgamma, pairw))
         end if
      end if
      status = GEE_OK
   end subroutine cluster_hg

   subroutine mean_prep(y, x, offset, waves, mean_links, variance_codes, beta, pr, d, mu, v, vmu, status)
      real(dp), intent(in) :: y(:) !! Cluster responses.
      real(dp), intent(in) :: x(:, :) !! Cluster mean design matrix.
      real(dp), intent(in) :: offset(:) !! Cluster mean offsets.
      integer, intent(in) :: waves(:) !! Cluster wave indices.
      integer, intent(in) :: mean_links(:) !! Mean-link codes by wave.
      integer, intent(in) :: variance_codes(:) !! Variance codes by wave.
      real(dp), intent(in) :: beta(:) !! Mean regression coefficients.
      real(dp), intent(out) :: pr(:) !! Standardized Pearson residuals excluding scale.
      real(dp), intent(out) :: d(:, :) !! Standardized derivative matrix dmu/dbeta divided by sqrt(V(mu)).
      real(dp), intent(out) :: mu(:) !! Fitted means.
      real(dp), intent(out) :: v(:) !! Variance-function values.
      real(dp), intent(out) :: vmu(:) !! Derivatives dV/dmu.
      integer, intent(out) :: status !! GEE_OK or GEE_ERR_INVALID_MEAN.
      real(dp) :: eta
      real(dp) :: dm
      integer :: i
      integer :: link_code
      integer :: var_code

      do i = 1, size(y)
         link_code = mean_links(waves(i))
         var_code = variance_codes(waves(i))
         eta = dot_product(x(i, :), beta) + offset(i)
         mu(i) = link_inverse(eta, link_code)
         if (.not. valid_mean(mu(i), var_code)) then
            status = GEE_ERR_INVALID_MEAN
            return
         end if
         v(i) = variance_function(mu(i), var_code)
         vmu(i) = variance_derivative(mu(i), var_code)
         if (v(i) <= 0.0_dp) then
            status = GEE_ERR_INVALID_MEAN
            return
         end if
         dm = link_derivative(eta, link_code)
         d(i, :) = dm * x(i, :) / sqrt(v(i))
         pr(i) = (y(i) - mu(i)) / sqrt(v(i))
      end do
      status = GEE_OK
   end subroutine mean_prep

   subroutine compute_pr(y, x, offset, waves, mean_links, variance_codes, beta, pr, status)
      real(dp), intent(in) :: y(:) !! Full response vector.
      real(dp), intent(in) :: x(:, :) !! Full mean design matrix.
      real(dp), intent(in) :: offset(:) !! Full mean offsets.
      integer, intent(in) :: waves(:) !! Full wave indices.
      integer, intent(in) :: mean_links(:) !! Mean-link codes by wave.
      integer, intent(in) :: variance_codes(:) !! Variance codes by wave.
      real(dp), intent(in) :: beta(:) !! Mean coefficients.
      real(dp), intent(out) :: pr(:) !! Standardized Pearson residuals excluding scale.
      integer, intent(out) :: status !! GEE_OK or GEE_ERR_INVALID_MEAN.
      real(dp) :: eta
      real(dp) :: mu
      real(dp) :: v
      integer :: i

      do i = 1, size(y)
         eta = dot_product(x(i, :), beta) + offset(i)
         mu = link_inverse(eta, mean_links(waves(i)))
         if (.not. valid_mean(mu, variance_codes(waves(i)))) then
            status = GEE_ERR_INVALID_MEAN
            return
         end if
         v = variance_function(mu, variance_codes(waves(i)))
         if (v <= 0.0_dp) then
            status = GEE_ERR_INVALID_MEAN
            return
         end if
         pr(i) = (y(i) - mu) / sqrt(v)
      end do
      status = GEE_OK
   end subroutine compute_pr

   subroutine compute_phi(scale_offset, zsca, waves, scale_links, gamma, spec, phi)
      real(dp), intent(in) :: scale_offset(:) !! Scale offsets.
      real(dp), intent(in) :: zsca(:, :) !! Scale design matrix.
      integer, intent(in) :: waves(:) !! Wave indices.
      integer, intent(in) :: scale_links(:) !! Scale-link codes by wave.
      real(dp), intent(in) :: gamma(:) !! Scale coefficients.
      type(gee_spec), intent(in) :: spec !! Model specification including fixed-scale settings.
      real(dp), intent(out) :: phi(:) !! Scale values.
      integer :: i
      real(dp) :: eta

      if (spec%scale_fixed) then
         phi = spec%scale_value
         return
      end if
      do i = 1, size(phi)
         eta = dot_product(zsca(i, :), gamma) + scale_offset(i)
         phi(i) = max(epsilon(1.0_dp), link_inverse(eta, scale_links(waves(i))))
      end do
   end subroutine compute_phi

   subroutine scale_prep(scale_offset, zsca, waves, scale_links, gamma, phi, d2)
      real(dp), intent(in) :: scale_offset(:) !! Cluster scale offsets.
      real(dp), intent(in) :: zsca(:, :) !! Cluster scale design matrix.
      integer, intent(in) :: waves(:) !! Cluster wave indices.
      integer, intent(in) :: scale_links(:) !! Scale-link codes by wave.
      real(dp), intent(in) :: gamma(:) !! Scale coefficients.
      real(dp), intent(out) :: phi(:) !! Fitted scale values.
      real(dp), intent(out) :: d2(:, :) !! Scale derivative matrix dphi/dgamma.
      integer :: i
      real(dp) :: eta
      real(dp) :: dphi

      do i = 1, size(phi)
         eta = dot_product(zsca(i, :), gamma) + scale_offset(i)
         phi(i) = max(epsilon(1.0_dp), link_inverse(eta, scale_links(waves(i))))
         dphi = link_derivative(eta, scale_links(waves(i)))
         d2(i, :) = dphi * zsca(i, :)
      end do
   end subroutine scale_prep

   subroutine cluster_correlation(spec, alpha, zcor, z0, wave, cor_param, npair, rmat, status)
      type(gee_spec), intent(in) :: spec !! Model specification including correlation structure and link.
      real(dp), intent(in) :: alpha(:) !! Current association coefficients.
      real(dp), intent(in) :: zcor(:, :) !! Association design matrix.
      integer, intent(in) :: z0 !! Number of association-design rows before this cluster.
      integer, intent(in) :: wave(:) !! Cluster wave indices used for link and unstructured-pair indexing.
      real(dp), intent(in) :: cor_param(:) !! Real-valued correlation coordinates; AR(1) uses their pairwise distances.
      integer, intent(in) :: npair !! Number of strict upper-triangle pairs in the cluster.
      real(dp), intent(out) :: rmat(:, :) !! Working correlation matrix for the cluster.
      integer, intent(out) :: status !! GEE_OK or an error code.
      real(dp), allocatable :: rho(:)
      real(dp) :: eta
      integer :: k

      select case (spec%corstr)
      case (COR_INDEPENDENCE)
         allocate(rho(0))
         call working_correlation(rho, wave, COR_INDEPENDENCE, rmat, status)
      case (COR_EXCHANGEABLE)
         allocate(rho(1))
         eta = dot_product(zcor(z0 + 1, :), alpha)
         rho(1) = link_inverse(eta, spec%corr_link)
         call working_correlation(rho, wave, COR_EXCHANGEABLE, rmat, status)
      case (COR_AR1)
         allocate(rho(1))
         eta = dot_product(zcor(z0 + 1, :), alpha)
         rho(1) = link_inverse(eta, spec%corr_link)
         call working_ar1_real(rho(1), cor_param, rmat, status)
      case (COR_FIXED)
         allocate(rho(npair))
         do k = 1, npair
            rho(k) = link_inverse(dot_product(zcor(z0 + k, :), alpha), spec%corr_link)
         end do
         call working_correlation(rho, wave, COR_USERDEFINED, rmat, status)
      case (COR_UNSTRUCTURED, COR_USERDEFINED)
         allocate(rho(npair))
         do k = 1, npair
            rho(k) = link_inverse(dot_product(zcor(z0 + k, :), alpha), spec%corr_link)
         end do
         call working_correlation(rho, wave, COR_USERDEFINED, rmat, status)
      case default
         status = GEE_ERR_ARGUMENT
      end select
   end subroutine cluster_correlation

   subroutine cluster_correlation_and_e(spec, alpha, zcor, z0, wave, cor_param, npair, rmat, e, status)
      type(gee_spec), intent(in) :: spec !! Model specification including correlation structure and link.
      real(dp), intent(in) :: alpha(:) !! Current association coefficients.
      real(dp), intent(in) :: zcor(:, :) !! Association design matrix.
      integer, intent(in) :: z0 !! Number of association-design rows before this cluster.
      integer, intent(in) :: wave(:) !! Cluster wave indices used for link and unstructured-pair indexing.
      real(dp), intent(in) :: cor_param(:) !! Real-valued correlation coordinates; AR(1) uses their pairwise distances.
      integer, intent(in) :: npair !! Number of strict upper-triangle pairs in the cluster.
      real(dp), intent(out) :: rmat(:, :) !! Working correlation matrix.
      real(dp), intent(out) :: e(:, :) !! Derivative of pairwise correlations with respect to alpha.
      integer, intent(out) :: status !! GEE_OK or an error code.
      real(dp), allocatable :: rho(:)
      real(dp), allocatable :: drho(:)
      real(dp), allocatable :: base_deriv(:, :)
      real(dp) :: eta
      integer :: k

      e = 0.0_dp
      select case (spec%corstr)
      case (COR_EXCHANGEABLE, COR_AR1)
         allocate(rho(1), drho(1), base_deriv(npair, 1))
         eta = dot_product(zcor(z0 + 1, :), alpha)
         rho(1) = link_inverse(eta, spec%corr_link)
         drho(1) = link_derivative(eta, spec%corr_link)
         if (spec%corstr == COR_EXCHANGEABLE) then
            call working_correlation(rho, wave, COR_EXCHANGEABLE, rmat, status)
            base_deriv = 1.0_dp
         else
            call working_ar1_real(rho(1), cor_param, rmat, status)
            call ar1_pair_derivative(rho(1), cor_param, base_deriv(:, 1))
         end if
         if (status /= GEE_OK) return
         do k = 1, npair
            e(k, :) = base_deriv(k, 1) * drho(1) * zcor(z0 + 1, :)
         end do
      case (COR_UNSTRUCTURED, COR_USERDEFINED)
         allocate(rho(npair), drho(npair))
         do k = 1, npair
            eta = dot_product(zcor(z0 + k, :), alpha)
            rho(k) = link_inverse(eta, spec%corr_link)
            drho(k) = link_derivative(eta, spec%corr_link)
            e(k, :) = drho(k) * zcor(z0 + k, :)
         end do
         call working_correlation(rho, wave, COR_USERDEFINED, rmat, status)
      case default
         call cluster_correlation(spec, alpha, zcor, z0, wave, cor_param, npair, rmat, status)
      end select
   end subroutine cluster_correlation_and_e

   pure subroutine working_ar1_real(rho, coordinate, rmat, status)
      real(dp), intent(in) :: rho !! AR(1) base correlation parameter.
      real(dp), intent(in) :: coordinate(:) !! Real-valued correlation coordinates for observations in one cluster.
      real(dp), intent(out) :: rmat(:, :) !! AR(1) working correlation matrix.
      integer, intent(out) :: status !! GEE_OK or GEE_ERR_ARGUMENT/GEE_ERR_CORRELATION.
      real(dp) :: lag
      integer :: i
      integer :: j

      status = GEE_OK
      if (size(rmat, 1) /= size(coordinate) .or. size(rmat, 2) /= size(coordinate)) then
         status = GEE_ERR_ARGUMENT
         return
      end if
      rmat = 0.0_dp
      do i = 1, size(coordinate)
         rmat(i, i) = 1.0_dp
      end do
      do i = 1, size(coordinate) - 1
         do j = i + 1, size(coordinate)
            lag = abs(coordinate(j) - coordinate(i))
            if (rho < 0.0_dp .and. abs(lag - real(nint(lag), dp)) > 32.0_dp * epsilon(1.0_dp)) then
               status = GEE_ERR_CORRELATION
               return
            end if
            rmat(i, j) = rho ** lag
            rmat(j, i) = rmat(i, j)
         end do
      end do
      if (any(abs(rmat) > 1.0_dp + 100.0_dp * epsilon(1.0_dp))) status = GEE_ERR_CORRELATION
   end subroutine working_ar1_real

   pure subroutine ar1_pair_derivative(rho, wave, derivative)
      real(dp), intent(in) :: rho !! AR(1) correlation parameter.
      real(dp), intent(in) :: wave(:) !! Real-valued correlation coordinates for one cluster.
      real(dp), intent(out) :: derivative(:) !! Pair derivatives d rho^lag / d rho.
      integer :: i
      integer :: j
      integer :: k
      real(dp) :: lag

      k = 0
      do i = 1, size(wave) - 1
         do j = i + 1, size(wave)
            k = k + 1
            lag = abs(wave(j) - wave(i))
            if (lag <= 32.0_dp * epsilon(1.0_dp)) then
               derivative(k) = 0.0_dp
            else if (abs(lag - 1.0_dp) <= 32.0_dp * epsilon(1.0_dp)) then
               derivative(k) = 1.0_dp
            else if (abs(rho) <= tiny(1.0_dp)) then
               derivative(k) = 0.0_dp
            else
               derivative(k) = lag * rho ** (lag - 1.0_dp)
            end if
         end do
      end do
   end subroutine ar1_pair_derivative

   pure subroutine get_s_beta(d, pr, v, vmu, sbeta)
      real(dp), intent(in) :: d(:, :) !! Standardized mean derivative matrix.
      real(dp), intent(in) :: pr(:) !! Standardized Pearson residuals.
      real(dp), intent(in) :: v(:) !! Mean variance-function values.
      real(dp), intent(in) :: vmu(:) !! Derivatives of the variance function with respect to mean.
      real(dp), intent(out) :: sbeta(:, :) !! Derivative of squared Pearson residuals with respect to beta.
      real(dp) :: factor
      integer :: i

      do i = 1, size(pr)
         factor = -2.0_dp * pr(i) / sqrt(v(i)) - pr(i) * pr(i) * vmu(i) / v(i)
         sbeta(i, :) = factor * d(i, :)
      end do
   end subroutine get_s_beta

   pure subroutine get_z_beta(d, pr, v, vmu, phi, zbeta)
      real(dp), intent(in) :: d(:, :) !! Standardized mean derivative matrix.
      real(dp), intent(in) :: pr(:) !! Standardized Pearson residuals.
      real(dp), intent(in) :: v(:) !! Mean variance-function values.
      real(dp), intent(in) :: vmu(:) !! Derivatives dV/dmu.
      real(dp), intent(in) :: phi(:) !! Scale values.
      real(dp), intent(out) :: zbeta(:, :) !! Derivative of pairwise standardized residual products with respect to beta.
      integer :: i
      integer :: j
      integer :: k
      real(dp) :: scale

      k = 0
      do i = 1, size(pr) - 1
         do j = i + 1, size(pr)
            k = k + 1
            zbeta(k, :) = -pr(i) * d(j, :) - pr(j) * d(i, :) - 0.5_dp * pr(i) * pr(j) * &
               (vmu(i) * d(i, :) / sqrt(v(i)) + vmu(j) * d(j, :) / sqrt(v(j)))
            scale = 1.0_dp / sqrt(phi(i) * phi(j))
            zbeta(k, :) = scale * zbeta(k, :)
         end do
      end do
   end subroutine get_z_beta

   pure subroutine get_z_gamma(d2, pr, phi, zi, zgamma)
      real(dp), intent(in) :: d2(:, :) !! Scale derivative matrix.
      real(dp), intent(in) :: pr(:) !! Standardized Pearson residuals.
      real(dp), intent(in) :: phi(:) !! Scale values.
      real(dp), intent(in) :: zi(:) !! Pair products of scale-standardized residuals.
      real(dp), intent(out) :: zgamma(:, :) !! Derivative of pair residual products with respect to gamma.
      integer :: i
      integer :: j
      integer :: k

      k = 0
      do i = 1, size(pr) - 1
         do j = i + 1, size(pr)
            k = k + 1
            zgamma(k, :) = -0.5_dp * zi(k) * (d2(i, :) / phi(i) + d2(j, :) / phi(j))
         end do
      end do
   end subroutine get_z_gamma

   pure subroutine row_scale_matrix(a, scale)
      real(dp), intent(inout) :: a(:, :) !! Matrix whose rows are scaled in place.
      real(dp), intent(in) :: scale(:) !! Multiplicative factor for each matrix row.
      integer :: i

      do i = 1, size(a, 1)
         a(i, :) = scale(i) * a(i, :)
      end do
   end subroutine row_scale_matrix

   pure function row_scaled_copy(a, scale) result(b)
      real(dp), intent(in) :: a(:, :) !! Matrix to scale by rows.
      real(dp), intent(in) :: scale(:) !! Multiplicative factor for each row.
      real(dp) :: b(size(a, 1), size(a, 2))
      integer :: i

      do i = 1, size(a, 1)
         b(i, :) = scale(i) * a(i, :)
      end do
   end function row_scaled_copy

   pure logical function all_valid_fitted(x, offset, waves, mean_links, variance_codes, beta) result(ok)
      real(dp), intent(in) :: x(:, :) !! Mean design matrix.
      real(dp), intent(in) :: offset(:) !! Mean offsets.
      integer, intent(in) :: waves(:) !! Wave indices.
      integer, intent(in) :: mean_links(:) !! Mean-link codes by wave.
      integer, intent(in) :: variance_codes(:) !! Variance codes by wave.
      real(dp), intent(in) :: beta(:) !! Candidate mean coefficients.
      real(dp) :: mu
      integer :: i

      ok = .true.
      do i = 1, size(offset)
         mu = link_inverse(dot_product(x(i, :), beta) + offset(i), mean_links(waves(i)))
         if (.not. valid_mean(mu, variance_codes(waves(i)))) then
            ok = .false.
            return
         end if
      end do
   end function all_valid_fitted

   pure integer function association_rows(corstr, cluster_size) result(nrow)
      integer, intent(in) :: corstr !! Correlation structure identifier.
      integer, intent(in) :: cluster_size !! Number of observations in the cluster.

      select case (corstr)
      case (COR_INDEPENDENCE)
         nrow = 0
      case (COR_EXCHANGEABLE, COR_AR1)
         nrow = 1
      case default
         nrow = cluster_size * (cluster_size - 1) / 2
      end select
   end function association_rows

   pure logical function association_active(corstr, q) result(active)
      integer, intent(in) :: corstr !! Correlation structure identifier.
      integer, intent(in) :: q !! Number of association coefficients.

      active = q > 0 .and. corstr /= COR_INDEPENDENCE .and. corstr /= COR_FIXED
   end function association_active

   subroutine regularize_inactive_blocks(h, p, r, q, spec)
      real(dp), intent(inout) :: h(:, :) !! Full lower-block sensitivity matrix.
      integer, intent(in) :: p !! Number of mean coefficients.
      integer, intent(in) :: r !! Number of scale coefficients.
      integer, intent(in) :: q !! Number of association coefficients.
      type(gee_spec), intent(in) :: spec !! Model specification determining active parameter blocks.
      integer :: i

      if (spec%scale_fixed) then
         do i = 1, r
            h(p + i, p + i) = 1.0_dp
         end do
      end if
      if (.not. association_active(spec%corstr, q)) then
         do i = 1, q
            h(p + r + i, p + r + i) = 1.0_dp
         end do
      end if
   end subroutine regularize_inactive_blocks

   subroutine approximate_jackknife(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
      mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, result, status)
      real(dp), intent(in) :: y(:) !! Response vector.
      real(dp), intent(in) :: x(:, :) !! Mean design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      real(dp), intent(in) :: offset(:) !! Mean offsets.
      real(dp), intent(in) :: scale_offset(:) !! Scale offsets.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: zsca(:, :) !! Scale design matrix.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix.
      real(dp), intent(in) :: cor_param(:) !! Correlation metadata retained for parity.
      integer, intent(in) :: mean_links(:) !! Mean-link codes.
      integer, intent(in) :: variance_codes(:) !! Variance codes.
      integer, intent(in) :: scale_links(:) !! Scale-link codes.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(in) :: beta(:) !! Fitted mean coefficients.
      real(dp), intent(in) :: gamma(:) !! Fitted scale coefficients.
      real(dp), intent(in) :: alpha(:) !! Fitted association coefficients.
      type(gee_result), intent(inout) :: result !! Result receiving AJS covariance blocks.
      integer, intent(out) :: status !! GEE_OK or a numerical error code.
      real(dp), allocatable :: htotal(:, :)
      real(dp), allocatable :: hi(:, :)
      real(dp), allocatable :: score(:)
      real(dp), allocatable :: invomit(:, :)
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: cov(:, :)
      integer :: p
      integer :: r
      integer :: q
      integer :: l
      integer :: gidx
      integer :: info
      real(dp) :: factor

      p = size(beta)
      r = size(gamma)
      q = size(alpha)
      l = p + r + q
      allocate(htotal(l, l), hi(l, l), score(l), invomit(l, l), delta(l), cov(l, l))
      htotal = 0.0_dp
      do gidx = 1, size(cluster_sizes)
         call cluster_hg(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
            mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, gidx, hi, score, status)
         if (status /= GEE_OK) return
         htotal = htotal + hi
      end do
      call regularize_inactive_blocks(htotal, p, r, q, spec)
      cov = 0.0_dp
      do gidx = 1, size(cluster_sizes)
         call cluster_hg(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
            mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, gidx, hi, score, status)
         if (status /= GEE_OK) return
         call regularize_omitted_inactive(htotal - hi, p, r, q, spec, invomit, info)
         if (info /= 0) then
            status = GEE_ERR_SINGULAR
            return
         end if
         delta = matmul(invomit, score)
         cov = cov + outer_product(delta, delta)
      end do
      factor = real(size(cluster_sizes) - p - q - r, dp) / real(size(cluster_sizes), dp)
      cov = factor * cov
      result%vbeta_ajs = cov(1:p, 1:p)
      if (r > 0) result%vgamma_ajs = cov(p + 1:p + r, p + 1:p + r)
      if (q > 0) result%valpha_ajs = cov(p + r + 1:l, p + r + 1:l)
      status = GEE_OK
   end subroutine approximate_jackknife

   subroutine regularize_omitted_inactive(h, p, r, q, spec, hinv, status)
      real(dp), intent(in) :: h(:, :) !! Omitted-cluster sensitivity matrix.
      integer, intent(in) :: p !! Number of mean coefficients.
      integer, intent(in) :: r !! Number of scale coefficients.
      integer, intent(in) :: q !! Number of association coefficients.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(out) :: hinv(:, :) !! Inverse after inactive-block regularization.
      integer, intent(out) :: status !! Zero on successful inversion.
      real(dp), allocatable :: work(:, :)

      allocate(work(size(h, 1), size(h, 2)))
      work = h
      call regularize_inactive_blocks(work, p, r, q, spec)
      call inverse_checked(work, hinv, status)
   end subroutine regularize_omitted_inactive

   subroutine refit_jackknife(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
      mean_links, variance_codes, scale_links, spec, beta, gamma, alpha, max_iterations, vbeta, vgamma, valpha, status)
      real(dp), intent(in) :: y(:) !! Response vector.
      real(dp), intent(in) :: x(:, :) !! Mean design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes.
      real(dp), intent(in) :: offset(:) !! Mean offsets.
      real(dp), intent(in) :: scale_offset(:) !! Scale offsets.
      real(dp), intent(in) :: weights(:) !! Observation weights.
      integer, intent(in) :: waves(:) !! Wave indices.
      real(dp), intent(in) :: zsca(:, :) !! Scale design matrix.
      real(dp), intent(in) :: zcor(:, :) !! Correlation design matrix.
      real(dp), intent(in) :: cor_param(:) !! Correlation metadata retained for parity.
      integer, intent(in) :: mean_links(:) !! Mean-link codes.
      integer, intent(in) :: variance_codes(:) !! Variance codes.
      integer, intent(in) :: scale_links(:) !! Scale-link codes.
      type(gee_spec), intent(in) :: spec !! Model specification.
      real(dp), intent(in) :: beta(:) !! Full-data fitted mean coefficients.
      real(dp), intent(in) :: gamma(:) !! Full-data fitted scale coefficients.
      real(dp), intent(in) :: alpha(:) !! Full-data fitted association coefficients.
      integer, intent(in) :: max_iterations !! Iteration count for each leave-one-cluster refit.
      real(dp), intent(out) :: vbeta(:, :) !! Jackknife covariance for mean coefficients.
      real(dp), intent(out) :: vgamma(:, :) !! Jackknife covariance for scale coefficients.
      real(dp), intent(out) :: valpha(:, :) !! Jackknife covariance for association coefficients.
      integer, intent(out) :: status !! GEE_OK or a numerical error code.
      real(dp), allocatable :: b(:)
      real(dp), allocatable :: gm(:)
      real(dp), allocatable :: a(:)
      real(dp) :: factor
      integer :: gidx
      integer :: iterations
      integer :: fit_status

      vbeta = 0.0_dp
      vgamma = 0.0_dp
      valpha = 0.0_dp
      do gidx = 1, size(cluster_sizes)
         allocate(b(size(beta)), gm(size(gamma)), a(size(alpha)))
         b = beta
         gm = gamma
         a = alpha
         call estimate_parameters(y, x, cluster_sizes, offset, scale_offset, weights, waves, zsca, zcor, cor_param, &
            mean_links, variance_codes, scale_links, spec, b, gm, a, gidx, max_iterations, iterations, fit_status)
         if (fit_status /= GEE_OK .and. fit_status /= GEE_ERR_MAXITER) then
            status = fit_status
            return
         end if
         vbeta = vbeta + outer_product(b - beta, b - beta)
         if (size(gamma) > 0) vgamma = vgamma + outer_product(gm - gamma, gm - gamma)
         if (size(alpha) > 0) valpha = valpha + outer_product(a - alpha, a - alpha)
         deallocate(b, gm, a)
      end do
      factor = real(size(cluster_sizes) - size(beta) - size(gamma) - size(alpha), dp) / real(size(cluster_sizes), dp)
      vbeta = factor * vbeta
      vgamma = factor * vgamma
      valpha = factor * valpha
      status = GEE_OK
   end subroutine refit_jackknife

end module geepack_gee
