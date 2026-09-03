! Modern Fortran translation of the computational core of R/ordgee.R and src/ordgee.cc.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_ordinal
   use r_kinds, only : dp
   use geepack_status, only : GEE_OK, GEE_ERR_SHAPE, GEE_ERR_ARGUMENT, GEE_ERR_SINGULAR, &
      GEE_ERR_MAXITER, GEE_ERR_INVALID_MEAN
   use geepack_links, only : LINK_LOGIT, link_function, link_inverse, link_derivative
   use geepack_correlations, only : COR_INDEPENDENCE, COR_EXCHANGEABLE, COR_UNSTRUCTURED, &
      COR_USERDEFINED
   use geepack_design, only : gen_zodds
   use geepack_matrix, only : inverse_checked, solve_checked, outer_product, max_abs
   implicit none
   private

   type, public :: ordinal_spec
      integer :: mean_link = LINK_LOGIT
      integer :: corstr = COR_INDEPENDENCE
      logical :: constant_intercepts = .true.
      logical :: reverse_coding = .false.
      real(dp) :: tolerance = 1.0e-4_dp
      integer :: max_iterations = 25
   end type ordinal_spec

   type, public :: ordinal_result
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: vbeta(:, :)
      real(dp), allocatable :: valpha(:, :)
      real(dp), allocatable :: vbeta_naive(:, :)
      real(dp), allocatable :: valpha_naive(:, :)
      real(dp), allocatable :: valpha_stable(:, :)
      real(dp), allocatable :: influence(:, :)
      integer :: iterations = 0
      integer :: error = GEE_OK
   end type ordinal_result

   public :: fit_ordgee
   public :: odds_to_p11, p11_odds_derivative, p11_mean_derivatives

contains

   subroutine fit_ordgee(y, x, cluster_sizes, nlevels, spec, result, beta_initial, alpha_initial, &
      offset, weights, waves, z, odds_offset)
      integer, intent(in) :: y(:) !! Ordered response levels coded 1,...,nlevels in contiguous cluster order.
      real(dp), intent(in) :: x(:, :) !! Covariate matrix excluding the threshold intercept columns.
      integer, intent(in) :: cluster_sizes(:) !! Numbers of ordinal observations in contiguous clusters.
      integer, intent(in) :: nlevels !! Number of ordered response levels; must be at least two.
      type(ordinal_spec), intent(in) :: spec !! Mean link, association structure, coding, and iteration controls.
      type(ordinal_result), intent(out) :: result !! Fitted coefficients, covariance estimates, influences, and status.
      real(dp), optional, intent(in) :: beta_initial(:) !! Initial threshold and covariate coefficients.
      real(dp), optional, intent(in) :: alpha_initial(:) !! Initial log-odds-ratio association coefficients.
      real(dp), optional, intent(in) :: offset(:) !! Original observation-level mean-model offset, subtracted as in ordgee.
      real(dp), optional, intent(in) :: weights(:) !! Observation-level fitting weights; defaults to one.
      integer, optional, intent(in) :: waves(:) !! One-based within-cluster wave identifiers; defaults to observation order.
      real(dp), optional, intent(in) :: z(:, :) !! Odds-ratio design matrix in pair-major, threshold-pair-minor order.
      real(dp), optional, intent(in) :: odds_offset(:) !! Log-odds-ratio offsets, one per row of z.
      real(dp), allocatable :: ye(:)
      real(dp), allocatable :: xe(:, :)
      real(dp), allocatable :: offe(:)
      real(dp), allocatable :: we(:)
      real(dp), allocatable :: zwork(:, :)
      real(dp), allocatable :: ooff(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: alpha(:)
      integer, allocatable :: wavework(:)
      integer :: n
      integer :: c
      integer :: p
      integer :: q
      integer :: status
      integer :: iterations

      call initialize_ordinal_result(result, 0, 0, size(cluster_sizes))
      n = size(y)
      c = nlevels - 1
      if (c < 1 .or. size(x, 1) /= n .or. sum(cluster_sizes) /= n .or. any(cluster_sizes < 1)) then
         result%error = GEE_ERR_SHAPE
         return
      end if
      if (any(y < 1) .or. any(y > nlevels)) then
         result%error = GEE_ERR_ARGUMENT
         return
      end if
      if (.not. valid_ordinal_link(spec%mean_link) .or. .not. valid_ordinal_corstr(spec%corstr)) then
         result%error = GEE_ERR_ARGUMENT
         return
      end if

      call setup_ordinal_waves(cluster_sizes, waves, wavework, status)
      if (status /= GEE_OK) then
         result%error = status
         return
      end if
      call expand_ordinal_data(y, x, c, spec, wavework, offset, weights, ye, xe, offe, we, status)
      if (status /= GEE_OK) then
         result%error = status
         return
      end if
      p = size(xe, 2)

      if (present(z)) then
         allocate(zwork(size(z, 1), size(z, 2)))
         zwork = z
      else
         call gen_zodds(cluster_sizes, wavework, spec%corstr, c, zwork, status)
         if (status /= GEE_OK) then
            result%error = status
            return
         end if
      end if
      q = size(zwork, 2)
      if (size(zwork, 1) /= expected_odds_rows(cluster_sizes, c, spec%corstr)) then
         result%error = GEE_ERR_SHAPE
         return
      end if

      allocate(ooff(size(zwork, 1)))
      if (present(odds_offset)) then
         if (size(odds_offset) /= size(ooff)) then
            result%error = GEE_ERR_SHAPE
            return
         end if
         ooff = odds_offset
      else
         ooff = 0.0_dp
      end if

      allocate(beta(p))
      if (present(beta_initial)) then
         if (size(beta_initial) /= p) then
            result%error = GEE_ERR_SHAPE
            return
         end if
         beta = beta_initial
      else
         call ordinal_starting_beta(y, x, c, wavework, spec, beta)
      end if
      allocate(alpha(q))
      if (present(alpha_initial)) then
         if (size(alpha_initial) /= q) then
            result%error = GEE_ERR_SHAPE
            return
         end if
         alpha = alpha_initial
      else
         alpha = 0.0_dp
      end if

      call ordinal_estimate(ye, xe, offe, we, cluster_sizes, c, spec, zwork, ooff, beta, alpha, &
         iterations, status)
      call initialize_ordinal_result(result, p, q, size(cluster_sizes))
      result%beta = beta
      result%alpha = alpha
      result%iterations = iterations
      result%error = status
      if (status /= GEE_OK .and. status /= GEE_ERR_MAXITER) return

      call ordinal_covariance(ye, xe, offe, cluster_sizes, c, spec, zwork, ooff, beta, alpha, result, status)
      if (status /= GEE_OK) result%error = status
   end subroutine fit_ordgee

   subroutine initialize_ordinal_result(result, p, q, ncluster)
      type(ordinal_result), intent(out) :: result !! Result object to allocate and initialize.
      integer, intent(in) :: p !! Number of mean parameters.
      integer, intent(in) :: q !! Number of association parameters.
      integer, intent(in) :: ncluster !! Number of independent clusters.
      integer :: l

      l = p + q
      allocate(result%beta(p), result%alpha(q))
      allocate(result%vbeta(p, p), result%valpha(q, q))
      allocate(result%vbeta_naive(p, p), result%valpha_naive(q, q))
      allocate(result%valpha_stable(q, q), result%influence(l, ncluster))
      result%beta = 0.0_dp
      result%alpha = 0.0_dp
      result%vbeta = 0.0_dp
      result%valpha = 0.0_dp
      result%vbeta_naive = 0.0_dp
      result%valpha_naive = 0.0_dp
      result%valpha_stable = 0.0_dp
      result%influence = 0.0_dp
      result%iterations = 0
      result%error = GEE_OK
   end subroutine initialize_ordinal_result

   subroutine setup_ordinal_waves(cluster_sizes, waves, wavework, status)
      integer, intent(in) :: cluster_sizes(:) !! Cluster sizes in data order.
      integer, optional, intent(in) :: waves(:) !! Optional one-based wave identifiers.
      integer, allocatable, intent(out) :: wavework(:) !! Validated or generated wave identifiers.
      integer, intent(out) :: status !! GEE_OK on success or an argument error.
      integer :: pos
      integer :: g
      integer :: j
      integer :: n

      n = sum(cluster_sizes)
      allocate(wavework(n))
      if (present(waves)) then
         if (size(waves) /= n .or. any(waves < 1)) then
            status = GEE_ERR_SHAPE
            return
         end if
         wavework = waves
      else
         pos = 0
         do g = 1, size(cluster_sizes)
            do j = 1, cluster_sizes(g)
               wavework(pos + j) = j
            end do
            pos = pos + cluster_sizes(g)
         end do
      end if
      status = GEE_OK
   end subroutine setup_ordinal_waves

   subroutine expand_ordinal_data(y, x, c, spec, waves, offset, weights, ye, xe, offe, we, status)
      integer, intent(in) :: y(:) !! Ordered response levels.
      real(dp), intent(in) :: x(:, :) !! Base covariate matrix without threshold intercepts.
      integer, intent(in) :: c !! Number of cumulative indicators, nlevels minus one.
      type(ordinal_spec), intent(in) :: spec !! Ordinal coding and intercept specification.
      integer, intent(in) :: waves(:) !! One-based wave identifiers.
      real(dp), optional, intent(in) :: offset(:) !! Original observation-level mean offset.
      real(dp), optional, intent(in) :: weights(:) !! Observation weights.
      real(dp), allocatable, intent(out) :: ye(:) !! Expanded binary cumulative responses.
      real(dp), allocatable, intent(out) :: xe(:, :) !! Expanded threshold-plus-covariate design matrix.
      real(dp), allocatable, intent(out) :: offe(:) !! Expanded mean offsets with ordgee sign convention.
      real(dp), allocatable, intent(out) :: we(:) !! Expanded fitting weights.
      integer, intent(out) :: status !! GEE_OK on success or a shape error.
      real(dp), allocatable :: off(:)
      real(dp), allocatable :: wt(:)
      integer :: n
      integer :: pc
      integer :: nthresh
      integer :: i
      integer :: k
      integer :: row
      integer :: col

      n = size(y)
      pc = size(x, 2)
      if (size(waves) /= n) then
         status = GEE_ERR_SHAPE
         allocate(ye(0), xe(0, 0), offe(0), we(0))
         return
      end if
      allocate(off(n), wt(n))
      off = 0.0_dp
      wt = 1.0_dp
      if (present(offset)) then
         if (size(offset) /= n) then
            status = GEE_ERR_SHAPE
            allocate(ye(0), xe(0, 0), offe(0), we(0))
            return
         end if
         off = offset
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            status = GEE_ERR_SHAPE
            allocate(ye(0), xe(0, 0), offe(0), we(0))
            return
         end if
         wt = weights
      end if
      if (spec%constant_intercepts) then
         nthresh = c
      else
         nthresh = maxval(waves) * c
      end if
      allocate(ye(n * c), xe(n * c, nthresh + pc), offe(n * c), we(n * c))
      xe = 0.0_dp
      do i = 1, n
         do k = 1, c
            row = (i - 1) * c + k
            if (spec%reverse_coding) then
               ye(row) = merge(1.0_dp, 0.0_dp, y(i) <= k)
            else
               ye(row) = merge(1.0_dp, 0.0_dp, y(i) > k)
            end if
            if (spec%constant_intercepts) then
               col = k
            else
               col = (waves(i) - 1) * c + k
            end if
            xe(row, col) = 1.0_dp
            if (pc > 0) xe(row, nthresh + 1:nthresh + pc) = x(i, :)
            offe(row) = -off(i)
            we(row) = wt(i)
         end do
      end do
      status = GEE_OK
   end subroutine expand_ordinal_data

   subroutine ordinal_starting_beta(y, x, c, waves, spec, beta)
      integer, intent(in) :: y(:) !! Ordered response levels.
      real(dp), intent(in) :: x(:, :) !! Base covariates, used only to determine coefficient count.
      integer, intent(in) :: c !! Number of cumulative indicators.
      integer, intent(in) :: waves(:) !! Wave identifiers.
      type(ordinal_spec), intent(in) :: spec !! Link and intercept specification.
      real(dp), intent(out) :: beta(:) !! Starting coefficient vector.
      real(dp) :: prob
      integer :: nthresh
      integer :: k
      integer :: w
      integer :: idx
      integer :: countw

      beta = 0.0_dp
      if (spec%constant_intercepts) then
         nthresh = c
         do k = 1, c
            if (spec%reverse_coding) then
               prob = real(count(y <= k), dp) / real(size(y), dp)
            else
               prob = real(count(y > k), dp) / real(size(y), dp)
            end if
            prob = min(1.0_dp - 1.0e-5_dp, max(1.0e-5_dp, prob))
            beta(k) = link_function(prob, spec%mean_link)
         end do
      else
         nthresh = maxval(waves) * c
         do w = 1, maxval(waves)
            countw = count(waves == w)
            do k = 1, c
               idx = (w - 1) * c + k
               if (countw == 0) cycle
               if (spec%reverse_coding) then
                  prob = real(count((waves == w) .and. (y <= k)), dp) / real(countw, dp)
               else
                  prob = real(count((waves == w) .and. (y > k)), dp) / real(countw, dp)
               end if
               prob = min(1.0_dp - 1.0e-5_dp, max(1.0e-5_dp, prob))
               beta(idx) = link_function(prob, spec%mean_link)
            end do
         end do
      end if
      if (size(beta) > nthresh + size(x, 2)) beta(nthresh + size(x, 2) + 1:) = 0.0_dp
   end subroutine ordinal_starting_beta

   subroutine ordinal_estimate(y, x, offset, weights, cluster_sizes, c, spec, z, odds_offset, beta, alpha, &
      iterations, status)
      real(dp), intent(in) :: y(:) !! Expanded cumulative binary responses.
      real(dp), intent(in) :: x(:, :) !! Expanded mean design matrix.
      real(dp), intent(in) :: offset(:) !! Expanded mean offsets.
      real(dp), intent(in) :: weights(:) !! Expanded fitting weights.
      integer, intent(in) :: cluster_sizes(:) !! Original ordinal cluster sizes.
      integer, intent(in) :: c !! Number of cumulative indicators per observation.
      type(ordinal_spec), intent(in) :: spec !! Ordinal fit controls and structures.
      real(dp), intent(in) :: z(:, :) !! Odds-ratio design matrix.
      real(dp), intent(in) :: odds_offset(:) !! Log-odds-ratio offsets.
      real(dp), intent(inout) :: beta(:) !! Mean coefficients updated in place.
      real(dp), intent(inout) :: alpha(:) !! Association coefficients updated in place.
      integer, intent(out) :: iterations !! Number of completed GEE iterations.
      integer, intent(out) :: status !! GEE_OK, GEE_ERR_MAXITER, or numerical error.
      real(dp) :: dbeta
      real(dp) :: dalpha
      integer :: iter

      status = GEE_OK
      iterations = 0
      do iter = 1, spec%max_iterations
         call ordinal_update_beta(y, x, offset, weights, cluster_sizes, c, spec, z, odds_offset, &
            beta, alpha, dbeta, status)
         if (status /= GEE_OK) return
         call ordinal_update_alpha(y, x, offset, cluster_sizes, c, spec, z, odds_offset, &
            beta, alpha, dalpha, status)
         if (status /= GEE_OK) return
         iterations = iter
         if (max(dbeta, dalpha) <= spec%tolerance) return
      end do
      status = GEE_ERR_MAXITER
   end subroutine ordinal_estimate

   subroutine ordinal_update_beta(y, x, offset, weights, cluster_sizes, c, spec, z, odds_offset, &
      beta, alpha, delta, status)
      real(dp), intent(in) :: y(:) !! Expanded cumulative binary responses.
      real(dp), intent(in) :: x(:, :) !! Expanded mean design matrix.
      real(dp), intent(in) :: offset(:) !! Expanded mean offsets.
      real(dp), intent(in) :: weights(:) !! Expanded fitting weights.
      integer, intent(in) :: cluster_sizes(:) !! Original ordinal cluster sizes.
      integer, intent(in) :: c !! Number of cumulative indicators per observation.
      type(ordinal_spec), intent(in) :: spec !! Ordinal model specification.
      real(dp), intent(in) :: z(:, :) !! Odds-ratio design matrix.
      real(dp), intent(in) :: odds_offset(:) !! Log-odds-ratio offsets.
      real(dp), intent(inout) :: beta(:) !! Mean coefficients updated by one Fisher-scoring step.
      real(dp), intent(in) :: alpha(:) !! Current log-odds-ratio coefficients.
      real(dp), intent(out) :: delta !! Maximum absolute mean-coefficient update.
      integer, intent(out) :: status !! GEE_OK or a numerical error.
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: gvec(:)
      real(dp), allocatable :: step(:)
      real(dp), allocatable :: pr(:)
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: vinv(:, :)
      real(dp), allocatable :: dw(:, :)
      real(dp), allocatable :: prw(:)
      integer :: p
      integer :: gi
      integer :: s
      integer :: obs0
      integer :: z0
      integer :: m
      integer :: zr
      integer :: info

      p = size(beta)
      allocate(h(p, p), gvec(p), step(p))
      h = 0.0_dp
      gvec = 0.0_dp
      obs0 = 0
      z0 = 0
      do gi = 1, size(cluster_sizes)
         s = cluster_sizes(gi)
         m = s * c
         zr = cluster_odds_rows(s, c, spec%corstr)
         allocate(pr(m), d(m, p), mu(m), v(m, m), vinv(m, m), dw(m, p), prw(m))
         call ordinal_cluster_beta(y(obs0 * c + 1:(obs0 + s) * c), &
            x(obs0 * c + 1:(obs0 + s) * c, :), offset(obs0 * c + 1:(obs0 + s) * c), &
            s, c, spec, z_segment(z, z0, zr), offset_segment(odds_offset, z0, zr), beta, alpha, &
            pr, d, mu, v, status)
         if (status /= GEE_OK) return
         call inverse_checked(v, vinv, info)
         if (info /= 0) then
            status = GEE_ERR_SINGULAR
            return
         end if
         dw = d
         prw = pr
         call apply_root_weights(dw, prw, weights(obs0 * c + 1:(obs0 + s) * c))
         h = h + matmul(transpose(dw), matmul(vinv, dw))
         gvec = gvec + matmul(transpose(dw), matmul(vinv, prw))
         deallocate(pr, d, mu, v, vinv, dw, prw)
         obs0 = obs0 + s
         z0 = z0 + zr
      end do
      call solve_checked(h, gvec, step, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      beta = beta + step
      delta = max_abs(step)
      status = GEE_OK
   end subroutine ordinal_update_beta

   subroutine ordinal_update_alpha(y, x, offset, cluster_sizes, c, spec, z, odds_offset, beta, alpha, delta, status)
      real(dp), intent(in) :: y(:) !! Expanded cumulative binary responses.
      real(dp), intent(in) :: x(:, :) !! Expanded mean design matrix.
      real(dp), intent(in) :: offset(:) !! Expanded mean offsets.
      integer, intent(in) :: cluster_sizes(:) !! Original ordinal cluster sizes.
      integer, intent(in) :: c !! Number of cumulative indicators per observation.
      type(ordinal_spec), intent(in) :: spec !! Ordinal model specification.
      real(dp), intent(in) :: z(:, :) !! Odds-ratio design matrix.
      real(dp), intent(in) :: odds_offset(:) !! Log-odds-ratio offsets.
      real(dp), intent(in) :: beta(:) !! Current mean coefficients.
      real(dp), intent(inout) :: alpha(:) !! Association coefficients updated by one scoring step.
      real(dp), intent(out) :: delta !! Maximum absolute association-coefficient update.
      integer, intent(out) :: status !! GEE_OK or a numerical error.
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: gvec(:)
      real(dp), allocatable :: step(:)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: pr(:)
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: u(:)
      real(dp), allocatable :: v2(:, :)
      real(dp), allocatable :: d2(:, :)
      real(dp), allocatable :: vinv(:, :)
      integer :: q
      integer :: gi
      integer :: s
      integer :: obs0
      integer :: z0
      integer :: zr
      integer :: j
      integer :: k
      integer :: pair0
      integer :: info

      q = size(alpha)
      if (q == 0 .or. spec%corstr == COR_INDEPENDENCE) then
         delta = 0.0_dp
         status = GEE_OK
         return
      end if
      allocate(h(q, q), gvec(q), step(q))
      h = 0.0_dp
      gvec = 0.0_dp
      obs0 = 0
      z0 = 0
      do gi = 1, size(cluster_sizes)
         s = cluster_sizes(gi)
         zr = cluster_odds_rows(s, c, spec%corstr)
         allocate(mu(s * c), pr(s * c), d(s * c, size(beta)))
         call ordinal_mean_only(y(obs0 * c + 1:(obs0 + s) * c), x(obs0 * c + 1:(obs0 + s) * c, :), &
            offset(obs0 * c + 1:(obs0 + s) * c), spec%mean_link, beta, mu, pr, d, status)
         if (status /= GEE_OK) return
         pair0 = 0
         do j = 1, s - 1
            do k = j + 1, s
               allocate(u(c * c), v2(c * c, c * c), d2(c * c, q), vinv(c * c, c * c))
               call ordinal_pair_alpha(pr((j - 1) * c + 1:j * c), pr((k - 1) * c + 1:k * c), &
                  mu((j - 1) * c + 1:j * c), mu((k - 1) * c + 1:k * c), &
                  z(z0 + pair0 * c * c + 1:z0 + (pair0 + 1) * c * c, :), &
                  odds_offset(z0 + pair0 * c * c + 1:z0 + (pair0 + 1) * c * c), &
                  spec%reverse_coding, alpha, u, v2, d2)
               call inverse_checked(v2, vinv, info)
               if (info /= 0) then
                  status = GEE_ERR_SINGULAR
                  return
               end if
               h = h + matmul(transpose(d2), matmul(vinv, d2))
               gvec = gvec + matmul(transpose(d2), matmul(vinv, u))
               deallocate(u, v2, d2, vinv)
               pair0 = pair0 + 1
            end do
         end do
         deallocate(mu, pr, d)
         obs0 = obs0 + s
         z0 = z0 + zr
      end do
      call solve_checked(h, gvec, step, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      alpha = alpha + step
      delta = max_abs(step)
      status = GEE_OK
   end subroutine ordinal_update_alpha

   subroutine ordinal_cluster_beta(y, x, offset, s, c, spec, z, odds_offset, beta, alpha, pr, d, mu, v, status)
      real(dp), intent(in) :: y(:) !! Expanded responses for one cluster.
      real(dp), intent(in) :: x(:, :) !! Expanded mean design for one cluster.
      real(dp), intent(in) :: offset(:) !! Expanded mean offsets for one cluster.
      integer, intent(in) :: s !! Number of ordinal observations in the cluster.
      integer, intent(in) :: c !! Number of cumulative indicators per observation.
      type(ordinal_spec), intent(in) :: spec !! Ordinal model specification.
      real(dp), intent(in) :: z(:, :) !! Cluster odds-ratio design rows.
      real(dp), intent(in) :: odds_offset(:) !! Cluster log-odds-ratio offsets.
      real(dp), intent(in) :: beta(:) !! Current mean coefficients.
      real(dp), intent(in) :: alpha(:) !! Current association coefficients.
      real(dp), intent(out) :: pr(:) !! Mean residuals y-mu for the cluster.
      real(dp), intent(out) :: d(:, :) !! Derivative of cluster means with respect to beta.
      real(dp), intent(out) :: mu(:) !! Fitted cumulative probabilities.
      real(dp), intent(out) :: v(:, :) !! Working covariance of expanded cumulative indicators.
      integer, intent(out) :: status !! GEE_OK or invalid-mean status.
      real(dp), allocatable :: psi(:)
      real(dp), allocatable :: cross(:, :)
      integer :: i
      integer :: j
      integer :: pair
      integer :: r1
      integer :: r2
      integer :: p0

      call ordinal_mean_only(y, x, offset, spec%mean_link, beta, mu, pr, d, status)
      if (status /= GEE_OK) return
      v = 0.0_dp
      do i = 1, s
         r1 = (i - 1) * c
         v(r1 + 1:r1 + c, r1 + 1:r1 + c) = cumulative_variance(mu(r1 + 1:r1 + c), spec%reverse_coding)
      end do
      if (size(alpha) == 0 .or. s == 1) return
      pair = 0
      do i = 1, s - 1
         r1 = (i - 1) * c
         do j = i + 1, s
            r2 = (j - 1) * c
            p0 = pair * c * c
            allocate(psi(c * c), cross(c, c))
            psi = exp(matmul(z(p0 + 1:p0 + c * c, :), alpha) + odds_offset(p0 + 1:p0 + c * c))
            call cross_covariance(mu(r1 + 1:r1 + c), mu(r2 + 1:r2 + c), psi, cross)
            v(r1 + 1:r1 + c, r2 + 1:r2 + c) = cross
            v(r2 + 1:r2 + c, r1 + 1:r1 + c) = transpose(cross)
            deallocate(psi, cross)
            pair = pair + 1
         end do
      end do
   end subroutine ordinal_cluster_beta

   subroutine ordinal_mean_only(y, x, offset, link_code, beta, mu, pr, d, status)
      real(dp), intent(in) :: y(:) !! Expanded cumulative responses.
      real(dp), intent(in) :: x(:, :) !! Expanded mean design matrix.
      real(dp), intent(in) :: offset(:) !! Expanded mean offsets.
      integer, intent(in) :: link_code !! Ordinal mean-link code, logit/probit/cloglog.
      real(dp), intent(in) :: beta(:) !! Current mean coefficients.
      real(dp), intent(out) :: mu(:) !! Fitted cumulative probabilities.
      real(dp), intent(out) :: pr(:) !! Residuals y-mu.
      real(dp), intent(out) :: d(:, :) !! Mean derivative with respect to beta.
      integer, intent(out) :: status !! GEE_OK or invalid-mean status.
      real(dp) :: eta
      real(dp) :: dm
      integer :: i

      do i = 1, size(y)
         eta = dot_product(x(i, :), beta) + offset(i)
         mu(i) = link_inverse(eta, link_code)
         if (mu(i) <= 0.0_dp .or. mu(i) >= 1.0_dp) then
            status = GEE_ERR_INVALID_MEAN
            return
         end if
         dm = link_derivative(eta, link_code)
         pr(i) = y(i) - mu(i)
         d(i, :) = dm * x(i, :)
      end do
      status = GEE_OK
   end subroutine ordinal_mean_only

   subroutine ordinal_pair_alpha(pr1, pr2, mu1, mu2, z, odds_offset, reverse, alpha, u, v2, d2)
      real(dp), intent(in) :: pr1(:) !! Residual vector for the first ordinal observation in a pair.
      real(dp), intent(in) :: pr2(:) !! Residual vector for the second ordinal observation in a pair.
      real(dp), intent(in) :: mu1(:) !! Cumulative means for the first observation.
      real(dp), intent(in) :: mu2(:) !! Cumulative means for the second observation.
      real(dp), intent(in) :: z(:, :) !! Pair-specific odds-ratio design with c squared rows.
      real(dp), intent(in) :: odds_offset(:) !! Pair-specific log-odds-ratio offsets.
      logical, intent(in) :: reverse !! Whether cumulative indicators use reverse coding.
      real(dp), intent(in) :: alpha(:) !! Current log-odds-ratio coefficients.
      real(dp), intent(out) :: u(:) !! Association residual S minus its mean.
      real(dp), intent(out) :: v2(:, :) !! Covariance of association residuals.
      real(dp), intent(out) :: d2(:, :) !! Derivative of expected pair products with respect to alpha.
      real(dp), allocatable :: psi(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: p11(:, :)
      real(dp), allocatable :: esst(:, :)
      integer :: c
      integer :: i
      integer :: j
      integer :: k

      c = size(mu1)
      allocate(psi(c * c), sigma(c * c), p11(c, c), esst(c * c, c * c))
      psi = exp(matmul(z, alpha) + odds_offset)
      k = 0
      do i = 1, c
         do j = 1, c
            k = k + 1
            p11(i, j) = odds_to_p11(psi(k), mu1(i), mu2(j))
            sigma(k) = p11(i, j) - mu1(i) * mu2(j)
            u(k) = pr1(i) * pr2(j) - sigma(k)
            d2(k, :) = p11_odds_derivative(psi(k), mu1(i), mu2(j)) * psi(k) * z(k, :)
         end do
      end do
      call esst_matrix(mu1, mu2, p11, reverse, esst)
      v2 = esst - outer_product(sigma, sigma)
   end subroutine ordinal_pair_alpha

   subroutine ordinal_covariance(y, x, offset, cluster_sizes, c, spec, z, odds_offset, beta, alpha, result, status)
      real(dp), intent(in) :: y(:) !! Expanded cumulative responses.
      real(dp), intent(in) :: x(:, :) !! Expanded mean design matrix.
      real(dp), intent(in) :: offset(:) !! Expanded mean offsets.
      integer, intent(in) :: cluster_sizes(:) !! Original ordinal cluster sizes.
      integer, intent(in) :: c !! Number of cumulative indicators per observation.
      type(ordinal_spec), intent(in) :: spec !! Ordinal model specification.
      real(dp), intent(in) :: z(:, :) !! Odds-ratio design matrix.
      real(dp), intent(in) :: odds_offset(:) !! Log-odds-ratio offsets.
      real(dp), intent(in) :: beta(:) !! Fitted mean coefficients.
      real(dp), intent(in) :: alpha(:) !! Fitted association coefficients.
      type(ordinal_result), intent(inout) :: result !! Result receiving covariance and influence estimates.
      integer, intent(out) :: status !! GEE_OK or a numerical error.
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: hi(:, :)
      real(dp), allocatable :: meat(:, :)
      real(dp), allocatable :: score(:)
      real(dp), allocatable :: hinv(:, :)
      real(dp), allocatable :: robust(:, :)
      real(dp), allocatable :: ftotal(:, :)
      real(dp), allocatable :: ameat(:, :)
      real(dp), allocatable :: astable(:, :)
      real(dp), allocatable :: finv(:, :)
      integer :: p
      integer :: q
      integer :: l
      integer :: gi
      integer :: obs0
      integer :: z0
      integer :: zr
      integer :: info

      p = size(beta)
      q = size(alpha)
      l = p + q
      allocate(h(l, l), hi(l, l), meat(l, l), score(l), hinv(l, l), robust(l, l))
      h = 0.0_dp
      meat = 0.0_dp
      if (q > 0) then
         allocate(ftotal(q, q), ameat(q, q))
         ftotal = 0.0_dp
         ameat = 0.0_dp
      end if
      obs0 = 0
      z0 = 0
      do gi = 1, size(cluster_sizes)
         zr = cluster_odds_rows(cluster_sizes(gi), c, spec%corstr)
         call ordinal_cluster_hg(y(obs0 * c + 1:(obs0 + cluster_sizes(gi)) * c), &
            x(obs0 * c + 1:(obs0 + cluster_sizes(gi)) * c, :), &
            offset(obs0 * c + 1:(obs0 + cluster_sizes(gi)) * c), cluster_sizes(gi), c, spec, &
            z_segment(z, z0, zr), offset_segment(odds_offset, z0, zr), beta, alpha, hi, score, status)
         if (status /= GEE_OK) return
         h = h + hi
         meat = meat + outer_product(score, score)
         result%influence(:, gi) = score
         if (q > 0) then
            ftotal = ftotal + hi(p + 1:l, p + 1:l)
            ameat = ameat + outer_product(score(p + 1:l), score(p + 1:l))
         end if
         obs0 = obs0 + cluster_sizes(gi)
         z0 = z0 + zr
      end do
      call inverse_checked(h, hinv, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      robust = matmul(hinv, matmul(meat, transpose(hinv)))
      result%influence = matmul(hinv, result%influence)
      result%vbeta = robust(1:p, 1:p)
      result%vbeta_naive = hinv(1:p, 1:p)
      if (q > 0) then
         result%valpha = robust(p + 1:l, p + 1:l)
         result%valpha_naive = hinv(p + 1:l, p + 1:l)
         allocate(finv(q, q), astable(q, q))
         call inverse_checked(ftotal, finv, info)
         if (info /= 0) then
            status = GEE_ERR_SINGULAR
            return
         end if
         astable = matmul(finv, matmul(ameat, transpose(finv)))
         result%valpha_stable = astable
      end if
      status = GEE_OK
   end subroutine ordinal_covariance

   subroutine ordinal_cluster_hg(y, x, offset, s, c, spec, z, odds_offset, beta, alpha, h, score, status)
      real(dp), intent(in) :: y(:) !! Expanded responses for one cluster.
      real(dp), intent(in) :: x(:, :) !! Expanded mean design for one cluster.
      real(dp), intent(in) :: offset(:) !! Expanded mean offsets for one cluster.
      integer, intent(in) :: s !! Number of ordinal observations in the cluster.
      integer, intent(in) :: c !! Number of cumulative indicators per observation.
      type(ordinal_spec), intent(in) :: spec !! Ordinal model specification.
      real(dp), intent(in) :: z(:, :) !! Cluster odds-ratio design rows.
      real(dp), intent(in) :: odds_offset(:) !! Cluster log-odds-ratio offsets.
      real(dp), intent(in) :: beta(:) !! Fitted mean coefficients.
      real(dp), intent(in) :: alpha(:) !! Fitted association coefficients.
      real(dp), intent(out) :: h(:, :) !! Cluster sensitivity matrix in beta/alpha block order.
      real(dp), intent(out) :: score(:) !! Cluster estimating-function vector.
      integer, intent(out) :: status !! GEE_OK or a numerical error.
      real(dp), allocatable :: pr(:)
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: vinv(:, :)
      real(dp), allocatable :: u(:)
      real(dp), allocatable :: v2(:, :)
      real(dp), allocatable :: d2(:, :)
      real(dp), allocatable :: v2inv(:, :)
      real(dp), allocatable :: ubeta(:, :)
      real(dp), allocatable :: psi(:)
      integer :: p
      integer :: q
      integer :: j
      integer :: k
      integer :: pair
      integer :: p0
      integer :: info

      p = size(beta)
      q = size(alpha)
      h = 0.0_dp
      score = 0.0_dp
      allocate(pr(s * c), d(s * c, p), mu(s * c), v(s * c, s * c), vinv(s * c, s * c))
      call ordinal_cluster_beta(y, x, offset, s, c, spec, z, odds_offset, beta, alpha, pr, d, mu, v, status)
      if (status /= GEE_OK) return
      call inverse_checked(v, vinv, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      h(1:p, 1:p) = matmul(transpose(d), matmul(vinv, d))
      score(1:p) = matmul(transpose(d), matmul(vinv, pr))
      if (q == 0 .or. s == 1) then
         status = GEE_OK
         return
      end if

      pair = 0
      do j = 1, s - 1
         do k = j + 1, s
            p0 = pair * c * c
            allocate(u(c * c), v2(c * c, c * c), d2(c * c, q), v2inv(c * c, c * c))
            allocate(ubeta(c * c, p), psi(c * c))
            call ordinal_pair_alpha(pr((j - 1) * c + 1:j * c), pr((k - 1) * c + 1:k * c), &
               mu((j - 1) * c + 1:j * c), mu((k - 1) * c + 1:k * c), &
               z(p0 + 1:p0 + c * c, :), odds_offset(p0 + 1:p0 + c * c), &
               spec%reverse_coding, alpha, u, v2, d2)
            call inverse_checked(v2, v2inv, info)
            if (info /= 0) then
               status = GEE_ERR_SINGULAR
               return
            end if
            h(p + 1:p + q, p + 1:p + q) = h(p + 1:p + q, p + 1:p + q) + &
               matmul(transpose(d2), matmul(v2inv, d2))
            score(p + 1:p + q) = score(p + 1:p + q) + matmul(transpose(d2), matmul(v2inv, u))
            psi = exp(matmul(z(p0 + 1:p0 + c * c, :), alpha) + odds_offset(p0 + 1:p0 + c * c))
            call association_beta_derivative(mu((j - 1) * c + 1:j * c), mu((k - 1) * c + 1:k * c), &
               psi, d((j - 1) * c + 1:j * c, :), d((k - 1) * c + 1:k * c, :), &
               pr((j - 1) * c + 1:j * c), pr((k - 1) * c + 1:k * c), ubeta)
            h(p + 1:p + q, 1:p) = h(p + 1:p + q, 1:p) - &
               matmul(transpose(d2), matmul(v2inv, ubeta))
            deallocate(u, v2, d2, v2inv, ubeta, psi)
            pair = pair + 1
         end do
      end do
      status = GEE_OK
   end subroutine ordinal_cluster_hg

   subroutine association_beta_derivative(mu1, mu2, psi, d1, d2, pr1, pr2, out)
      real(dp), intent(in) :: mu1(:) !! Cumulative means of the first ordinal observation.
      real(dp), intent(in) :: mu2(:) !! Cumulative means of the second ordinal observation.
      real(dp), intent(in) :: psi(:) !! Pairwise cumulative odds ratios in row-major threshold order.
      real(dp), intent(in) :: d1(:, :) !! Mean derivatives for the first observation.
      real(dp), intent(in) :: d2(:, :) !! Mean derivatives for the second observation.
      real(dp), intent(in) :: pr1(:) !! Mean residuals for the first observation.
      real(dp), intent(in) :: pr2(:) !! Mean residuals for the second observation.
      real(dp), intent(out) :: out(:, :) !! Negative derivative contribution used by the association equation.
      real(dp) :: deriv(2)
      integer :: c
      integer :: i
      integer :: j
      integer :: k

      c = size(mu1)
      k = 0
      do i = 1, c
         do j = 1, c
            k = k + 1
            call p11_mean_derivatives(psi(k), mu1(i), mu2(j), deriv)
            deriv(1) = deriv(1) - mu2(j)
            deriv(2) = deriv(2) - mu1(i)
            out(k, :) = (-pr2(j) - deriv(1)) * d1(i, :) + (-pr1(i) - deriv(2)) * d2(j, :)
         end do
      end do
   end subroutine association_beta_derivative

   pure real(dp) function odds_to_p11(psi, mu1, mu2) result(p11)
      real(dp), intent(in) :: psi !! Positive cumulative odds ratio.
      real(dp), intent(in) :: mu1 !! First marginal cumulative probability in [0,1].
      real(dp), intent(in) :: mu2 !! Second marginal cumulative probability in [0,1].
      real(dp) :: exp1
      real(dp) :: disc

      if (abs(psi - 1.0_dp) < 0.001_dp) then
         p11 = mu1 * mu2
      else
         exp1 = 1.0_dp + (mu1 + mu2) * (psi - 1.0_dp)
         disc = exp1 * exp1 + 4.0_dp * psi * (1.0_dp - psi) * mu1 * mu2
         p11 = 0.5_dp / (psi - 1.0_dp) * (exp1 - sqrt(max(0.0_dp, disc)))
      end if
   end function odds_to_p11

   pure real(dp) function p11_odds_derivative(psi, mu1, mu2) result(ans)
      real(dp), intent(in) :: psi !! Positive cumulative odds ratio.
      real(dp), intent(in) :: mu1 !! First marginal cumulative probability.
      real(dp), intent(in) :: mu2 !! Second marginal cumulative probability.
      real(dp) :: e1
      real(dp) :: e2
      real(dp) :: e3
      real(dp) :: e5
      real(dp) :: e7
      real(dp) :: e8
      real(dp) :: e9
      real(dp) :: e10
      real(dp) :: e12
      real(dp) :: e14
      real(dp) :: e23

      if (abs(psi - 1.0_dp) < 0.001_dp) then
         ans = mu1 * mu2 * (-(mu1 + mu2) + mu1 * mu2 + 1.0_dp)
         return
      end if
      e1 = psi - 1.0_dp
      e2 = 0.5_dp / e1
      e3 = mu1 + mu2
      e5 = 1.0_dp + e3 * e1
      e7 = 4.0_dp * psi
      e8 = 1.0_dp - psi
      e9 = e7 * e8
      e10 = e9 * mu1
      e12 = e5 * e5 + e10 * mu2
      e14 = e5 - sqrt(e12)
      e23 = 1.0_dp / sqrt(e12)
      ans = e2 * (e3 - 0.5_dp * ((2.0_dp * e3 * e5 + (4.0_dp * e8 - e7) * mu1 * mu2) * e23)) - &
         0.5_dp / (e1 * e1) * e14
   end function p11_odds_derivative

   pure subroutine p11_mean_derivatives(psi, mu1, mu2, deriv)
      real(dp), intent(in) :: psi !! Positive cumulative odds ratio.
      real(dp), intent(in) :: mu1 !! First marginal cumulative probability.
      real(dp), intent(in) :: mu2 !! Second marginal cumulative probability.
      real(dp), intent(out) :: deriv(2) !! Derivatives of p11 with respect to mu1 and mu2.
      real(dp) :: e1
      real(dp) :: e2
      real(dp) :: e3
      real(dp) :: e5
      real(dp) :: e7
      real(dp) :: e8
      real(dp) :: e9
      real(dp) :: e10
      real(dp) :: e12
      real(dp) :: e23
      real(dp) :: e33

      if (abs(psi - 1.0_dp) < 0.001_dp) then
         deriv(1) = mu2
         deriv(2) = mu1
         return
      end if
      e1 = psi - 1.0_dp
      e2 = 0.5_dp / e1
      e3 = mu1 + mu2
      e5 = 1.0_dp + e3 * e1
      e7 = 4.0_dp * psi
      e8 = 1.0_dp - psi
      e9 = e7 * e8
      e10 = e9 * mu1
      e12 = e5 * e5 + e10 * mu2
      e23 = 1.0_dp / sqrt(e12)
      e33 = 2.0_dp * e1 * e5
      deriv(1) = e2 * (e1 - 0.5_dp * ((e33 + e9 * mu2) * e23))
      deriv(2) = e2 * (e1 - 0.5_dp * ((e33 + e10) * e23))
   end subroutine p11_mean_derivatives

   pure function cumulative_variance(mu, reverse) result(v)
      real(dp), intent(in) :: mu(:) !! Cumulative probabilities for one ordinal observation.
      logical, intent(in) :: reverse !! Whether cumulative indicators use reverse coding.
      real(dp) :: v(size(mu), size(mu))
      integer :: i
      integer :: j
      integer :: ij

      do i = 1, size(mu)
         do j = 1, size(mu)
            if (reverse) then
               ij = max(i, j)
            else
               ij = min(i, j)
            end if
            v(i, j) = mu(ij) - mu(i) * mu(j)
         end do
      end do
   end function cumulative_variance

   subroutine cross_covariance(mu1, mu2, psi, v)
      real(dp), intent(in) :: mu1(:) !! Cumulative means for the first observation.
      real(dp), intent(in) :: mu2(:) !! Cumulative means for the second observation.
      real(dp), intent(in) :: psi(:) !! Cumulative odds ratios in row-major threshold order.
      real(dp), intent(out) :: v(:, :) !! Cross-covariance of cumulative indicators.
      integer :: i
      integer :: j
      integer :: k

      k = 0
      do i = 1, size(mu1)
         do j = 1, size(mu2)
            k = k + 1
            v(i, j) = odds_to_p11(psi(k), mu1(i), mu2(j)) - mu1(i) * mu2(j)
         end do
      end do
   end subroutine cross_covariance

   subroutine esst_matrix(mu1, mu2, p11, reverse, esst)
      real(dp), intent(in) :: mu1(:) !! Cumulative means for the first observation.
      real(dp), intent(in) :: mu2(:) !! Cumulative means for the second observation.
      real(dp), intent(in) :: p11(:, :) !! Joint cumulative probabilities for the observation pair.
      logical, intent(in) :: reverse !! Whether cumulative indicators use reverse coding.
      real(dp), intent(out) :: esst(:, :) !! Expected outer product of centered pair-products.
      real(dp), allocatable :: block(:, :)
      integer :: c
      integer :: c1
      integer :: c3
      integer :: r1
      integer :: r3

      c = size(mu1)
      esst = 0.0_dp
      allocate(block(c, c))
      do c1 = 1, c
         r1 = (c1 - 1) * c
         do c3 = c1, c
            r3 = (c3 - 1) * c
            call esst_block(mu1, mu2, p11, c1, c3, reverse, block)
            esst(r1 + 1:r1 + c, r3 + 1:r3 + c) = block
            if (c3 > c1) esst(r3 + 1:r3 + c, r1 + 1:r1 + c) = transpose(block)
         end do
      end do
   end subroutine esst_matrix

   subroutine esst_block(mu1, mu2, p11, c1, c3, reverse, block)
      real(dp), intent(in) :: mu1(:) !! Cumulative means for the first observation.
      real(dp), intent(in) :: mu2(:) !! Cumulative means for the second observation.
      real(dp), intent(in) :: p11(:, :) !! Joint cumulative probabilities for the observation pair.
      integer, intent(in) :: c1 !! First cumulative-index component of the pair product.
      integer, intent(in) :: c3 !! Second cumulative-index component of the pair product.
      logical, intent(in) :: reverse !! Whether cumulative indicators use reverse coding.
      real(dp), intent(out) :: block(:, :) !! c-by-c ESST block for indices c1 and c3.
      integer :: c2
      integer :: c4
      integer :: c13
      integer :: c24

      if (reverse) then
         c13 = max(c1, c3)
      else
         c13 = min(c1, c3)
      end if
      do c2 = 1, size(mu2)
         do c4 = c2, size(mu2)
            if (reverse) then
               c24 = max(c2, c4)
            else
               c24 = min(c2, c4)
            end if
            block(c2, c4) = p11(c13, c24) - p11(c13, c2) * mu2(c4) - &
               p11(c13, c4) * mu2(c2) + mu1(c13) * mu2(c2) * mu2(c4) - &
               p11(c1, c24) * mu1(c3) + p11(c1, c2) * mu1(c3) * mu2(c4) + &
               p11(c1, c4) * mu1(c3) * mu2(c2) - &
               3.0_dp * mu1(c1) * mu1(c3) * mu2(c2) * mu2(c4) - &
               p11(c3, c24) * mu1(c1) + p11(c3, c2) * mu1(c1) * mu2(c4) + &
               p11(c3, c4) * mu1(c1) * mu2(c2) + mu1(c1) * mu1(c3) * mu2(c24)
            if (c4 > c2) block(c4, c2) = block(c2, c4)
         end do
      end do
   end subroutine esst_block

   subroutine apply_root_weights(d, pr, weights)
      real(dp), intent(inout) :: d(:, :) !! Mean-derivative matrix to row-scale by sqrt(weights).
      real(dp), intent(inout) :: pr(:) !! Residual vector to scale by sqrt(weights).
      real(dp), intent(in) :: weights(:) !! Nonnegative expanded observation weights.
      real(dp) :: sw
      integer :: i

      do i = 1, size(weights)
         sw = sqrt(weights(i))
         d(i, :) = sw * d(i, :)
         pr(i) = sw * pr(i)
      end do
   end subroutine apply_root_weights

   pure integer function cluster_odds_rows(cluster_size, c, corstr) result(nrow)
      integer, intent(in) :: cluster_size !! Number of ordinal observations in one cluster.
      integer, intent(in) :: c !! Number of cumulative indicators per ordinal observation.
      integer, intent(in) :: corstr !! Ordinal association structure code.

      if (corstr == COR_INDEPENDENCE) then
         nrow = 0
      else
         nrow = cluster_size * (cluster_size - 1) / 2 * c * c
      end if
   end function cluster_odds_rows

   pure integer function expected_odds_rows(cluster_sizes, c, corstr) result(nrow)
      integer, intent(in) :: cluster_sizes(:) !! Original ordinal cluster sizes.
      integer, intent(in) :: c !! Number of cumulative indicators per observation.
      integer, intent(in) :: corstr !! Ordinal association structure code.

      if (corstr == COR_INDEPENDENCE) then
         nrow = 0
      else
         nrow = sum(cluster_sizes * (cluster_sizes - 1) / 2) * c * c
      end if
   end function expected_odds_rows

   pure logical function valid_ordinal_link(link_code) result(ok)
      integer, intent(in) :: link_code !! Candidate mean-link code.

      ok = link_code == LINK_LOGIT .or. link_code == 3 .or. link_code == 4
   end function valid_ordinal_link

   pure logical function valid_ordinal_corstr(corstr) result(ok)
      integer, intent(in) :: corstr !! Candidate ordinal association structure code.

      ok = corstr == COR_INDEPENDENCE .or. corstr == COR_EXCHANGEABLE .or. &
         corstr == COR_UNSTRUCTURED .or. corstr == COR_USERDEFINED
   end function valid_ordinal_corstr

   function z_segment(z, start, count_rows) result(seg)
      real(dp), intent(in) :: z(:, :) !! Full association design matrix.
      integer, intent(in) :: start !! Number of association rows preceding the requested segment.
      integer, intent(in) :: count_rows !! Number of rows in the requested segment.
      real(dp), allocatable :: seg(:, :)

      if (count_rows == 0) then
         allocate(seg(0, size(z, 2)))
      else
         allocate(seg(count_rows, size(z, 2)))
         seg = z(start + 1:start + count_rows, :)
      end if
   end function z_segment

   function offset_segment(offset, start, count_rows) result(seg)
      real(dp), intent(in) :: offset(:) !! Full log-odds-ratio offset vector.
      integer, intent(in) :: start !! Number of offset elements preceding the requested segment.
      integer, intent(in) :: count_rows !! Number of offset elements in the requested segment.
      real(dp), allocatable :: seg(:)

      if (count_rows == 0) then
         allocate(seg(0))
      else
         allocate(seg(count_rows))
         seg = offset(start + 1:start + count_rows)
      end if
   end function offset_segment

end module geepack_ordinal
