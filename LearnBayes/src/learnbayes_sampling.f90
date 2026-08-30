module learnbayes_sampling
   use learnbayes_distributions, only: dmnorm, dmt, rmnorm, rmt
   use learnbayes_kinds, only: dp
   use learnbayes_linalg, only: cholesky_lower, inverse_matrix
   use learnbayes_rng, only: rng_discrete, rng_normal, rng_state, rng_uniform
   use learnbayes_types, only: importance_result, laplace_result, log_density_callback, mcmc_result
   implicit none
   private

   real(dp), parameter :: pi_dp = acos(-1.0_dp)

   public :: eval_callback
   public :: gibbs
   public :: impsampling
   public :: indepmetrop
   public :: laplace
   public :: rejectsampling
   public :: rwmetrop
   public :: simcontour
   public :: sir

contains

   function eval_callback(callback, theta) result(value)
      type(log_density_callback), intent(in) :: callback !! Bound user callback plus its numeric data and parameter context.
      real(dp), intent(in) :: theta(:) !! Parameter vector passed to the callback.
      real(dp) :: value
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: params(:)

      if (.not. associated(callback%eval)) then
         value = -huge(1.0_dp)
         return
      end if
      if (allocated(callback%data)) then
         allocate(data(size(callback%data, 1), size(callback%data, 2)))
         data = callback%data
      else
         allocate(data(0, 0))
      end if
      if (allocated(callback%params)) then
         allocate(params(size(callback%params)))
         params = callback%params
      else
         allocate(params(0))
      end if
      value = callback%eval(theta, data, params)
   end function eval_callback

   subroutine gibbs(rng, logpost, start, m, scale, result)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for all coordinate proposals and accept/reject draws.
      type(log_density_callback), intent(in) :: logpost !! User log posterior for the Metropolis-within-Gibbs sampler.
      real(dp), intent(in) :: start(:) !! Initial parameter vector.
      integer, intent(in) :: m !! Number of complete Gibbs sweeps to simulate.
      real(dp), intent(in) :: scale(:) !! Coordinate-wise Gaussian random-walk standard deviations.
      type(mcmc_result), intent(out) :: result !! Simulated parameter matrix and coordinate acceptance rates.
      real(dp), allocatable :: theta0(:)
      real(dp), allocatable :: theta1(:)
      real(dp) :: f0
      real(dp) :: f1
      real(dp) :: logu
      integer :: i
      integer :: j
      integer :: p

      p = size(start)
      allocate(result%par(m, p), result%accept_by_parameter(p), theta0(p), theta1(p))
      result%accept_by_parameter = 0.0_dp
      theta0 = start
      f0 = eval_callback(logpost, theta0)
      do i = 1, m
         do j = 1, p
            theta1 = theta0
            theta1(j) = theta0(j) + rng_normal(rng)*scale(j)
            f1 = eval_callback(logpost, theta1)
            logu = log(rng_uniform(rng))
            if (logu < min(0.0_dp, f1 - f0)) then
               theta0 = theta1
               f0 = f1
               result%accept_by_parameter(j) = result%accept_by_parameter(j) + 1.0_dp
            end if
            result%par(i, :) = theta0
         end do
      end do
      result%accept_by_parameter = result%accept_by_parameter/real(m, dp)
      result%accept_rate = sum(result%accept_by_parameter)/real(p, dp)
   end subroutine gibbs

   subroutine rwmetrop(rng, logpost, proposal_var, proposal_scale, start, m, result, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for multivariate proposals and accept/reject draws.
      type(log_density_callback), intent(in) :: logpost !! User-defined log posterior for the random-walk Metropolis sampler.
      real(dp), intent(in) :: proposal_var(:, :) !! Positive-definite proposal covariance matrix before scalar scaling.
      real(dp), intent(in) :: proposal_scale !! Nonnegative multiplier applied to proposal standard deviations.
      real(dp), intent(in) :: start(:) !! Initial parameter vector.
      integer, intent(in) :: m !! Number of Metropolis iterations to simulate.
      type(mcmc_result), intent(out) :: result !! Simulated parameter matrix and overall acceptance fraction.
      integer, intent(out) :: info !! Zero on success; nonzero if proposal dimensions or Cholesky factorization fail.
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: candidate(:)
      real(dp), allocatable :: z(:)
      real(dp) :: f0
      real(dp) :: f1
      integer :: i
      integer :: j
      integer :: p
      integer :: accepted

      p = size(start)
      info = 0
      if (size(proposal_var, 1) /= p .or. size(proposal_var, 2) /= p) then
         info = -1
         return
      end if
      allocate(l(p, p), theta(p), candidate(p), z(p), result%par(m, p))
      call cholesky_lower(proposal_var, l, info)
      if (info /= 0) return
      theta = start
      f0 = eval_callback(logpost, theta)
      accepted = 0
      do i = 1, m
         do j = 1, p
            z(j) = rng_normal(rng)
         end do
         candidate = theta + proposal_scale*matmul(l, z)
         f1 = eval_callback(logpost, candidate)
         if (log(rng_uniform(rng)) < min(0.0_dp, f1 - f0)) then
            theta = candidate
            f0 = f1
            accepted = accepted + 1
         end if
         result%par(i, :) = theta
      end do
      result%accept_rate = real(accepted, dp)/real(m, dp)
      allocate(result%accept_by_parameter(0))
   end subroutine rwmetrop

   subroutine indepmetrop(rng, logpost, proposal_mu, proposal_var, start, m, result, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for independent Gaussian proposals.
      type(log_density_callback), intent(in) :: logpost !! User-defined log posterior for the independence Metropolis sampler.
      real(dp), intent(in) :: proposal_mu(:) !! Mean vector of the independent multivariate-normal proposal.
      real(dp), intent(in) :: proposal_var(:, :) !! Positive-definite covariance of the independent proposal.
      real(dp), intent(in) :: start(:) !! Initial parameter vector.
      integer, intent(in) :: m !! Number of Metropolis iterations to simulate.
      type(mcmc_result), intent(out) :: result !! Simulated parameter matrix and overall acceptance fraction.
      integer, intent(out) :: info !! Zero on success; nonzero if proposal dimensions or factorization fail.
      real(dp), allocatable :: draw(:, :)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: candidate(:)
      real(dp) :: f0
      real(dp) :: f1
      real(dp) :: logratio
      integer :: i
      integer :: p
      integer :: accepted

      p = size(start)
      info = 0
      if (size(proposal_mu) /= p .or. size(proposal_var, 1) /= p .or. size(proposal_var, 2) /= p) then
         info = -1
         return
      end if
      allocate(draw(1, p), theta(p), candidate(p), result%par(m, p))
      theta = start
      f0 = eval_callback(logpost, theta)
      accepted = 0
      do i = 1, m
         call rmnorm(rng, 1, proposal_mu, proposal_var, draw, info)
         if (info /= 0) return
         candidate = draw(1, :)
         f1 = eval_callback(logpost, candidate)
         logratio = dmnorm(theta, proposal_mu, proposal_var, .true.) - &
            dmnorm(candidate, proposal_mu, proposal_var, .true.) + f1 - f0
         if (log(rng_uniform(rng)) < min(0.0_dp, logratio)) then
            theta = candidate
            f0 = f1
            accepted = accepted + 1
         end if
         result%par(i, :) = theta
      end do
      result%accept_rate = real(accepted, dp)/real(m, dp)
      allocate(result%accept_by_parameter(0))
   end subroutine indepmetrop

   subroutine laplace(logpost, start, result, max_iter, tol)
      type(log_density_callback), intent(in) :: logpost !! User-defined log posterior or log integrand to maximize.
      real(dp), intent(in) :: start(:) !! Starting parameter vector for Nelder-Mead optimization.
      type(laplace_result), intent(out) :: result !! Posterior mode, covariance approximation, log integral, and convergence flag.
      integer, intent(in), optional :: max_iter !! Optional maximum Nelder-Mead iterations; defaults to 1000.
      real(dp), intent(in), optional :: tol !! Optional simplex convergence tolerance; defaults to 1e-8.
      real(dp), allocatable :: simplex(:, :)
      real(dp), allocatable :: f(:)
      real(dp), allocatable :: centroid(:)
      real(dp), allocatable :: xr(:)
      real(dp), allocatable :: xe(:)
      real(dp), allocatable :: xc(:)
      real(dp), allocatable :: hess(:, :)
      real(dp), allocatable :: invh(:, :)
      real(dp) :: alpha
      real(dp) :: gamma
      real(dp) :: rho
      real(dp) :: sigma
      real(dp) :: fr
      real(dp) :: fe
      real(dp) :: fc
      real(dp) :: tolerance
      real(dp) :: spread
      real(dp) :: simplex_spread
      real(dp) :: step
      real(dp) :: logdet
      integer :: i
      integer :: iter
      integer :: iter_max
      integer :: n
      integer :: info

      n = size(start)
      iter_max = 1000
      if (present(max_iter)) iter_max = max_iter
      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      alpha = 1.0_dp
      gamma = 2.0_dp
      rho = 0.5_dp
      sigma = 0.5_dp
      allocate(simplex(n + 1, n), f(n + 1), centroid(n), xr(n), xe(n), xc(n))
      simplex(1, :) = start
      do i = 1, n
         simplex(i + 1, :) = start
         step = 0.05_dp*max(1.0_dp, abs(start(i)))
         simplex(i + 1, i) = simplex(i + 1, i) + step
      end do
      do i = 1, n + 1
         f(i) = eval_callback(logpost, simplex(i, :))
      end do

      result%converged = .false.
      iter = 0
      do iter = 1, iter_max
         call sort_simplex_desc(simplex, f)
         spread = maxval(abs(f - f(1)))
         simplex_spread = maxval(abs(simplex(2:n + 1, :) - spread_rows(simplex(1, :), n)))
         if (spread <= tolerance*(1.0_dp + abs(f(1))) .and. simplex_spread <= sqrt(tolerance)) then
            result%converged = .true.
            exit
         end if
         centroid = sum(simplex(1:n, :), dim=1)/real(n, dp)
         xr = centroid + alpha*(centroid - simplex(n + 1, :))
         fr = eval_callback(logpost, xr)
         if (fr > f(1)) then
            xe = centroid + gamma*(xr - centroid)
            fe = eval_callback(logpost, xe)
            if (fe > fr) then
               simplex(n + 1, :) = xe
               f(n + 1) = fe
            else
               simplex(n + 1, :) = xr
               f(n + 1) = fr
            end if
         else if (fr > f(n)) then
            simplex(n + 1, :) = xr
            f(n + 1) = fr
         else
            if (fr > f(n + 1)) then
               xc = centroid + rho*(xr - centroid)
            else
               xc = centroid + rho*(simplex(n + 1, :) - centroid)
            end if
            fc = eval_callback(logpost, xc)
            if (fc > f(n + 1)) then
               simplex(n + 1, :) = xc
               f(n + 1) = fc
            else
               do i = 2, n + 1
                  simplex(i, :) = simplex(1, :) + sigma*(simplex(i, :) - simplex(1, :))
                  f(i) = eval_callback(logpost, simplex(i, :))
               end do
            end if
         end if
      end do
      call sort_simplex_desc(simplex, f)
      allocate(result%mode(n), result%var(n, n), hess(n, n), invh(n, n))
      result%mode = simplex(1, :)
      result%iterations = min(iter, iter_max)
      call numerical_hessian(logpost, result%mode, hess)
      call inverse_matrix(-hess, invh, info)
      if (info /= 0) then
         result%var = 0.0_dp
         result%log_integral = -huge(1.0_dp)
         result%converged = .false.
         return
      end if
      result%var = invh
      call logdet_spd(result%var, logdet, info)
      if (info /= 0) then
         result%log_integral = -huge(1.0_dp)
         result%converged = .false.
         return
      end if
      result%log_integral = 0.5_dp*real(n, dp)*log(2.0_dp*pi_dp) + 0.5_dp*logdet + f(1)
   end subroutine laplace

   subroutine numerical_hessian(logpost, x, hess)
      type(log_density_callback), intent(in) :: logpost !! User log density whose Hessian is approximated by central differences.
      real(dp), intent(in) :: x(:) !! Point at which the Hessian is evaluated.
      real(dp), intent(out) :: hess(:, :) !! Symmetric numerical Hessian matrix.
      real(dp), allocatable :: xp(:)
      real(dp), allocatable :: xm(:)
      real(dp), allocatable :: xpp(:)
      real(dp), allocatable :: xpm(:)
      real(dp), allocatable :: xmp(:)
      real(dp), allocatable :: xmm(:)
      real(dp), allocatable :: step(:)
      real(dp) :: f0
      integer :: i
      integer :: j
      integer :: n

      n = size(x)
      allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n), step(n))
      do i = 1, n
         step(i) = 1.0e-4_dp*max(1.0_dp, abs(x(i)))
      end do
      f0 = eval_callback(logpost, x)
      hess = 0.0_dp
      do i = 1, n
         xp = x
         xm = x
         xp(i) = xp(i) + step(i)
         xm(i) = xm(i) - step(i)
         hess(i, i) = (eval_callback(logpost, xp) - 2.0_dp*f0 + eval_callback(logpost, xm))/(step(i)*step(i))
         do j = i + 1, n
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = xpp(i) + step(i)
            xpp(j) = xpp(j) + step(j)
            xpm(i) = xpm(i) + step(i)
            xpm(j) = xpm(j) - step(j)
            xmp(i) = xmp(i) - step(i)
            xmp(j) = xmp(j) + step(j)
            xmm(i) = xmm(i) - step(i)
            xmm(j) = xmm(j) - step(j)
            hess(i, j) = (eval_callback(logpost, xpp) - eval_callback(logpost, xpm) - &
               eval_callback(logpost, xmp) + eval_callback(logpost, xmm))/(4.0_dp*step(i)*step(j))
            hess(j, i) = hess(i, j)
         end do
      end do
   end subroutine numerical_hessian

   subroutine logdet_spd(a, logdet, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix.
      real(dp), intent(out) :: logdet !! Natural logarithm of det(a) on success.
      integer, intent(out) :: info !! Zero on success; nonzero if Cholesky factorization fails.
      real(dp), allocatable :: l(:, :)
      integer :: i
      integer :: n

      n = size(a, 1)
      allocate(l(n, n))
      call cholesky_lower(a, l, info)
      if (info /= 0) then
         logdet = -huge(1.0_dp)
         return
      end if
      logdet = 0.0_dp
      do i = 1, n
         logdet = logdet + 2.0_dp*log(l(i, i))
      end do
   end subroutine logdet_spd

   subroutine impsampling(rng, logf, statistic, proposal_mean, proposal_var, proposal_df, n, result, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used to draw the multivariate-t proposal sample.
      type(log_density_callback), intent(in) :: logf !! Target log density known up to a proportionality constant.
      type(log_density_callback), intent(in) :: statistic !! Scalar function h(theta) whose posterior expectation is estimated.
      real(dp), intent(in) :: proposal_mean(:) !! Location vector of the multivariate-t importance proposal.
      real(dp), intent(in) :: proposal_var(:, :) !! Scale matrix of the multivariate-t importance proposal.
      real(dp), intent(in) :: proposal_df !! Positive degrees of freedom for the multivariate-t proposal.
      integer, intent(in) :: n !! Number of proposal draws.
      type(importance_result), intent(out) :: result !! Importance estimate, standard error, draws, and stabilized weights.
      integer, intent(out) :: info !! Zero on success; nonzero if the proposal draw fails.
      real(dp), allocatable :: lf(:)
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: lp(:)
      real(dp) :: md
      real(dp) :: denom
      integer :: i
      integer :: d

      d = size(proposal_mean)
      allocate(result%theta(n, d), result%weight(n), lf(n), h(n), lp(n))
      call rmt(rng, n, proposal_mean, proposal_var, proposal_df, result%theta, info)
      if (info /= 0) return
      do i = 1, n
         lf(i) = eval_callback(logf, result%theta(i, :))
         h(i) = eval_callback(statistic, result%theta(i, :))
         lp(i) = dmt(result%theta(i, :), proposal_mean, proposal_var, proposal_df, .true.)
      end do
      md = maxval(lf - lp)
      result%weight = exp(lf - lp - md)
      denom = sum(result%weight)
      result%estimate = sum(result%weight*h)/denom
      result%se = sqrt(sum((h - result%estimate)**2*result%weight**2))/denom
   end subroutine impsampling

   subroutine rejectsampling(rng, logf, proposal_mean, proposal_var, proposal_df, dmax, n, theta, n_accept, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for proposal and uniform acceptance draws.
      type(log_density_callback), intent(in) :: logf !! Target log density known up to proportionality.
      real(dp), intent(in) :: proposal_mean(:) !! Location vector of the multivariate-t rejection proposal.
      real(dp), intent(in) :: proposal_var(:, :) !! Scale matrix of the multivariate-t rejection proposal.
      real(dp), intent(in) :: proposal_df !! Positive proposal degrees of freedom.
      real(dp), intent(in) :: dmax !! Log envelope constant such that log target - log proposal <= dmax.
      integer, intent(in) :: n !! Number of proposal draws to attempt, matching the R routine's batch semantics.
      real(dp), allocatable, intent(out) :: theta(:, :) !! Accepted draws; allocated with n_accept rows.
      integer, intent(out) :: n_accept !! Number of proposals accepted from the batch.
      integer, intent(out) :: info !! Zero on success; nonzero if multivariate-t proposal generation fails.
      real(dp), allocatable :: candidate(:, :)
      logical, allocatable :: keep(:)
      real(dp) :: logprob
      integer :: i
      integer :: j
      integer :: d

      d = size(proposal_mean)
      allocate(candidate(n, d), keep(n))
      call rmt(rng, n, proposal_mean, proposal_var, proposal_df, candidate, info)
      if (info /= 0) then
         allocate(theta(0, d))
         n_accept = 0
         return
      end if
      keep = .false.
      do i = 1, n
         logprob = eval_callback(logf, candidate(i, :)) - &
            dmt(candidate(i, :), proposal_mean, proposal_var, proposal_df, .true.) - dmax
         keep(i) = log(rng_uniform(rng)) < min(0.0_dp, logprob)
      end do
      n_accept = count(keep)
      allocate(theta(n_accept, d))
      j = 0
      do i = 1, n
         if (keep(i)) then
            j = j + 1
            theta(j, :) = candidate(i, :)
         end if
      end do
   end subroutine rejectsampling

   subroutine sir(rng, logf, proposal_mean, proposal_var, proposal_df, n, theta, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for proposal draws and weighted resampling.
      type(log_density_callback), intent(in) :: logf !! Target log density known up to proportionality.
      real(dp), intent(in) :: proposal_mean(:) !! Location vector of the multivariate-t proposal.
      real(dp), intent(in) :: proposal_var(:, :) !! Scale matrix of the multivariate-t proposal.
      real(dp), intent(in) :: proposal_df !! Positive degrees of freedom of the multivariate-t proposal.
      integer, intent(in) :: n !! Number of proposal draws and resampled output draws.
      real(dp), intent(out) :: theta(:, :) !! Sampling-importance-resampling output matrix shaped (n,d).
      integer, intent(out) :: info !! Zero on success; nonzero if multivariate-t proposal generation fails.
      real(dp), allocatable :: proposal(:, :)
      real(dp), allocatable :: logw(:)
      real(dp), allocatable :: weight(:)
      real(dp) :: md
      integer :: i
      integer :: index
      integer :: d

      d = size(proposal_mean)
      allocate(proposal(n, d), logw(n), weight(n))
      call rmt(rng, n, proposal_mean, proposal_var, proposal_df, proposal, info)
      if (info /= 0) return
      do i = 1, n
         logw(i) = eval_callback(logf, proposal(i, :)) - &
            dmt(proposal(i, :), proposal_mean, proposal_var, proposal_df, .true.)
      end do
      md = maxval(logw)
      weight = exp(logw - md)
      do i = 1, n
         index = rng_discrete(rng, weight)
         theta(i, :) = proposal(index, :)
      end do
   end subroutine sir

   subroutine simcontour(rng, logf, limits, m, x, y)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for grid resampling and cell jitter.
      type(log_density_callback), intent(in) :: logf !! Two-parameter target log density evaluated on a 50-by-50 grid.
      real(dp), intent(in) :: limits(4) !! Rectangle xmin, xmax, ymin, ymax defining the contour simulation grid.
      integer, intent(in) :: m !! Number of approximate draws returned from the gridded distribution.
      real(dp), intent(out) :: x(:) !! Simulated first-coordinate values, with length m.
      real(dp), intent(out) :: y(:) !! Simulated second-coordinate values, with length m.
      real(dp), allocatable :: weight(:)
      real(dp) :: theta(2)
      real(dp) :: dx
      real(dp) :: dy
      real(dp) :: mx
      integer :: i
      integer :: ix
      integer :: iy
      integer :: index
      integer, parameter :: ng = 50

      allocate(weight(ng*ng))
      dx = (limits(2) - limits(1))/real(ng - 1, dp)
      dy = (limits(4) - limits(3))/real(ng - 1, dp)
      index = 0
      do iy = 1, ng
         theta(2) = limits(3) + real(iy - 1, dp)*dy
         do ix = 1, ng
            theta(1) = limits(1) + real(ix - 1, dp)*dx
            index = index + 1
            weight(index) = eval_callback(logf, theta)
         end do
      end do
      mx = maxval(weight)
      weight = exp(weight - mx)
      do i = 1, m
         index = rng_discrete(rng, weight)
         iy = (index - 1)/ng + 1
         ix = index - (iy - 1)*ng
         x(i) = limits(1) + real(ix - 1, dp)*dx + (rng_uniform(rng) - 0.5_dp)*dx
         y(i) = limits(3) + real(iy - 1, dp)*dy + (rng_uniform(rng) - 0.5_dp)*dy
      end do
   end subroutine simcontour

   pure function spread_rows(row, nrow) result(out)
      real(dp), intent(in) :: row(:) !! Reference simplex row replicated for coordinate-spread calculation.
      integer, intent(in) :: nrow !! Number of replicated rows required.
      real(dp) :: out(nrow, size(row))
      integer :: i

      do i = 1, nrow
         out(i, :) = row
      end do
   end function spread_rows

   subroutine sort_simplex_desc(simplex, f)
      real(dp), intent(inout) :: simplex(:, :) !! Nelder-Mead simplex rows reordered from highest to lowest objective.
      real(dp), intent(inout) :: f(:) !! Objective values reordered consistently with simplex.
      real(dp), allocatable :: row(:)
      real(dp) :: key
      integer :: i
      integer :: j
      integer :: n

      n = size(simplex, 2)
      allocate(row(n))
      do i = 2, size(f)
         key = f(i)
         row = simplex(i, :)
         j = i - 1
         do while (j >= 1)
            if (f(j) >= key) exit
            f(j + 1) = f(j)
            simplex(j + 1, :) = simplex(j, :)
            j = j - 1
         end do
         f(j + 1) = key
         simplex(j + 1, :) = row
      end do
   end subroutine sort_simplex_desc

end module learnbayes_sampling
