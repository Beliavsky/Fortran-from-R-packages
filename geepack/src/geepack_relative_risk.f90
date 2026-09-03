! Translation of the computational COPY construction in R/relative-risk-regression.R.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_relative_risk
   use r_kinds, only : dp
   use geepack_status, only : GEE_OK, GEE_ERR_SHAPE, GEE_ERR_ARGUMENT
   use geepack_links, only : LINK_LOG, LINK_IDENTITY, VAR_BINOMIAL
   use geepack_gee, only : gee_spec, gee_result, fit_geese
   implicit none
   private

   public :: make_relative_risk_copy, fit_relative_risk

contains

   subroutine make_relative_risk_copy(y, x, cluster_sizes, ncopy, ycopy, xcopy, cluster_sizes_copy, &
      weights_copy, status, offset, offset_copy, waves, waves_copy)
      real(dp), intent(in) :: y(:) !! Binary responses coded zero or one.
      real(dp), intent(in) :: x(:, :) !! Mean-model design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Original contiguous cluster sizes.
      integer, intent(in) :: ncopy !! COPY tuning count; must exceed one.
      real(dp), allocatable, intent(out) :: ycopy(:) !! Original responses followed by their complements.
      real(dp), allocatable, intent(out) :: xcopy(:, :) !! Design matrix duplicated for the complement copy.
      integer, allocatable, intent(out) :: cluster_sizes_copy(:) !! Original cluster-size sequence repeated twice.
      real(dp), allocatable, intent(out) :: weights_copy(:) !! COPY weights 1-1/ncopy and 1/ncopy.
      integer, intent(out) :: status !! GEE_OK or an argument/shape error.
      real(dp), optional, intent(in) :: offset(:) !! Optional original mean-model offsets.
      real(dp), allocatable, optional, intent(out) :: offset_copy(:) !! Duplicated offsets when requested.
      integer, optional, intent(in) :: waves(:) !! Optional original wave identifiers.
      integer, allocatable, optional, intent(out) :: waves_copy(:) !! Duplicated wave identifiers when requested.
      integer :: n
      integer :: g
      integer :: pos
      integer :: j

      n = size(y)
      if (ncopy <= 1 .or. size(x, 1) /= n .or. sum(cluster_sizes) /= n .or. &
          any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
         status = GEE_ERR_ARGUMENT
         allocate(ycopy(0), xcopy(0, 0), cluster_sizes_copy(0), weights_copy(0))
         return
      end if
      allocate(ycopy(2 * n), xcopy(2 * n, size(x, 2)), cluster_sizes_copy(2 * size(cluster_sizes)))
      allocate(weights_copy(2 * n))
      ycopy(1:n) = y
      ycopy(n + 1:2 * n) = 1.0_dp - y
      xcopy(1:n, :) = x
      xcopy(n + 1:2 * n, :) = x
      cluster_sizes_copy(1:size(cluster_sizes)) = cluster_sizes
      cluster_sizes_copy(size(cluster_sizes) + 1:) = cluster_sizes
      weights_copy(1:n) = 1.0_dp - 1.0_dp / real(ncopy, dp)
      weights_copy(n + 1:2 * n) = 1.0_dp / real(ncopy, dp)

      if (present(offset_copy)) then
         allocate(offset_copy(2 * n))
         if (present(offset)) then
            if (size(offset) /= n) then
               status = GEE_ERR_SHAPE
               return
            end if
            offset_copy(1:n) = offset
            offset_copy(n + 1:) = offset
         else
            offset_copy = 0.0_dp
         end if
      end if
      if (present(waves_copy)) then
         allocate(waves_copy(2 * n))
         if (present(waves)) then
            if (size(waves) /= n) then
               status = GEE_ERR_SHAPE
               return
            end if
            waves_copy(1:n) = waves
            waves_copy(n + 1:) = waves
         else
            pos = 0
            do g = 1, size(cluster_sizes)
               do j = 1, cluster_sizes(g)
                  waves_copy(pos + j) = j
                  waves_copy(n + pos + j) = j
               end do
               pos = pos + cluster_sizes(g)
            end do
         end if
      end if
      status = GEE_OK
   end subroutine make_relative_risk_copy

   subroutine fit_relative_risk(y, x, cluster_sizes, ncopy, corstr, result, beta_initial, alpha_initial, &
      offset, waves, tolerance, max_iterations)
      real(dp), intent(in) :: y(:) !! Binary responses coded zero or one.
      real(dp), intent(in) :: x(:, :) !! Relative-risk mean-model design matrix.
      integer, intent(in) :: cluster_sizes(:) !! Original contiguous cluster sizes.
      integer, intent(in) :: ncopy !! COPY tuning count; upstream default is 1000.
      integer, intent(in) :: corstr !! Working-correlation structure code COR_*.
      type(gee_result), intent(out) :: result !! Fitted log-relative-risk GEE result.
      real(dp), optional, intent(in) :: beta_initial(:) !! Optional starting log-risk coefficients.
      real(dp), optional, intent(in) :: alpha_initial(:) !! Optional starting working-correlation coefficients.
      real(dp), optional, intent(in) :: offset(:) !! Optional observation-level log-risk offsets.
      integer, optional, intent(in) :: waves(:) !! Optional wave identifiers.
      real(dp), optional, intent(in) :: tolerance !! Optional GEE convergence tolerance.
      integer, optional, intent(in) :: max_iterations !! Optional maximum GEE iteration count.
      real(dp), allocatable :: yc(:)
      real(dp), allocatable :: xc(:, :)
      real(dp), allocatable :: wc(:)
      real(dp), allocatable :: oc(:)
      real(dp), allocatable :: b0(:)
      integer, allocatable :: cc(:)
      integer, allocatable :: wavc(:)
      type(gee_spec) :: spec
      integer :: status

      call make_relative_risk_copy(y, x, cluster_sizes, ncopy, yc, xc, cc, wc, status, &
         offset=offset, offset_copy=oc, waves=waves, waves_copy=wavc)
      if (status /= GEE_OK) then
         result%error = status
         return
      end if
      allocate(b0(size(x, 2)))
      if (present(beta_initial)) then
         if (size(beta_initial) /= size(b0)) then
            result%error = GEE_ERR_SHAPE
            return
         end if
         b0 = beta_initial
      else
         b0 = 0.0_dp
         if (size(b0) > 0) then
            if (all(abs(x(:, 1) - 1.0_dp) <= 10.0_dp * epsilon(1.0_dp))) b0(1) = log(0.5_dp)
         end if
      end if
      spec%corstr = corstr
      spec%scale_fixed = .true.
      spec%scale_value = 1.0_dp
      allocate(spec%mean_links(1), spec%variance_codes(1), spec%scale_links(1))
      spec%mean_links = LINK_LOG
      spec%variance_codes = VAR_BINOMIAL
      spec%scale_links = LINK_IDENTITY
      if (present(tolerance)) spec%tolerance = tolerance
      if (present(max_iterations)) spec%max_iterations = max_iterations
      if (present(alpha_initial)) then
         call fit_geese(yc, xc, cc, b0, spec, result, alpha_initial=alpha_initial, offset=oc, &
            weights=wc, waves=wavc)
      else
         call fit_geese(yc, xc, cc, b0, spec, result, offset=oc, weights=wc, waves=wavc)
      end if
   end subroutine fit_relative_risk

end module geepack_relative_risk
