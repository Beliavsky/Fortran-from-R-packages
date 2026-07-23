! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_cgarch
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : normal_quantile
   use rmgarch_types, only : copula_garch_fit_result, dist_gaussian
   use rmgarch_univariate, only : fit_marginal_garch11, filter_garch11
   use rmgarch_copula, only : fit_copula, parametric_uniform_transform, &
      copula_score_transform, simulate_static_copula, simulate_dynamic_copula
   use rmgarch_dcc, only : dcc_filter
   implicit none
   private

   public :: fit_copula_garch11, filter_copula_garch11
   public :: simulate_fitted_copula_garch11

contains

   function fit_copula_garch11(data, copula_distribution, time_varying, p, q, g, &
      shape, estimate_shape, max_iterations) result(fit)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: copula_distribution, p, q, g, max_iterations
      logical, intent(in), optional :: time_varying, estimate_shape
      real(dp), intent(in), optional :: shape
      type(copula_garch_fit_result) :: fit
      integer, allocatable :: marginal_distributions(:)
      real(dp), allocatable :: marginal_shapes(:)
      integer :: n, m, dist, pp, qq, gg, maxit, j
      real(dp) :: nu
      logical :: dynamic, fit_shape

      n = size(data,1)
      m = size(data,2)
      dist = dist_gaussian
      if (present(copula_distribution)) dist = copula_distribution
      pp = 1
      qq = 1
      gg = 0
      if (present(p)) pp = max(0,p)
      if (present(q)) qq = max(0,q)
      if (present(g)) gg = max(0,g)
      dynamic = .true.
      if (present(time_varying)) dynamic = time_varying
      nu = 8.0_dp
      if (present(shape)) nu = shape
      fit_shape = .false.
      if (present(estimate_shape)) fit_shape = estimate_shape
      maxit = 500
      if (present(max_iterations)) maxit = max(1,max_iterations)
      if (n <= 2 .or. m < 1) then
         fit%status = 2
         return
      end if

      allocate(fit%margins(m),fit%sigma(n,m),fit%standardized(n,m), &
         fit%uniforms(n,m),fit%scores(n,m))
      allocate(marginal_distributions(m),marginal_shapes(m))
      call fit_marginal_garch11(data,fit%margins,fit%standardized,fit%sigma,maxit)
      marginal_distributions = dist_gaussian
      marginal_shapes = 0.0_dp
      call parametric_uniform_transform(fit%standardized,marginal_distributions, &
         marginal_shapes,fit%uniforms)
      call copula_score_transform(fit%uniforms,dist,fit%scores,nu)
      fit%copula = fit_copula(fit%scores,distribution=dist,time_varying=dynamic, &
         p=pp,q=qq,g=gg,shape=nu,estimate_shape=fit_shape,max_iterations=maxit)
      fit%status = fit%copula%status
      do j = 1, m
         fit%status = max(fit%status,fit%margins(j)%status)
      end do
   end function fit_copula_garch11

   subroutine filter_copula_garch11(data, fit, sigma, standardized, uniforms, scores, q, r, valid)
      real(dp), intent(in) :: data(:,:)
      type(copula_garch_fit_result), intent(in) :: fit
      real(dp), intent(out) :: sigma(size(data,1),size(data,2))
      real(dp), intent(out) :: standardized(size(data,1),size(data,2))
      real(dp), intent(out) :: uniforms(size(data,1),size(data,2))
      real(dp), intent(out) :: scores(size(data,1),size(data,2))
      real(dp), intent(out) :: q(size(data,2),size(data,2),size(data,1))
      real(dp), intent(out) :: r(size(data,2),size(data,2),size(data,1))
      logical, intent(out), optional :: valid
      integer, allocatable :: marginal_distributions(:)
      real(dp), allocatable :: marginal_shapes(:), residuals(:), ll(:)
      integer :: j, t, m, n
      logical :: ok, margin_ok
      real(dp) :: dummy_ll

      n = size(data,1)
      m = size(data,2)
      sigma = 0.0_dp
      standardized = 0.0_dp
      uniforms = 0.5_dp
      scores = 0.0_dp
      q = 0.0_dp
      r = 0.0_dp
      ok = allocated(fit%margins) .and. size(fit%margins) == m
      if (.not. ok) then
         if (present(valid)) valid = .false.
         return
      end if
      allocate(marginal_distributions(m),marginal_shapes(m),residuals(n),ll(n))
      do j = 1, m
         call filter_garch11(data(:,j),fit%margins(j)%mean,fit%margins(j)%omega, &
            fit%margins(j)%alpha,fit%margins(j)%beta,residuals,sigma(:,j), &
            standardized(:,j),dummy_ll,margin_ok)
         ok = ok .and. margin_ok
      end do
      marginal_distributions = dist_gaussian
      marginal_shapes = 0.0_dp
      call parametric_uniform_transform(standardized,marginal_distributions,marginal_shapes,uniforms)
      call copula_score_transform(uniforms,fit%copula%distribution,scores,fit%copula%shape)
      if (fit%copula%time_varying .and. allocated(fit%copula%dcc%spec%alpha)) then
         call dcc_filter(scores,fit%copula%dcc%spec,q,r,ll,valid=margin_ok)
         ok = ok .and. margin_ok
      else if (allocated(fit%copula%correlation)) then
         do t = 1, n
            q(:,:,t) = fit%copula%correlation
            r(:,:,t) = fit%copula%correlation
         end do
      else
         ok = .false.
      end if
      if (present(valid)) valid = ok
   end subroutine filter_copula_garch11

   subroutine simulate_fitted_copula_garch11(nobs, fit, returns, burn, valid)
      integer, intent(in) :: nobs
      type(copula_garch_fit_result), intent(in) :: fit
      real(dp), intent(out) :: returns(nobs,size(fit%margins))
      integer, intent(in), optional :: burn
      logical, intent(out), optional :: valid
      real(dp), allocatable :: scores(:,:), uniforms(:,:), q(:,:,:), r(:,:,:), h(:), residual(:)
      real(dp) :: innovation
      integer :: nburn, ntotal, t, j, m
      logical :: ok

      returns = 0.0_dp
      ok = nobs > 0 .and. allocated(fit%margins) .and. allocated(fit%copula%correlation)
      if (.not. ok) then
         if (present(valid)) valid = .false.
         return
      end if
      nburn = 300
      if (present(burn)) nburn = max(0,burn)
      ntotal = nobs+nburn
      m = size(fit%margins)
      allocate(scores(ntotal,m),uniforms(ntotal,m),q(m,m,ntotal),r(m,m,ntotal),h(m),residual(m))
      if (fit%copula%time_varying .and. allocated(fit%copula%dcc%qbar)) then
         call simulate_dynamic_copula(ntotal,fit%copula%dcc%spec,fit%copula%dcc%qbar, &
            scores,uniforms,q,r,burn=0)
      else
         call simulate_static_copula(ntotal,fit%copula%correlation,fit%copula%distribution, &
            scores,uniforms,fit%copula%shape,ok)
      end if
      do j = 1, m
         h(j) = fit%margins(j)%omega/max(1.0_dp-fit%margins(j)%alpha-fit%margins(j)%beta,1.0e-6_dp)
         residual(j) = 0.0_dp
      end do
      do t = 1, ntotal
         do j = 1, m
            if (t > 1) h(j) = fit%margins(j)%omega+fit%margins(j)%alpha* &
               residual(j)**2+fit%margins(j)%beta*h(j)
            innovation = normal_quantile(uniforms(t,j))
            residual(j) = sqrt(max(h(j),1.0e-12_dp))*innovation
            if (t > nburn) returns(t-nburn,j) = fit%margins(j)%mean+residual(j)
         end do
      end do
      if (present(valid)) valid = ok
   end subroutine simulate_fitted_copula_garch11

end module rmgarch_cgarch
