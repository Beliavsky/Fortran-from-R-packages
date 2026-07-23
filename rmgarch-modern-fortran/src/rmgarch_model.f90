! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_model
   use rmgarch_kinds, only : dp
   use rmgarch_types, only : multivariate_garch_fit, dist_gaussian
   use rmgarch_univariate, only : fit_marginal_garch11
   use rmgarch_dcc, only : fit_dcc, fit_dcc11
   implicit none
   private

   public :: fit_two_step_dcc, fit_two_step_dcc_general, conditional_covariance

contains

   function fit_two_step_dcc(data, asymmetric, distribution, shape, max_iterations) result(fit)
      real(dp), intent(in) :: data(:,:)
      logical, intent(in), optional :: asymmetric
      integer, intent(in), optional :: distribution, max_iterations
      real(dp), intent(in), optional :: shape
      type(multivariate_garch_fit) :: fit
      integer :: n, m, t, i, j, dist
      logical :: asym
      real(dp) :: nu

      n = size(data,1); m = size(data,2)
      allocate(fit%margins(m),fit%sigma(n,m),fit%standardized(n,m))
      call fit_marginal_garch11(data,fit%margins,fit%standardized,fit%sigma,max_iterations)
      asym = .false.; if (present(asymmetric)) asym = asymmetric
      dist = dist_gaussian; if (present(distribution)) dist = distribution
      nu = 8.0_dp; if (present(shape)) nu = shape
      fit%dcc = fit_dcc11(fit%standardized,asym,dist,nu,max_iterations)
      allocate(fit%covariance(m,m,n))
      do t = 1, n
         do j = 1, m
            do i = 1, m
               fit%covariance(i,j,t) = fit%sigma(t,i)*fit%dcc%r(i,j,t)*fit%sigma(t,j)
            end do
         end do
      end do
      fit%status = fit%dcc%status
      do i = 1, m
         if (fit%margins(i)%status /= 0) fit%status = max(fit%status,fit%margins(i)%status)
      end do
   end function fit_two_step_dcc

   subroutine conditional_covariance(sigma, correlation, covariance)
      real(dp), intent(in) :: sigma(:,:), correlation(:,:,:)
      real(dp), intent(out) :: covariance(size(sigma,2),size(sigma,2),size(sigma,1))
      integer :: t, i, j
      do t = 1, size(sigma,1)
         do j = 1, size(sigma,2)
            do i = 1, size(sigma,2)
               covariance(i,j,t) = sigma(t,i)*correlation(i,j,t)*sigma(t,j)
            end do
         end do
      end do
   end subroutine conditional_covariance

   function fit_two_step_dcc_general(data, p, q, g, distribution, shape, estimate_shape, &
      max_iterations) result(fit)
      !! Two-step Gaussian marginal GARCH(1,1) plus general DCC(p,q,g).
      real(dp), intent(in) :: data(:,:)
      integer, intent(in) :: p, q, g
      integer, intent(in), optional :: distribution, max_iterations
      real(dp), intent(in), optional :: shape
      logical, intent(in), optional :: estimate_shape
      type(multivariate_garch_fit) :: fit
      integer :: n, m, t, i, j, dist, maxit
      real(dp) :: nu
      logical :: fit_shape

      n = size(data,1)
      m = size(data,2)
      dist = dist_gaussian
      if (present(distribution)) dist = distribution
      nu = 8.0_dp
      if (present(shape)) nu = shape
      fit_shape = .false.
      if (present(estimate_shape)) fit_shape = estimate_shape
      maxit = 500
      if (present(max_iterations)) maxit = max(1,max_iterations)
      allocate(fit%margins(m),fit%sigma(n,m),fit%standardized(n,m))
      call fit_marginal_garch11(data,fit%margins,fit%standardized,fit%sigma,maxit)
      fit%dcc = fit_dcc(fit%standardized,p=p,q=q,g=g,distribution=dist,shape=nu, &
         estimate_shape=fit_shape,max_iterations=maxit)
      allocate(fit%covariance(m,m,n))
      if (allocated(fit%dcc%r)) then
         do t = 1, n
            do j = 1, m
               do i = 1, m
                  fit%covariance(i,j,t) = fit%sigma(t,i)*fit%dcc%r(i,j,t)*fit%sigma(t,j)
               end do
            end do
         end do
      else
         fit%covariance = 0.0_dp
      end if
      fit%status = fit%dcc%status
      do i = 1, m
         fit%status = max(fit%status,fit%margins(i)%status)
      end do
   end function fit_two_step_dcc_general

end module rmgarch_model
