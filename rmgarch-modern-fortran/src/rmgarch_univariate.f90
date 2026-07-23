! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_univariate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : pi, sample_mean, sample_variance
   use rmgarch_optimizer, only : optimizer_result, nelder_mead
   use rmgarch_types, only : univariate_garch_fit
   implicit none
   private

   type :: garch_context
      real(dp), allocatable :: y(:)
   end type garch_context

   public :: filter_garch11, fit_garch11, fit_marginal_garch11

contains

   subroutine filter_garch11(y, mean, omega, alpha, beta, residuals, sigma, standardized, log_likelihood, valid)
      real(dp), intent(in) :: y(:), mean, omega, alpha, beta
      real(dp), intent(out) :: residuals(size(y)), sigma(size(y)), standardized(size(y))
      real(dp), intent(out), optional :: log_likelihood
      logical, intent(out), optional :: valid
      real(dp) :: h(size(y)), v0, ll
      integer :: t, n
      logical :: ok

      n = size(y)
      ok = n >= 2 .and. omega > 0.0_dp .and. alpha >= 0.0_dp .and. beta >= 0.0_dp .and. alpha+beta < 1.0_dp
      if (.not. ok) then
         residuals = 0.0_dp; sigma = huge(1.0_dp); standardized = 0.0_dp
         if (present(log_likelihood)) log_likelihood = -huge(1.0_dp)
         if (present(valid)) valid = .false.
         return
      end if
      residuals = y-mean
      v0 = max(sample_variance(y),omega/max(1.0_dp-alpha-beta,1.0e-6_dp))
      h(1) = max(v0,1.0e-12_dp)
      do t = 2, n
         h(t) = omega+alpha*residuals(t-1)**2+beta*h(t-1)
         h(t) = max(h(t),1.0e-12_dp)
      end do
      sigma = sqrt(h)
      standardized = residuals/sigma
      ll = -0.5_dp*sum(log(2.0_dp*pi)+log(h)+standardized**2)
      if (present(log_likelihood)) log_likelihood = ll
      if (present(valid)) valid = ieee_is_finite(ll)
   end subroutine filter_garch11

   function fit_garch11(y, max_iterations) result(fit)
      real(dp), intent(in) :: y(:)
      integer, intent(in), optional :: max_iterations
      type(univariate_garch_fit) :: fit
      type(garch_context) :: context
      type(optimizer_result) :: opt
      real(dp) :: x0(4), v, ea, eb, denom
      integer :: maxit
      logical :: ok

      maxit = 500
      if (present(max_iterations)) maxit = max_iterations
      allocate(context%y(size(y))); context%y = y
      v = max(sample_variance(y),1.0e-8_dp)
      x0(1) = sample_mean(y)
      x0(2) = log(max(0.05_dp*v,1.0e-10_dp))
      x0(3) = log(0.05_dp/0.05_dp)
      x0(4) = log(0.90_dp/0.05_dp)
      opt = nelder_mead(garch_objective,x0,context,step=0.15_dp,tolerance=1.0e-8_dp,max_iterations=maxit)
      ea = exp(max(-30.0_dp,min(30.0_dp,opt%x(3))))
      eb = exp(max(-30.0_dp,min(30.0_dp,opt%x(4))))
      denom = 1.0_dp+ea+eb
      fit%mean = opt%x(1)
      fit%omega = exp(max(-40.0_dp,min(20.0_dp,opt%x(2))))
      fit%alpha = 0.999_dp*ea/denom
      fit%beta = 0.999_dp*eb/denom
      fit%iterations = opt%iterations
      fit%status = opt%status
      allocate(fit%residuals(size(y)),fit%sigma(size(y)),fit%standardized(size(y)))
      call filter_garch11(y,fit%mean,fit%omega,fit%alpha,fit%beta,fit%residuals,fit%sigma, &
         fit%standardized,fit%log_likelihood,ok)
      if (.not. ok) fit%status = 2
   end function fit_garch11

   subroutine fit_marginal_garch11(data, fits, standardized, sigma, max_iterations)
      real(dp), intent(in) :: data(:,:)
      type(univariate_garch_fit), intent(out) :: fits(size(data,2))
      real(dp), intent(out) :: standardized(size(data,1),size(data,2)), sigma(size(data,1),size(data,2))
      integer, intent(in), optional :: max_iterations
      integer :: j
      do j = 1, size(data,2)
         fits(j) = fit_garch11(data(:,j),max_iterations)
         standardized(:,j) = fits(j)%standardized
         sigma(:,j) = fits(j)%sigma
      end do
   end subroutine fit_marginal_garch11

   function garch_objective(x, generic_context) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value, omega, alpha, beta, ea, eb, denom, ll
      real(dp), allocatable :: residuals(:), sigma(:), z(:)
      logical :: ok
      select type (context => generic_context)
      type is (garch_context)
         ea = exp(max(-30.0_dp,min(30.0_dp,x(3))))
         eb = exp(max(-30.0_dp,min(30.0_dp,x(4))))
         denom = 1.0_dp+ea+eb
         omega = exp(max(-40.0_dp,min(20.0_dp,x(2))))
         alpha = 0.999_dp*ea/denom
         beta = 0.999_dp*eb/denom
         allocate(residuals(size(context%y)),sigma(size(context%y)),z(size(context%y)))
         call filter_garch11(context%y,x(1),omega,alpha,beta,residuals,sigma,z,ll,ok)
         if (ok) then
            value = -ll
         else
            value = huge(1.0_dp)/100.0_dp
         end if
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function garch_objective

end module rmgarch_univariate
