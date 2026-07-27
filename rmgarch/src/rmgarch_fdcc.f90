! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_fdcc
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : covariance_matrix, outer_product, normalize_covariance, &
      make_positive_definite
   use rmgarch_distributions, only : multivariate_normal_logpdf, random_multivariate_normal
   use rmgarch_optimizer, only : optimizer_result, nelder_mead
   use rmgarch_types, only : fdcc_fit_result
   implicit none
   private

   type :: fdcc_context
      real(dp), allocatable :: z(:,:)
      integer, allocatable :: group_index(:)
      integer :: ngroups = 0
   end type fdcc_context

   public :: fdcc_filter, fdcc_log_likelihood, fit_fdcc11
   public :: fdcc_forecast, simulate_fdcc, fdcc_parameter_matrices

contains

   subroutine fdcc_filter(z, a, b, c, q, r, loglikelihoods, qbar_out, valid)
      real(dp), intent(in) :: z(:,:)
      real(dp), intent(in) :: a(:,:), b(:,:), c(:,:)
      real(dp), intent(out) :: q(size(z,2),size(z,2),size(z,1))
      real(dp), intent(out) :: r(size(z,2),size(z,2),size(z,1))
      real(dp), intent(out) :: loglikelihoods(size(z,1))
      real(dp), intent(out), optional :: qbar_out(size(z,2),size(z,2))
      logical, intent(out), optional :: valid
      real(dp) :: qbar(size(z,2),size(z,2)), intercept(size(z,2),size(z,2))
      real(dp) :: aa(size(z,2),size(z,2)), bb(size(z,2),size(z,2))
      real(dp) :: zero_mean(size(z,2))
      integer :: n, m, p, qq, maxlag, t, i
      logical :: ok, step_ok

      n = size(z,1)
      m = size(z,2)
      p = size(a,2)
      qq = size(b,2)
      ok = size(a,1) == m .and. size(b,1) == m .and. size(c,1) == m .and. &
         size(c,2) == m .and. n > max(1,max(p,qq))
      ok = ok .and. all(a >= 0.0_dp) .and. all(b >= 0.0_dp)
      q = 0.0_dp
      r = 0.0_dp
      loglikelihoods = -huge(1.0_dp)
      if (.not. ok) then
         if (present(valid)) valid = .false.
         return
      end if
      qbar = make_positive_definite(covariance_matrix(z),1.0e-10_dp)
      intercept = c*qbar
      zero_mean = 0.0_dp
      maxlag = max(1,max(p,qq))
      do t = 1, min(maxlag,n)
         q(:,:,t) = qbar
         r(:,:,t) = normalize_covariance(qbar)
         loglikelihoods(t) = multivariate_normal_logpdf(z(t,:),zero_mean,r(:,:,t),step_ok)
         ok = ok .and. step_ok
      end do
      do t = maxlag+1, n
         q(:,:,t) = intercept
         do i = 1, p
            aa = outer_product(a(:,i),a(:,i))
            q(:,:,t) = q(:,:,t)+aa*outer_product(z(t-i,:),z(t-i,:))
         end do
         do i = 1, qq
            bb = outer_product(b(:,i),b(:,i))
            q(:,:,t) = q(:,:,t)+bb*q(:,:,t-i)
         end do
         q(:,:,t) = make_positive_definite(q(:,:,t),1.0e-10_dp)
         r(:,:,t) = normalize_covariance(q(:,:,t))
         loglikelihoods(t) = multivariate_normal_logpdf(z(t,:),zero_mean,r(:,:,t),step_ok)
         ok = ok .and. step_ok
      end do
      if (present(qbar_out)) qbar_out = qbar
      if (present(valid)) valid = ok
   end subroutine fdcc_filter

   function fdcc_log_likelihood(z, a, b, c, valid) result(value)
      real(dp), intent(in) :: z(:,:), a(:,:), b(:,:), c(:,:)
      logical, intent(out), optional :: valid
      real(dp) :: value
      real(dp), allocatable :: q(:,:,:), r(:,:,:), ll(:)
      logical :: ok

      allocate(q(size(z,2),size(z,2),size(z,1)),r(size(z,2),size(z,2),size(z,1)),ll(size(z,1)))
      call fdcc_filter(z,a,b,c,q,r,ll,valid=ok)
      if (ok) then
         value = sum(ll)
      else
         value = -huge(1.0_dp)
      end if
      if (present(valid)) valid = ok
   end function fdcc_log_likelihood

   function fit_fdcc11(z, group_index, max_iterations) result(fit)
      real(dp), intent(in) :: z(:,:)
      integer, intent(in) :: group_index(:)
      integer, intent(in), optional :: max_iterations
      type(fdcc_fit_result) :: fit
      type(fdcc_context) :: context
      type(optimizer_result) :: opt
      real(dp), allocatable :: x0(:), ga(:), gb(:)
      integer :: ngroups, g, maxit, k
      logical :: ok

      ok = size(group_index) == size(z,2) .and. size(group_index) > 0
      if (ok) ok = minval(group_index) == 1
      if (ok) then
         ngroups = maxval(group_index)
         do g = 1, ngroups
            ok = ok .and. any(group_index == g)
         end do
      else
         ngroups = 0
      end if
      if (.not. ok) then
         fit%status = 3
         fit%message = 'group_index must be contiguous from one'
         return
      end if

      context%ngroups = ngroups
      allocate(context%z(size(z,1),size(z,2)),context%group_index(size(group_index)))
      context%z = z
      context%group_index = group_index
      allocate(x0(2*ngroups))
      do g = 1, ngroups
         x0(2*g-1) = log(0.01_dp/0.18_dp)
         x0(2*g) = log(0.81_dp/0.18_dp)
      end do
      maxit = 800
      if (present(max_iterations)) maxit = max_iterations
      opt = nelder_mead(fdcc_objective,x0,context,step=0.15_dp,tolerance=1.0e-8_dp, &
         max_iterations=maxit)
      call decode_fdcc_parameters(opt%x,ga,gb)
      allocate(fit%group_alpha(ngroups),fit%group_beta(ngroups),fit%group_index(size(group_index)))
      fit%group_alpha = ga
      fit%group_beta = gb
      fit%group_index = group_index
      call fdcc_parameter_matrices(group_index,ga,gb,fit%a,fit%b,fit%c,ok)
      allocate(fit%qbar(size(z,2),size(z,2)))
      allocate(fit%q(size(z,2),size(z,2),size(z,1)),fit%r(size(z,2),size(z,2),size(z,1)))
      allocate(fit%loglikelihoods(size(z,1)))
      if (ok) call fdcc_filter(z,fit%a,fit%b,fit%c,fit%q,fit%r,fit%loglikelihoods,fit%qbar,ok)
      if (ok) fit%log_likelihood = sum(fit%loglikelihoods)
      fit%iterations = opt%iterations
      fit%status = merge(0,2,ok)
      if (opt%status /= 0 .and. fit%status == 0) fit%status = opt%status
      if (fit%status == 0) then
         fit%message = 'converged'
      else if (ok) then
         fit%message = 'valid fit; optimizer stopped before convergence'
      else
         fit%message = 'invalid fitted recursion'
      end if
      k = size(opt%x)
      if (ok) then
         fit%aic = -2.0_dp*fit%log_likelihood+2.0_dp*real(k,dp)
         fit%bic = -2.0_dp*fit%log_likelihood+log(real(size(z,1),dp))*real(k,dp)
      end if
   end function fit_fdcc11

   subroutine fdcc_parameter_matrices(group_index, group_alpha, group_beta, a, b, c, valid)
      integer, intent(in) :: group_index(:)
      real(dp), intent(in) :: group_alpha(:), group_beta(:)
      real(dp), allocatable, intent(out) :: a(:,:), b(:,:), c(:,:)
      logical, intent(out), optional :: valid
      integer :: i, m
      real(dp) :: aa(size(group_index),size(group_index)), bb(size(group_index),size(group_index))
      logical :: ok

      m = size(group_index)
      ok = size(group_alpha) == size(group_beta) .and. m > 0
      if (ok) ok = minval(group_index) >= 1 .and. maxval(group_index) <= size(group_alpha)
      allocate(a(m,1),b(m,1),c(m,m))
      if (ok) then
         do i = 1, m
            a(i,1) = group_alpha(group_index(i))
            b(i,1) = group_beta(group_index(i))
         end do
         aa = outer_product(a(:,1),a(:,1))
         bb = outer_product(b(:,1),b(:,1))
         c = 1.0_dp-aa-bb
         ok = minval(c) > 0.0_dp
      else
         a = 0.0_dp
         b = 0.0_dp
         c = 0.0_dp
      end if
      if (present(valid)) valid = ok
   end subroutine fdcc_parameter_matrices

   subroutine fdcc_forecast(a, b, c, qbar, last_q, last_z, horizons, q_forecast, r_forecast)
      real(dp), intent(in) :: a(:,:), b(:,:), c(:,:), qbar(:,:), last_q(:,:), last_z(:)
      integer, intent(in) :: horizons
      real(dp), intent(out) :: q_forecast(size(qbar,1),size(qbar,2),horizons)
      real(dp), intent(out) :: r_forecast(size(qbar,1),size(qbar,2),horizons)
      real(dp) :: aa(size(qbar,1),size(qbar,2)), bb(size(qbar,1),size(qbar,2))
      integer :: h

      if (horizons <= 0) return
      aa = outer_product(a(:,1),a(:,1))
      bb = outer_product(b(:,1),b(:,1))
      q_forecast(:,:,1) = c*qbar+aa*outer_product(last_z,last_z)+bb*last_q
      q_forecast(:,:,1) = make_positive_definite(q_forecast(:,:,1),1.0e-10_dp)
      r_forecast(:,:,1) = normalize_covariance(q_forecast(:,:,1))
      do h = 2, horizons
         q_forecast(:,:,h) = c*qbar+aa*qbar+bb*q_forecast(:,:,h-1)
         q_forecast(:,:,h) = make_positive_definite(q_forecast(:,:,h),1.0e-10_dp)
         r_forecast(:,:,h) = normalize_covariance(q_forecast(:,:,h))
      end do
   end subroutine fdcc_forecast

   subroutine simulate_fdcc(nobs, a, b, c, qbar, z, q, r, burn)
      integer, intent(in) :: nobs
      real(dp), intent(in) :: a(:,:), b(:,:), c(:,:), qbar(:,:)
      real(dp), intent(out) :: z(nobs,size(qbar,1))
      real(dp), intent(out) :: q(size(qbar,1),size(qbar,2),nobs)
      real(dp), intent(out) :: r(size(qbar,1),size(qbar,2),nobs)
      integer, intent(in), optional :: burn
      real(dp), allocatable :: zall(:,:), qall(:,:,:), rall(:,:,:)
      real(dp) :: aa(size(qbar,1),size(qbar,2)), bb(size(qbar,1),size(qbar,2))
      real(dp) :: zero_mean(size(qbar,1)), draw(size(qbar,1))
      integer :: nburn, ntotal, t, m
      logical :: ok

      nburn = 300
      if (present(burn)) nburn = max(0,burn)
      ntotal = nobs+nburn
      m = size(qbar,1)
      allocate(zall(ntotal,m),qall(m,m,ntotal),rall(m,m,ntotal))
      aa = outer_product(a(:,1),a(:,1))
      bb = outer_product(b(:,1),b(:,1))
      zero_mean = 0.0_dp
      qall(:,:,1) = make_positive_definite(qbar)
      rall(:,:,1) = normalize_covariance(qall(:,:,1))
      do t = 1, ntotal
         if (t > 1) then
            qall(:,:,t) = c*qbar+aa*outer_product(zall(t-1,:),zall(t-1,:))+bb*qall(:,:,t-1)
            qall(:,:,t) = make_positive_definite(qall(:,:,t),1.0e-10_dp)
            rall(:,:,t) = normalize_covariance(qall(:,:,t))
         end if
         call random_multivariate_normal(zero_mean,rall(:,:,t),draw,ok)
         if (ok) then
            zall(t,:) = draw
         else
            zall(t,:) = 0.0_dp
         end if
      end do
      z = zall(nburn+1:ntotal,:)
      q = qall(:,:,nburn+1:ntotal)
      r = rall(:,:,nburn+1:ntotal)
   end subroutine simulate_fdcc

   function fdcc_objective(x, generic_context) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value, ll
      real(dp), allocatable :: ga(:), gb(:), a(:,:), b(:,:), c(:,:)
      logical :: ok

      select type (context => generic_context)
      type is (fdcc_context)
         call decode_fdcc_parameters(x,ga,gb)
         call fdcc_parameter_matrices(context%group_index,ga,gb,a,b,c,ok)
         if (ok) ll = fdcc_log_likelihood(context%z,a,b,c,ok)
         if (ok) then
            value = -ll
         else
            value = huge(1.0_dp)/100.0_dp
         end if
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function fdcc_objective

   subroutine decode_fdcc_parameters(x, group_alpha, group_beta)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: group_alpha(:), group_beta(:)
      real(dp) :: exa, exb, denom
      integer :: g, ngroups

      ngroups = size(x)/2
      allocate(group_alpha(ngroups),group_beta(ngroups))
      do g = 1, ngroups
         exa = exp(max(-30.0_dp,min(30.0_dp,x(2*g-1))))
         exb = exp(max(-30.0_dp,min(30.0_dp,x(2*g))))
         denom = 1.0_dp+exa+exb
         group_alpha(g) = sqrt(0.999_dp*exa/denom)
         group_beta(g) = sqrt(0.999_dp*exb/denom)
      end do
   end subroutine decode_fdcc_parameters

end module rmgarch_fdcc
