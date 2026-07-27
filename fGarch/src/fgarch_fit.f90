! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

module fgarch_fit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fgarch_kinds, only : dp
   use fgarch_distributions, only : dist_norm, dist_snorm, dist_std, dist_sstd, &
      dist_ged, dist_sged, distribution_pdf
   use fgarch_optimizer, only : optimizer_result, nelder_mead
   use fgarch_types, only : garch_spec, garch_fit_result, distribution_fit_result, &
      make_garch_spec, model_garch, model_aparch
   use fgarch_models, only : garch_log_likelihood
   implicit none
   private

   type :: dist_fit_context
      real(dp), allocatable :: y(:)
      integer :: kind = dist_norm
      logical :: fit_shape = .false.
      logical :: fit_skew = .false.
   end type dist_fit_context

   type :: garch_fit_context
      real(dp), allocatable :: y(:)
      integer :: model = model_garch
      integer :: cond_dist = dist_norm
      logical :: fit_mean = .true.
      logical :: fit_shape = .false.
      logical :: fit_skew = .false.
   end type garch_fit_context

   public :: fit_distribution, fit_garch11, fit_aparch11

contains

   function fit_distribution(y, kind, fit_shape, fit_skew, max_iterations) result(fit)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: kind
      logical, intent(in), optional :: fit_shape, fit_skew
      integer, intent(in), optional :: max_iterations
      type(distribution_fit_result) :: fit
      type(dist_fit_context) :: context
      type(optimizer_result) :: opt
      real(dp), allocatable :: x0(:)
      real(dp) :: mean0, sd0
      logical :: use_shape, use_skew
      integer :: npar, idx, maxit

      use_shape = kind == dist_std .or. kind == dist_sstd .or. kind == dist_ged .or. kind == dist_sged
      if (present(fit_shape)) use_shape = fit_shape .and. use_shape
      use_skew = kind == dist_snorm .or. kind == dist_sstd .or. kind == dist_sged
      if (present(fit_skew)) use_skew = fit_skew .and. use_skew
      maxit = 1500
      if (present(max_iterations)) maxit = max_iterations

      context%y = y
      context%kind = kind
      context%fit_shape = use_shape
      context%fit_skew = use_skew
      mean0 = sum(y)/real(size(y),dp)
      sd0 = sqrt(max(sum((y-mean0)**2)/real(max(1,size(y)-1),dp),1.0e-12_dp))
      npar = 2+merge(1,0,use_shape)+merge(1,0,use_skew)
      allocate(x0(npar))
      x0(1) = mean0
      x0(2) = log(sd0)
      idx = 3
      if (use_shape) then
         if (kind == dist_std .or. kind == dist_sstd) then
            x0(idx) = log(5.0_dp-2.05_dp)
         else
            x0(idx) = log(2.0_dp-0.2_dp)
         end if
         idx = idx+1
      end if
      if (use_skew) x0(idx) = 0.0_dp

      opt = nelder_mead(distribution_objective,x0,context,step=0.15_dp, &
                        tolerance=1.0e-8_dp,max_iterations=maxit)
      call decode_distribution_parameters(opt%x,context,fit%mean,fit%sd,fit%shape,fit%skew)
      fit%log_likelihood = -opt%objective
      fit%iterations = opt%iterations
      fit%status = opt%status
   end function fit_distribution

   function distribution_objective(x, generic_context) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value
      real(dp) :: mean, sd, shape, skew, density
      integer :: i

      select type (context => generic_context)
      type is (dist_fit_context)
         call decode_distribution_parameters(x,context,mean,sd,shape,skew)
         value = 0.0_dp
         do i = 1, size(context%y)
            density = distribution_pdf((context%y(i)-mean)/sd,context%kind,shape,skew)/sd
            if (density <= tiny(1.0_dp) .or. .not. finite_value(density)) then
               value = huge(1.0_dp)/100.0_dp
               return
            end if
            value = value-log(density)
         end do
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function distribution_objective

   subroutine decode_distribution_parameters(x, context, mean, sd, shape, skew)
      real(dp), intent(in) :: x(:)
      type(dist_fit_context), intent(in) :: context
      real(dp), intent(out) :: mean, sd, shape, skew
      integer :: idx

      mean = x(1)
      sd = exp(max(-30.0_dp,min(30.0_dp,x(2))))
      shape = 5.0_dp
      skew = 1.0_dp
      idx = 3
      if (context%fit_shape) then
         if (context%kind == dist_std .or. context%kind == dist_sstd) then
            shape = 2.05_dp+exp(max(-20.0_dp,min(10.0_dp,x(idx))))
         else
            shape = 0.2_dp+exp(max(-20.0_dp,min(5.0_dp,x(idx))))
         end if
         idx = idx+1
      end if
      if (context%fit_skew) skew = exp(max(-5.0_dp,min(5.0_dp,x(idx))))
   end subroutine decode_distribution_parameters

   function fit_garch11(y, cond_dist, fit_mean, fit_shape, fit_skew, max_iterations) result(fit)
      real(dp), intent(in) :: y(:)
      integer, intent(in), optional :: cond_dist
      logical, intent(in), optional :: fit_mean, fit_shape, fit_skew
      integer, intent(in), optional :: max_iterations
      type(garch_fit_result) :: fit

      fit = fit_garch_core(y,model_garch,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations)
   end function fit_garch11

   function fit_aparch11(y, cond_dist, fit_mean, fit_shape, fit_skew, max_iterations) result(fit)
      real(dp), intent(in) :: y(:)
      integer, intent(in), optional :: cond_dist
      logical, intent(in), optional :: fit_mean, fit_shape, fit_skew
      integer, intent(in), optional :: max_iterations
      type(garch_fit_result) :: fit

      fit = fit_garch_core(y,model_aparch,cond_dist,fit_mean,fit_shape,fit_skew,max_iterations)
   end function fit_aparch11

   function fit_garch_core(y, model, cond_dist, fit_mean, fit_shape, fit_skew, max_iterations) result(fit)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: model
      integer, intent(in), optional :: cond_dist
      logical, intent(in), optional :: fit_mean, fit_shape, fit_skew
      integer, intent(in), optional :: max_iterations
      type(garch_fit_result) :: fit
      type(garch_fit_context) :: context
      type(optimizer_result) :: opt
      real(dp), allocatable :: x0(:)
      real(dp) :: mean0, var0
      integer :: kind, npar, idx, maxit, k
      logical :: use_mean, use_shape, use_skew

      kind = dist_norm
      if (present(cond_dist)) kind = cond_dist
      use_mean = .true.
      if (present(fit_mean)) use_mean = fit_mean
      use_shape = kind == dist_std .or. kind == dist_sstd .or. kind == dist_ged .or. kind == dist_sged
      if (present(fit_shape)) use_shape = fit_shape .and. use_shape
      use_skew = kind == dist_snorm .or. kind == dist_sstd .or. kind == dist_sged
      if (present(fit_skew)) use_skew = fit_skew .and. use_skew
      maxit = 2500
      if (present(max_iterations)) maxit = max_iterations

      context%y = y
      context%model = model
      context%cond_dist = kind
      context%fit_mean = use_mean
      context%fit_shape = use_shape
      context%fit_skew = use_skew
      mean0 = sum(y)/real(size(y),dp)
      var0 = max(sum((y-mean0)**2)/real(max(1,size(y)-1),dp),1.0e-10_dp)

      npar = merge(1,0,use_mean)+3+merge(2,0,model==model_aparch)+ &
             merge(1,0,use_shape)+merge(1,0,use_skew)
      allocate(x0(npar))
      idx = 1
      if (use_mean) then
         x0(idx) = mean0
         idx = idx+1
      end if
      x0(idx) = log(max(0.05_dp*var0,1.0e-12_dp)); idx = idx+1
      x0(idx) = log(0.08_dp/(1.0_dp-0.08_dp-0.88_dp)); idx = idx+1
      x0(idx) = log(0.88_dp/(1.0_dp-0.08_dp-0.88_dp)); idx = idx+1
      if (model == model_aparch) then
         x0(idx) = 0.0_dp; idx = idx+1
         x0(idx) = 0.0_dp; idx = idx+1
      end if
      if (use_shape) then
         if (kind == dist_std .or. kind == dist_sstd) then
            x0(idx) = log(5.0_dp-2.05_dp)
         else
            x0(idx) = log(2.0_dp-0.2_dp)
         end if
         idx = idx+1
      end if
      if (use_skew) x0(idx) = 0.0_dp

      opt = nelder_mead(garch_objective,x0,context,step=0.12_dp, &
                        tolerance=2.0e-7_dp,max_iterations=maxit)
      call decode_garch_parameters(opt%x,context,fit%spec)
      allocate(fit%residuals(size(y)),fit%sigma(size(y)))
      fit%log_likelihood = garch_log_likelihood(y,fit%spec,fit%residuals,fit%sigma)
      k = size(x0)
      fit%aic = -2.0_dp*fit%log_likelihood+2.0_dp*real(k,dp)
      fit%bic = -2.0_dp*fit%log_likelihood+log(real(size(y),dp))*real(k,dp)
      fit%iterations = opt%iterations
      fit%evaluations = opt%evaluations
      fit%status = opt%status
      if (opt%status == 0) then
         fit%message = 'converged'
      else
         fit%message = 'iteration limit reached; inspect estimates'
      end if
   end function fit_garch_core

   function garch_objective(x, generic_context) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value, llh
      type(garch_spec) :: spec

      select type (context => generic_context)
      type is (garch_fit_context)
         call decode_garch_parameters(x,context,spec)
         llh = garch_log_likelihood(context%y,spec)
         if (.not. finite_value(llh) .or. llh <= -huge(1.0_dp)/10.0_dp) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = -llh
         end if
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function garch_objective

   subroutine decode_garch_parameters(x, context, spec)
      real(dp), intent(in) :: x(:)
      type(garch_fit_context), intent(in) :: context
      type(garch_spec), intent(out) :: spec
      real(dp) :: ea, eb, denom
      integer :: idx

      spec = make_garch_spec(1,1,context%model,context%cond_dist)
      idx = 1
      if (context%fit_mean) then
         spec%mean = x(idx)
         idx = idx+1
      else
         spec%mean = 0.0_dp
      end if
      spec%omega = exp(max(-35.0_dp,min(20.0_dp,x(idx)))); idx = idx+1
      ea = exp(max(-20.0_dp,min(20.0_dp,x(idx)))); idx = idx+1
      eb = exp(max(-20.0_dp,min(20.0_dp,x(idx)))); idx = idx+1
      denom = 1.0_dp+ea+eb
      spec%alpha(1) = 0.998_dp*ea/denom
      spec%beta(1) = 0.998_dp*eb/denom
      spec%gamma(1) = 0.0_dp
      spec%delta = 2.0_dp
      if (context%model == model_aparch) then
         spec%gamma(1) = 0.999_dp*tanh(x(idx)); idx = idx+1
         spec%delta = 0.2_dp+3.8_dp/(1.0_dp+exp(-max(-30.0_dp,min(30.0_dp,x(idx)))))
         idx = idx+1
      end if
      spec%shape = 5.0_dp
      if (context%fit_shape) then
         if (context%cond_dist == dist_std .or. context%cond_dist == dist_sstd) then
            spec%shape = 2.05_dp+exp(max(-20.0_dp,min(10.0_dp,x(idx))))
         else
            spec%shape = 0.2_dp+exp(max(-20.0_dp,min(5.0_dp,x(idx))))
         end if
         idx = idx+1
      end if
      spec%skew = 1.0_dp
      if (context%fit_skew) spec%skew = exp(max(-5.0_dp,min(5.0_dp,x(idx))))
   end subroutine decode_garch_parameters

   pure elemental logical function finite_value(x) result(value)
      real(dp), intent(in) :: x
      value = ieee_is_finite(x)
   end function finite_value

end module fgarch_fit
