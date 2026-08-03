! SPDX-License-Identifier: GPL-2.0-only
module glmnet_relax
   use glmnet_kinds, only : dp
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument
   use glmnet_types, only : glmnet_control_type, glmnet_path_result, &
      glmnet_family_gaussian, glmnet_family_binomial, glmnet_family_poisson, &
      glmnet_family_multinomial, glmnet_family_mgaussian, glmnet_family_cox
   use glmnet_gaussian, only : fit_gaussian_path, fit_mgaussian_path
   use glmnet_glm, only : fit_binomial_path, fit_poisson_path
   use glmnet_multinomial, only : fit_multinomial_path
   use glmnet_cox, only : fit_cox_path
   implicit none
   private
   public :: relax_glmnet, relax_multinomial, relax_mgaussian, relax_cox
contains
   subroutine relax_glmnet(fit, x, y, relaxed, control, weights, offset, maxp)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:), y(:)
      type(glmnet_path_result), intent(out) :: relaxed
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights(:), offset(:)
      integer, intent(in), optional :: maxp
      type(glmnet_control_type) :: ctl
      type(glmnet_path_result) :: subfit = glmnet_path_result()
      real(dp), allocatable :: xs(:,:), w(:), o(:), lambda_zero(:)
      integer, allocatable :: active(:)
      integer :: l, j, count_active, maximum_active
      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      maximum_active = huge(1)
      if (present(maxp)) maximum_active = maxp
      call copy_structure(fit, relaxed)
      if (size(x, 1) /= size(y) .or. size(x, 2) /= fit%nvars) then
         relaxed%status = glmnet_invalid_argument
         return
      end if
      allocate(w(size(y)), o(size(y)), lambda_zero(1))
      w = 1.0_dp; o = 0.0_dp; lambda_zero = 1.0e-10_dp
      if (present(weights)) w = weights
      if (present(offset)) o = offset
      ctl%alpha = 0.0_dp
      ctl%nlambda = 1
      ctl%deviance_max = 1.0_dp
      ctl%fractional_deviance = 0.0_dp
      do l = 1, fit%nlambda
         count_active = count(maxval(abs(fit%beta(:, :, l)), dim=2) > 1.0e-12_dp)
         if (count_active == 0 .or. count_active > maximum_active) cycle
         allocate(active(count_active), xs(size(x, 1), count_active))
         count_active = 0
         do j = 1, fit%nvars
            if (maxval(abs(fit%beta(j, :, l))) > 1.0e-12_dp) then
               count_active = count_active + 1
               active(count_active) = j
               xs(:, count_active) = x(:, j)
            end if
         end do
         select case (fit%family_code)
         case (glmnet_family_gaussian)
            call fit_gaussian_path(xs, y, subfit, ctl, weights_in=w, offset_in=o, &
               lambda_in=lambda_zero)
         case (glmnet_family_binomial)
            call fit_binomial_path(xs, y, subfit, ctl, weights_in=w, offset_in=o, &
               lambda_in=lambda_zero)
         case (glmnet_family_poisson)
            call fit_poisson_path(xs, y, subfit, ctl, weights_in=w, offset_in=o, &
               lambda_in=lambda_zero)
         case default
            relaxed%status = glmnet_invalid_argument
            return
         end select
         if (subfit%nlambda == 1) then
            relaxed%beta(:, :, l) = 0.0_dp
            do j = 1, size(active)
               relaxed%beta(active(j), 1, l) = subfit%beta(j, 1, 1)
            end do
            relaxed%intercept(:, l) = subfit%intercept(:, 1)
            relaxed%dev_ratio(l) = subfit%dev_ratio(1)
            relaxed%objective(l) = subfit%objective(1)
            relaxed%df(l) = size(active)
            relaxed%converged(l) = subfit%converged(1)
         end if
         deallocate(active, xs)
      end do
   end subroutine relax_glmnet

   subroutine relax_multinomial(fit, x, class_id, relaxed, control, weights, offset, maxp)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: class_id(:)
      type(glmnet_path_result), intent(out) :: relaxed
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights(:), offset(:,:)
      integer, intent(in), optional :: maxp
      type(glmnet_control_type) :: ctl
      type(glmnet_path_result) :: subfit = glmnet_path_result()
      real(dp), allocatable :: xs(:,:), w(:), o(:,:), lambda_zero(:)
      integer, allocatable :: active(:)
      integer :: l, j, count_active, maximum_active
      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      maximum_active = huge(1); if (present(maxp)) maximum_active = maxp
      call copy_structure(fit, relaxed)
      if (fit%family_code /= glmnet_family_multinomial .or. &
          size(x, 1) /= size(class_id) .or. size(x, 2) /= fit%nvars) then
         relaxed%status = glmnet_invalid_argument
         return
      end if
      allocate(w(size(class_id)), o(size(class_id), fit%nout), lambda_zero(1))
      w = 1.0_dp; o = 0.0_dp; lambda_zero = 1.0e-10_dp
      if (present(weights)) w = weights
      if (present(offset)) o = offset
      ctl%alpha = 0.0_dp; ctl%nlambda = 1
      do l = 1, fit%nlambda
         count_active = count(maxval(abs(fit%beta(:, :, l)), dim=2) > 1.0e-12_dp)
         if (count_active == 0 .or. count_active > maximum_active) cycle
         allocate(active(count_active), xs(size(x, 1), count_active))
         count_active = 0
         do j = 1, fit%nvars
            if (maxval(abs(fit%beta(j, :, l))) > 1.0e-12_dp) then
               count_active = count_active + 1
               active(count_active) = j
               xs(:, count_active) = x(:, j)
            end if
         end do
         call fit_multinomial_path(xs, class_id, subfit, ctl, weights_in=w, &
            offset_in=o, lambda_in=lambda_zero, nclass=fit%nout)
         if (subfit%nlambda == 1) then
            relaxed%beta(:, :, l) = 0.0_dp
            do j = 1, size(active)
               relaxed%beta(active(j), :, l) = subfit%beta(j, :, 1)
            end do
            relaxed%intercept(:, l) = subfit%intercept(:, 1)
            relaxed%dev_ratio(l) = subfit%dev_ratio(1)
            relaxed%objective(l) = subfit%objective(1)
            relaxed%df(l) = size(active)
            relaxed%converged(l) = subfit%converged(1)
         end if
         deallocate(active, xs)
      end do
   end subroutine relax_multinomial

   subroutine relax_mgaussian(fit, x, y, relaxed, control, weights, offset, maxp)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:), y(:,:)
      type(glmnet_path_result), intent(out) :: relaxed
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights(:), offset(:,:)
      integer, intent(in), optional :: maxp
      type(glmnet_control_type) :: ctl
      type(glmnet_path_result) :: subfit = glmnet_path_result()
      real(dp), allocatable :: xs(:,:), w(:), o(:,:), lambda_zero(:)
      integer, allocatable :: active(:)
      integer :: l, j, count_active, maximum_active
      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      maximum_active = huge(1); if (present(maxp)) maximum_active = maxp
      call copy_structure(fit, relaxed)
      if (fit%family_code /= glmnet_family_mgaussian .or. size(y, 1) /= size(x, 1) .or. &
          size(y, 2) /= fit%nout .or. size(x, 2) /= fit%nvars) then
         relaxed%status = glmnet_invalid_argument
         return
      end if
      allocate(w(size(x, 1)), o(size(x, 1), fit%nout), lambda_zero(1))
      w = 1.0_dp; o = 0.0_dp; lambda_zero = 1.0e-10_dp
      if (present(weights)) w = weights
      if (present(offset)) o = offset
      ctl%alpha = 0.0_dp; ctl%nlambda = 1
      do l = 1, fit%nlambda
         count_active = count(maxval(abs(fit%beta(:, :, l)), dim=2) > 1.0e-12_dp)
         if (count_active == 0 .or. count_active > maximum_active) cycle
         allocate(active(count_active), xs(size(x, 1), count_active))
         count_active = 0
         do j = 1, fit%nvars
            if (maxval(abs(fit%beta(j, :, l))) > 1.0e-12_dp) then
               count_active = count_active + 1
               active(count_active) = j
               xs(:, count_active) = x(:, j)
            end if
         end do
         call fit_mgaussian_path(xs, y, subfit, ctl, weights_in=w, offset_in=o, &
            lambda_in=lambda_zero)
         if (subfit%nlambda == 1) then
            relaxed%beta(:, :, l) = 0.0_dp
            do j = 1, size(active)
               relaxed%beta(active(j), :, l) = subfit%beta(j, :, 1)
            end do
            relaxed%intercept(:, l) = subfit%intercept(:, 1)
            relaxed%dev_ratio(l) = subfit%dev_ratio(1)
            relaxed%objective(l) = subfit%objective(1)
            relaxed%df(l) = size(active)
            relaxed%converged(l) = subfit%converged(1)
         end if
         deallocate(active, xs)
      end do
   end subroutine relax_mgaussian


   subroutine relax_cox(fit, x, start_time, stop_time, event, relaxed, control, &
      weights, offset, strata, efron, maxp)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:), start_time(:), stop_time(:)
      integer, intent(in) :: event(:)
      type(glmnet_path_result), intent(out) :: relaxed
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights(:), offset(:)
      integer, intent(in), optional :: strata(:), maxp
      logical, intent(in), optional :: efron
      type(glmnet_control_type) :: ctl
      type(glmnet_path_result) :: subfit = glmnet_path_result()
      real(dp), allocatable :: xs(:,:), w(:), o(:), lambda_zero(:)
      integer, allocatable :: active(:), s(:)
      integer :: l, j, count_active, maximum_active
      logical :: ef
      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      maximum_active = huge(1)
      if (present(maxp)) maximum_active = maxp
      call copy_structure(fit, relaxed)
      if (fit%family_code /= glmnet_family_cox .or. size(x, 1) /= size(start_time) .or. &
          size(stop_time) /= size(start_time) .or. size(event) /= size(start_time) .or. &
          size(x, 2) /= fit%nvars) then
         relaxed%status = glmnet_invalid_argument
         return
      end if
      allocate(w(size(start_time)), o(size(start_time)), s(size(start_time)), lambda_zero(1))
      w = 1.0_dp
      o = 0.0_dp
      s = 1
      lambda_zero = 1.0e-10_dp
      if (present(weights)) w = weights
      if (present(offset)) o = offset
      if (present(strata)) s = strata
      ef = fit%efron
      if (present(efron)) ef = efron
      ctl%alpha = 0.0_dp
      ctl%nlambda = 1
      do l = 1, fit%nlambda
         count_active = count(maxval(abs(fit%beta(:, :, l)), dim=2) > 1.0e-12_dp)
         if (count_active == 0 .or. count_active > maximum_active) cycle
         allocate(active(count_active), xs(size(x, 1), count_active))
         count_active = 0
         do j = 1, fit%nvars
            if (maxval(abs(fit%beta(j, :, l))) > 1.0e-12_dp) then
               count_active = count_active + 1
               active(count_active) = j
               xs(:, count_active) = x(:, j)
            end if
         end do
         call fit_cox_path(xs, start_time, stop_time, event, subfit, ctl, weights_in=w, &
            offset_in=o, strata_in=s, lambda_in=lambda_zero, efron=ef)
         if (subfit%nlambda == 1) then
            relaxed%beta(:, :, l) = 0.0_dp
            do j = 1, size(active)
               relaxed%beta(active(j), 1, l) = subfit%beta(j, 1, 1)
            end do
            relaxed%intercept(:, l) = 0.0_dp
            relaxed%dev_ratio(l) = subfit%dev_ratio(1)
            relaxed%objective(l) = subfit%objective(1)
            relaxed%df(l) = size(active)
            relaxed%converged(l) = subfit%converged(1)
         end if
         deallocate(active, xs)
      end do
   end subroutine relax_cox

   subroutine copy_structure(source, destination)
      type(glmnet_path_result), intent(in) :: source
      type(glmnet_path_result), intent(out) :: destination
      destination = source
      destination%status = glmnet_success
   end subroutine copy_structure
end module glmnet_relax
