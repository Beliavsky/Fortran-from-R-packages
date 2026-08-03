! SPDX-License-Identifier: GPL-2.0-only
module glmnet_control
   use glmnet_kinds, only : dp
   use glmnet_types, only : glmnet_control_type
   implicit none
   private
   public :: default_glmnet_control, update_glmnet_control
contains
   pure function default_glmnet_control() result(control)
      type(glmnet_control_type) :: control
      control = glmnet_control_type()
   end function default_glmnet_control

   subroutine update_glmnet_control(control, alpha, nlambda, lambda_min_ratio, &
      threshold, max_iterations, max_outer_iterations, standardize, intercept, &
      grouped, trace)
      type(glmnet_control_type), intent(inout) :: control
      real(dp), intent(in), optional :: alpha, lambda_min_ratio, threshold
      integer, intent(in), optional :: nlambda, max_iterations, max_outer_iterations
      logical, intent(in), optional :: standardize, intercept, grouped, trace
      if (present(alpha)) control%alpha = min(max(alpha, 0.0_dp), 1.0_dp)
      if (present(nlambda)) control%nlambda = max(nlambda, 1)
      if (present(lambda_min_ratio)) control%lambda_min_ratio = lambda_min_ratio
      if (present(threshold)) control%threshold = max(threshold, epsilon(1.0_dp))
      if (present(max_iterations)) control%max_iterations = max(max_iterations, 1)
      if (present(max_outer_iterations)) &
         control%max_outer_iterations = max(max_outer_iterations, 1)
      if (present(standardize)) control%standardize = standardize
      if (present(intercept)) control%intercept = intercept
      if (present(grouped)) control%grouped = grouped
      if (present(trace)) control%trace = trace
   end subroutine update_glmnet_control
end module glmnet_control
