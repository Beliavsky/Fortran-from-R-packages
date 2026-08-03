! SPDX-License-Identifier: GPL-2.0-only
module glmnet_families
   use glmnet_kinds, only : dp
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument
   implicit none
   private
   public :: gaussian_identity_working
contains
   subroutine gaussian_identity_working(y, eta, base_weight, working, irls_weight, &
      deviance, status)
      real(dp), intent(in) :: y(:), eta(:), base_weight(:)
      real(dp), intent(out) :: working(:), irls_weight(:)
      real(dp), intent(out) :: deviance
      integer, intent(out) :: status
      if (size(eta) /= size(y) .or. size(base_weight) /= size(y) .or. &
          size(working) /= size(y) .or. size(irls_weight) /= size(y)) then
         status = glmnet_invalid_argument
         deviance = 0.0_dp
         return
      end if
      working = y
      irls_weight = base_weight
      deviance = sum(base_weight * (y - eta) ** 2)
      status = glmnet_success
   end subroutine gaussian_identity_working
end module glmnet_families
