! SPDX-License-Identifier: GPL-3.0-only
module garchito
   use garchito_kinds, only : dp
   use garchito_types, only : garchito_control, garchito_result, &
      garchito_success, garchito_max_iterations, garchito_invalid_input, &
      garchito_numerical_failure
   use garchito_models, only : unified_est, realized_est, realized_est_option
   implicit none
   private

   public :: dp
   public :: garchito_control, garchito_result
   public :: garchito_success, garchito_max_iterations
   public :: garchito_invalid_input, garchito_numerical_failure
   public :: unified_est, realized_est, realized_est_option
end module garchito
