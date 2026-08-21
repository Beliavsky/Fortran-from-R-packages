! SPDX-License-Identifier: GPL-2.0-or-later
module numderiv
   use numderiv_kinds, only : dp
   use numderiv_types, only : deriv_options, gend_result, first_deriv_options, &
      hessian_options, nd_success, nd_invalid_argument, nd_nonfinite_value, &
      nd_shape_mismatch
   use numderiv_callbacks, only : scalar_real_function, vector_real_function, &
      scalar_complex_function, vector_complex_function
   use numderiv_core, only : grad, grad_elementwise, grad_complex, &
      grad_elementwise_complex, jacobian, jacobian_complex, hessian, &
      hessian_complex, gend
   implicit none
   public
end module numderiv
