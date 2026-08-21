! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_tensor
  !! Public bridge to the bundled tensorA translation.  This closes the
  !! dynamic named-index/high-rank tensor gap without changing the simpler
  !! explicit-array APIs used by the rest of compositions-fortran.
  use tensora, only: tensor_t, tensor, to_tensor, tensor_from_real, tensor_from_complex, &
    tensor_zeros, tensor_ones, scalar_tensor, tensor_is_real, real_data, &
    reorder_tensor, reorder_tensor_pos, reorder_tensor_names, pos_tensor, &
    rename_axis, rename_first_axes, mark_tensor, mul_tensor, mul_tensor_pos, mul_tensor_names, &
    trace_tensor, margin_tensor, diagmul_tensor, delta_tensor, diag_tensor, tripledelta_tensor, one_tensor, &
    add_tensor, sub_tensor, elem_mul_tensor, elem_div_tensor, scale_tensor, repeat_tensor, slice_tensor, &
    undrop_tensor, untensor_tensor, bind_tensor, einstein_pair, riemann_pair, contraname_tensor, &
    is_covariate_tensor, is_contravariate_tensor, inv_tensor, solve_tensor, svd_tensor, chol_tensor, &
    power_tensor, opnorm_tensor, opnorm_by_tensor, to_matrix_tensor, norm_tensor, norm_tensor_along, &
    mean_tensor, var_tensor, cov_tensor, drag_tensor, positions_by_name, contraname, &
    is_covariate_name, is_contravariate_name, as_covariate_name, as_contravariate_name
  implicit none
  private
  public :: tensor_t, tensor, to_tensor, tensor_from_real, tensor_from_complex
  public :: tensor_zeros, tensor_ones, scalar_tensor, tensor_is_real, real_data
  public :: reorder_tensor, reorder_tensor_pos, reorder_tensor_names, pos_tensor
  public :: rename_axis, rename_first_axes, mark_tensor, mul_tensor, mul_tensor_pos, mul_tensor_names
  public :: trace_tensor, margin_tensor, diagmul_tensor, delta_tensor, diag_tensor, tripledelta_tensor, one_tensor
  public :: add_tensor, sub_tensor, elem_mul_tensor, elem_div_tensor, scale_tensor, repeat_tensor, slice_tensor
  public :: undrop_tensor, untensor_tensor, bind_tensor, einstein_pair, riemann_pair, contraname_tensor
  public :: is_covariate_tensor, is_contravariate_tensor, inv_tensor, solve_tensor, svd_tensor, chol_tensor
  public :: power_tensor, opnorm_tensor, opnorm_by_tensor, to_matrix_tensor
  public :: norm_tensor, norm_tensor_along, mean_tensor, var_tensor, cov_tensor, drag_tensor
  public :: positions_by_name, contraname, is_covariate_name, is_contravariate_name
  public :: as_covariate_name, as_contravariate_name
end module compositions_tensor
