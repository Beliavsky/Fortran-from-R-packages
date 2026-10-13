! SPDX-License-Identifier: BSD-3-Clause
module transform_test_callbacks
   use r_kinds, only: dp
   implicit none
   real(dp) :: visited_sum = 0.0_dp
   integer :: visited_count = 0
contains
   pure function add_values(values) result(value)
      real(dp), intent(in) :: values(:) !! All selected numeric components of a draw.
      real(dp) :: value
      value = sum(values)
   end function add_values
   subroutine accumulate(values, draw_id)
      real(dp), intent(in) :: values(:) !! Variables in the current draw.
      integer, intent(in) :: draw_id !! Expected sequential draw index.
      if (draw_id /= visited_count + 1) error stop 'Draw iteration order failed'
      visited_sum = visited_sum + sum(values)
      visited_count = visited_count + 1
   end subroutine accumulate
end module transform_test_callbacks

program test_transforms
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use posterior_transforms
   use transform_test_callbacks
   implicit none
   real(dp) :: x(3, 2, 2), matrix(2, 2, 3), nan, inf, cast(5), values(2, 3)
   real(dp), allocatable :: modified(:, :, :), applied(:, :), diagonal(:, :), constructed(:, :, :)
   integer :: i, status, matches(2, 3)
   procedure(draw_scalar_function), pointer :: callback
   procedure(draw_action), pointer :: action
   callback => add_values
   action => accumulate
   x = reshape([(real(i, dp), i = 1, 12)], shape(x))
   call mutate_variables(x, callback, modified, status = status)
   if (status /= 0 .or. size(modified, 3) /= 3) error stop 'Mutation append failed'
   if (any(modified(:, :, 3) /= x(:, :, 1) + x(:, :, 2))) error stop 'Mutation values failed'
   call mutate_variables(x, callback, modified, variable = 1, status = status)
   if (any(modified(:, :, 1) /= x(:, :, 1) + x(:, :, 2))) error stop 'Mutation replace failed'
   if (any(modified(:, :, 2) /= x(:, :, 2))) error stop 'Mutation overwrote other variable'
   call mutate_variables(x, callback, modified, variable = 0, status = status)
   if (status == 0 .or. size(modified) /= 0) error stop 'Invalid replacement not rejected'
   call for_each_draw(x, action, status)
   if (status /= 0 .or. visited_count /= 6) error stop 'Draw callback count failed'
   if (visited_sum /= sum(x)) error stop 'Draw callback values failed'
   matrix = reshape([(real(i, dp), i = 1, 12)], shape(matrix))
   call rvar_apply_margin(matrix, 1, callback, applied, status)
   if (status /= 0 .or. any(applied /= sum(matrix, dim = 3))) error stop 'Row margin apply failed'
   call rvar_apply_margin(matrix, 2, callback, applied, status)
   if (status /= 0 .or. any(applied /= sum(matrix, dim = 2))) error stop 'Column margin apply failed'
   call rvar_apply_margin(matrix, 0, callback, applied, status)
   if (status == 0 .or. size(applied) /= 0) error stop 'Invalid margin not rejected'
   nullify(callback)
   call mutate_variables(x, callback, modified, status = status)
   if (status == 0) error stop 'Null callback not rejected'
   diagonal = rvar_diagonal(matrix)
   if (any(diagonal(:, 1) /= matrix(:, 1, 1))) error stop 'First diagonal failed'
   if (any(diagonal(:, 2) /= matrix(:, 2, 2))) error stop 'Second diagonal failed'
   constructed = rvar_diagonal_matrix(diagonal, 3, 4)
   if (any(constructed(:, 3, 3) /= diagonal(:, 1))) error stop 'Diagonal recycling failed'
   if (any(constructed(:, 1, 3) /= 0.0_dp)) error stop 'Off diagonal values failed'
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   values = reshape([1.0_dp, 3.0_dp, nan, inf, 3.0_dp, 4.0_dp], shape(values))
   matches = rvar_match(values, [3.0_dp, nan, inf, 3.0_dp])
   if (any(reshape(matches, [6]) /= [0, 1, 2, 3, 1, 0])) error stop 'Numeric/NaN match failed'
   if (count(rvar_in(values, [3.0_dp, inf])) /= 3) error stop 'Membership failed'
   cast = rvar_as_integer([1.9_dp, -1.9_dp, nan, inf, 2147483648.0_dp])
   if (any(cast(:2) /= [1.0_dp, -1.0_dp])) error stop 'Integer truncation failed'
   if (.not. all(ieee_is_nan(cast(3:)))) error stop 'Invalid integer casts failed'
   cast = rvar_as_logical([0.0_dp, -1.9_dp, nan, inf, 5.0_dp])
   if (any(cast([1, 2, 4, 5]) /= [0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])) error stop 'Logical cast failed'
   if (.not. ieee_is_nan(cast(3))) error stop 'Logical NaN preservation failed'
   print *, 'posterior numeric transforms: property and regression checks passed'
end program test_transforms
