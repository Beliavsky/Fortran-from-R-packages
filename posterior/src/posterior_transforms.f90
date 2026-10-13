! SPDX-License-Identifier: BSD-3-Clause
! Numeric specializations of posterior draw/rvar transformations.
module posterior_transforms
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only: dp
   implicit none
   private
   public :: rvar_match, rvar_in, rvar_as_integer, rvar_as_logical
   public :: rvar_diagonal, rvar_diagonal_matrix, mutate_variables, rvar_apply_margin
   public :: for_each_draw, draw_scalar_function, draw_action
   abstract interface
      pure function draw_scalar_function(values) result(value)
         import dp
         real(dp), intent(in) :: values(:) !! Values of all selected variables for a single draw.
         real(dp) :: value
      end function draw_scalar_function
      subroutine draw_action(values, draw_id)
         import dp
         real(dp), intent(in) :: values(:) !! Values of all variables for this draw.
         integer, intent(in) :: draw_id !! One-based draw index, with iterations varying before chains.
      end subroutine draw_action
   end interface
contains
   pure function rvar_match(x, table, nomatch) result(values)
      real(dp), intent(in) :: x(:, :) !! Draws by random-variable elements, matched independently at each position.
      real(dp), intent(in) :: table(:) !! Numeric lookup values; NaN matches the first NaN, with NA/NaN not distinguished.
      integer, intent(in), optional :: nomatch !! Integer for unmatched values, default zero instead of R NA_integer_.
      integer :: values(size(x, 1), size(x, 2)), i, j, k, missing
      logical :: equal
      missing = 0
      if (present(nomatch)) missing = nomatch
      values = missing
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            do k = 1, size(table)
               equal = .false.
               if (ieee_is_nan(x(i, j))) then
                  equal = ieee_is_nan(table(k))
               else if (.not. ieee_is_nan(table(k))) then
                  equal = x(i, j) == table(k)
               end if
               if (.not. equal) cycle
               values(i, j) = k
               exit
            end do
         end do
      end do
   end function rvar_match

   pure function rvar_in(x, table) result(values)
      real(dp), intent(in) :: x(:, :) !! Draws by elements to test for lookup membership.
      real(dp), intent(in) :: table(:) !! Numeric lookup values, including an optional matching NaN.
      logical :: values(size(x, 1), size(x, 2))
      values = rvar_match(x, table, nomatch = 0) > 0
   end function rvar_in

   pure elemental function rvar_as_integer(x) result(value)
      real(dp), intent(in) :: x !! Numeric draw, truncated toward zero; R 32-bit integer overflow or NaN gives NaN.
      real(dp) :: value, truncated
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. ieee_is_finite(x)) return
      truncated = aint(x)
      if (truncated > 2147483647.0_dp .or. truncated < -2147483647.0_dp) return
      value = truncated
   end function rvar_as_integer

   pure elemental function rvar_as_logical(x) result(value)
      real(dp), intent(in) :: x !! Numeric draw converted to 0/1, with NaN preserved; infinities are true.
      real(dp) :: value
      value = x
      if (.not. ieee_is_nan(x)) value = merge(1.0_dp, 0.0_dp, x /= 0.0_dp)
   end function rvar_as_logical

   pure function rvar_diagonal(x) result(values)
      real(dp), intent(in) :: x(:, :, :) !! Draws by rows by columns of a random matrix.
      real(dp) :: values(size(x, 1), min(size(x, 2), size(x, 3)))
      integer :: j
      do j = 1, size(values, 2)
         values(:, j) = x(:, j, j)
      end do
   end function rvar_diagonal

   pure function rvar_diagonal_matrix(x, nrow, ncol) result(values)
      real(dp), intent(in) :: x(:, :) !! Draws by diagonal-vector elements, recycled to the requested diagonal length.
      integer, intent(in) :: nrow !! Nonnegative number of output matrix rows.
      integer, intent(in) :: ncol !! Nonnegative number of output matrix columns.
      real(dp) :: values(size(x, 1), max(0, nrow), max(0, ncol))
      integer :: j
      values = 0.0_dp
      if (size(x, 2) == 0) return
      do j = 1, min(nrow, ncol)
         values(:, j, j) = x(:, mod(j - 1, size(x, 2)) + 1)
      end do
   end function rvar_diagonal_matrix

   pure subroutine mutate_variables(x, transform, values, variable, status)
      real(dp), intent(in) :: x(:, :, :) !! Iterations by chains by variables whose per-draw values enter the transform.
      procedure(draw_scalar_function), pointer, intent(in) :: transform !! Associated pure callback returning one derived scalar per draw.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Draws with one derived variable appended or an existing variable replaced.
      integer, intent(in), optional :: variable !! Existing variable index to replace; absent appends at the end.
      integer, intent(out), optional :: status !! Zero success, one null callback or invalid replacement index.
      integer :: i, j, target, nvars
      if (present(status)) status = 1
      allocate(values(0, 0, 0))
      if (.not. associated(transform)) return
      nvars = size(x, 3)
      target = nvars + 1
      if (present(variable)) then
         if (variable < 1 .or. variable > nvars) return
         target = variable
      else
         nvars = nvars + 1
      end if
      deallocate(values)
      allocate(values(size(x, 1), size(x, 2), nvars))
      values(:, :, :size(x, 3)) = x
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            values(i, j, target) = transform(x(i, j, :))
         end do
      end do
      if (present(status)) status = 0
   end subroutine mutate_variables

   pure subroutine rvar_apply_margin(x, margin, transform, values, status)
      real(dp), intent(in) :: x(:, :, :) !! Draws by rows by columns of a random matrix.
      integer, intent(in) :: margin !! Retained matrix margin: 1 retains rows, 2 retains columns.
      procedure(draw_scalar_function), pointer, intent(in) :: transform !! Pure scalar reduction applied separately within each draw.
      real(dp), allocatable, intent(out) :: values(:, :) !! Draws by retained matrix margin of reduction results.
      integer, intent(out), optional :: status !! Zero success, one unsupported margin or unassociated callback.
      integer :: i, j, n
      if (present(status)) status = 1
      allocate(values(0, 0))
      if (.not. associated(transform)) return
      if (margin < 1 .or. margin > 2) return
      n = size(x, margin + 1)
      deallocate(values)
      allocate(values(size(x, 1), n))
      do j = 1, n
         do i = 1, size(x, 1)
            if (margin == 1) then
               values(i, j) = transform(x(i, j, :))
            else
               values(i, j) = transform(x(i, :, j))
            end if
         end do
      end do
      if (present(status)) status = 0
   end subroutine rvar_apply_margin

   subroutine for_each_draw(x, action, status)
      real(dp), intent(in) :: x(:, :, :) !! Iterations by chains by variables, visited in chain-major order.
      procedure(draw_action), pointer, intent(in) :: action !! Associated callback receiving variables and one-based draw index.
      integer, intent(out), optional :: status !! Zero success, one unassociated callback.
      integer :: i, j
      if (present(status)) status = 1
      if (.not. associated(action)) return
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            call action(x(i, j, :), (j - 1) * size(x, 1) + i)
         end do
      end do
      if (present(status)) status = 0
   end subroutine for_each_draw
end module posterior_transforms
