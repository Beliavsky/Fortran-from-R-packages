! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_utils
   use mgcv_kinds, only : dp
   use mgcv_linalg, only : kronecker_product
   implicit none
   private
   public :: not_exp, not_log, not_exp2, not_log2
   public :: null_space_dimension, unique_rows, exclude_too_far
   public :: tensor_product_model_matrix, tensor_product_penalties
   public :: difference_matrix, diagonal_extract, diagonal_set
   public :: center_columns, standardize_vector, log_expm1, log1p_stable

contains

   elemental function not_exp(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y, z
      if (x > 1.0_dp) then
         y = exp(1.0_dp) * (x * x + 1.0_dp) / 2.0_dp
      else if (x > -1.0_dp) then
         y = exp(x)
      else
         z = -x
         y = 2.0_dp / (exp(1.0_dp) * (z * z + 1.0_dp))
      end if
   end function not_exp

   elemental function not_log(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y, z
      if (x <= 0.0_dp) then
         y = -huge(1.0_dp)
      else if (x > exp(1.0_dp)) then
         y = sqrt(2.0_dp * x / exp(1.0_dp) - 1.0_dp)
      else if (x > exp(-1.0_dp)) then
         y = log(x)
      else
         z = 1.0_dp / x
         y = -sqrt(2.0_dp * z / exp(1.0_dp) - 1.0_dp)
      end if
   end function not_log

   elemental function not_exp2(x, d) result(y)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: d
      real(dp) :: y, dd
      dd = 20.0_dp; if (present(d)) dd = d
      y = exp(dd * sin(x / dd))
   end function not_exp2

   elemental function not_log2(x, d) result(y)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: d
      real(dp) :: y, dd, z
      dd = 20.0_dp; if (present(d)) dd = d
      if (x <= 0.0_dp) then
         y = -0.5_dp * acos(-1.0_dp) * dd
      else
         z = max(-1.0_dp, min(1.0_dp, log(x) / dd))
         y = dd * asin(z)
      end if
   end function not_log2

   elemental function log1p_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (abs(x) < 1.0e-4_dp) then
         y = x * (1.0_dp + x * (-0.5_dp + x * (1.0_dp / 3.0_dp + &
             x * (-0.25_dp + 0.2_dp * x))))
      else
         y = log(1.0_dp + x)
      end if
   end function log1p_stable

   elemental function log_expm1(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (x > 50.0_dp) then
         y = x + log1p_stable(-exp(-x))
      else if (x < 1.0e-5_dp) then
         y = log(x) + log1p_stable(0.5_dp * x + x * x / 6.0_dp)
      else
         y = log(exp(x) - 1.0_dp)
      end if
   end function log_expm1

   integer function null_space_dimension(d, m) result(dim)
      integer, intent(in) :: d, m
      integer :: i
      real(dp) :: value
      if (d < 1 .or. m < 0) then
         dim = 0; return
      end if
      value = 1.0_dp
      do i = 1, d
         value = value * real(m + i - 1, dp) / real(i, dp)
      end do
      dim = nint(value)
   end function null_space_dimension

   subroutine unique_rows(x, unique, index, counts, tol)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: unique(:, :)
      integer, allocatable, intent(out) :: index(:), counts(:)
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: work(:, :)
      real(dp) :: eps
      integer :: i, j, n, p, nu
      logical :: found

      n = size(x, 1); p = size(x, 2)
      eps = 0.0_dp; if (present(tol)) eps = max(0.0_dp, tol)
      allocate(work(n, p), index(n), counts(n))
      work = 0.0_dp; counts = 0; nu = 0
      do i = 1, n
         found = .false.
         do j = 1, nu
            if (maxval(abs(x(i, :) - work(j, :))) <= eps) then
               index(i) = j; counts(j) = counts(j) + 1; found = .true.; exit
            end if
         end do
         if (.not. found) then
            nu = nu + 1; work(nu, :) = x(i, :); index(i) = nu; counts(nu) = 1
         end if
      end do
      allocate(unique(nu, p)); unique = work(1:nu, :)
      counts = counts(1:nu)
   end subroutine unique_rows

   subroutine exclude_too_far(grid_x, grid_y, data_x, data_y, distance, excluded)
      real(dp), intent(in) :: grid_x(:), grid_y(:), data_x(:), data_y(:), distance
      logical, allocatable, intent(out) :: excluded(:)
      real(dp) :: xr, yr, gx, gy, dx, dy, best
      integer :: i, j
      if (size(grid_x) /= size(grid_y) .or. size(data_x) /= size(data_y)) then
         allocate(excluded(0)); return
      end if
      xr = maxval(data_x) - minval(data_x); yr = maxval(data_y) - minval(data_y)
      if (xr <= 0.0_dp) xr = 1.0_dp
      if (yr <= 0.0_dp) yr = 1.0_dp
      allocate(excluded(size(grid_x)))
      do i = 1, size(grid_x)
         gx = (grid_x(i) - minval(data_x)) / xr
         gy = (grid_y(i) - minval(data_y)) / yr
         best = huge(1.0_dp)
         do j = 1, size(data_x)
            dx = gx - (data_x(j) - minval(data_x)) / xr
            dy = gy - (data_y(j) - minval(data_y)) / yr
            best = min(best, sqrt(dx * dx + dy * dy))
         end do
         excluded(i) = best > distance
      end do
   end subroutine exclude_too_far

   function tensor_product_model_matrix(a, b) result(x)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable :: x(:, :)
      integer :: i, j, k, col
      if (size(a, 1) /= size(b, 1)) then
         allocate(x(0, 0)); return
      end if
      allocate(x(size(a, 1), size(a, 2) * size(b, 2)))
      do i = 1, size(a, 1)
         col = 0
         do j = 1, size(a, 2)
            do k = 1, size(b, 2)
               col = col + 1
               x(i, col) = a(i, j) * b(i, k)
            end do
         end do
      end do
   end function tensor_product_model_matrix

   subroutine tensor_product_penalties(s1, s2, penalties)
      real(dp), intent(in) :: s1(:, :), s2(:, :)
      real(dp), allocatable, intent(out) :: penalties(:, :, :)
      real(dp), allocatable :: i1(:, :), i2(:, :)
      integer :: i
      allocate(i1(size(s1, 1), size(s1, 1)), i2(size(s2, 1), size(s2, 1)))
      i1 = 0.0_dp; i2 = 0.0_dp
      do i = 1, size(i1, 1); i1(i, i) = 1.0_dp; end do
      do i = 1, size(i2, 1); i2(i, i) = 1.0_dp; end do
      allocate(penalties(size(s1, 1) * size(s2, 1), &
                         size(s1, 2) * size(s2, 2), 2))
      penalties(:, :, 1) = kronecker_product(s1, i2)
      penalties(:, :, 2) = kronecker_product(i1, s2)
   end subroutine tensor_product_penalties

   function difference_matrix(n, order) result(d)
      integer, intent(in) :: n, order
      real(dp), allocatable :: d(:, :), current(:, :), next(:, :)
      integer :: i, r, m
      if (n <= 0 .or. order < 0 .or. order >= n) then
         allocate(d(0, 0)); return
      end if
      allocate(current(n, n)); current = 0.0_dp
      do i = 1, n; current(i, i) = 1.0_dp; end do
      m = n
      do r = 1, order
         allocate(next(m - 1, n)); next = current(2:m, :) - current(1:m - 1, :)
         call move_alloc(next, current); m = m - 1
      end do
      allocate(d(m, n)); d = current
   end function difference_matrix

   function diagonal_extract(a, k) result(v)
      real(dp), intent(in) :: a(:, :)
      integer, intent(in), optional :: k
      real(dp), allocatable :: v(:)
      integer :: kk, n, i, j, len
      kk = 0; if (present(k)) kk = k
      if (kk >= 0) then
         len = min(size(a, 1), size(a, 2) - kk)
         if (len < 0) len = 0
         allocate(v(len))
         do n = 1, len; i = n; j = n + kk; v(n) = a(i, j); end do
      else
         len = min(size(a, 1) + kk, size(a, 2))
         if (len < 0) len = 0
         allocate(v(len))
         do n = 1, len; i = n - kk; j = n; v(n) = a(i, j); end do
      end if
   end function diagonal_extract

   subroutine diagonal_set(a, value, k)
      real(dp), intent(inout) :: a(:, :)
      real(dp), intent(in) :: value(:)
      integer, intent(in), optional :: k
      integer :: kk, n, i, j, len
      kk = 0; if (present(k)) kk = k
      if (kk >= 0) then
         len = min(size(a, 1), size(a, 2) - kk)
         do n = 1, min(len, size(value)); i = n; j = n + kk; a(i, j) = value(n); end do
      else
         len = min(size(a, 1) + kk, size(a, 2))
         do n = 1, min(len, size(value)); i = n - kk; j = n; a(i, j) = value(n); end do
      end if
   end subroutine diagonal_set

   subroutine center_columns(x, means)
      real(dp), intent(inout) :: x(:, :)
      real(dp), allocatable, intent(out) :: means(:)
      integer :: j
      allocate(means(size(x, 2)))
      do j = 1, size(x, 2)
         means(j) = sum(x(:, j)) / real(max(1, size(x, 1)), dp)
         x(:, j) = x(:, j) - means(j)
      end do
   end subroutine center_columns

   subroutine standardize_vector(x, mean_x, scale_x, z)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: mean_x, scale_x
      real(dp), allocatable, intent(out) :: z(:)
      mean_x = sum(x) / real(max(1, size(x)), dp)
      scale_x = sqrt(sum((x - mean_x)**2) / real(max(1, size(x) - 1), dp))
      if (scale_x <= sqrt(epsilon(1.0_dp))) scale_x = 1.0_dp
      allocate(z(size(x))); z = (x - mean_x) / scale_x
   end subroutine standardize_vector

end module mgcv_utils
