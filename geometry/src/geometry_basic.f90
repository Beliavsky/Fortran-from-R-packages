module geometry_basic
   use geometry_kinds, only : dp
   use geometry_linalg, only : cross3
   implicit none
   private

   public :: cart2pol_points, pol2cart_points, cart2sph_points, sph2cart_points
   public :: dot_columns, dot_rows, extprod3d_rows, polyarea_vector
   public :: matmax_rows, matmin_rows, matsort_rows, matorder_rows, unique_rows
   public :: entry_value_real, entry_set_real
   public :: mesh_dcircle, mesh_drectangle, mesh_dsphere, mesh_hunif
   public :: mesh_diff_values, mesh_union_values, mesh_intersect_values
   public :: rbox_points

contains

   pure subroutine cart2pol_points(cart, polar)
      real(dp), intent(in) :: cart(:, :) !! Cartesian points with two or three columns.
      real(dp), allocatable, intent(out) :: polar(:, :) !! Polar or cylindrical points with matching row count.
      integer :: n, d

      n = size(cart, 1)
      d = size(cart, 2)
      if (d == 2) then
         allocate(polar(n, 2))
      else if (d == 3) then
         allocate(polar(n, 3))
      else
         allocate(polar(0, 0))
         return
      end if
      polar(:, 1) = atan2(cart(:, 2), cart(:, 1))
      polar(:, 2) = sqrt(cart(:, 1)**2 + cart(:, 2)**2)
      if (d == 3) polar(:, 3) = cart(:, 3)
   end subroutine cart2pol_points

   pure subroutine pol2cart_points(polar, cart)
      real(dp), intent(in) :: polar(:, :) !! Polar or cylindrical points with two or three columns.
      real(dp), allocatable, intent(out) :: cart(:, :) !! Cartesian points with matching row count.
      integer :: n, d

      n = size(polar, 1)
      d = size(polar, 2)
      if (d == 2) then
         allocate(cart(n, 2))
      else if (d == 3) then
         allocate(cart(n, 3))
      else
         allocate(cart(0, 0))
         return
      end if
      cart(:, 1) = polar(:, 2)*cos(polar(:, 1))
      cart(:, 2) = polar(:, 2)*sin(polar(:, 1))
      if (d == 3) cart(:, 3) = polar(:, 3)
   end subroutine pol2cart_points

   pure subroutine cart2sph_points(cart, spherical)
      real(dp), intent(in) :: cart(:, :) !! Cartesian points with exactly three columns.
      real(dp), allocatable, intent(out) :: spherical(:, :) !! Spherical points theta, phi, radius.
      integer :: n

      if (size(cart, 2) /= 3) then
         allocate(spherical(0, 0))
         return
      end if
      n = size(cart, 1)
      allocate(spherical(n, 3))
      spherical(:, 1) = atan2(cart(:, 2), cart(:, 1))
      spherical(:, 2) = atan2(cart(:, 3), sqrt(cart(:, 1)**2 + cart(:, 2)**2))
      spherical(:, 3) = sqrt(sum(cart*cart, dim=2))
   end subroutine cart2sph_points

   pure subroutine sph2cart_points(spherical, cart)
      real(dp), intent(in) :: spherical(:, :) !! Spherical points theta, phi, radius.
      real(dp), allocatable, intent(out) :: cart(:, :) !! Cartesian points with three columns.
      integer :: n

      if (size(spherical, 2) /= 3) then
         allocate(cart(0, 0))
         return
      end if
      n = size(spherical, 1)
      allocate(cart(n, 3))
      cart(:, 1) = spherical(:, 3)*cos(spherical(:, 2))*cos(spherical(:, 1))
      cart(:, 2) = spherical(:, 3)*cos(spherical(:, 2))*sin(spherical(:, 1))
      cart(:, 3) = spherical(:, 3)*sin(spherical(:, 2))
   end subroutine sph2cart_points

   pure subroutine dot_columns(x, y, values)
      real(dp), intent(in) :: x(:, :) !! First matrix of vectors stored by columns.
      real(dp), intent(in) :: y(:, :) !! Second matrix with the same shape as x.
      real(dp), allocatable, intent(out) :: values(:) !! Column-wise dot products.
      if (any(shape(x) /= shape(y))) then
         allocate(values(0))
         return
      end if
      allocate(values(size(x, 2)))
      values = sum(x*y, dim=1)
   end subroutine dot_columns

   pure subroutine dot_rows(x, y, values)
      real(dp), intent(in) :: x(:, :) !! First matrix of vectors stored by rows.
      real(dp), intent(in) :: y(:, :) !! Second matrix with the same shape as x.
      real(dp), allocatable, intent(out) :: values(:) !! Row-wise dot products.
      if (any(shape(x) /= shape(y))) then
         allocate(values(0))
         return
      end if
      allocate(values(size(x, 1)))
      values = sum(x*y, dim=2)
   end subroutine dot_rows

   pure subroutine extprod3d_rows(x, y, z)
      real(dp), intent(in) :: x(:, :) !! First set of row-wise three-dimensional vectors.
      real(dp), intent(in) :: y(:, :) !! Second set of row-wise three-dimensional vectors.
      real(dp), allocatable, intent(out) :: z(:, :) !! Row-wise cross products.
      integer :: i, n

      if (size(x, 2) /= 3 .or. any(shape(x) /= shape(y))) then
         allocate(z(0, 0))
         return
      end if
      n = size(x, 1)
      allocate(z(n, 3))
      do i = 1, n
         z(i, :) = cross3(x(i, :), y(i, :))
      end do
   end subroutine extprod3d_rows

   pure function polyarea_vector(x, y) result(area)
      real(dp), intent(in) :: x(:) !! Polygon x-coordinates in boundary order.
      real(dp), intent(in) :: y(:) !! Polygon y-coordinates in boundary order.
      real(dp) :: area
      integer :: i, j, n

      area = 0.0_dp
      if (size(x) /= size(y) .or. size(x) < 3) return
      n = size(x)
      j = n
      do i = 1, n
         area = area + x(j)*y(i) - x(i)*y(j)
         j = i
      end do
      area = 0.5_dp*abs(area)
   end function polyarea_vector

   pure subroutine matmax_rows(x, values)
      real(dp), intent(in) :: x(:, :) !! Matrix whose row maxima are required.
      real(dp), allocatable, intent(out) :: values(:) !! Maximum value from each matrix row.
      allocate(values(size(x, 1)))
      if (size(x, 2) == 0) then
         values = -huge(1.0_dp)
      else
         values = maxval(x, dim=2)
      end if
   end subroutine matmax_rows

   pure subroutine matmin_rows(x, values)
      real(dp), intent(in) :: x(:, :) !! Matrix whose row minima are required.
      real(dp), allocatable, intent(out) :: values(:) !! Minimum value from each matrix row.
      allocate(values(size(x, 1)))
      if (size(x, 2) == 0) then
         values = huge(1.0_dp)
      else
         values = minval(x, dim=2)
      end if
   end subroutine matmin_rows

   pure subroutine matsort_rows(x, sorted)
      real(dp), intent(in) :: x(:, :) !! Matrix whose rows are sorted independently.
      real(dp), allocatable, intent(out) :: sorted(:, :) !! Matrix with each input row sorted ascending.
      real(dp) :: key
      integer :: i, j, k

      allocate(sorted(size(x, 1), size(x, 2)))
      sorted = x
      do i = 1, size(sorted, 1)
         do j = 2, size(sorted, 2)
            key = sorted(i, j)
            k = j - 1
            do while (k >= 1)
               if (sorted(i, k) <= key) exit
               sorted(i, k + 1) = sorted(i, k)
               k = k - 1
            end do
            sorted(i, k + 1) = key
         end do
      end do
   end subroutine matsort_rows

   pure subroutine matorder_rows(x, order_index)
      real(dp), intent(in) :: x(:, :) !! Matrix whose rows are ordered lexicographically by columns.
      integer, allocatable, intent(out) :: order_index(:) !! Row permutation giving ascending lexicographic order.
      integer :: i, j, key

      allocate(order_index(size(x, 1)))
      order_index = [(i, i=1, size(x, 1))]
      do i = 2, size(order_index)
         key = order_index(i)
         j = i - 1
         do while (j >= 1)
            if (.not. row_less(x(key, :), x(order_index(j), :))) exit
            order_index(j + 1) = order_index(j)
            j = j - 1
         end do
         order_index(j + 1) = key
      end do
   end subroutine matorder_rows

   pure subroutine unique_rows(x, rows_are_sets, unique_x)
      real(dp), intent(in) :: x(:, :) !! Numeric matrix from which duplicate rows are removed.
      logical, intent(in) :: rows_are_sets !! Treat entries within each row as an unordered set when true.
      real(dp), allocatable, intent(out) :: unique_x(:, :) !! Sorted unique rows.
      real(dp), allocatable :: work(:, :), tmp(:, :)
      integer, allocatable :: ord(:)
      integer :: i, count

      if (rows_are_sets) then
         call matsort_rows(x, work)
      else
         allocate(work(size(x, 1), size(x, 2)))
         work = x
      end if
      call matorder_rows(work, ord)
      allocate(tmp(size(work, 1), size(work, 2)))
      count = 0
      do i = 1, size(ord)
         if (count == 0) then
            count = 1
            tmp(count, :) = work(ord(i), :)
         else if (any(abs(work(ord(i), :) - tmp(count, :)) > 0.0_dp)) then
            count = count + 1
            tmp(count, :) = work(ord(i), :)
         end if
      end do
      allocate(unique_x(count, size(work, 2)))
      if (count > 0) unique_x = tmp(1:count, :)
   end subroutine unique_rows

   pure subroutine entry_value_real(data, dims, idx, values, ok)
      real(dp), intent(in) :: data(:) !! Flat column-major storage of the source array.
      integer, intent(in) :: dims(:) !! Dimensions of the source array.
      integer, intent(in) :: idx(:, :) !! One-based subscripts, one array entry per row.
      real(dp), allocatable, intent(out) :: values(:) !! Values selected by idx.
      logical, intent(out) :: ok !! True when all dimensions and indices are valid.
      integer :: i, j, stride, offset

      ok = .false.
      allocate(values(size(idx, 1)))
      values = 0.0_dp
      if (size(idx, 2) /= size(dims)) return
      if (product(dims) /= size(data)) return
      do i = 1, size(idx, 1)
         stride = 1
         offset = 1
         do j = 1, size(dims)
            if (idx(i, j) < 1 .or. idx(i, j) > dims(j)) return
            offset = offset + (idx(i, j) - 1)*stride
            stride = stride*dims(j)
         end do
         values(i) = data(offset)
      end do
      ok = .true.
   end subroutine entry_value_real

   pure subroutine entry_set_real(data, dims, idx, values, ok)
      real(dp), intent(inout) :: data(:) !! Flat column-major array storage modified at selected entries.
      integer, intent(in) :: dims(:) !! Dimensions of the array represented by data.
      integer, intent(in) :: idx(:, :) !! One-based subscripts, one assignment per row.
      real(dp), intent(in) :: values(:) !! Replacement values, length one or number of selected entries.
      logical, intent(out) :: ok !! True when the replacement completed with valid subscripts.
      integer :: i, j, stride, offset, iv

      ok = .false.
      if (size(idx, 2) /= size(dims)) return
      if (product(dims) /= size(data)) return
      if (size(values) /= 1 .and. size(values) /= size(idx, 1)) return
      do i = 1, size(idx, 1)
         stride = 1
         offset = 1
         do j = 1, size(dims)
            if (idx(i, j) < 1 .or. idx(i, j) > dims(j)) return
            offset = offset + (idx(i, j) - 1)*stride
            stride = stride*dims(j)
         end do
         iv = min(i, size(values))
         data(offset) = values(iv)
      end do
      ok = .true.
   end subroutine entry_set_real

   pure subroutine mesh_dcircle(p, radius, distance)
      real(dp), intent(in) :: p(:, :) !! Two-dimensional points at which signed distance is evaluated.
      real(dp), intent(in) :: radius !! Circle radius, positive for a nondegenerate circle.
      real(dp), allocatable, intent(out) :: distance(:) !! Signed distance, negative inside the circle.
      allocate(distance(size(p, 1)))
      if (size(p, 2) /= 2) then
         distance = huge(1.0_dp)
         return
      end if
      distance = sqrt(sum(p*p, dim=2)) - radius
   end subroutine mesh_dcircle

   pure subroutine mesh_dsphere(p, radius, distance)
      real(dp), intent(in) :: p(:, :) !! Points in any dimension at which sphere distance is evaluated.
      real(dp), intent(in) :: radius !! Sphere radius, positive for a nondegenerate sphere.
      real(dp), allocatable, intent(out) :: distance(:) !! Signed distance, negative inside the sphere.
      allocate(distance(size(p, 1)))
      distance = sqrt(sum(p*p, dim=2)) - radius
   end subroutine mesh_dsphere

   pure subroutine mesh_drectangle(p, x1, y1, x2, y2, distance)
      real(dp), intent(in) :: p(:, :) !! Two-dimensional points at which rectangle distance is evaluated.
      real(dp), intent(in) :: x1 !! Lower x coordinate of the rectangle.
      real(dp), intent(in) :: y1 !! Lower y coordinate of the rectangle.
      real(dp), intent(in) :: x2 !! Upper x coordinate of the rectangle.
      real(dp), intent(in) :: y2 !! Upper y coordinate of the rectangle.
      real(dp), allocatable, intent(out) :: distance(:) !! Signed Euclidean distance, negative inside.
      real(dp) :: dx, dy, inside
      integer :: i

      allocate(distance(size(p, 1)))
      if (size(p, 2) /= 2) then
         distance = huge(1.0_dp)
         return
      end if
      do i = 1, size(p, 1)
         dx = max(max(x1 - p(i, 1), 0.0_dp), p(i, 1) - x2)
         dy = max(max(y1 - p(i, 2), 0.0_dp), p(i, 2) - y2)
         if (dx > 0.0_dp .or. dy > 0.0_dp) then
            distance(i) = sqrt(dx*dx + dy*dy)
         else
            inside = min(min(p(i, 1) - x1, x2 - p(i, 1)), min(p(i, 2) - y1, y2 - p(i, 2)))
            distance(i) = -inside
         end if
      end do
   end subroutine mesh_drectangle

   pure subroutine mesh_hunif(p, desired_length)
      real(dp), intent(in) :: p(:, :) !! Point matrix whose number of rows determines output length.
      real(dp), allocatable, intent(out) :: desired_length(:) !! Unit desired edge length at every point.
      allocate(desired_length(size(p, 1)))
      desired_length = 1.0_dp
   end subroutine mesh_hunif

   pure subroutine mesh_diff_values(distance_a, distance_b, distance)
      real(dp), intent(in) :: distance_a(:) !! Signed distances to region A.
      real(dp), intent(in) :: distance_b(:) !! Signed distances to region B.
      real(dp), allocatable, intent(out) :: distance(:) !! Signed distances for A minus B.
      allocate(distance(min(size(distance_a), size(distance_b))))
      distance = max(distance_a(1:size(distance)), -distance_b(1:size(distance)))
   end subroutine mesh_diff_values

   pure subroutine mesh_union_values(distance_a, distance_b, distance)
      real(dp), intent(in) :: distance_a(:) !! Signed distances to region A.
      real(dp), intent(in) :: distance_b(:) !! Signed distances to region B.
      real(dp), allocatable, intent(out) :: distance(:) !! Signed distances for the union of A and B.
      allocate(distance(min(size(distance_a), size(distance_b))))
      distance = min(distance_a(1:size(distance)), distance_b(1:size(distance)))
   end subroutine mesh_union_values

   pure subroutine mesh_intersect_values(distance_a, distance_b, distance)
      real(dp), intent(in) :: distance_a(:) !! Signed distances to region A.
      real(dp), intent(in) :: distance_b(:) !! Signed distances to region B.
      real(dp), allocatable, intent(out) :: distance(:) !! Signed distances for intersection of A and B.
      allocate(distance(min(size(distance_a), size(distance_b))))
      distance = max(distance_a(1:size(distance)), distance_b(1:size(distance)))
   end subroutine mesh_intersect_values

   subroutine rbox_points(n, d, bound, corner, points)
      integer, intent(in) :: n !! Number of random points uniformly generated in the hypercube.
      integer, intent(in) :: d !! Number of spatial dimensions.
      real(dp), intent(in) :: bound !! Half-width of the random-point hypercube.
      real(dp), intent(in) :: corner !! Half-width of deterministic corner cube; negative disables corners.
      real(dp), allocatable, intent(out) :: points(:, :) !! Generated corner points followed by random points.
      integer :: ncorners, i, j, bit
      real(dp), allocatable :: u(:, :)

      if (d < 1 .or. n < 0) then
         allocate(points(0, 0))
         return
      end if
      ncorners = 0
      if (corner >= 0.0_dp .and. d < bit_size(ncorners) - 1) ncorners = 2**d
      allocate(points(ncorners + n, d))
      do i = 1, ncorners
         do j = 1, d
            bit = ibits(i - 1, j - 1, 1)
            if (bit == 0) then
               points(i, j) = -corner
            else
               points(i, j) = corner
            end if
         end do
      end do
      if (n > 0) then
         allocate(u(n, d))
         call random_number(u)
         points(ncorners + 1:, :) = -bound + 2.0_dp*bound*u
      end if
   end subroutine rbox_points

   pure logical function row_less(a, b) result(value)
      real(dp), intent(in) :: a(:) !! First numeric row in a lexicographic comparison.
      real(dp), intent(in) :: b(:) !! Second numeric row in a lexicographic comparison.
      integer :: j
      value = .false.
      do j = 1, min(size(a), size(b))
         if (a(j) < b(j)) then
            value = .true.
            return
         else if (a(j) > b(j)) then
            return
         end if
      end do
   end function row_less

end module geometry_basic
