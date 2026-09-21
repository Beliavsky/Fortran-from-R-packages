module geometry_linalg
   use geometry_kinds, only : dp
   implicit none
   private

   public :: cross3, determinant, solve_linear, hyperplane_normal
   public :: simplex_measure, simplex_facet_measure, next_combination
   public :: combination_count, sort_integer, sort_real_with_index
   public :: lexicographic_less, row_norm

contains

   pure function row_norm(x) result(value)
      real(dp), intent(in) :: x(:) !! Vector whose Euclidean norm is required.
      real(dp) :: value
      value = sqrt(sum(x*x))
   end function row_norm

   pure function cross3(x, y) result(z)
      real(dp), intent(in) :: x(3) !! First three-dimensional vector.
      real(dp), intent(in) :: y(3) !! Second three-dimensional vector.
      real(dp) :: z(3)
      z(1) = x(2)*y(3) - x(3)*y(2)
      z(2) = x(3)*y(1) - x(1)*y(3)
      z(3) = x(1)*y(2) - x(2)*y(1)
   end function cross3

   pure function determinant(a) result(det)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose determinant is required.
      real(dp) :: det
      real(dp), allocatable :: work(:, :), tmp(:)
      real(dp) :: pivot, factor
      integer :: n, i, j, k, ipiv

      n = size(a, 1)
      if (n /= size(a, 2)) then
         det = 0.0_dp
         return
      end if
      if (n == 0) then
         det = 1.0_dp
         return
      end if
      allocate(work(n, n), tmp(n))
      work = a
      det = 1.0_dp
      do k = 1, n
         ipiv = k
         pivot = abs(work(k, k))
         do i = k + 1, n
            if (abs(work(i, k)) > pivot) then
               pivot = abs(work(i, k))
               ipiv = i
            end if
         end do
         if (pivot <= epsilon(1.0_dp)*max(1.0_dp, maxval(abs(work)))) then
            det = 0.0_dp
            return
         end if
         if (ipiv /= k) then
            tmp = work(k, :)
            work(k, :) = work(ipiv, :)
            work(ipiv, :) = tmp
            det = -det
         end if
         det = det*work(k, k)
         do i = k + 1, n
            factor = work(i, k)/work(k, k)
            work(i, k) = 0.0_dp
            do j = k + 1, n
               work(i, j) = work(i, j) - factor*work(k, j)
            end do
         end do
      end do
   end function determinant

   pure subroutine solve_linear(a, b, x, ok)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      real(dp), intent(out) :: x(:) !! Solution vector when the system is nonsingular.
      logical, intent(out) :: ok !! True when a stable nonsingular solve was completed.
      real(dp), allocatable :: work(:, :), rhs(:), tmp_row(:)
      real(dp) :: pivot, factor, rhs_tmp, scale
      integer :: n, i, j, k, ipiv

      n = size(a, 1)
      ok = .false.
      x = 0.0_dp
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) return
      allocate(work(n, n), rhs(n), tmp_row(n))
      work = a
      rhs = b
      scale = max(1.0_dp, maxval(abs(work)))
      do k = 1, n
         ipiv = k
         pivot = abs(work(k, k))
         do i = k + 1, n
            if (abs(work(i, k)) > pivot) then
               pivot = abs(work(i, k))
               ipiv = i
            end if
         end do
         if (pivot <= 100.0_dp*epsilon(1.0_dp)*scale) return
         if (ipiv /= k) then
            tmp_row = work(k, :)
            work(k, :) = work(ipiv, :)
            work(ipiv, :) = tmp_row
            rhs_tmp = rhs(k)
            rhs(k) = rhs(ipiv)
            rhs(ipiv) = rhs_tmp
         end if
         do i = k + 1, n
            factor = work(i, k)/work(k, k)
            work(i, k) = 0.0_dp
            do j = k + 1, n
               work(i, j) = work(i, j) - factor*work(k, j)
            end do
            rhs(i) = rhs(i) - factor*rhs(k)
         end do
      end do
      do i = n, 1, -1
         x(i) = rhs(i)
         if (i < n) x(i) = x(i) - dot_product(work(i, i + 1:n), x(i + 1:n))
         x(i) = x(i)/work(i, i)
      end do
      ok = .true.
   end subroutine solve_linear

   pure subroutine hyperplane_normal(points, normal, offset, ok)
      real(dp), intent(in) :: points(:, :) !! D-by-(D-1 simplex) rows defining a candidate hyperplane.
      real(dp), intent(out) :: normal(:) !! Unit normal vector to the candidate hyperplane.
      real(dp), intent(out) :: offset !! Hyperplane offset in normal dot x plus offset equals zero.
      logical, intent(out) :: ok !! True when the supplied points are affinely independent.
      real(dp), allocatable :: edges(:, :), minor(:, :)
      real(dp) :: nrm
      integer :: d, j, jj, col

      d = size(points, 2)
      normal = 0.0_dp
      offset = 0.0_dp
      ok = .false.
      if (size(points, 1) /= d .or. size(normal) /= d) return
      if (d == 1) then
         normal(1) = 1.0_dp
         offset = -points(1, 1)
         ok = .true.
         return
      end if
      allocate(edges(d - 1, d), minor(d - 1, d - 1))
      do jj = 1, d - 1
         edges(jj, :) = points(jj + 1, :) - points(1, :)
      end do
      do j = 1, d
         col = 0
         do jj = 1, d
            if (jj == j) cycle
            col = col + 1
            minor(:, col) = edges(:, jj)
         end do
         normal(j) = determinant(minor)
         if (mod(j + 1, 2) /= 0) normal(j) = -normal(j)
      end do
      nrm = sqrt(sum(normal*normal))
      if (nrm <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(points)))) return
      normal = normal/nrm
      offset = -dot_product(normal, points(1, :))
      ok = .true.
   end subroutine hyperplane_normal

   pure function simplex_measure(vertices) result(value)
      real(dp), intent(in) :: vertices(:, :) !! (D+1)-by-D simplex vertex coordinates.
      real(dp) :: value
      real(dp), allocatable :: edges(:, :)
      integer :: d, i

      d = size(vertices, 2)
      value = 0.0_dp
      if (size(vertices, 1) /= d + 1) return
      allocate(edges(d, d))
      do i = 1, d
         edges(i, :) = vertices(i + 1, :) - vertices(1, :)
      end do
      value = abs(determinant(edges))/real(factorial_int(d), dp)
   end function simplex_measure

   pure function simplex_facet_measure(vertices) result(value)
      real(dp), intent(in) :: vertices(:, :) !! D-by-D coordinates of a (D-1)-simplex embedded in D-space.
      real(dp) :: value
      real(dp), allocatable :: edges(:, :), gram(:, :)
      integer :: d, i, j

      d = size(vertices, 2)
      value = 0.0_dp
      if (size(vertices, 1) /= d) return
      if (d == 1) then
         value = 1.0_dp
         return
      end if
      allocate(edges(d - 1, d), gram(d - 1, d - 1))
      do i = 1, d - 1
         edges(i, :) = vertices(i + 1, :) - vertices(1, :)
      end do
      do i = 1, d - 1
         do j = 1, d - 1
            gram(i, j) = dot_product(edges(i, :), edges(j, :))
         end do
      end do
      value = sqrt(max(0.0_dp, determinant(gram)))/real(factorial_int(d - 1), dp)
   end function simplex_facet_measure

   pure integer function factorial_int(n) result(value)
      integer, intent(in) :: n !! Nonnegative integer factorial argument.
      integer :: i
      value = 1
      do i = 2, n
         value = value*i
      end do
   end function factorial_int

   pure integer function combination_count(n, k, cap) result(value)
      integer, intent(in) :: n !! Number of available objects.
      integer, intent(in) :: k !! Number of objects selected per combination.
      integer, intent(in) :: cap !! Maximum count returned before saturation.
      integer(kind=8) :: acc
      integer :: i, kk

      if (k < 0 .or. k > n) then
         value = 0
         return
      end if
      kk = min(k, n - k)
      acc = 1_8
      do i = 1, kk
         acc = acc*int(n - kk + i, 8)/int(i, 8)
         if (acc >= int(cap, 8)) then
            value = cap
            return
         end if
      end do
      value = int(acc)
   end function combination_count

   pure subroutine next_combination(c, n, has_next)
      integer, intent(inout) :: c(:) !! Current increasing combination, advanced in place.
      integer, intent(in) :: n !! Largest admissible element of the combination.
      logical, intent(out) :: has_next !! True when a successor combination exists.
      integer :: i, j, k

      k = size(c)
      has_next = .false.
      do i = k, 1, -1
         if (c(i) < n - k + i) then
            c(i) = c(i) + 1
            if (i < k) c(i + 1:k) = [(c(i) + j, j = 1, k - i)]
            has_next = .true.
            return
         end if
      end do
   end subroutine next_combination

   pure subroutine sort_integer(x)
      integer, intent(inout) :: x(:) !! Integer vector sorted into ascending order in place.
      integer :: i, j, key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_integer

   pure subroutine sort_real_with_index(x, idx)
      real(dp), intent(inout) :: x(:) !! Real keys sorted into ascending order in place.
      integer, intent(inout) :: idx(:) !! Permutation moved consistently with the real keys.
      real(dp) :: key
      integer :: ikey, i, j
      do i = 2, size(x)
         key = x(i)
         ikey = idx(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         x(j + 1) = key
         idx(j + 1) = ikey
      end do
   end subroutine sort_real_with_index

   pure logical function lexicographic_less(a, b) result(value)
      real(dp), intent(in) :: a(:) !! First row compared lexicographically.
      real(dp), intent(in) :: b(:) !! Second row compared lexicographically.
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
      value = size(a) < size(b)
   end function lexicographic_less

end module geometry_linalg
