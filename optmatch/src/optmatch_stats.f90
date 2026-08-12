! SPDX-License-Identifier: MIT
module optmatch_stats
   use optmatch_kinds, only : dp
   implicit none
   private
   public :: mean_value, sample_variance, covariance_matrix, rank_columns
   public :: symmetric_pseudoinverse, median_value, mad_value

contains

pure real(dp) function mean_value(x) result(ans)
   real(dp), intent(in) :: x(:)
   if (size(x) == 0) error stop 'optmatch: mean of empty vector'
   ans = sum(x) / real(size(x), dp)
end function mean_value

pure real(dp) function sample_variance(x) result(ans)
   real(dp), intent(in) :: x(:)
   real(dp) :: m
   if (size(x) < 2) then
      ans = 0.0_dp
      return
   end if
   m = mean_value(x)
   ans = sum((x - m)**2) / real(size(x) - 1, dp)
end function sample_variance

function median_value(x) result(ans)
   real(dp), intent(in) :: x(:)
   real(dp) :: ans
   real(dp), allocatable :: y(:)
   integer :: n
   n = size(x)
   if (n == 0) error stop 'optmatch: median of empty vector'
   y = x
   call insertion_sort(y)
   if (mod(n, 2) == 1) then
      ans = y((n + 1) / 2)
   else
      ans = 0.5_dp * (y(n / 2) + y(n / 2 + 1))
   end if
end function median_value

function mad_value(x) result(ans)
   real(dp), intent(in) :: x(:)
   real(dp) :: ans, med
   real(dp), allocatable :: dev(:)
   med = median_value(x)
   allocate(dev(size(x)))
   dev = abs(x - med)
   ans = 1.482602218505602_dp * median_value(dev)
end function mad_value

subroutine covariance_matrix(x, cov)
   real(dp), intent(in) :: x(:, :)
   real(dp), allocatable, intent(out) :: cov(:, :)
   real(dp), allocatable :: center(:)
   integer :: n, p, i, j
   n = size(x, 1)
   p = size(x, 2)
   allocate(cov(p, p), center(p))
   if (n < 2) then
      cov = 0.0_dp
      return
   end if
   do j = 1, p
      center(j) = sum(x(:, j)) / real(n, dp)
   end do
   cov = 0.0_dp
   do i = 1, p
      do j = i, p
         cov(i, j) = sum((x(:, i) - center(i)) * (x(:, j) - center(j))) / real(n - 1, dp)
         cov(j, i) = cov(i, j)
      end do
   end do
end subroutine covariance_matrix

subroutine rank_columns(x, ranks, any_ties)
   real(dp), intent(in) :: x(:, :)
   real(dp), allocatable, intent(out) :: ranks(:, :)
   logical, intent(out) :: any_ties
   integer :: n, p, j
   n = size(x, 1)
   p = size(x, 2)
   allocate(ranks(n, p))
   any_ties = .false.
   do j = 1, p
      call rank_vector(x(:, j), ranks(:, j), any_ties)
   end do
end subroutine rank_columns

subroutine rank_vector(x, ranks, any_ties)
   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: ranks(:)
   logical, intent(inout) :: any_ties
   integer, allocatable :: idx(:)
   integer :: n, i, j, k, tmp
   real(dp) :: r
   n = size(x)
   allocate(idx(n))
   idx = [(i, i = 1, n)]
   do i = 2, n
      tmp = idx(i)
      j = i - 1
      do while (j >= 1)
         if (x(idx(j)) <= x(tmp)) exit
         idx(j + 1) = idx(j)
         j = j - 1
      end do
      idx(j + 1) = tmp
   end do
   i = 1
   do while (i <= n)
      k = i
      do while (k < n)
         if (abs(x(idx(k + 1)) - x(idx(i))) > 0.0_dp) exit
         k = k + 1
      end do
      if (k > i) any_ties = .true.
      r = 0.5_dp * real(i + k, dp)
      do j = i, k
         ranks(idx(j)) = r
      end do
      i = k + 1
   end do
end subroutine rank_vector

subroutine symmetric_pseudoinverse(a, ainv, tol)
   real(dp), intent(in) :: a(:, :)
   real(dp), allocatable, intent(out) :: ainv(:, :)
   real(dp), intent(in), optional :: tol
   real(dp), allocatable :: v(:, :), d(:), work(:, :)
   real(dp) :: threshold, dmax
   integer :: n, i, j
   if (size(a, 1) /= size(a, 2)) error stop 'optmatch: pseudoinverse requires square matrix'
   n = size(a, 1)
   allocate(v(n, n), d(n), work(n, n), ainv(n, n))
   work = 0.5_dp * (a + transpose(a))
   call jacobi_symmetric(work, d, v)
   dmax = maxval(abs(d))
   if (present(tol)) then
      threshold = tol * max(1.0_dp, dmax)
   else
      threshold = 1.0e-10_dp * max(1.0_dp, dmax)
   end if
   ainv = 0.0_dp
   do i = 1, n
      if (abs(d(i)) > threshold) then
         do j = 1, n
            ainv(:, j) = ainv(:, j) + v(:, i) * v(j, i) / d(i)
         end do
      end if
   end do
end subroutine symmetric_pseudoinverse

subroutine jacobi_symmetric(a, d, v)
   real(dp), intent(inout) :: a(:, :)
   real(dp), intent(out) :: d(:)
   real(dp), intent(out) :: v(:, :)
   integer :: n, iter, p, q, i
   real(dp) :: maxoff, app, aqq, apq, phi, c, s, aip, aiq, vip, viq
   n = size(a, 1)
   v = 0.0_dp
   do i = 1, n
      v(i, i) = 1.0_dp
   end do
   if (n <= 1) then
      if (n == 1) d(1) = a(1, 1)
      return
   end if
   do iter = 1, max(50, 50 * n * n)
      maxoff = 0.0_dp
      p = 1
      q = 2
      do i = 1, n - 1
         call largest_in_row(a, i, maxoff, p, q)
      end do
      if (maxoff <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) exit
      app = a(p, p)
      aqq = a(q, q)
      apq = a(p, q)
      phi = 0.5_dp * atan2(2.0_dp * apq, aqq - app)
      c = cos(phi)
      s = sin(phi)
      do i = 1, n
         if (i /= p .and. i /= q) then
            aip = a(i, p)
            aiq = a(i, q)
            a(i, p) = c * aip - s * aiq
            a(p, i) = a(i, p)
            a(i, q) = s * aip + c * aiq
            a(q, i) = a(i, q)
         end if
      end do
      a(p, p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
      a(q, q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
      a(p, q) = 0.0_dp
      a(q, p) = 0.0_dp
      do i = 1, n
         vip = v(i, p)
         viq = v(i, q)
         v(i, p) = c * vip - s * viq
         v(i, q) = s * vip + c * viq
      end do
   end do
   do i = 1, n
      d(i) = a(i, i)
   end do
end subroutine jacobi_symmetric

subroutine largest_in_row(a, irow, maxoff, p, q)
   real(dp), intent(in) :: a(:, :)
   integer, intent(in) :: irow
   real(dp), intent(inout) :: maxoff
   integer, intent(inout) :: p, q
   integer :: j
   do j = irow + 1, size(a, 2)
      if (abs(a(irow, j)) > maxoff) then
         maxoff = abs(a(irow, j))
         p = irow
         q = j
      end if
   end do
end subroutine largest_in_row

subroutine insertion_sort(x)
   real(dp), intent(inout) :: x(:)
   integer :: i, j
   real(dp) :: key
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
end subroutine insertion_sort

end module optmatch_stats
