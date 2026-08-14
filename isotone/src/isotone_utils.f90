module isotone_utils
   use isotone_kinds, only : dp
   use isotone_linalg, only : least_squares
   implicit none
   private
   public :: weighted_mean_value, weighted_median_value, weighted_fractile_value
   public :: weighted_midrange_value, constraint_values, component_labels
   public :: lagrange_multipliers, transpose_constraint_product, kkt_values
   public :: sort_real_indices, same_real
contains
   logical function same_real(a, b) result(same)
      real(dp), intent(in) :: a, b
      same = (.not. (a < b)) .and. (.not. (a > b))
   end function same_real

   real(dp) function weighted_mean_value(x, w) result(v)
      real(dp), intent(in) :: x(:), w(:)
      real(dp) :: sw
      sw = sum(w)
      if (abs(sw) <= tiny(1.0_dp)) then
         v = sum(x) / real(max(1,size(x)),dp)
      else
         v = dot_product(x,w) / sw
      end if
   end function weighted_mean_value

   subroutine sort_real_indices(x, idx, secondary)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: idx(:)
      real(dp), intent(in), optional :: secondary(:)
      integer :: i, j, key
      logical :: less
      idx = [(i, i=1,size(x))]
      do i = 2, size(idx)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            less = x(key) < x(idx(j))
            if (present(secondary)) then
               if (same_real(x(key),x(idx(j)))) less = secondary(key) < secondary(idx(j))
            end if
            if (.not. less) exit
            idx(j+1) = idx(j)
            j = j - 1
         end do
         idx(j+1) = key
      end do
   end subroutine sort_real_indices

   real(dp) function weighted_median_value(x, w) result(v)
      real(dp), intent(in) :: x(:), w(:)
      integer, allocatable :: idx(:)
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: low, total, dprev, dcur
      integer :: n, k
      n = size(x)
      if (n == 0) then
         v = 0.0_dp; return
      end if
      allocate(idx(n), xs(n), ws(n))
      call sort_real_indices(x, idx)
      xs = x(idx); ws = w(idx)
      total = sum(ws)
      low = 0.0_dp
      dprev = -total
      do k = 1, n
         low = low + ws(k)
         dcur = 2.0_dp * low - total
         if (dcur >= 0.0_dp) then
            if (abs(dcur) <= 10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(total)) .and. k < n) then
               v = (ws(k+1)*xs(k+1) + ws(k)*xs(k)) / max(tiny(1.0_dp), ws(k+1)+ws(k))
            else if (dprev < 0.0_dp) then
               v = xs(k)
            else
               v = xs(max(1,k-1))
            end if
            return
         end if
         dprev = dcur
      end do
      v = xs(n)
   end function weighted_median_value

   real(dp) function weighted_fractile_value(x, w, p) result(v)
      real(dp), intent(in) :: x(:), w(:), p
      integer, allocatable :: idx(:)
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: target, cum, total, tol
      integer :: n, k
      n = size(x)
      if (n == 0) then
         v = 0.0_dp; return
      end if
      allocate(idx(n), xs(n), ws(n))
      call sort_real_indices(x, idx)
      xs = x(idx); ws = w(idx)
      total = sum(ws)
      target = min(1.0_dp,max(0.0_dp,p)) * total
      tol = 10.0_dp * epsilon(1.0_dp) * max(1.0_dp,abs(total))
      cum = 0.0_dp
      do k = 1, n
         cum = cum + ws(k)
         if (cum >= target - tol) then
            if (abs(cum-target) <= tol .and. k < n) then
               v = (ws(k+1)*xs(k+1) + ws(k)*xs(k)) / max(tiny(1.0_dp), ws(k+1)+ws(k))
            else
               v = xs(k)
            end if
            return
         end if
      end do
      v = xs(n)
   end function weighted_fractile_value

   real(dp) function weighted_midrange_value(x, w) result(v)
      real(dp), intent(in) :: x(:), w(:)
      real(dp) :: s, t
      integer :: n, i, j, i0, j0
      n = size(x)
      if (n <= 1) then
         if (n == 1) then
            v = x(1)
         else
            v = 0.0_dp
         end if
         return
      end if
      s = -1.0_dp; i0 = 1; j0 = 2
      do i = 1, n - 1
         do j = i + 1, n
            t = w(i)*w(j)*abs(x(i)-x(j))/max(tiny(1.0_dp),w(i)+w(j))
            if (t > s) then
               s = t; i0 = i; j0 = j
            end if
         end do
      end do
      v = (w(i0)*x(i0)+w(j0)*x(j0))/max(tiny(1.0_dp),w(i0)+w(j0))
   end function weighted_midrange_value

   subroutine constraint_values(isomat, x, ax)
      integer, intent(in) :: isomat(:,:)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: ax(:)
      integer :: k
      do k = 1, size(isomat,1)
         ax(k) = x(isomat(k,2)) - x(isomat(k,1))
      end do
   end subroutine constraint_values

   subroutine component_labels(n, isomat, active, labels, ncomp)
      integer, intent(in) :: n, isomat(:,:), active(:)
      integer, intent(out) :: labels(:)
      integer, intent(out) :: ncomp
      integer, allocatable :: parent(:), rootmap(:)
      integer :: i, k, r1, r2, r
      parent = [(i, i=1,n)]
      do k = 1, size(active)
         if (active(k) <= 0) cycle
         r1 = find_root(parent,isomat(active(k),1))
         r2 = find_root(parent,isomat(active(k),2))
         if (r1 /= r2) parent(r2) = r1
      end do
      allocate(rootmap(n)); rootmap = 0
      ncomp = 0
      do i = 1, n
         r = find_root(parent,i)
         if (rootmap(r) == 0) then
            ncomp = ncomp + 1
            rootmap(r) = ncomp
         end if
         labels(i) = rootmap(r)
      end do
   contains
      integer function find_root(p, x) result(r0)
         integer, intent(inout) :: p(:)
         integer, intent(in) :: x
         integer :: q, nxt
         q = x
         do while (p(q) /= q)
            q = p(q)
         end do
         r0 = q
         q = x
         do while (p(q) /= q)
            nxt = p(q); p(q) = r0; q = nxt
         end do
      end function find_root
   end subroutine component_labels

   subroutine lagrange_multipliers(isomat, active, grad, lambda, ok)
      integer, intent(in) :: isomat(:,:), active(:)
      real(dp), intent(in) :: grad(:)
      real(dp), allocatable, intent(out) :: lambda(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: at(:,:)
      integer :: n, m, k, e
      n = size(grad); m = size(active)
      allocate(lambda(m))
      if (m == 0) then
         ok = .true.; return
      end if
      allocate(at(n,m)); at = 0.0_dp
      do k = 1, m
         e = active(k)
         at(isomat(e,2),k) = 1.0_dp
         at(isomat(e,1),k) = -1.0_dp
      end do
      call least_squares(at, grad, lambda, ok)
   end subroutine lagrange_multipliers

   subroutine transpose_constraint_product(isomat, lambda, h)
      integer, intent(in) :: isomat(:,:)
      real(dp), intent(in) :: lambda(:)
      real(dp), intent(out) :: h(:)
      integer :: k
      h = 0.0_dp
      do k = 1, size(lambda)
         h(isomat(k,2)) = h(isomat(k,2)) + lambda(k)
         h(isomat(k,1)) = h(isomat(k,1)) - lambda(k)
      end do
   end subroutine transpose_constraint_product

   subroutine kkt_values(grad, ax, lambda, alambda, values)
      real(dp), intent(in) :: grad(:), ax(:), lambda(:), alambda(:)
      real(dp), intent(out) :: values(4)
      values(1) = maxval(abs(grad-alambda))
      values(2) = minval(ax)
      values(3) = minval(lambda)
      values(4) = dot_product(ax,lambda)
   end subroutine kkt_values
end module isotone_utils
