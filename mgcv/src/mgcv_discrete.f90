! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_discrete
   use mgcv_kinds, only : dp
   implicit none
   private
   public :: xwxd, xwyd, xbd, diag_xvxd, ij_xvxd

contains

   function xwxd(x, w) result(value)
      real(dp), intent(in) :: x(:, :), w(:)
      real(dp), allocatable :: value(:, :), wx(:, :)
      if (size(x, 1) /= size(w)) then; allocate(value(0, 0)); return; end if
      allocate(wx(size(x, 1), size(x, 2)))
      wx = x * spread(w, 2, size(x, 2))
      value = matmul(transpose(x), wx)
   end function xwxd

   function xwyd(x, w, y) result(value)
      real(dp), intent(in) :: x(:, :), w(:), y(:)
      real(dp), allocatable :: value(:)
      if (size(x, 1) /= size(w) .or. size(y) /= size(w)) then; allocate(value(0)); return; end if
      value = matmul(transpose(x), w * y)
   end function xwyd

   function xbd(x, beta) result(value)
      real(dp), intent(in) :: x(:, :), beta(:)
      real(dp), allocatable :: value(:)
      if (size(x, 2) /= size(beta)) then; allocate(value(0)); return; end if
      value = matmul(x, beta)
   end function xbd

   function diag_xvxd(x, v, d) result(value)
      real(dp), intent(in) :: x(:, :), v(:, :)
      real(dp), intent(in), optional :: d(:)
      real(dp), allocatable :: value(:), xv(:, :)
      integer :: i
      if (size(v, 1) /= size(x, 2) .or. size(v, 2) /= size(x, 2)) then
         allocate(value(0)); return
      end if
      xv = matmul(x, v)
      allocate(value(size(x, 1)))
      do i = 1, size(x, 1)
         value(i) = dot_product(xv(i, :), x(i, :))
      end do
      if (present(d)) then
         if (size(d) == size(value)) value = value * d
      end if
   end function diag_xvxd

   function ij_xvxd(xi, xj, v, d) result(value)
      real(dp), intent(in) :: xi(:, :), xj(:, :), v(:, :)
      real(dp), intent(in), optional :: d(:)
      real(dp), allocatable :: value(:), xiv(:, :)
      integer :: i
      if (size(xi, 1) /= size(xj, 1) .or. size(xi, 2) /= size(v, 1) .or. &
          size(xj, 2) /= size(v, 2)) then
         allocate(value(0)); return
      end if
      xiv = matmul(xi, v)
      allocate(value(size(xi, 1)))
      do i = 1, size(xi, 1)
         value(i) = dot_product(xiv(i, :), xj(i, :))
      end do
      if (present(d)) then
         if (size(d) == size(value)) value = value * d
      end if
   end function ij_xvxd

end module mgcv_discrete
