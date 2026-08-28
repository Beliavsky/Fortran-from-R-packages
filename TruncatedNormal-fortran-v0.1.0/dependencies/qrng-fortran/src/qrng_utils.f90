! Computational utilities translated from qrng 0.0-11 R/auxiliaries.R.
! Upstream qrng license: GPL-2 | GPL-3.
module qrng_utils_mod
use r_mod, only: dp
implicit none
private
public :: to_array_matrix, to_array_3d

contains

function to_array_matrix(x, f) result(y)
real(dp), intent(in) :: x(:,:)
integer, intent(in) :: f
real(dp), allocatable :: y(:,:)
integer :: n, d, i, k, j
if (f < 1) error stop "to_array_matrix: f must be at least 1"
n = size(x,1)
if (mod(size(x,2),f) /= 0) error stop "to_array_matrix: f must divide size(x,2)"
d = size(x,2) / f
allocate(y(n*f,d))
do i = 1, n
   do k = 1, f
      do j = 1, d
         y((i-1)*f+k,j) = x(i,(k-1)*d+j)
      end do
   end do
end do
end function to_array_matrix

function to_array_3d(x, f) result(y)
real(dp), intent(in) :: x(:,:)
integer, intent(in) :: f
real(dp), allocatable :: y(:,:,:)
integer :: n, d, i, k, j
if (f < 1) error stop "to_array_3d: f must be at least 1"
n = size(x,1)
if (mod(size(x,2),f) /= 0) error stop "to_array_3d: f must divide size(x,2)"
d = size(x,2) / f
allocate(y(n,f,d))
do j = 1, d
   do k = 1, f
      do i = 1, n
         y(i,k,j) = x(i,(j-1)*f+k)
      end do
   end do
end do
end function to_array_3d

end module qrng_utils_mod
