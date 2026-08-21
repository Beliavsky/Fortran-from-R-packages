! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_matrix
   use mgcv_kinds, only : dp
   use mgcv_linalg, only : symmetric_root, jacobi_eigen, matrix_rank
   implicit none
   private
   public :: mroot, mini_root, rrank, slanczos

contains

   subroutine mroot(a, root, rank, status, tol)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: root(:, :)
      integer, intent(out) :: rank, status
      real(dp), intent(in), optional :: tol
      call symmetric_root(a, root, rank, status, tol)
   end subroutine mroot

   subroutine mini_root(s, root, rank, status, tol)
      real(dp), intent(in) :: s(:, :)
      real(dp), allocatable, intent(out) :: root(:, :)
      integer, intent(out) :: rank, status
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: transposed(:, :)
      call symmetric_root(s, transposed, rank, status, tol)
      if (status == 0) then
         allocate(root(rank, size(s, 1))); root = transpose(transposed)
      else
         allocate(root(0, 0))
      end if
   end subroutine mini_root

   integer function rrank(a, tol) result(rank)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in), optional :: tol
      rank = matrix_rank(a, tol)
   end function rrank

   subroutine slanczos(a, k, values, vectors, status, largest)
      real(dp), intent(in) :: a(:, :)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: values(:), vectors(:, :)
      integer, intent(out) :: status
      logical, intent(in), optional :: largest
      real(dp), allocatable :: all_values(:), all_vectors(:, :)
      logical :: high
      integer :: kk, n
      high = .true.; if (present(largest)) high = largest
      call jacobi_eigen(a, all_values, all_vectors, status)
      if (status /= 0) then; allocate(values(0), vectors(0, 0)); return; end if
      n = size(all_values); kk = min(max(0, k), n)
      allocate(values(kk), vectors(n, kk))
      if (high) then
         values = all_values(1:kk); vectors = all_vectors(:, 1:kk)
      else
         values = all_values(n:n - kk + 1:-1); vectors = all_vectors(:, n:n - kk + 1:-1)
      end if
   end subroutine slanczos

end module mgcv_matrix
