! SPDX-License-Identifier: GPL-3.0-or-later
module cccp_linalg
   use cccp_kinds, only : dp
   use r_linalg, only : shared_solve_system => solve_system
   use r_linalg, only : shared_spd_inverse_logdet => spd_inverse_logdet
   use r_linalg, only : shared_symmetric_eigen => symmetric_eigen
   use r_linalg, only : shared_symmetrize => symmetrize
   implicit none
   private
   public :: solve_system, symmetric_eigenvalues, spd_inverse_logdet
   public :: equality_particular, vector_norm2, symmetrize

contains

   pure function vector_norm2(x) result(ans)
      real(dp), intent(in) :: x(:)
      real(dp) :: ans
      intrinsic :: norm2

      ans = norm2(x)
   end function vector_norm2

   pure function symmetrize(a) result(s)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: s(size(a,1), size(a,2))
      s = shared_symmetrize(a)
   end function symmetrize

   subroutine solve_system(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:,:)
      integer :: n, i

      n = size(b)
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      allocate(aa(n,n))
      aa = a
      do i = 1, n
         aa(i,i) = aa(i,i) + 1.0e-12_dp * max(1.0_dp, abs(aa(i,i)))
      end do
      call shared_solve_system(aa, b, x, info)
   end subroutine solve_system

   subroutine symmetric_eigenvalues(a, w, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info
      real(dp), allocatable :: values(:), vectors(:,:)
      integer :: n

      n = size(a,1)
      if (size(a,2) /= n .or. size(w) /= n) then
         info = -1
         return
      end if
      call shared_symmetric_eigen(a, values, vectors, info)
      if (info == 0) w = values
   end subroutine symmetric_eigenvalues

   subroutine spd_inverse_logdet(a, ainv, logdet, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(:,:)
      real(dp), intent(out) :: logdet
      integer, intent(out) :: info
      real(dp), allocatable :: shared_inverse(:,:)
      integer :: n

      n = size(a,1)
      if (size(a,2) /= n .or. any(shape(ainv) /= [n,n])) then
         info = -1
         logdet = huge(1.0_dp)
         return
      end if
      call shared_spd_inverse_logdet(a, shared_inverse, logdet, info)
      if (info == 0) ainv = shared_inverse
   end subroutine spd_inverse_logdet

   subroutine equality_particular(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aat(:,:), y(:)
      integer :: p, n

      p = size(a,1)
      n = size(a,2)
      if (size(b) /= p .or. size(x) /= n) then
         info = -1
         return
      end if
      if (p == 0) then
         x = 0.0_dp
         info = 0
         return
      end if
      allocate(aat(p,p), y(p))
      aat = matmul(a, transpose(a))
      call solve_system(aat, b, y, info)
      if (info == 0) x = matmul(transpose(a), y)
   end subroutine equality_particular

end module cccp_linalg
