! SPDX-License-Identifier: GPL-3.0-only
module scs_sparse
   use scs_kinds, only : dp, i4
   use scs_types, only : scs_matrix
   implicit none
   private
   public :: dense_to_csc, dense_upper_to_csc, csc_to_dense
   public :: accum_by_a, accum_by_atrans, accum_by_p, validate_matrix
contains

   subroutine dense_to_csc(a, s, tol)
      real(dp), intent(in) :: a(:,:)
      type(scs_matrix), intent(out) :: s
      real(dp), intent(in), optional :: tol
      real(dp) :: ztol
      integer :: r, c, k, nnz
      ztol = 0.0_dp
      if (present(tol)) ztol = max(0.0_dp, tol)
      s%m = int(size(a,1), i4)
      s%n = int(size(a,2), i4)
      nnz = count(abs(a) > ztol)
      allocate(s%x(nnz), s%i(nnz), s%p(s%n+1))
      k = 0
      s%p(1) = 1_i4
      do c = 1, size(a,2)
         do r = 1, size(a,1)
            if (abs(a(r,c)) > ztol) then
               k = k + 1
               s%x(k) = a(r,c)
               s%i(k) = int(r, i4)
            end if
         end do
         s%p(c+1) = int(k+1, i4)
      end do
   end subroutine dense_to_csc

   subroutine dense_upper_to_csc(a, s, tol)
      real(dp), intent(in) :: a(:,:)
      type(scs_matrix), intent(out) :: s
      real(dp), intent(in), optional :: tol
      real(dp) :: ztol
      integer :: r, c, k, nnz
      if (size(a,1) /= size(a,2)) error stop 'dense_upper_to_csc: matrix must be square'
      ztol = 0.0_dp
      if (present(tol)) ztol = max(0.0_dp, tol)
      nnz = 0
      do c = 1, size(a,2)
         do r = 1, c
            if (abs(a(r,c)) > ztol) nnz = nnz + 1
         end do
      end do
      s%m = int(size(a,1), i4)
      s%n = int(size(a,2), i4)
      allocate(s%x(nnz), s%i(nnz), s%p(s%n+1))
      k = 0
      s%p(1) = 1_i4
      do c = 1, size(a,2)
         do r = 1, c
            if (abs(a(r,c)) > ztol) then
               k = k + 1
               s%x(k) = a(r,c)
               s%i(k) = int(r, i4)
            end if
         end do
         s%p(c+1) = int(k+1, i4)
      end do
   end subroutine dense_upper_to_csc

   function csc_to_dense(s, symmetric_upper) result(a)
      type(scs_matrix), intent(in) :: s
      logical, intent(in), optional :: symmetric_upper
      real(dp), allocatable :: a(:,:)
      logical :: sym
      integer :: c, k, r
      sym = .false.
      if (present(symmetric_upper)) sym = symmetric_upper
      allocate(a(s%m,s%n), source=0.0_dp)
      do c = 1, s%n
         do k = s%p(c), s%p(c+1)-1
            r = s%i(k)
            a(r,c) = a(r,c) + s%x(k)
            if (sym .and. r /= c) a(c,r) = a(c,r) + s%x(k)
         end do
      end do
   end function csc_to_dense

   subroutine accum_by_a(a, x, y)
      type(scs_matrix), intent(in) :: a
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:)
      integer :: c, k
      do c = 1, a%n
         do k = a%p(c), a%p(c+1)-1
            y(a%i(k)) = y(a%i(k)) + a%x(k) * x(c)
         end do
      end do
   end subroutine accum_by_a

   subroutine accum_by_atrans(a, x, y)
      type(scs_matrix), intent(in) :: a
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:)
      integer :: c, k
      do c = 1, a%n
         do k = a%p(c), a%p(c+1)-1
            y(c) = y(c) + a%x(k) * x(a%i(k))
         end do
      end do
   end subroutine accum_by_atrans

   subroutine accum_by_p(p, x, y)
      type(scs_matrix), intent(in) :: p
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:)
      integer :: c, k, r
      do c = 1, p%n
         do k = p%p(c), p%p(c+1)-1
            r = p%i(k)
            y(c) = y(c) + p%x(k) * x(r)
            if (r /= c) y(r) = y(r) + p%x(k) * x(c)
         end do
      end do
   end subroutine accum_by_p

   logical function validate_matrix(a, symmetric_upper) result(ok)
      type(scs_matrix), intent(in) :: a
      logical, intent(in), optional :: symmetric_upper
      logical :: sym
      integer :: c, k
      ok = .false.
      sym = .false.
      if (present(symmetric_upper)) sym = symmetric_upper
      if (.not. allocated(a%x) .or. .not. allocated(a%i) .or. .not. allocated(a%p)) return
      if (size(a%p) /= a%n + 1) return
      if (a%p(1) /= 1 .or. a%p(a%n+1) /= size(a%x)+1) return
      if (size(a%i) /= size(a%x)) return
      if (any(a%i < 1) .or. any(a%i > a%m)) return
      if (any(a%p(2:) < a%p(:a%n))) return
      if (sym) then
         if (a%m /= a%n) return
         do c = 1, a%n
            do k = a%p(c), a%p(c+1)-1
               if (a%i(k) > c) return
            end do
         end do
      end if
      ok = .true.
   end function validate_matrix
end module scs_sparse
