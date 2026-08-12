module qap_core
   use qap_kinds, only : dp
   implicit none
   private
   public :: qap_obj, qap_swap_delta, qap_is_permutation, qap_validate

contains

   function qap_obj(A, B, perm) result(obj)
      real(dp), intent(in) :: A(:,:), B(:,:)
      integer, intent(in) :: perm(:)
      real(dp) :: obj
      integer :: i, j, n

      n = size(perm)
      if (size(A,1) /= n .or. size(A,2) /= n) error stop "qap_obj: A has wrong shape"
      if (size(B,1) /= n .or. size(B,2) /= n) error stop "qap_obj: B has wrong shape"
      if (.not. qap_is_permutation(perm)) error stop "qap_obj: perm is not a permutation"

      obj = 0.0_dp
      do j = 1, n
         do i = 1, n
            obj = obj + A(i,j) * B(perm(i), perm(j))
         end do
      end do
   end function qap_obj

   function qap_swap_delta(A, B, perm, i1, i2) result(delta)
      real(dp), intent(in) :: A(:,:), B(:,:)
      integer, intent(in) :: perm(:)
      integer, intent(in) :: i1, i2
      real(dp) :: delta
      integer :: j, n, ibild, jbild, kbild

      n = size(perm)
      if (i1 < 1 .or. i1 > n .or. i2 < 1 .or. i2 > n) then
         error stop "qap_swap_delta: swap index out of range"
      end if
      if (i1 == i2) then
         delta = 0.0_dp
         return
      end if

      ibild = perm(i1)
      jbild = perm(i2)
      delta = 0.0_dp
      do j = 1, n
         if (j /= i1 .and. j /= i2) then
            kbild = perm(j)
            delta = delta - (A(i1,j) - A(i2,j)) * &
               (B(ibild,kbild) - B(jbild,kbild))
         end if
      end do
      delta = 2.0_dp * delta - (A(i1,i1) - A(i2,i2)) * &
         (B(ibild,ibild) - B(jbild,jbild))
   end function qap_swap_delta

   function qap_is_permutation(perm) result(ok)
      integer, intent(in) :: perm(:)
      logical :: ok
      logical, allocatable :: seen(:)
      integer :: i, n

      n = size(perm)
      allocate(seen(n), source=.false.)
      ok = .true.
      do i = 1, n
         if (perm(i) < 1 .or. perm(i) > n) then
            ok = .false.
            return
         end if
         if (seen(perm(i))) then
            ok = .false.
            return
         end if
         seen(perm(i)) = .true.
      end do
   end function qap_is_permutation

   subroutine qap_validate(A, B)
      real(dp), intent(in) :: A(:,:), B(:,:)
      integer :: n
      real(dp) :: scale_a, scale_b, tol_a, tol_b

      n = size(A,1)
      if (n < 2) error stop "qap: matrix dimension must be at least 2"
      if (size(A,2) /= n) error stop "qap: A must be square"
      if (size(B,1) /= n .or. size(B,2) /= n) error stop "qap: A and B do not conform"
      if (any(A < 0.0_dp) .or. any(B < 0.0_dp)) then
         error stop "qap: A and B must contain nonnegative values"
      end if

      scale_a = max(1.0_dp, maxval(abs(A)))
      scale_b = max(1.0_dp, maxval(abs(B)))
      tol_a = 100.0_dp * epsilon(1.0_dp) * scale_a
      tol_b = 100.0_dp * epsilon(1.0_dp) * scale_b
      if (maxval(abs(A - transpose(A))) > tol_a) error stop "qap: A must be symmetric"
      if (maxval(abs(B - transpose(B))) > tol_b) error stop "qap: B must be symmetric"
   end subroutine qap_validate

end module qap_core
