! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
module fda_fd
   use r_kinds, only : dp
   use fda_basis, only : basis_type, basis_breaks, eval_basis
   use fda_numeric, only : quadset
   implicit none
   private

   type, public :: fd_type
      type(basis_type) :: basis
      real(dp), allocatable :: coefs(:, :)
   end type fd_type

   public :: make_fd
   public :: eval_fd
   public :: mean_fd
   public :: center_fd
   public :: inprod_basis
   public :: inprod_fd

contains

   pure subroutine make_fd(coefs, basis, fdobj, info)
      real(dp), intent(in) :: coefs(:, :) !! Basis coefficient matrix with basis functions in rows and replications in columns.
      type(basis_type), intent(in) :: basis !! Basis associated with the coefficient rows.
      type(fd_type), intent(out) :: fdobj !! Constructed functional-data object.
      integer, intent(out) :: info !! Zero on success; nonzero when coefficient and basis dimensions disagree.

      info = 0
      if (size(coefs, 1) /= basis%nbasis) then
         info = 1
         return
      end if
      fdobj%basis = basis
      allocate(fdobj%coefs(size(coefs, 1), size(coefs, 2)))
      fdobj%coefs = coefs
   end subroutine make_fd

   pure subroutine eval_fd(x, fdobj, nderiv, values, info)
      real(dp), intent(in) :: x(:) !! Argument values at which all functional-data replications are evaluated.
      type(fd_type), intent(in) :: fdobj !! Functional-data coefficients and basis to evaluate.
      integer, intent(in) :: nderiv !! Nonnegative derivative order applied to the basis before evaluation.
      real(dp), allocatable, intent(out) :: values(:, :) !! Evaluations: argument values by replication.
      integer, intent(out) :: info !! Zero on success; otherwise a basis evaluation or dimension error code.
      real(dp), allocatable :: basismat(:, :)

      if (.not. allocated(fdobj%coefs)) then
         allocate(values(0, 0))
         info = 1
         return
      end if
      call eval_basis(x, fdobj%basis, nderiv, basismat, info)
      if (info /= 0) then
         allocate(values(0, 0))
         return
      end if
      allocate(values(size(x), size(fdobj%coefs, 2)))
      values = matmul(basismat, fdobj%coefs)
   end subroutine eval_fd

   pure subroutine mean_fd(fdobj, meanobj, info)
      type(fd_type), intent(in) :: fdobj !! Functional-data object whose replications are averaged coefficientwise.
      type(fd_type), intent(out) :: meanobj !! Single-replication functional-data object containing the mean function.
      integer, intent(out) :: info !! Zero on success; nonzero when the input has no allocated replications.
      real(dp), allocatable :: coefs(:, :)

      if (.not. allocated(fdobj%coefs) .or. size(fdobj%coefs, 2) < 1) then
         info = 1
         return
      end if
      allocate(coefs(size(fdobj%coefs, 1), 1))
      coefs(:, 1) = sum(fdobj%coefs, dim=2) / real(size(fdobj%coefs, 2), dp)
      call make_fd(coefs, fdobj%basis, meanobj, info)
   end subroutine mean_fd

   pure subroutine center_fd(fdobj, centered, meanobj, info)
      type(fd_type), intent(in) :: fdobj !! Functional-data object whose replication mean is to be removed.
      type(fd_type), intent(out) :: centered !! Functional data after subtracting the coefficientwise replication mean.
      type(fd_type), intent(out), optional :: meanobj !! Optional single-replication mean function removed from `fdobj`.
      integer, intent(out) :: info !! Zero on success; nonzero when the input has no allocated replications.
      type(fd_type) :: local_mean
      real(dp), allocatable :: coefs(:, :)

      call mean_fd(fdobj, local_mean, info)
      if (info /= 0) return
      allocate(coefs(size(fdobj%coefs, 1), size(fdobj%coefs, 2)))
      coefs = fdobj%coefs - spread(local_mean%coefs(:, 1), 2, size(fdobj%coefs, 2))
      call make_fd(coefs, fdobj%basis, centered, info)
      if (present(meanobj)) meanobj = local_mean
   end subroutine center_fd

   subroutine inprod_basis(basis1, basis2, nderiv1, nderiv2, result, info, nquad)
      type(basis_type), intent(in) :: basis1 !! First basis in the pairwise L2 inner products.
      type(basis_type), intent(in) :: basis2 !! Second basis in the pairwise L2 inner products.
      integer, intent(in) :: nderiv1 !! Derivative order applied to the first basis before integration.
      integer, intent(in) :: nderiv2 !! Derivative order applied to the second basis before integration.
      real(dp), allocatable, intent(out) :: result(:, :) !! Allocated matrix with shape `(basis1%nbasis, basis2%nbasis)`.
      integer, intent(out) :: info !! Zero on success; nonzero for disjoint ranges or evaluation/quadrature failure.
      integer, intent(in), optional :: nquad !! Odd Simpson point count per merged natural interval; defaults to nine.
      real(dp), allocatable :: b1(:), b2(:), breaks(:), points(:), values1(:, :), values2(:, :), weights(:)
      integer :: nq

      info = 0
      if (abs(basis1%rangeval(1) - basis2%rangeval(1)) > 64.0_dp * epsilon(1.0_dp) .or. &
          abs(basis1%rangeval(2) - basis2%rangeval(2)) > 64.0_dp * epsilon(1.0_dp)) then
         allocate(result(0, 0))
         info = 1
         return
      end if
      call basis_breaks(basis1, b1)
      call basis_breaks(basis2, b2)
      call merge_break_vectors(b1, b2, breaks)
      nq = 9
      if (present(nquad)) nq = nquad
      call quadset(breaks, nq, points, weights, info)
      if (info /= 0) then
         allocate(result(0, 0))
         return
      end if
      call eval_basis(points, basis1, nderiv1, values1, info)
      if (info /= 0) then
         allocate(result(0, 0))
         return
      end if
      call eval_basis(points, basis2, nderiv2, values2, info)
      if (info /= 0) then
         allocate(result(0, 0))
         return
      end if
      allocate(result(basis1%nbasis, basis2%nbasis))
      result = matmul(transpose(values1), values2 * spread(weights, 2, basis2%nbasis))
   end subroutine inprod_basis

   subroutine inprod_fd(fdobj1, fdobj2, nderiv1, nderiv2, result, info, nquad)
      type(fd_type), intent(in) :: fdobj1 !! First functional-data object in the replication-by-replication inner products.
      type(fd_type), intent(in) :: fdobj2 !! Second functional-data object in the replication-by-replication inner products.
      integer, intent(in) :: nderiv1 !! Derivative order applied to `fdobj1` before integration.
      integer, intent(in) :: nderiv2 !! Derivative order applied to `fdobj2` before integration.
      real(dp), allocatable, intent(out) :: result(:, :) !! Pairwise inner products: `fdobj1` reps by `fdobj2` reps.
      integer, intent(out) :: info !! Zero on success; otherwise an inner-product or object-dimension error code.
      integer, intent(in), optional :: nquad !! Odd Simpson point count per natural interval; defaults to nine.
      real(dp), allocatable :: gram(:, :)

      if (.not. allocated(fdobj1%coefs) .or. .not. allocated(fdobj2%coefs)) then
         allocate(result(0, 0))
         info = 1
         return
      end if
      call inprod_basis(fdobj1%basis, fdobj2%basis, nderiv1, nderiv2, gram, info, nquad)
      if (info /= 0) then
         allocate(result(0, 0))
         return
      end if
      allocate(result(size(fdobj1%coefs, 2), size(fdobj2%coefs, 2)))
      result = matmul(transpose(fdobj1%coefs), matmul(gram, fdobj2%coefs))
   end subroutine inprod_fd

   pure subroutine merge_break_vectors(a, b, merged)
      real(dp), intent(in) :: a(:) !! First sorted break-point vector.
      real(dp), intent(in) :: b(:) !! Second sorted break-point vector.
      real(dp), allocatable, intent(out) :: merged(:) !! Sorted union of both vectors with near-duplicates removed.
      real(dp), allocatable :: temp(:)
      real(dp) :: value
      integer :: i, j, k

      allocate(temp(size(a) + size(b)))
      i = 1
      j = 1
      k = 0
      do while (i <= size(a) .or. j <= size(b))
         if (j > size(b)) then
            value = a(i)
            i = i + 1
         else if (i > size(a)) then
            value = b(j)
            j = j + 1
         else if (a(i) <= b(j)) then
            value = a(i)
            i = i + 1
         else
            value = b(j)
            j = j + 1
         end if
         if (k == 0) then
            k = 1
            temp(k) = value
         else if (abs(value - temp(k)) > 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value))) then
            k = k + 1
            temp(k) = value
         end if
      end do
      allocate(merged(k))
      merged = temp(1:k)
   end subroutine merge_break_vectors

end module fda_fd
