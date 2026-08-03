! This file is part of multiAssetOptions-fortran.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module mao_sparse
   use mao_kinds, only: dp
   use mao_status, only: status_type, clear_status, set_status, &
      mao_invalid_argument, mao_allocation_error, mao_solver_failure
   implicit none
   private

   type, public :: csr_matrix
      integer :: nrow = 0
      integer :: ncol = 0
      integer, allocatable :: row_ptr(:)
      integer, allocatable :: col_ind(:)
      real(dp), allocatable :: value(:)
   end type csr_matrix

   public :: csr_matvec, csr_to_dense, csr_shifted_identity
   public :: bicgstab_solve, csr_diagonal

contains

   subroutine csr_matvec(a, x, y)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(:)
      integer :: i, k

      y = 0.0_dp
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i+1) - 1
            y(i) = y(i) + a%value(k) * x(a%col_ind(k))
         end do
      end do
   end subroutine csr_matvec

   subroutine csr_diagonal(a, diagonal, status)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable, intent(out) :: diagonal(:)
      type(status_type), intent(out) :: status
      integer :: i, k, stat
      logical :: found

      call clear_status(status)
      if (a%nrow /= a%ncol) then
         call set_status(status, mao_invalid_argument, &
            'diagonal requires a square matrix')
         return
      end if
      allocate(diagonal(a%nrow), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate diagonal')
         return
      end if
      diagonal = 0.0_dp
      do i = 1, a%nrow
         found = .false.
         do k = a%row_ptr(i), a%row_ptr(i+1) - 1
            if (a%col_ind(k) == i) then
               diagonal(i) = a%value(k)
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            call set_status(status, mao_invalid_argument, &
               'CSR matrix is missing an explicit diagonal entry')
            return
         end if
      end do
   end subroutine csr_diagonal

   subroutine csr_to_dense(a, dense, status)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable, intent(out) :: dense(:,:)
      type(status_type), intent(out) :: status
      integer :: i, k, stat

      call clear_status(status)
      allocate(dense(a%nrow,a%ncol), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate dense matrix')
         return
      end if
      dense = 0.0_dp
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i+1) - 1
            dense(i,a%col_ind(k)) = dense(i,a%col_ind(k)) + a%value(k)
         end do
      end do
   end subroutine csr_to_dense

   subroutine csr_shifted_identity(operator, scale, extra_diagonal, system, status)
      type(csr_matrix), intent(in) :: operator
      real(dp), intent(in) :: scale
      real(dp), intent(in) :: extra_diagonal(:)
      type(csr_matrix), intent(out) :: system
      type(status_type), intent(out) :: status
      integer :: i, k, stat
      logical :: found

      call clear_status(status)
      if (operator%nrow /= operator%ncol .or. &
          size(extra_diagonal) /= operator%nrow) then
         call set_status(status, mao_invalid_argument, &
            'invalid shifted-identity dimensions')
         return
      end if

      system%nrow = operator%nrow
      system%ncol = operator%ncol
      allocate(system%row_ptr(size(operator%row_ptr)), &
         system%col_ind(size(operator%col_ind)), &
         system%value(size(operator%value)), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate system matrix')
         return
      end if
      system%row_ptr = operator%row_ptr
      system%col_ind = operator%col_ind
      system%value = scale * operator%value

      do i = 1, system%nrow
         found = .false.
         do k = system%row_ptr(i), system%row_ptr(i+1) - 1
            if (system%col_ind(k) == i) then
               system%value(k) = system%value(k) + 1.0_dp + extra_diagonal(i)
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            call set_status(status, mao_invalid_argument, &
               'operator is missing an explicit diagonal entry')
            return
         end if
      end do
   end subroutine csr_shifted_identity

   subroutine bicgstab_solve(a, b, x, tolerance, max_iterations, &
      iterations, residual_norm, status)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: b(:)
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: tolerance
      integer, intent(in) :: max_iterations
      integer, intent(out) :: iterations
      real(dp), intent(out) :: residual_norm
      type(status_type), intent(out) :: status

      real(dp), allocatable :: r(:), rhat(:), p(:), v(:), s(:), t(:)
      real(dp), allocatable :: phat(:), shat(:), ax(:), diag(:)
      real(dp) :: alpha, beta, omega, rho, rho_old, denom
      real(dp) :: norm_b, tiny_value
      integer :: i, stat
      type(status_type) :: local_status

      call clear_status(status)
      iterations = 0
      residual_norm = huge(1.0_dp)
      if (a%nrow /= a%ncol .or. size(b) /= a%nrow .or. &
          size(x) /= a%ncol .or. tolerance <= 0.0_dp .or. &
          max_iterations < 1) then
         call set_status(status, mao_invalid_argument, &
            'invalid BiCGSTAB dimensions or controls')
         return
      end if

      allocate(r(a%nrow), rhat(a%nrow), p(a%nrow), v(a%nrow), &
         s(a%nrow), t(a%nrow), phat(a%nrow), shat(a%nrow), &
         ax(a%nrow), stat=stat)
      if (stat /= 0) then
         call set_status(status, mao_allocation_error, &
            'unable to allocate BiCGSTAB workspace')
         return
      end if
      call csr_diagonal(a, diag, local_status)
      if (.not. local_status%ok()) then
         status = local_status
         return
      end if

      tiny_value = 100.0_dp * tiny(1.0_dp)
      do i = 1, size(diag)
         if (abs(diag(i)) <= tiny_value) diag(i) = 1.0_dp
      end do

      call csr_matvec(a, x, ax)
      r = b - ax
      rhat = r
      p = 0.0_dp
      v = 0.0_dp
      alpha = 1.0_dp
      omega = 1.0_dp
      rho_old = 1.0_dp
      norm_b = sqrt(max(dot_product(b/diag,b/diag), tiny_value))
      residual_norm = sqrt(dot_product(r/diag,r/diag)) / norm_b
      if (residual_norm <= tolerance) return

      do iterations = 1, max_iterations
         rho = dot_product(rhat,r)
         if (abs(rho) <= tiny_value) then
            call set_status(status, mao_solver_failure, &
               'BiCGSTAB breakdown: rho is zero')
            return
         end if
         if (iterations == 1) then
            p = r
         else
            if (abs(omega) <= tiny_value) then
               call set_status(status, mao_solver_failure, &
                  'BiCGSTAB breakdown: omega is zero')
               return
            end if
            beta = (rho / rho_old) * (alpha / omega)
            p = r + beta * (p - omega * v)
         end if

         phat = p / diag
         call csr_matvec(a, phat, v)
         denom = dot_product(rhat,v)
         if (abs(denom) <= tiny_value) then
            call set_status(status, mao_solver_failure, &
               'BiCGSTAB breakdown: alpha denominator is zero')
            return
         end if
         alpha = rho / denom
         s = r - alpha * v
         residual_norm = sqrt(dot_product(s/diag,s/diag)) / norm_b
         if (residual_norm <= tolerance) then
            x = x + alpha * phat
            return
         end if

         shat = s / diag
         call csr_matvec(a, shat, t)
         denom = dot_product(t,t)
         if (denom <= tiny_value) then
            call set_status(status, mao_solver_failure, &
               'BiCGSTAB breakdown: omega denominator is zero')
            return
         end if
         omega = dot_product(t,s) / denom
         x = x + alpha * phat + omega * shat
         r = s - omega * t
         residual_norm = sqrt(dot_product(r/diag,r/diag)) / norm_b
         if (residual_norm <= tolerance) return
         rho_old = rho
      end do

      iterations = max_iterations
      call set_status(status, mao_solver_failure, &
         'BiCGSTAB did not converge within max_iterations')
   end subroutine bicgstab_solve

end module mao_sparse
