! SPDX-License-Identifier: GPL-3.0-only
! Sparse quasi-definite KKT backend for SCS. The numerical LDL factorization
! is provided by the Apache-2.0 scs_qdldl module translated from QDLDL.
module scs_ldlt
   use scs_kinds, only : dp, i4
   use scs_types, only : scs_matrix
   use scs_qdldl, only : qdldl_etree, qdldl_factor, qdldl_solve
   implicit none
   private
   public :: scs_ldlt_factor

   type, public :: scs_ldlt_factor
      integer(i4) :: n = 0_i4
      integer(i4) :: n_primal = 0_i4
      integer(i4) :: n_dual = 0_i4
      integer(i4) :: kkt_nnz = 0_i4
      integer(i4) :: factor_nnz = 0_i4
      integer(i4) :: factorizations = 0_i4
      integer(i4) :: symbolic_analyses = 0_i4
      logical :: analyzed = .false.
      type(scs_matrix) :: kkt
      integer(i4), allocatable :: diag_idx(:), etree(:), lnz(:), lp(:), li(:)
      integer(i4), allocatable :: etree_work(:), y_idx(:), elim(:), next_space(:)
      real(dp), allocatable :: p_diag(:), lx(:), d(:), dinv(:), y_vals(:)
      logical, allocatable :: markers(:)
   contains
      procedure :: factorize => ldlt_factorize
      procedure :: solve => ldlt_solve
   end type scs_ldlt_factor
contains

   subroutine ldlt_factorize(self, a, p, has_p, diag_r, ok)
      class(scs_ldlt_factor), intent(inout) :: self
      type(scs_matrix), intent(in) :: a, p
      logical, intent(in) :: has_p
      real(dp), intent(in) :: diag_r(:)
      logical, intent(out) :: ok
      integer(i4) :: status

      ok = .false.
      if (.not. self%analyzed) then
         call build_kkt_upper(a, p, has_p, diag_r, self%kkt, self%diag_idx, self%p_diag)
         self%n_primal = a%n
         self%n_dual = a%m
         self%n = a%n + a%m
         self%kkt_nnz = int(size(self%kkt%x), i4)
         allocate(self%etree(self%n), self%lnz(self%n), self%etree_work(self%n))
         status = qdldl_etree(self%n, self%kkt%p, self%kkt%i, self%etree_work, self%lnz, self%etree)
         if (status < 0_i4) return
         self%factor_nnz = status
         allocate(self%lp(self%n+1), self%li(self%factor_nnz), self%lx(self%factor_nnz))
         allocate(self%d(self%n), self%dinv(self%n), self%markers(self%n), self%y_idx(self%n))
         allocate(self%elim(self%n), self%next_space(self%n), self%y_vals(self%n))
         self%analyzed = .true.
         self%symbolic_analyses = self%symbolic_analyses + 1_i4
      else
         if (a%n /= self%n_primal .or. a%m /= self%n_dual) return
         call update_kkt_diagonal(self, diag_r)
      end if

      status = qdldl_factor(self%n, self%kkt%p, self%kkt%i, self%kkt%x, &
                            self%lp, self%li, self%lx, self%d, self%dinv, self%lnz, self%etree, &
                            self%markers, self%y_idx, self%elim, self%next_space, self%y_vals)
      if (status < self%n_primal) return
      self%factorizations = self%factorizations + 1_i4
      ok = .true.
   end subroutine ldlt_factorize

   subroutine ldlt_solve(self, b)
      class(scs_ldlt_factor), intent(in) :: self
      real(dp), intent(inout) :: b(:)
      if (size(b) < self%n) error stop 'scs_ldlt: right-hand side too short'
      call qdldl_solve(self%n, self%lp, self%li, self%lx, self%dinv, b)
   end subroutine ldlt_solve

   subroutine update_kkt_diagonal(self, diag_r)
      class(scs_ldlt_factor), intent(inout) :: self
      real(dp), intent(in) :: diag_r(:)
      integer(i4) :: j
      do j = 1_i4, self%n_primal
         self%kkt%x(self%diag_idx(j)) = self%p_diag(j) + diag_r(j)
      end do
      do j = 1_i4, self%n_dual
         self%kkt%x(self%diag_idx(self%n_primal+j)) = -diag_r(self%n_primal+j)
      end do
   end subroutine update_kkt_diagonal

   subroutine build_kkt_upper(a, p, has_p, diag_r, kkt, diag_idx, p_diag)
      type(scs_matrix), intent(in) :: a, p
      logical, intent(in) :: has_p
      real(dp), intent(in) :: diag_r(:)
      type(scs_matrix), intent(out) :: kkt
      integer(i4), allocatable, intent(out) :: diag_idx(:)
      real(dp), allocatable, intent(out) :: p_diag(:)
      integer(i4), allocatable :: counts(:), next(:)
      integer(i4) :: n, m, nm, c, r, q, pos, nnz, pk

      n = a%n
      m = a%m
      nm = n + m
      allocate(counts(nm), source=1_i4)
      allocate(p_diag(n), source=0.0_dp)

      if (has_p) then
         do c = 1_i4, n
            do pk = p%p(c), p%p(c+1_i4)-1_i4
               r = p%i(pk)
               if (r > c) cycle
               if (r == c) then
                  p_diag(c) = p_diag(c) + p%x(pk)
               else
                  counts(c) = counts(c) + 1_i4
               end if
            end do
         end do
      end if
      do c = 1_i4, n
         do q = a%p(c), a%p(c+1_i4)-1_i4
            r = a%i(q)
            counts(n+r) = counts(n+r) + 1_i4
         end do
      end do

      nnz = sum(counts)
      kkt%m = nm
      kkt%n = nm
      allocate(kkt%x(nnz), kkt%i(nnz), kkt%p(nm+1), diag_idx(nm), next(nm))
      kkt%p(1) = 1_i4
      do c = 1_i4, nm
         kkt%p(c+1_i4) = kkt%p(c) + counts(c)
      end do
      next = kkt%p(1:nm)

      do c = 1_i4, n
         if (has_p) then
            do pk = p%p(c), p%p(c+1_i4)-1_i4
               r = p%i(pk)
               if (r >= c) cycle
               pos = next(c)
               kkt%i(pos) = r
               kkt%x(pos) = p%x(pk)
               next(c) = pos + 1_i4
            end do
         end if
         pos = next(c)
         kkt%i(pos) = c
         kkt%x(pos) = p_diag(c) + diag_r(c)
         diag_idx(c) = pos
         next(c) = pos + 1_i4
      end do

      do c = 1_i4, n
         do q = a%p(c), a%p(c+1_i4)-1_i4
            r = a%i(q)
            pos = next(n+r)
            kkt%i(pos) = c
            kkt%x(pos) = a%x(q)
            next(n+r) = pos + 1_i4
         end do
      end do
      do r = 1_i4, m
         pos = next(n+r)
         kkt%i(pos) = n+r
         kkt%x(pos) = -diag_r(n+r)
         diag_idx(n+r) = pos
         next(n+r) = pos + 1_i4
      end do
   end subroutine build_kkt_upper

end module scs_ldlt
