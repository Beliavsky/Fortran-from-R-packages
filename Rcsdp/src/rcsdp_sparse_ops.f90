! Sparse constraint and Schur-complement kernels translated from CSDP 6.1.1.
! Original software is distributed under the Common Public License 1.0.
module rcsdp_sparse_ops
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_problem, csdp_block_matrix, csdp_sparse_block, csdp_matrix, csdp_diag
   use rcsdp_block_ops, only : zero_mat, mat_mult
   implicit none
   private

   type, public :: csdp_block_index
      integer, allocatable :: icon(:)
      integer, allocatable :: pos(:)
      logical, allocatable :: is_sparse(:)
   end type csdp_block_index

   type, public :: csdp_sparse_workspace
      type(csdp_block_index), allocatable :: byblock(:)
      integer :: constraint_nnz = 0
      integer :: sparse_blocks = 0
      integer :: dense_blocks = 0
   end type csdp_sparse_workspace

   public :: build_sparse_workspace, op_a_sparse, op_at_sparse
   public :: op_o_sparse, apply_o_sparse

contains

   subroutine build_sparse_workspace(prob, work)
      type(csdp_problem), intent(in) :: prob
      type(csdp_sparse_workspace), intent(out) :: work
      integer, allocatable :: count(:), next(:)
      integer :: i, p, ib, q, nb, ne, n

      nb = size(prob%c%block)
      allocate(work%byblock(nb), count(nb), next(nb))
      count = 0
      work%constraint_nnz = 0
      work%sparse_blocks = 0
      work%dense_blocks = 0

      do i = 1, size(prob%a)
         do p = 1, size(prob%a(i)%block)
            ib = prob%a(i)%block(p)%blocknum
            count(ib) = count(ib) + 1
            work%constraint_nnz = work%constraint_nnz + prob%a(i)%block(p)%nnz()
         end do
      end do

      do ib = 1, nb
         allocate(work%byblock(ib)%icon(count(ib)), work%byblock(ib)%pos(count(ib)), &
            work%byblock(ib)%is_sparse(count(ib)))
      end do

      next = 0
      do i = 1, size(prob%a)
         do p = 1, size(prob%a(i)%block)
            ib = prob%a(i)%block(p)%blocknum
            next(ib) = next(ib) + 1
            q = next(ib)
            work%byblock(ib)%icon(q) = i
            work%byblock(ib)%pos(q) = p
            ne = prob%a(i)%block(p)%nnz()
            n = prob%a(i)%block(p)%n
            ! CSDP's original classification from easysdp.c/readprob.c:
            ! dense iff numentries > 0.25*blocksize and numentries > 15.
            work%byblock(ib)%is_sparse(q) = prob%c%block(ib)%category == csdp_diag .or. &
               .not. (real(ne,dp) > 0.25_dp*real(n,dp) .and. ne > 15)
            if (work%byblock(ib)%is_sparse(q)) then
               work%sparse_blocks = work%sparse_blocks + 1
            else
               work%dense_blocks = work%dense_blocks + 1
            end if
         end do
      end do
   end subroutine build_sparse_workspace

   subroutine op_a_sparse(prob, x, result)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), intent(in) :: x
      real(dp), intent(out) :: result(:)
      integer :: icon, p, e, ib, ii, jj
      real(dp) :: ent, s

      if (size(result) /= size(prob%a)) error stop 'op_a_sparse: result has wrong size'
      result = 0.0_dp
      do icon = 1, size(prob%a)
         s = 0.0_dp
         do p = 1, size(prob%a(icon)%block)
            ib = prob%a(icon)%block(p)%blocknum
            if (prob%c%block(ib)%category == csdp_diag) then
               do e = 1, prob%a(icon)%block(p)%nnz()
                  ii = prob%a(icon)%block(p)%i(e)
                  s = s + prob%a(icon)%block(p)%v(e)*x%block(ib)%diag(ii)
               end do
            else
               do e = 1, prob%a(icon)%block(p)%nnz()
                  ii = prob%a(icon)%block(p)%i(e)
                  jj = prob%a(icon)%block(p)%j(e)
                  ent = prob%a(icon)%block(p)%v(e)
                  if (ii == jj) then
                     s = s + ent*x%block(ib)%mat(ii,jj)
                  else
                     s = s + ent*(x%block(ib)%mat(ii,jj) + x%block(ib)%mat(jj,ii))
                  end if
               end do
            end if
         end do
         result(icon) = s
      end do
   end subroutine op_a_sparse

   subroutine op_at_sparse(prob, y, result)
      type(csdp_problem), intent(in) :: prob
      real(dp), intent(in) :: y(:)
      type(csdp_block_matrix), intent(out) :: result
      integer :: icon, p, e, ib, ii, jj
      real(dp) :: val

      if (size(y) /= size(prob%a)) error stop 'op_at_sparse: y has wrong size'
      result = prob%c
      call zero_mat(result)
      do icon = 1, size(prob%a)
         if (abs(y(icon)) <= tiny(1.0_dp)) cycle
         do p = 1, size(prob%a(icon)%block)
            ib = prob%a(icon)%block(p)%blocknum
            if (prob%c%block(ib)%category == csdp_diag) then
               do e = 1, prob%a(icon)%block(p)%nnz()
                  ii = prob%a(icon)%block(p)%i(e)
                  result%block(ib)%diag(ii) = result%block(ib)%diag(ii) + &
                     y(icon)*prob%a(icon)%block(p)%v(e)
               end do
            else
               do e = 1, prob%a(icon)%block(p)%nnz()
                  ii = prob%a(icon)%block(p)%i(e)
                  jj = prob%a(icon)%block(p)%j(e)
                  val = y(icon)*prob%a(icon)%block(p)%v(e)
                  result%block(ib)%mat(ii,jj) = result%block(ib)%mat(ii,jj) + val
                  if (ii /= jj) result%block(ib)%mat(jj,ii) = result%block(ib)%mat(jj,ii) + val
               end do
            end if
         end do
      end do
   end subroutine op_at_sparse

   subroutine op_o_sparse(prob, zi, x, o, work, sparse_pairs, dense_products)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), intent(in) :: zi, x
      real(dp), intent(out) :: o(:,:)
      type(csdp_sparse_workspace), intent(in) :: work
      integer, intent(out), optional :: sparse_pairs, dense_products
      real(dp), allocatable :: ai(:,:), tmp(:,:), transformed(:,:)
      real(dp) :: contrib
      integer :: ib, a, b, i, j, nr, n, npair, ndense

      if (size(o,1) /= size(prob%a) .or. size(o,2) /= size(prob%a)) error stop 'op_o_sparse: wrong shape'
      o = 0.0_dp
      npair = 0
      ndense = 0

      do ib = 1, size(work%byblock)
         nr = size(work%byblock(ib)%icon)
         if (nr == 0) cycle

         if (prob%c%block(ib)%category == csdp_diag) then
            do a = 1, nr
               i = work%byblock(ib)%icon(a)
               do b = a, nr
                  j = work%byblock(ib)%icon(b)
                  contrib = diag_pair(prob%a(i)%block(work%byblock(ib)%pos(a)), &
                     prob%a(j)%block(work%byblock(ib)%pos(b)), zi%block(ib)%diag, x%block(ib)%diag)
                  o(i,j) = o(i,j) + contrib
                  npair = npair + 1
               end do
            end do
         else
            n = prob%c%block(ib)%n
            allocate(ai(n,n), tmp(n,n), transformed(n,n))
            do a = 1, nr
               i = work%byblock(ib)%icon(a)
               if (.not. work%byblock(ib)%is_sparse(a)) then
                  call sparse_block_to_dense(prob%a(i)%block(work%byblock(ib)%pos(a)), ai)
                  tmp = matmul(zi%block(ib)%mat, ai)
                  transformed = matmul(tmp, x%block(ib)%mat)
                  ndense = ndense + 2
                  do b = a, nr
                     j = work%byblock(ib)%icon(b)
                        o(i,j) = o(i,j) + sparse_trace(prob%a(j)%block(work%byblock(ib)%pos(b)), transformed)
                  end do
               else
                  do b = a, nr
                     j = work%byblock(ib)%icon(b)
                        o(i,j) = o(i,j) + matrix_pair(prob%a(i)%block(work%byblock(ib)%pos(a)), &
                        prob%a(j)%block(work%byblock(ib)%pos(b)), zi%block(ib)%mat, x%block(ib)%mat)
                     npair = npair + 1
                  end do
               end if
            end do
            deallocate(ai, tmp, transformed)
         end if
      end do

      ! The CSDP Schur matrix is used as a symmetric positive-definite system.
      ! The block cross-index is ordered by constraint number, so contributions
      ! above were accumulated in the upper triangle only.
      do j = 1, size(o,2)
         do i = j + 1, size(o,1)
            o(i,j) = o(j,i)
         end do
      end do

      if (present(sparse_pairs)) sparse_pairs = npair
      if (present(dense_products)) dense_products = ndense
   end subroutine op_o_sparse

   subroutine apply_o_sparse(prob, zi, x, y, result)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), intent(in) :: zi, x
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: result(:)
      type(csdp_block_matrix) :: aty, tmp, prod

      if (size(y) /= size(prob%a) .or. size(result) /= size(prob%a)) error stop 'apply_o_sparse: wrong size'
      call op_at_sparse(prob, y, aty)
      tmp = prob%c
      prod = prob%c
      call mat_mult(1.0_dp, 0.0_dp, zi, aty, tmp)
      call mat_mult(1.0_dp, 0.0_dp, tmp, x, prod)
      call op_a_sparse(prob, prod, result)
   end subroutine apply_o_sparse

   real(dp) function diag_pair(a, b, zi, x) result(s)
      type(csdp_sparse_block), intent(in) :: a, b
      real(dp), intent(in) :: zi(:), x(:)
      integer :: p, q
      s = 0.0_dp
      do p = 1, a%nnz()
         do q = 1, b%nnz()
            if (a%i(p) == b%i(q)) s = s + a%v(p)*b%v(q)*zi(a%i(p))*x(a%i(p))
         end do
      end do
   end function diag_pair

   real(dp) function matrix_pair(a, b, zi, x) result(s)
      type(csdp_sparse_block), intent(in) :: a, b
      real(dp), intent(in) :: zi(:,:), x(:,:)
      integer :: ia, ib, p, q, r, t
      real(dp) :: ea, eb
      s = 0.0_dp
      do ia = 1, a%nnz()
         ea = a%v(ia)
         p = a%i(ia)
         q = a%j(ia)
         if (p == q) then
            do ib = 1, b%nnz()
               eb = b%v(ib)
               r = b%i(ib)
               t = b%j(ib)
               if (r == t) then
                  s = s + ea*eb*zi(r,q)*x(t,p)
               else
                  s = s + ea*eb*(zi(r,q)*x(t,p) + zi(t,q)*x(r,p))
               end if
            end do
         else
            do ib = 1, b%nnz()
               eb = b%v(ib)
               r = b%i(ib)
               t = b%j(ib)
               if (r == t) then
                  s = s + ea*eb*(zi(r,q)*x(t,p) + zi(r,p)*x(t,q))
               else
                  s = s + ea*eb*(zi(r,q)*x(t,p) + zi(r,p)*x(t,q) + &
                     zi(t,q)*x(r,p) + zi(t,p)*x(r,q))
               end if
            end do
         end if
      end do
   end function matrix_pair

   real(dp) function sparse_trace(a, w) result(s)
      type(csdp_sparse_block), intent(in) :: a
      real(dp), intent(in) :: w(:,:)
      integer :: e, i, j
      s = 0.0_dp
      do e = 1, a%nnz()
         i = a%i(e)
         j = a%j(e)
         if (i == j) then
            s = s + a%v(e)*w(i,j)
         else
            s = s + a%v(e)*(w(i,j) + w(j,i))
         end if
      end do
   end function sparse_trace

   subroutine sparse_block_to_dense(a, d)
      type(csdp_sparse_block), intent(in) :: a
      real(dp), intent(out) :: d(:,:)
      integer :: e, i, j
      d = 0.0_dp
      do e = 1, a%nnz()
         i = a%i(e)
         j = a%j(e)
         d(i,j) = d(i,j) + a%v(e)
         if (i /= j) d(j,i) = d(j,i) + a%v(e)
      end do
   end subroutine sparse_block_to_dense

end module rcsdp_sparse_ops
