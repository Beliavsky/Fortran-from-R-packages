program test_sparse_schur
   use rcsdp
   use rcsdp_block_ops, only : chol, inverse_from_chol
   implicit none
   type(csdp_problem) :: p
   type(csdp_block_matrix), allocatable :: adense(:)
   type(csdp_block_matrix) :: x, z, zchol, zi, rinv, at1, at2
   type(csdp_sparse_workspace) :: work
   real(dp), allocatable :: y(:), avec1(:), avec2(:), o1(:,:), o2(:,:)
   real(dp) :: scale
   integer :: info, i, spairs, dprods

   call read_sdpa_sparse('data/theta1.dat-s', p, info)
   if (info /= 0) error stop 'test_sparse_schur: read failed'
   call build_dense_constraints(p, adense)
   call build_sparse_workspace(p, work)
   call initsoln(p, x, y, z, adense)

   allocate(avec1(size(p%b)), avec2(size(p%b)))
   call op_a(p, x, avec1, adense)
   call op_a_sparse(p, x, avec2)
   if (maxval(abs(avec1-avec2)) > 1.0e-12_dp) error stop 'sparse A(X) mismatch'

   do i = 1, size(y)
      y(i) = sin(real(i,dp))*0.01_dp
   end do
   call op_at(p, y, at1, adense)
   call op_at_sparse(p, y, at2)
   if (block_maxdiff(at1,at2) > 1.0e-13_dp) error stop 'sparse A''(y) mismatch'

   zchol = z
   info = chol(zchol)
   if (info /= 0) error stop 'initial Z Cholesky failed'
   call inverse_from_chol(zchol, zi, rinv, info)
   if (info /= 0) error stop 'initial Z inverse failed'

   allocate(o1(size(p%b),size(p%b)), o2(size(p%b),size(p%b)))
   call op_o(p, zi, x, o1, adense)
   call op_o_sparse(p, zi, x, o2, work, spairs, dprods)
   scale = max(1.0_dp,maxval(abs(o1)))
   if (maxval(abs(o1-o2)) > 5.0e-13_dp*scale) error stop 'sparse Schur assembly mismatch'
   if (spairs <= 0) error stop 'sparse pair kernel was not exercised'
   if (work%constraint_nnz <= 0) error stop 'constraint nnz statistics missing'
   if (work%sparse_blocks <= 0) error stop 'sparse block classification missing'

   print *, 'test_sparse_schur: PASS'

contains

   real(dp) function block_maxdiff(a,b) result(v)
      type(csdp_block_matrix), intent(in) :: a,b
      integer :: k
      v = 0.0_dp
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            v = max(v,maxval(abs(a%block(k)%diag-b%block(k)%diag)))
         else
            v = max(v,maxval(abs(a%block(k)%mat-b%block(k)%mat)))
         end if
      end do
   end function block_maxdiff

end program test_sparse_schur
