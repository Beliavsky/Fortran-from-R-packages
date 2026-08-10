program test_fill_products
   use rcsdp
   use rcsdp_block_ops, only : mat_mult, chol, inverse_from_chol
   implicit none
   type(csdp_problem) :: p
   type(csdp_fill_workspace) :: fill
   type(csdp_block_matrix) :: x, z, zchol, zi, rinv, c1, c2
   real(dp), allocatable :: y(:), v1(:), v2(:)
   real(dp) :: err, scale
   integer :: info, ib, e, ns, nd

   call read_sdpa_sparse('data/theta1.dat-s',p,info)
   if (info /= 0) error stop 'test_fill_products: read failed'
   call build_fill_workspace(p,fill)
   call initsoln(p,x,y,z)
   if (fill%fill_nnz <= 0) error stop 'fill pattern was empty'

   zchol=z
   info=chol(zchol)
   if (info /= 0) error stop 'Z Cholesky failed'
   call inverse_from_chol(zchol,zi,rinv,info)
   if (info /= 0) error stop 'Z inverse failed'

   ! spB: C is itself supported on makefill by construction.
   c1=p%c; c2=p%c
   call mat_mult(1.0_dp,0.0_dp,zi,p%c,c1)
   ns=0; nd=0
   call mat_multspb(1.0_dp,0.0_dp,zi,p%c,c2,fill,1.10_dp,ns,nd)
   err=block_maxdiff(c1,c2)
   if (err > 1.0e-11_dp) error stop 'mat_multspb mismatch'

   ! spA: C is supported on fill; compare the complete product.
   call mat_mult(1.0_dp,0.0_dp,p%c,x,c1)
   call mat_multspa(1.0_dp,0.0_dp,p%c,x,c2,fill,1.10_dp,ns,nd)
   err=block_maxdiff(c1,c2)
   if (err > 1.0e-11_dp) error stop 'mat_multspa mismatch'

   ! spC only promises entries in fill.
   call mat_mult(1.0_dp,0.0_dp,zi,x,c1)
   call mat_multspc(1.0_dp,0.0_dp,zi,x,c2,fill,1.10_dp,ns,nd)
   err=0.0_dp
   scale=1.0_dp
   do ib=1,size(fill%block)
      if (p%c%block(ib)%category == csdp_diag) then
         err=max(err,maxval(abs(c1%block(ib)%diag-c2%block(ib)%diag)))
         scale=max(scale,maxval(abs(c1%block(ib)%diag)))
      else
         do e=1,fill%block(ib)%nnz()
            err=max(err,abs(c1%block(ib)%mat(fill%block(ib)%i(e),fill%block(ib)%j(e))- &
               c2%block(ib)%mat(fill%block(ib)%i(e),fill%block(ib)%j(e))))
         end do
         scale=max(scale,maxval(abs(c1%block(ib)%mat)))
      end if
   end do
   if (err > 1.0e-12_dp*scale) error stop 'mat_multspc mismatch on fill'

   allocate(v1(size(p%b)),v2(size(p%b)))
   y=[(sin(real(e,dp))*0.01_dp,e=1,size(y))]
   call apply_o_sparse(p,zi,x,y,v1)
   call apply_o_fill(p,zi,x,y,v2,fill,1.10_dp,ns,nd)
   if (maxval(abs(v1-v2)) > 1.0e-10_dp*max(1.0_dp,maxval(abs(v1)))) &
      error stop 'apply_o_fill mismatch'

   if (ns <= 0) error stop 'sparse fill products were not exercised'
   print *, 'test_fill_products: PASS'

contains

   real(dp) function block_maxdiff(a,b) result(v)
      type(csdp_block_matrix), intent(in) :: a,b
      integer :: k
      v=0.0_dp
      do k=1,size(a%block)
         if (a%block(k)%category == csdp_diag) then
            v=max(v,maxval(abs(a%block(k)%diag-b%block(k)%diag)))
         else
            v=max(v,maxval(abs(a%block(k)%mat-b%block(k)%mat)))
         end if
      end do
   end function block_maxdiff
end program test_fill_products
