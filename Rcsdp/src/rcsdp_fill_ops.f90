! Fill-restricted block products translated from CSDP 6.1.1 makefill.c and
! mat_multsp.c.  See LICENSE (CPL-1.0).
module rcsdp_fill_ops
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_problem, csdp_block_matrix, csdp_matrix, csdp_diag
   use rcsdp_block_ops, only : zero_mat
   use rcsdp_sparse_ops, only : op_at_sparse, op_a_sparse
   implicit none
   private

   type, public :: csdp_fill_block
      integer :: n = 0
      integer, allocatable :: i(:)
      integer, allocatable :: j(:)
      integer, allocatable :: row_ptr(:)
      integer, allocatable :: col_ptr(:)
      integer, allocatable :: col_entry(:)
   contains
      procedure :: nnz => fill_block_nnz
      procedure :: density => fill_block_density
   end type csdp_fill_block

   type, public :: csdp_fill_workspace
      type(csdp_fill_block), allocatable :: block(:)
      integer :: fill_nnz = 0
      integer :: full_entries = 0
   end type csdp_fill_workspace

   public :: build_fill_workspace
   public :: mat_multspa, mat_multspb, mat_multspc
   public :: apply_o_fill

contains

   pure integer function fill_block_nnz(this) result(n)
      class(csdp_fill_block), intent(in) :: this
      if (allocated(this%i)) then
         n = size(this%i)
      else
         n = 0
      end if
   end function fill_block_nnz

   pure real(dp) function fill_block_density(this) result(v)
      class(csdp_fill_block), intent(in) :: this
      if (this%n <= 0) then
         v = 0.0_dp
      else
         v = real(this%nnz(),dp)/real(this%n*this%n,dp)
      end if
   end function fill_block_density

   subroutine build_fill_workspace(prob, fill)
      type(csdp_problem), intent(in) :: prob
      type(csdp_fill_workspace), intent(out) :: fill
      logical, allocatable :: mask(:,:)
      integer :: ib, n, i, j, icon, p, e, k, cnt

      allocate(fill%block(size(prob%c%block)))
      fill%fill_nnz = 0
      fill%full_entries = 0

      do ib = 1, size(prob%c%block)
         n = prob%c%block(ib)%n
         fill%block(ib)%n = n
         if (prob%c%block(ib)%category == csdp_diag) then
            fill%full_entries = fill%full_entries + n
            allocate(fill%block(ib)%i(n), fill%block(ib)%j(n), fill%block(ib)%row_ptr(n+1), &
               fill%block(ib)%col_ptr(n+1), fill%block(ib)%col_entry(n))
            do i = 1, n
               fill%block(ib)%i(i) = i
               fill%block(ib)%j(i) = i
               fill%block(ib)%row_ptr(i) = i
               fill%block(ib)%col_ptr(i) = i
               fill%block(ib)%col_entry(i) = i
            end do
            fill%block(ib)%row_ptr(n+1) = n+1
            fill%block(ib)%col_ptr(n+1) = n+1
            fill%fill_nnz = fill%fill_nnz + n
            cycle
         end if
         fill%full_entries = fill%full_entries + n*n

         allocate(mask(n,n))
         mask = .false.
         do i = 1, n
            mask(i,i) = .true.
         end do

         ! CSDP makefill includes every nonzero of C as well as the diagonal.
         do j = 1, n
            do i = 1, n
               if (abs(prob%c%block(ib)%mat(i,j)) > tiny(1.0_dp)) mask(i,j) = .true.
            end do
         end do

         ! Add all symmetric positions touched by all constraint blocks.
         do icon = 1, size(prob%a)
            do p = 1, size(prob%a(icon)%block)
               if (prob%a(icon)%block(p)%blocknum /= ib) cycle
               do e = 1, prob%a(icon)%block(p)%nnz()
                  i = prob%a(icon)%block(p)%i(e)
                  j = prob%a(icon)%block(p)%j(e)
                  mask(i,j) = .true.
                  mask(j,i) = .true.
               end do
            end do
         end do

         cnt = count(mask)
         allocate(fill%block(ib)%i(cnt), fill%block(ib)%j(cnt))
         k = 0
         ! Match the original makefill row-major traversal; ordering is not
         ! mathematically significant but deterministic ordering helps tests.
         do i = 1, n
            do j = 1, n
               if (.not. mask(i,j)) cycle
               k = k + 1
               fill%block(ib)%i(k) = i
               fill%block(ib)%j(k) = j
            end do
         end do
         call build_fill_indices(fill%block(ib))
         fill%fill_nnz = fill%fill_nnz + cnt
         deallocate(mask)
      end do
   end subroutine build_fill_workspace

   subroutine build_fill_indices(fb)
      type(csdp_fill_block), intent(inout) :: fb
      integer, allocatable :: count_col(:), next_col(:)
      integer :: e, p, q, n

      n=fb%n
      allocate(fb%row_ptr(n+1),fb%col_ptr(n+1),fb%col_entry(fb%nnz()))
      fb%row_ptr=fb%nnz()+1
      do e=1,fb%nnz()
         p=fb%i(e)
         if (fb%row_ptr(p)==fb%nnz()+1) fb%row_ptr(p)=e
      end do
      fb%row_ptr(n+1)=fb%nnz()+1
      do p=n,1,-1
         if (fb%row_ptr(p)==fb%nnz()+1) fb%row_ptr(p)=fb%row_ptr(p+1)
      end do

      allocate(count_col(n),next_col(n))
      count_col=0
      do e=1,fb%nnz()
         count_col(fb%j(e))=count_col(fb%j(e))+1
      end do
      fb%col_ptr(1)=1
      do q=1,n
         fb%col_ptr(q+1)=fb%col_ptr(q)+count_col(q)
      end do
      next_col=fb%col_ptr(1:n)
      do e=1,fb%nnz()
         q=fb%j(e)
         fb%col_entry(next_col(q))=e
         next_col(q)=next_col(q)+1
      end do
   end subroutine build_fill_indices

   subroutine scale_output(scale2, c)
      real(dp), intent(in) :: scale2
      type(csdp_block_matrix), intent(inout) :: c
      integer :: ib
      if (abs(scale2) <= tiny(1.0_dp)) then
         call zero_mat(c)
      else if (abs(scale2-1.0_dp) > tiny(1.0_dp)) then
         do ib = 1, size(c%block)
            if (c%block(ib)%category == csdp_diag) then
               c%block(ib)%diag = scale2*c%block(ib)%diag
            else
               c%block(ib)%mat = scale2*c%block(ib)%mat
            end if
         end do
      end if
   end subroutine scale_output

   ! C = scale1*A*B + scale2*C, exploiting that B is nonzero only on fill.
   subroutine mat_multspb(scale1, scale2, a, b, c, fill, density_limit, sparse_count, dense_count)
      real(dp), intent(in) :: scale1, scale2
      type(csdp_block_matrix), intent(in) :: a, b
      type(csdp_block_matrix), intent(inout) :: c
      type(csdp_fill_workspace), intent(in) :: fill
      real(dp), intent(in), optional :: density_limit
      integer, intent(inout), optional :: sparse_count, dense_count
      real(dp) :: lim, temp
      integer :: ib, e, p, q

      lim = 0.01_dp
      if (present(density_limit)) lim = density_limit
      call scale_output(scale2,c)
      if (abs(scale1) <= tiny(1.0_dp)) return

      do ib = 1, size(a%block)
         if (a%block(ib)%category == csdp_diag) then
            c%block(ib)%diag = c%block(ib)%diag + scale1*a%block(ib)%diag*b%block(ib)%diag
         else if (fill%block(ib)%density() > lim) then
            c%block(ib)%mat = c%block(ib)%mat + scale1*matmul(a%block(ib)%mat,b%block(ib)%mat)
            if (present(dense_count)) dense_count = dense_count + 1
         else
            !$omp parallel do default(shared) private(q,e,p,temp) schedule(static)
            do q = 1, fill%block(ib)%n
               do e = fill%block(ib)%col_ptr(q), fill%block(ib)%col_ptr(q+1)-1
                  p = fill%block(ib)%i(fill%block(ib)%col_entry(e))
                  temp = scale1*b%block(ib)%mat(p,q)
                  c%block(ib)%mat(:,q) = c%block(ib)%mat(:,q) + temp*a%block(ib)%mat(:,p)
               end do
            end do
            !$omp end parallel do
            if (present(sparse_count)) sparse_count = sparse_count + 1
         end if
      end do
   end subroutine mat_multspb

   ! C = scale1*A*B + scale2*C, exploiting that A is nonzero only on fill.
   subroutine mat_multspa(scale1, scale2, a, b, c, fill, density_limit, sparse_count, dense_count)
      real(dp), intent(in) :: scale1, scale2
      type(csdp_block_matrix), intent(in) :: a, b
      type(csdp_block_matrix), intent(inout) :: c
      type(csdp_fill_workspace), intent(in) :: fill
      real(dp), intent(in), optional :: density_limit
      integer, intent(inout), optional :: sparse_count, dense_count
      real(dp) :: lim, temp
      integer :: ib, e, p, q

      lim = 0.01_dp
      if (present(density_limit)) lim = density_limit
      call scale_output(scale2,c)
      if (abs(scale1) <= tiny(1.0_dp)) return

      do ib = 1, size(a%block)
         if (a%block(ib)%category == csdp_diag) then
            c%block(ib)%diag = c%block(ib)%diag + scale1*a%block(ib)%diag*b%block(ib)%diag
         else if (fill%block(ib)%density() > lim) then
            c%block(ib)%mat = c%block(ib)%mat + scale1*matmul(a%block(ib)%mat,b%block(ib)%mat)
            if (present(dense_count)) dense_count = dense_count + 1
         else
            !$omp parallel do default(shared) private(p,e,q,temp) schedule(static)
            do p = 1, fill%block(ib)%n
               do e = fill%block(ib)%row_ptr(p), fill%block(ib)%row_ptr(p+1)-1
                  q = fill%block(ib)%j(e)
                  temp = scale1*a%block(ib)%mat(p,q)
                  c%block(ib)%mat(p,:) = c%block(ib)%mat(p,:) + temp*b%block(ib)%mat(q,:)
               end do
            end do
            !$omp end parallel do
            if (present(sparse_count)) sparse_count = sparse_count + 1
         end if
      end do
   end subroutine mat_multspa

   ! C = scale1*A*B + scale2*C, but generate only entries in the fill set.
   subroutine mat_multspc(scale1, scale2, a, b, c, fill, density_limit, sparse_count, dense_count)
      real(dp), intent(in) :: scale1, scale2
      type(csdp_block_matrix), intent(in) :: a, b
      type(csdp_block_matrix), intent(inout) :: c
      type(csdp_fill_workspace), intent(in) :: fill
      real(dp), intent(in), optional :: density_limit
      integer, intent(inout), optional :: sparse_count, dense_count
      real(dp) :: lim, temp
      integer :: ib, e, p, q

      lim = 0.01_dp
      if (present(density_limit)) lim = density_limit
      call scale_output(scale2,c)
      if (abs(scale1) <= tiny(1.0_dp)) return

      do ib = 1, size(a%block)
         if (a%block(ib)%category == csdp_diag) then
            c%block(ib)%diag = c%block(ib)%diag + scale1*a%block(ib)%diag*b%block(ib)%diag
         else if (fill%block(ib)%density() > lim) then
            c%block(ib)%mat = c%block(ib)%mat + scale1*matmul(a%block(ib)%mat,b%block(ib)%mat)
            if (present(dense_count)) dense_count = dense_count + 1
         else
            !$omp parallel do default(shared) private(e,p,q,temp) schedule(static)
            do e = 1, fill%block(ib)%nnz()
               p = fill%block(ib)%i(e)
               q = fill%block(ib)%j(e)
               temp = dot_product(a%block(ib)%mat(p,:),b%block(ib)%mat(:,q))
               c%block(ib)%mat(p,q) = c%block(ib)%mat(p,q) + scale1*temp
            end do
            !$omp end parallel do
            if (present(sparse_count)) sparse_count = sparse_count + 1
         end if
      end do
   end subroutine mat_multspc

   subroutine apply_o_fill(prob, zi, x, y, result, fill, density_limit, sparse_count, dense_count)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), intent(in) :: zi, x
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: result(:)
      type(csdp_fill_workspace), intent(in) :: fill
      real(dp), intent(in), optional :: density_limit
      integer, intent(inout), optional :: sparse_count, dense_count
      type(csdp_block_matrix) :: aty, tmp, prod

      if (size(y) /= size(prob%a) .or. size(result) /= size(prob%a)) error stop 'apply_o_fill: wrong size'
      call op_at_sparse(prob,y,aty)
      tmp = prob%c
      prod = prob%c
      call mat_multspa(1.0_dp,0.0_dp,aty,x,tmp,fill,density_limit,sparse_count,dense_count)
      call mat_multspc(1.0_dp,0.0_dp,zi,tmp,prod,fill,density_limit,sparse_count,dense_count)
      call op_a_sparse(prob,prod,result)
   end subroutine apply_o_fill

end module rcsdp_fill_ops
