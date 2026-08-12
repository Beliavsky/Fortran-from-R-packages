! DSDP-style abstract symmetric data matrices: dense, sparse, and low rank.
! DSDP copyright/license: see licenses/DSDP-LICENSE.
module rdsdp_data
   use rdsdp_kinds, only : dp
   use rdsdp_types, only : dsdp_block, dsdp_sparse_sym, dsdp_lowrank_sym, &
      dsdp_data_dense, dsdp_data_sparse, dsdp_data_lowrank
   implicit none
   private
   public :: sparse_from_dense, sparse_from_upper, sparse_to_dense
   public :: get_data_dense, data_add_scaled, data_dot, data_fnorm2, data_nnz
   public :: set_constraint_lowrank, set_objective_lowrank
   public :: data_schur_pair, data_storage, compress_sdp_block, densify_sdp_block

contains

   subroutine sparse_from_dense(a,sp,tol)
      real(dp), intent(in) :: a(:,:)
      type(dsdp_sparse_sym), intent(out) :: sp
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      integer :: i,j,k,n,nnz
      n=size(a,1)
      if (size(a,2)/=n) error stop 'sparse_from_dense: matrix must be square'
      eps=0.0_dp; if (present(tol)) eps=max(0.0_dp,tol)
      nnz=0
      do j=1,n
         do i=1,n
            if (abs(a(i,j))>eps) nnz=nnz+1
         end do
      end do
      sp%n=n; sp%nnz=nnz
      allocate(sp%row(nnz),sp%col(nnz),sp%val(nnz))
      k=0
      do j=1,n
         do i=1,n
            if (abs(a(i,j))>eps) then
               k=k+1; sp%row(k)=i; sp%col(k)=j; sp%val(k)=a(i,j)
            end if
         end do
      end do
   end subroutine sparse_from_dense

   subroutine sparse_from_upper(n,row,col,val,sp)
      integer, intent(in) :: n,row(:),col(:)
      real(dp), intent(in) :: val(:)
      type(dsdp_sparse_sym), intent(out) :: sp
      integer :: k,p,nnz
      if (size(row)/=size(col) .or. size(row)/=size(val)) error stop 'sparse_from_upper: size mismatch'
      nnz=0
      do k=1,size(val)
         if (abs(val(k))<=tiny(1.0_dp)) cycle
         nnz=nnz+1
         if (row(k)/=col(k)) nnz=nnz+1
      end do
      sp%n=n; sp%nnz=nnz
      allocate(sp%row(nnz),sp%col(nnz),sp%val(nnz))
      p=0
      do k=1,size(val)
         if (abs(val(k))<=tiny(1.0_dp)) cycle
         p=p+1; sp%row(p)=row(k); sp%col(p)=col(k); sp%val(p)=val(k)
         if (row(k)/=col(k)) then
            p=p+1; sp%row(p)=col(k); sp%col(p)=row(k); sp%val(p)=val(k)
         end if
      end do
   end subroutine sparse_from_upper

   subroutine sparse_to_dense(sp,a)
      type(dsdp_sparse_sym), intent(in) :: sp
      real(dp), allocatable, intent(out) :: a(:,:)
      integer :: k
      allocate(a(sp%n,sp%n)); a=0.0_dp
      do k=1,sp%nnz
         a(sp%row(k),sp%col(k))=a(sp%row(k),sp%col(k))+sp%val(k)
      end do
   end subroutine sparse_to_dense

   integer function data_storage(block,idx) result(storage)
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: idx
      if (idx==0) then
         storage=block%c_storage
      else
         if (allocated(block%a_storage)) then
            storage=block%a_storage(idx)
         else
            storage=dsdp_data_dense
         end if
      end if
   end function data_storage

   subroutine get_data_dense(block,idx,a)
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: idx
      real(dp), allocatable, intent(out) :: a(:,:)
      integer :: storage,k
      storage=data_storage(block,idx)
      select case(storage)
      case(dsdp_data_dense)
         allocate(a(block%n,block%n))
         if (idx==0) then
            if (.not.allocated(block%c)) error stop 'get_data_dense: dense objective unavailable'
            a=block%c
         else
            if (.not.allocated(block%a)) error stop 'get_data_dense: dense constraints unavailable'
            a=block%a(:,:,idx)
         end if
      case(dsdp_data_sparse)
         if (idx==0) then
            call sparse_to_dense(block%c_sparse,a)
         else
            call sparse_to_dense(block%a_sparse(idx),a)
         end if
      case(dsdp_data_lowrank)
         allocate(a(block%n,block%n)); a=0.0_dp
         if (idx==0) then
            do k=1,block%c_lowrank%rank
               a=a+block%c_lowrank%coeff(k)*outer(block%c_lowrank%vec(:,k),block%c_lowrank%vec(:,k))
            end do
         else
            do k=1,block%a_lowrank(idx)%rank
               a=a+block%a_lowrank(idx)%coeff(k)*outer(block%a_lowrank(idx)%vec(:,k),block%a_lowrank(idx)%vec(:,k))
            end do
         end if
      case default
         error stop 'get_data_dense: unknown storage'
      end select
   end subroutine get_data_dense

   subroutine data_add_scaled(block,idx,alpha,x)
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: idx
      real(dp), intent(in) :: alpha
      real(dp), intent(inout) :: x(:,:)
      integer :: storage
      real(dp), allocatable :: a(:,:)
      if (abs(alpha)<=tiny(1.0_dp)) return
      storage=data_storage(block,idx)
      select case(storage)
      case(dsdp_data_dense)
         if (idx==0) then
            x=x+alpha*block%c
         else
            x=x+alpha*block%a(:,:,idx)
         end if
      case(dsdp_data_sparse)
         if (idx==0) then
            call add_sparse(block%c_sparse,alpha,x)
         else
            call add_sparse(block%a_sparse(idx),alpha,x)
         end if
      case(dsdp_data_lowrank)
         if (idx==0) then
            call add_lowrank(block%c_lowrank,alpha,x)
         else
            call add_lowrank(block%a_lowrank(idx),alpha,x)
         end if
      case default
         call get_data_dense(block,idx,a); x=x+alpha*a
      end select
   end subroutine data_add_scaled

   real(dp) function data_dot(block,idx,x) result(v)
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: idx
      real(dp), intent(in) :: x(:,:)
      integer :: storage
      v=0.0_dp; storage=data_storage(block,idx)
      select case(storage)
      case(dsdp_data_dense)
         if (idx==0) then; v=sum(block%c*x); else; v=sum(block%a(:,:,idx)*x); end if
      case(dsdp_data_sparse)
         if (idx==0) then; v=dot_sparse(block%c_sparse,x); else; v=dot_sparse(block%a_sparse(idx),x); end if
      case(dsdp_data_lowrank)
         if (idx==0) then; v=dot_lowrank(block%c_lowrank,x); else; v=dot_lowrank(block%a_lowrank(idx),x); end if
      end select
   end function data_dot

   subroutine add_sparse(sp,alpha,x)
      type(dsdp_sparse_sym), intent(in) :: sp
      real(dp), intent(in) :: alpha
      real(dp), intent(inout) :: x(:,:)
      integer :: k
      do k=1,sp%nnz
         x(sp%row(k),sp%col(k))=x(sp%row(k),sp%col(k))+alpha*sp%val(k)
      end do
   end subroutine add_sparse

   subroutine add_lowrank(lr,alpha,x)
      type(dsdp_lowrank_sym), intent(in) :: lr
      real(dp), intent(in) :: alpha
      real(dp), intent(inout) :: x(:,:)
      integer :: k
      do k=1,lr%rank
         x=x+alpha*lr%coeff(k)*outer(lr%vec(:,k),lr%vec(:,k))
      end do
   end subroutine add_lowrank

   real(dp) function dot_sparse(sp,x) result(v)
      type(dsdp_sparse_sym), intent(in) :: sp
      real(dp), intent(in) :: x(:,:)
      integer :: k
      v=0.0_dp
      do k=1,sp%nnz; v=v+sp%val(k)*x(sp%row(k),sp%col(k)); end do
   end function dot_sparse

   real(dp) function dot_lowrank(lr,x) result(v)
      type(dsdp_lowrank_sym), intent(in) :: lr
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable :: tmp(:)
      integer :: k
      allocate(tmp(lr%n)); v=0.0_dp
      do k=1,lr%rank
         tmp=matmul(x,lr%vec(:,k))
         v=v+lr%coeff(k)*dot_product(lr%vec(:,k),tmp)
      end do
   end function dot_lowrank

   real(dp) function data_fnorm2(block,idx) result(v)
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: idx
      real(dp), allocatable :: a(:,:)
      call get_data_dense(block,idx,a)
      v=sum(a*a)
   end function data_fnorm2

   integer function data_nnz(block,idx) result(nnz)
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: idx
      integer :: storage
      storage=data_storage(block,idx)
      select case(storage)
      case(dsdp_data_sparse)
         if (idx==0) then; nnz=block%c_sparse%nnz; else; nnz=block%a_sparse(idx)%nnz; end if
      case(dsdp_data_lowrank)
         if (idx==0) then; nnz=block%c_lowrank%n*block%c_lowrank%n; else; nnz=block%a_lowrank(idx)%n**2; end if
      case default
         if (idx==0) then
            if (.not.allocated(block%c)) then; nnz=0; return; end if
            nnz=count(abs(block%c)>tiny(1.0_dp))
         else
            if (.not.allocated(block%a)) then; nnz=0; return; end if
            nnz=count(abs(block%a(:,:,idx))>tiny(1.0_dp))
         end if
      end select
   end function data_nnz

   subroutine compress_sdp_block(block,m,threshold,release_dense)
      type(dsdp_block), intent(inout) :: block
      integer, intent(in) :: m
      real(dp), intent(in) :: threshold
      logical, intent(in), optional :: release_dense
      logical :: rel
      integer :: i,n2,nnz
      rel=.true.; if (present(release_dense)) rel=release_dense
      if (.not.allocated(block%a_storage)) then
         allocate(block%a_storage(m)); block%a_storage=dsdp_data_dense
      end if
      if (.not.allocated(block%a_sparse)) allocate(block%a_sparse(m))
      if (.not.allocated(block%a_lowrank)) allocate(block%a_lowrank(m))
      n2=block%n*block%n
      if (allocated(block%c) .and. block%c_storage==dsdp_data_dense) then
         nnz=count(abs(block%c)>tiny(1.0_dp))
         if (real(nnz,dp)<=threshold*real(n2,dp)) then
            call sparse_from_dense(block%c,block%c_sparse)
            block%c_storage=dsdp_data_sparse
            if (rel) deallocate(block%c)
         else
            block%c_storage=dsdp_data_dense
         end if
      end if
      if (allocated(block%a)) then
         do i=1,m
            if (block%a_storage(i)/=dsdp_data_dense) cycle
            nnz=count(abs(block%a(:,:,i))>tiny(1.0_dp))
            if (real(nnz,dp)<=threshold*real(n2,dp)) then
               call sparse_from_dense(block%a(:,:,i),block%a_sparse(i))
               block%a_storage(i)=dsdp_data_sparse
            end if
         end do
         if (rel .and. all(block%a_storage/=dsdp_data_dense)) deallocate(block%a)
      end if
   end subroutine compress_sdp_block

   subroutine densify_sdp_block(block,m)
      type(dsdp_block), intent(inout) :: block
      integer, intent(in) :: m
      real(dp), allocatable :: tmp(:,:), anew(:,:,:)
      integer :: i
      if (block%c_storage/=dsdp_data_dense) then
         call get_data_dense(block,0,tmp)
         if (allocated(block%c)) deallocate(block%c)
         allocate(block%c(block%n,block%n)); block%c=tmp
         block%c_storage=dsdp_data_dense
      end if
      if (.not.allocated(block%a_storage)) then
         allocate(block%a_storage(m)); block%a_storage=dsdp_data_dense
         return
      end if
      if (any(block%a_storage/=dsdp_data_dense)) then
         allocate(anew(block%n,block%n,m))
         do i=1,m
            call get_data_dense(block,i,tmp); anew(:,:,i)=tmp
         end do
         if (allocated(block%a)) deallocate(block%a)
         call move_alloc(anew,block%a)
         block%a_storage=dsdp_data_dense
      end if
   end subroutine densify_sdp_block

   subroutine set_constraint_lowrank(block,idx,coeff,vec)
      type(dsdp_block), intent(inout) :: block
      integer, intent(in) :: idx
      real(dp), intent(in) :: coeff(:),vec(:,:)
      integer :: m
      if (size(vec,1)/=block%n .or. size(vec,2)/=size(coeff)) error stop 'set_constraint_lowrank: size mismatch'
      if (.not.allocated(block%a_storage)) error stop 'set_constraint_lowrank: a_storage not allocated'
      m=size(block%a_storage)
      if (idx<1 .or. idx>m) error stop 'set_constraint_lowrank: index out of range'
      if (.not.allocated(block%a_lowrank)) allocate(block%a_lowrank(m))
      block%a_lowrank(idx)%n=block%n; block%a_lowrank(idx)%rank=size(coeff)
      allocate(block%a_lowrank(idx)%coeff(size(coeff)),block%a_lowrank(idx)%vec(size(vec,1),size(vec,2)))
      block%a_lowrank(idx)%coeff=coeff; block%a_lowrank(idx)%vec=vec
      block%a_storage(idx)=dsdp_data_lowrank
   end subroutine set_constraint_lowrank

   subroutine set_objective_lowrank(block,coeff,vec)
      type(dsdp_block), intent(inout) :: block
      real(dp), intent(in) :: coeff(:),vec(:,:)
      if (size(vec,1)/=block%n .or. size(vec,2)/=size(coeff)) error stop 'set_objective_lowrank: size mismatch'
      block%c_lowrank%n=block%n; block%c_lowrank%rank=size(coeff)
      allocate(block%c_lowrank%coeff(size(coeff)),block%c_lowrank%vec(size(vec,1),size(vec,2)))
      block%c_lowrank%coeff=coeff; block%c_lowrank%vec=vec; block%c_storage=dsdp_data_lowrank
      if (allocated(block%c)) deallocate(block%c)
   end subroutine set_objective_lowrank

   real(dp) function data_schur_pair(block,i,j,sinv) result(v)
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: i,j
      real(dp), intent(in) :: sinv(:,:)
      integer :: si,sj
      real(dp), allocatable :: ai(:,:),bj(:,:)
      si=data_storage(block,i); sj=data_storage(block,j)
      if (si==dsdp_data_sparse .and. sj==dsdp_data_sparse) then
         v=sparse_sparse_pair(block%a_sparse(i),block%a_sparse(j),sinv)
      else if (si==dsdp_data_lowrank .and. sj==dsdp_data_lowrank) then
         v=lowrank_lowrank_pair(block%a_lowrank(i),block%a_lowrank(j),sinv)
      else if (si==dsdp_data_sparse .and. sj==dsdp_data_lowrank) then
         v=sparse_lowrank_pair(block%a_sparse(i),block%a_lowrank(j),sinv)
      else if (si==dsdp_data_lowrank .and. sj==dsdp_data_sparse) then
         v=sparse_lowrank_pair(block%a_sparse(j),block%a_lowrank(i),sinv)
      else
         ! For any pair involving a dense matrix, form S*A_j*S once locally.
         call get_data_dense(block,j,ai)
         bj=matmul(matmul(sinv,ai),sinv)
         v=data_dot(block,i,bj)
      end if
   end function data_schur_pair

   real(dp) function sparse_sparse_pair(a,b,s) result(v)
      type(dsdp_sparse_sym), intent(in) :: a,b
      real(dp), intent(in) :: s(:,:)
      integer :: p,q
      v=0.0_dp
      do q=1,b%nnz
         do p=1,a%nnz
            v=v+a%val(p)*b%val(q)*s(b%col(q),a%row(p))*s(a%col(p),b%row(q))
         end do
      end do
   end function sparse_sparse_pair

   real(dp) function lowrank_lowrank_pair(a,b,s) result(v)
      type(dsdp_lowrank_sym), intent(in) :: a,b
      real(dp), intent(in) :: s(:,:)
      real(dp), allocatable :: sv(:,:)
      real(dp) :: t
      integer :: p,q
      allocate(sv(size(s,1),b%rank)); sv=matmul(s,b%vec)
      v=0.0_dp
      do q=1,b%rank
         do p=1,a%rank
            t=dot_product(a%vec(:,p),sv(:,q))
            v=v+a%coeff(p)*b%coeff(q)*t*t
         end do
      end do
   end function lowrank_lowrank_pair

   real(dp) function sparse_lowrank_pair(a,b,s) result(v)
      type(dsdp_sparse_sym), intent(in) :: a
      type(dsdp_lowrank_sym), intent(in) :: b
      real(dp), intent(in) :: s(:,:)
      real(dp), allocatable :: sv(:,:)
      integer :: e,q
      allocate(sv(size(s,1),b%rank)); sv=matmul(s,b%vec)
      v=0.0_dp
      do q=1,b%rank
         do e=1,a%nnz
            v=v+a%val(e)*b%coeff(q)*sv(a%row(e),q)*sv(a%col(e),q)
         end do
      end do
   end function sparse_lowrank_pair

   pure function outer(x,y) result(a)
      real(dp), intent(in) :: x(:),y(:)
      real(dp) :: a(size(x),size(y))
      integer :: j
      do j=1,size(y); a(:,j)=x*y(j); end do
   end function outer

end module rdsdp_data
