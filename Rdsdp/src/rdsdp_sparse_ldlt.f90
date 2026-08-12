! Sparse symmetric LDL^T factorization for Schur systems.
!
! The numerical factorization is a modern Fortran translation of the core
! QDLDL elimination-tree/factor/solve algorithm (Apache-2.0).  RCM ordering
! and the cache/reuse layer are original to this translation.  See
! licenses/QDLDL-APACHE-2.0.txt and TRANSLATION_NOTES.md.
module rdsdp_sparse_ldlt
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rdsdp_kinds, only : dp
   implicit none
   private

   type, public :: sparse_ldlt_cache
      integer :: n = 0
      logical :: analyzed = .false.
      integer, allocatable :: perm(:)
      integer, allocatable :: ap(:), ai(:)
      integer, allocatable :: lnz(:), etree(:)
      integer, allocatable :: lp(:), li(:)
      real(dp), allocatable :: lx(:), d(:), dinv(:)
      integer :: symbolic_analyses = 0
      integer :: numeric_factorizations = 0
      integer :: matrix_nnz = 0
      integer :: factor_nnz = 0
   end type sparse_ldlt_cache

   public :: sparse_ldlt_solve, reset_sparse_ldlt_cache

contains

   subroutine reset_sparse_ldlt_cache(cache)
      type(sparse_ldlt_cache), intent(inout) :: cache
      cache%n=0; cache%analyzed=.false.
      if (allocated(cache%perm)) deallocate(cache%perm)
      if (allocated(cache%ap)) deallocate(cache%ap)
      if (allocated(cache%ai)) deallocate(cache%ai)
      if (allocated(cache%lnz)) deallocate(cache%lnz)
      if (allocated(cache%etree)) deallocate(cache%etree)
      if (allocated(cache%lp)) deallocate(cache%lp)
      if (allocated(cache%li)) deallocate(cache%li)
      if (allocated(cache%lx)) deallocate(cache%lx)
      if (allocated(cache%d)) deallocate(cache%d)
      if (allocated(cache%dinv)) deallocate(cache%dinv)
      cache%matrix_nnz=0; cache%factor_nnz=0
      cache%symbolic_analyses=0; cache%numeric_factorizations=0
   end subroutine reset_sparse_ldlt_cache

   subroutine sparse_ldlt_solve(a,b,x,reg,drop_tol,cache,ok)
      real(dp), intent(in) :: a(:,:),b(:),reg,drop_tol
      real(dp), intent(out) :: x(:)
      type(sparse_ldlt_cache), intent(inout) :: cache
      logical, intent(out) :: ok
      integer, allocatable :: perm(:),ap(:),ai(:),work(:)
      integer, allocatable :: y_idx(:),elim(:),next_space(:)
      logical, allocatable :: markers(:)
      real(dp), allocatable :: ax(:),rhs(:),y_vals(:)
      integer :: n,sum_lnz,npositive,i,old_symbolic,old_numeric
      logical :: same_pattern

      n=size(b)
      if (size(a,1)/=n .or. size(a,2)/=n .or. size(x)/=n) then
         ok=.false.; return
      end if

      call rcm_order(a,drop_tol,perm)
      call dense_upper_csc(a,perm,reg,drop_tol,ap,ai,ax)

      same_pattern=cache%analyzed .and. cache%n==n
      if (same_pattern) then
         same_pattern=allocated(cache%perm) .and. allocated(cache%ap) .and. allocated(cache%ai)
      end if
      if (same_pattern) then
         same_pattern=all(cache%perm==perm) .and. size(cache%ap)==size(ap) .and. size(cache%ai)==size(ai)
      end if
      if (same_pattern) then
         same_pattern=all(cache%ap==ap) .and. all(cache%ai==ai)
      end if

      if (.not.same_pattern) then
         old_symbolic=cache%symbolic_analyses; old_numeric=cache%numeric_factorizations
         call reset_sparse_ldlt_cache(cache)
         cache%symbolic_analyses=old_symbolic; cache%numeric_factorizations=old_numeric
         cache%n=n
         call move_alloc(perm,cache%perm)
         call move_alloc(ap,cache%ap)
         call move_alloc(ai,cache%ai)
         allocate(cache%lnz(n),cache%etree(n),work(n))
         sum_lnz=qdldl_etree(n,cache%ap,cache%ai,work,cache%lnz,cache%etree)
         if (sum_lnz<0) then
            ok=.false.; return
         end if
         allocate(cache%lp(n+1),cache%li(sum_lnz),cache%lx(sum_lnz),cache%d(n),cache%dinv(n))
         cache%factor_nnz=sum_lnz
         cache%matrix_nnz=size(cache%ai)
         cache%analyzed=.true.
         cache%symbolic_analyses=cache%symbolic_analyses+1
      end if

      ! If move_alloc occurred above, AX still corresponds to cache pattern.  If
      ! the pattern was reused, PERM/AP/AI are temporary and can simply vanish.
      allocate(markers(n),y_idx(n),elim(n),next_space(n),y_vals(n))
      markers=.false.; y_idx=0; elim=0; next_space=0; y_vals=0.0_dp
      npositive=qdldl_factor(n,cache%ap,cache%ai,ax,cache%lp,cache%li,cache%lx, &
         cache%d,cache%dinv,cache%lnz,cache%etree,markers,y_idx,elim,next_space,y_vals)
      cache%numeric_factorizations=cache%numeric_factorizations+1
      if (npositive/=n) then
         ok=.false.; return
      end if

      allocate(rhs(n)); rhs=b(cache%perm)
      call qdldl_solve(n,cache%lp,cache%li,cache%lx,cache%dinv,rhs)
      x=0.0_dp
      do i=1,n
         x(cache%perm(i))=rhs(i)
      end do
      ok=all(ieee_is_finite(x))
   end subroutine sparse_ldlt_solve

   subroutine dense_upper_csc(a,perm,reg,drop_tol,ap,ai,ax)
      real(dp), intent(in) :: a(:,:),reg,drop_tol
      integer, intent(in) :: perm(:)
      integer, allocatable, intent(out) :: ap(:),ai(:)
      real(dp), allocatable, intent(out) :: ax(:)
      integer :: n,i,j,k,nnz,oi,oj
      real(dp) :: v,scale,tol
      n=size(perm); nnz=0
      do j=1,n
         oj=perm(j)
         do i=1,j
            oi=perm(i)
            v=a(oi,oj)
            if (i==j) v=v+max(0.0_dp,reg)
            scale=sqrt(max(abs(a(oi,oi))*abs(a(oj,oj)),1.0_dp))
            tol=max(0.0_dp,drop_tol)*scale
            if (i==j .or. abs(v)>tol) nnz=nnz+1
         end do
      end do
      allocate(ap(n+1),ai(nnz),ax(nnz)); k=1; ap(1)=1
      do j=1,n
         oj=perm(j)
         do i=1,j
            oi=perm(i)
            v=a(oi,oj)
            if (i==j) v=v+max(0.0_dp,reg)
            scale=sqrt(max(abs(a(oi,oi))*abs(a(oj,oj)),1.0_dp))
            tol=max(0.0_dp,drop_tol)*scale
            if (i==j .or. abs(v)>tol) then
               ai(k)=i; ax(k)=v; k=k+1
            end if
         end do
         ap(j+1)=k
      end do
   end subroutine dense_upper_csc

   subroutine rcm_order(a,drop_tol,perm)
      real(dp), intent(in) :: a(:,:),drop_tol
      integer, allocatable, intent(out) :: perm(:)
      integer, allocatable :: degree(:),queue(:),nbr(:),order(:)
      logical, allocatable :: seen(:)
      integer :: n,i,j,k,start,head,tail,v,nnbr,tmp,p
      real(dp) :: scale,tol
      n=size(a,1)
      allocate(perm(n),degree(n),queue(n),nbr(n),order(n),seen(n))
      degree=0; seen=.false.; k=0
      do i=1,n
         do j=1,n
            if (i==j) cycle
            scale=sqrt(max(abs(a(i,i))*abs(a(j,j)),1.0_dp)); tol=max(0.0_dp,drop_tol)*scale
            if (abs(a(i,j))>tol) degree(i)=degree(i)+1
         end do
      end do
      do while(k<n)
         start=0
         do i=1,n
            if (.not.seen(i)) then
               if (start==0) then
                  start=i
               else if (degree(i)<degree(start)) then
                  start=i
               end if
            end if
         end do
         head=1; tail=1; queue(1)=start; seen(start)=.true.
         do while(head<=tail)
            v=queue(head); head=head+1; k=k+1; order(k)=v
            nnbr=0
            do j=1,n
               if (j==v .or. seen(j)) cycle
               scale=sqrt(max(abs(a(v,v))*abs(a(j,j)),1.0_dp)); tol=max(0.0_dp,drop_tol)*scale
               if (abs(a(v,j))>tol) then
                  nnbr=nnbr+1; nbr(nnbr)=j
               end if
            end do
            ! Stable insertion sort by degree, then index.
            do i=2,nnbr
               tmp=nbr(i); p=i-1
               do while(p>=1)
                  if (degree(nbr(p))<degree(tmp)) exit
                  if (degree(nbr(p))==degree(tmp) .and. nbr(p)<tmp) exit
                  nbr(p+1)=nbr(p); p=p-1
               end do
               nbr(p+1)=tmp
            end do
            do i=1,nnbr
               if (.not.seen(nbr(i))) then
                  tail=tail+1; queue(tail)=nbr(i); seen(nbr(i))=.true.
               end if
            end do
         end do
      end do
      do i=1,n
         perm(i)=order(n-i+1)
      end do
   end subroutine rcm_order

   integer function qdldl_etree(n,ap,ai,work,lnz,etree) result(sum_lnz)
      integer, intent(in) :: n,ap(:),ai(:)
      integer, intent(out) :: work(:),lnz(:),etree(:)
      integer :: i,j,p
      sum_lnz=0; work=0; lnz=0; etree=0
      do i=1,n
         if (ap(i)==ap(i+1)) then; sum_lnz=-1; return; end if
      end do
      do j=1,n
         work(j)=j
         do p=ap(j),ap(j+1)-1
            i=ai(p)
            if (i>j .or. i<1) then; sum_lnz=-1; return; end if
            do while(work(i)/=j)
               if (etree(i)==0) etree(i)=j
               lnz(i)=lnz(i)+1; work(i)=j; i=etree(i)
               if (i==0) exit
            end do
         end do
      end do
      sum_lnz=sum(lnz)
   end function qdldl_etree

   integer function qdldl_factor(n,ap,ai,ax,lp,li,lx,d,dinv,lnz,etree, &
      markers,y_idx,elim,next_space,y_vals) result(npositive)
      integer, intent(in) :: n,ap(:),ai(:),lnz(:),etree(:)
      real(dp), intent(in) :: ax(:)
      integer, intent(out) :: lp(:),li(:)
      real(dp), intent(out) :: lx(:),d(:),dinv(:)
      logical, intent(inout) :: markers(:)
      integer, intent(inout) :: y_idx(:),elim(:),next_space(:)
      real(dp), intent(inout) :: y_vals(:)
      integer :: i,j,k,p,nnz_y,bidx,cidx,next_idx,nnz_e,tmp_idx
      real(dp) :: y_c
      npositive=0; lp(1)=1
      do i=1,n
         lp(i+1)=lp(i)+lnz(i); markers(i)=.false.; y_vals(i)=0.0_dp
         d(i)=0.0_dp; next_space(i)=lp(i)
      end do
      if (size(lx)>0) lx=0.0_dp
      if (size(li)>0) li=0
      do k=1,n
         nnz_y=0
         do p=ap(k),ap(k+1)-1
            bidx=ai(p)
            if (bidx==k) then; d(k)=d(k)+ax(p); cycle; end if
            if (bidx>k .or. bidx<1) then; npositive=-1; return; end if
            y_vals(bidx)=y_vals(bidx)+ax(p); next_idx=bidx
            if (.not.markers(next_idx)) then
               markers(next_idx)=.true.; elim(1)=next_idx; nnz_e=1; next_idx=etree(bidx)
               do while(next_idx/=0 .and. next_idx<k)
                  if (markers(next_idx)) exit
                  markers(next_idx)=.true.; nnz_e=nnz_e+1; elim(nnz_e)=next_idx; next_idx=etree(next_idx)
               end do
               do while(nnz_e>0)
                  nnz_y=nnz_y+1; y_idx(nnz_y)=elim(nnz_e); nnz_e=nnz_e-1
               end do
            end if
         end do
         do i=nnz_y,1,-1
            cidx=y_idx(i); tmp_idx=next_space(cidx); y_c=y_vals(cidx)
            do j=lp(cidx),tmp_idx-1
               y_vals(li(j))=y_vals(li(j))-lx(j)*y_c
            end do
            li(tmp_idx)=k; lx(tmp_idx)=y_c*dinv(cidx); d(k)=d(k)-y_c*lx(tmp_idx)
            next_space(cidx)=tmp_idx+1; y_vals(cidx)=0.0_dp; markers(cidx)=.false.
         end do
         if (abs(d(k))<=tiny(1.0_dp) .or. .not.ieee_is_finite(d(k))) then; npositive=-1; return; end if
         if (d(k)>0.0_dp) npositive=npositive+1
         dinv(k)=1.0_dp/d(k)
      end do
   end function qdldl_factor

   subroutine qdldl_solve(n,lp,li,lx,dinv,x)
      integer, intent(in) :: n,lp(:),li(:)
      real(dp), intent(in) :: lx(:),dinv(:)
      real(dp), intent(inout) :: x(:)
      integer :: i,j
      real(dp) :: val
      do i=1,n
         val=x(i)
         do j=lp(i),lp(i+1)-1
            x(li(j))=x(li(j))-lx(j)*val
         end do
      end do
      x(1:n)=x(1:n)*dinv(1:n)
      do i=n,1,-1
         val=x(i)
         do j=lp(i),lp(i+1)-1
            val=val-lx(j)*x(li(j))
         end do
         x(i)=val
      end do
   end subroutine qdldl_solve

end module rdsdp_sparse_ldlt
