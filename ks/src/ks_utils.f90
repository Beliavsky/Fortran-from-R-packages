! SPDX-License-Identifier: GPL-2.0-only
module ks_utils
   use ks_kinds, only: dp
   use ks_linalg, only: covariance_matrix, matrix_sqrt, spd_inverse
   implicit none
   private
   public :: vec, vech, invvec, invvech, pre_scale, pre_sphere
   public :: block_indices, matrix_power_int, kron, kron_power, row_kron_power
   public :: symmetrizer_apply, symmetrizer_matrix, lp_grid_diff
contains
   function vec(a,byrow) result(v)
      real(dp), intent(in) :: a(:,:)
      logical, intent(in), optional :: byrow
      real(dp), allocatable :: v(:)
      logical :: br
      integer :: i,j,k
      br=.false.; if(present(byrow)) br=byrow
      allocate(v(size(a)))
      k=0
      if (br) then
         do i=1,size(a,1); do j=1,size(a,2); k=k+1; v(k)=a(i,j); end do; end do
      else
         do j=1,size(a,2); do i=1,size(a,1); k=k+1; v(k)=a(i,j); end do; end do
      end if
   end function vec

   function vech(a) result(v)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: v(:)
      integer :: d,j,i,k
      d=size(a,1)
      if (size(a,2)/=d) then; allocate(v(0)); return; end if
      allocate(v(d*(d+1)/2)); k=0
      do j=1,d
         do i=j,d
            k=k+1; v(k)=a(i,j)
         end do
      end do
   end function vech

   subroutine invvec(v,nrow,ncol,a,byrow)
      real(dp), intent(in) :: v(:)
      integer, intent(in) :: nrow,ncol
      real(dp), intent(out) :: a(nrow,ncol)
      logical, intent(in), optional :: byrow
      logical :: br
      integer :: i,j,k
      br=.false.; if(present(byrow)) br=byrow
      if(size(v)/=nrow*ncol) then; a=0.0_dp; return; end if
      k=0
      if(br) then
         do i=1,nrow; do j=1,ncol; k=k+1; a(i,j)=v(k); end do; end do
      else
         do j=1,ncol; do i=1,nrow; k=k+1; a(i,j)=v(k); end do; end do
      end if
   end subroutine invvec

   subroutine invvech(v,a,info)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable, intent(out) :: a(:,:)
      integer, intent(out), optional :: info
      integer :: d,i,j,k
      d=nint((-1.0_dp+sqrt(1.0_dp+8.0_dp*real(size(v),dp)))/2.0_dp)
      if(d*(d+1)/2/=size(v)) then
         allocate(a(0,0)); if(present(info)) info=-1; return
      end if
      allocate(a(d,d)); a=0.0_dp; k=0
      do j=1,d
         do i=j,d
            k=k+1; a(i,j)=v(k); a(j,i)=v(k)
         end do
      end do
      if(present(info)) info=0
   end subroutine invvech

   subroutine pre_scale(x,xstar,center,scale)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: xstar(size(x,1),size(x,2))
      real(dp), intent(out), optional :: center(size(x,2)),scale(size(x,2))
      real(dp) :: mu(size(x,2)),s(size(x,2)),z(size(x,1),size(x,2))
      integer :: j,n
      n=size(x,1); mu=sum(x,dim=1)/real(max(n,1),dp); z=x-spread(mu,1,n)
      do j=1,size(x,2)
         if(n>1) then
            s(j)=sqrt(sum(z(:,j)**2)/real(n-1,dp))
         else
            s(j)=1.0_dp
         end if
         if(s(j)<=sqrt(tiny(1.0_dp))) s(j)=1.0_dp
         xstar(:,j)=z(:,j)/s(j)
      end do
      if(present(center)) center=mu
      if(present(scale)) scale=s
   end subroutine pre_scale

   subroutine pre_sphere(x,xstar,center,root,info)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: xstar(size(x,1),size(x,2))
      real(dp), intent(out), optional :: center(size(x,2)),root(size(x,2),size(x,2))
      integer, intent(out), optional :: info
      real(dp) :: mu(size(x,2)),s(size(x,2),size(x,2)),r(size(x,2),size(x,2)),rinv(size(x,2),size(x,2))
      integer :: ierr
      call covariance_matrix(x,s,mu)
      call matrix_sqrt(s,r,ierr)
      if(ierr==0) call spd_inverse(r,rinv,ierr)
      if(ierr==0) then
         xstar=matmul(x-spread(mu,1,size(x,1)),rinv)
      else
         xstar=0.0_dp
      end if
      if(present(center)) center=mu
      if(present(root)) root=r
      if(present(info)) info=ierr
   end subroutine pre_sphere

   subroutine block_indices(nx,ny,d,r,diff,block_limit,bounds)
      integer, intent(in) :: nx,ny,d,r
      logical, intent(in) :: diff
      integer, intent(in), optional :: block_limit
      integer, allocatable, intent(out) :: bounds(:,:)
      integer :: lim,npg,nb,k,lo,hi
      lim=1000000; if(present(block_limit)) lim=block_limit
      if(diff) then
         npg=max(lim/max(nx*d**r,1),1)
      else
         npg=max(lim/max(nx,1),1)
      end if
      nb=(ny+npg-1)/npg
      allocate(bounds(nb,2))
      do k=1,nb
         lo=1+(k-1)*npg; hi=min(ny,k*npg)
         bounds(k,:)=[lo,hi]
      end do
   end subroutine block_indices

   function matrix_power_int(a,n) result(b)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: n
      real(dp) :: b(size(a,1),size(a,2)),base(size(a,1),size(a,2)),ainv(size(a,1),size(a,2))
      integer :: p,info,i
      b=0.0_dp; do i=1,size(a,1); b(i,i)=1.0_dp; end do
      if(n==0) return
      if(n<0) then
         call spd_inverse(a,ainv,info); if(info/=0) then; b=0.0_dp; return; end if
         base=ainv; p=-n
      else
         base=a; p=n
      end if
      do while(p>0)
         if(mod(p,2)==1) b=matmul(b,base)
         p=p/2
         if(p>0) base=matmul(base,base)
      end do
   end function matrix_power_int

   function kron(a,b) result(c)
      real(dp), intent(in) :: a(:,:),b(:,:)
      real(dp), allocatable :: c(:,:)
      integer :: i,j,m,n,p,q
      m=size(a,1); n=size(a,2); p=size(b,1); q=size(b,2)
      allocate(c(m*p,n*q))
      do j=1,n; do i=1,m
         c((i-1)*p+1:i*p,(j-1)*q+1:j*q)=a(i,j)*b
      end do; end do
   end function kron

   function kron_power(a,power) result(c)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: power
      real(dp), allocatable :: c(:,:),tmp(:,:)
      integer :: k
      if(power==0) then; allocate(c(1,1)); c=1.0_dp; return; end if
      c=a
      do k=2,power
         tmp=kron(c,a); call move_alloc(tmp,c)
      end do
   end function kron_power

   function row_kron_power(a,power) result(out)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: power
      real(dp), allocatable :: out(:,:)
      integer :: n,d,m,i,code,k,c
      n=size(a,1); d=size(a,2); m=d**power
      if(power==0) then; allocate(out(n,1)); out=1.0_dp; return; end if
      allocate(out(n,m))
      do i=1,n
         do code=0,m-1
            c=code; out(i,code+1)=1.0_dp
            do k=1,power
               out(i,code+1)=out(i,code+1)*a(i,mod(c,d)+1); c=c/d
            end do
         end do
      end do
   end function row_kron_power

   recursive subroutine permute_positions(a,l,perms,count)
      integer, intent(inout) :: a(:)
      integer, intent(in) :: l
      integer, intent(inout) :: perms(:,:),count
      integer :: i,t
      if(l>size(a)) then
         count=count+1; perms(count,:)=a; return
      end if
      do i=l,size(a)
         t=a(l); a(l)=a(i); a(i)=t
         call permute_positions(a,l+1,perms,count)
         t=a(l); a(l)=a(i); a(i)=t
      end do
   end subroutine permute_positions

   integer function factorial_int(n) result(v)
      integer,intent(in)::n
      integer::i
      v=1; do i=2,n; v=v*i; end do
   end function factorial_int

   subroutine decode_code(code,d,r,idx)
      integer,intent(in)::code,d,r
      integer,intent(out)::idx(r)
      integer::c,k
      c=code
      do k=1,r; idx(k)=mod(c,d)+1; c=c/d; end do
   end subroutine decode_code

   integer function encode_idx(idx,d) result(code)
      integer,intent(in)::idx(:),d
      integer::k,p
      code=0; p=1
      do k=1,size(idx); code=code+(idx(k)-1)*p; p=p*d; end do
   end function encode_idx

   subroutine symmetrizer_apply(d,r,v,sv)
      integer,intent(in)::d,r
      real(dp),intent(in)::v(:)
      real(dp),intent(out)::sv(size(v))
      integer :: nperm,count,code,p,k
      integer, allocatable :: perms(:,:),base(:),idx(:),idxp(:)
      if(r<=1) then; sv=v; return; end if
      if(size(v)/=d**r) then; sv=0.0_dp; return; end if
      nperm=factorial_int(r)
      allocate(perms(nperm,r),base(r),idx(r),idxp(r))
      base=[(k,k=1,r)]; count=0; call permute_positions(base,1,perms,count)
      do code=0,d**r-1
         call decode_code(code,d,r,idx); sv(code+1)=0.0_dp
         do p=1,nperm
            do k=1,r; idxp(k)=idx(perms(p,k)); end do
            sv(code+1)=sv(code+1)+v(encode_idx(idxp,d)+1)
         end do
         sv(code+1)=sv(code+1)/real(nperm,dp)
      end do
   end subroutine symmetrizer_apply

   subroutine symmetrizer_matrix(d,r,s)
      integer,intent(in)::d,r
      real(dp),allocatable,intent(out)::s(:,:)
      real(dp),allocatable::e(:),sv(:)
      integer::m,j
      m=d**r; allocate(s(m,m),e(m),sv(m)); s=0.0_dp
      do j=1,m
         e=0.0_dp; e(j)=1.0_dp; call symmetrizer_apply(d,r,e,sv); s(:,j)=sv
      end do
   end subroutine symmetrizer_matrix

   function lp_grid_diff(f1,f2,spacing,p) result(val)
      real(dp),intent(in)::f1(:),f2(:),spacing(:)
      real(dp),intent(in),optional::p
      real(dp)::val,pp,cell
      pp=2.0_dp; if(present(p)) pp=p
      cell=product(spacing)
      val=sum(abs(f1-f2)**pp)*cell
   end function lp_grid_diff
end module ks_utils
