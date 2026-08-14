module isotone_mregnn
   use isotone_kinds, only : dp
   use isotone_linalg, only : qr_basis
   use nnls, only : nnls_result, nnls_fit, NNLS_SUCCESS
   implicit none
   private
   type, public :: mregnn_result
      real(dp), allocatable :: xb(:)
      real(dp), allocatable :: lambda(:)
      real(dp) :: f = 0.0_dp
      integer :: status = 0
      integer :: rank = 0
   end type mregnn_result
   public :: mregnn, mregnn_monotone, mregnn_positive
contains
   subroutine mregnn(x,y,a,result)
      real(dp),intent(in)::x(:,:),y(:),a(:,:)
      type(mregnn_result),intent(out)::result
      real(dp),allocatable::q(:,:),u(:),v(:,:),rhs(:),tmp(:)
      type(nnls_result)::nr
      integer::n,r,m,rank
      n=size(x,1);m=size(a,1)
      if(size(y)/=n .or. size(a,2)/=n) then;result%status=1;return;end if
      call qr_basis(x,q,rank);r=rank;result%rank=rank
      if(r==0) then
         allocate(result%xb(n),result%lambda(m));result%xb=0.0_dp;result%lambda=0.0_dp
         result%f=dot_product(y,y);result%status=0;return
      end if
      allocate(u(r),v(r,m),rhs(r),tmp(r))
      u=matmul(transpose(q),y)
      v=-matmul(transpose(q),transpose(a))
      call nnls_fit(v,u,nr)
      if(nr%mode/=NNLS_SUCCESS) then;result%status=2;return;end if
      rhs=u-matmul(v,nr%x)
      allocate(result%xb(n),result%lambda(m))
      result%xb=matmul(q,rhs);result%lambda=nr%x
      result%f=sum((y-result%xb)**2);result%status=0
   end subroutine mregnn

   subroutine mregnn_monotone(x,y,result)
      real(dp),intent(in)::x(:,:),y(:)
      type(mregnn_result),intent(out)::result
      real(dp),allocatable::a(:,:)
      integer::n,i
      n=size(x,1);allocate(a(max(0,n-1),n));a=0.0_dp
      do i=1,n-1;a(i,i)=-1.0_dp;a(i,i+1)=1.0_dp;end do
      call mregnn(x,y,a,result)
   end subroutine mregnn_monotone

   subroutine mregnn_positive(x,y,result)
      real(dp),intent(in)::x(:,:),y(:)
      type(mregnn_result),intent(out)::result
      real(dp),allocatable::a(:,:)
      integer::n,i
      n=size(x,1);allocate(a(n,n));a=0.0_dp
      do i=1,n;a(i,i)=1.0_dp;end do
      call mregnn(x,y,a,result)
   end subroutine mregnn_positive
end module isotone_mregnn
