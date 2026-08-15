module ld_linalg
use ld_kinds, only: dp
implicit none
private
public :: chol_lower, solve_lower, solve_upper, solve_spd, inverse_spd, logdet_spd
public :: cov_to_cor, symmetrize, make_positive_definite, outer_product, sample_covariance
public :: weighted_mean, weighted_covariance, trace_matrix
contains
subroutine chol_lower(a, l, info)
   real(dp), intent(in) :: a(:,:)
   real(dp), intent(out) :: l(:,:)
   integer, intent(out) :: info
   integer :: i, j, k, n
   real(dp) :: s
   n = size(a,1); l = 0.0_dp; info = 0
   if (size(a,2) /= n .or. size(l,1) /= n .or. size(l,2) /= n) then
      info = -1; return
   end if
   do i=1,n
      do j=1,i
         s = a(i,j)
         do k=1,j-1
            s = s - l(i,k)*l(j,k)
         end do
         if (i == j) then
            if (s <= 0.0_dp .or. s /= s) then
               info = i; return
            end if
            l(i,j) = sqrt(s)
         else
            l(i,j) = s/l(j,j)
         end if
      end do
   end do
end subroutine chol_lower

subroutine solve_lower(l, b, x)
   real(dp), intent(in) :: l(:,:), b(:)
   real(dp), intent(out) :: x(:)
   integer :: i, k, n
   n=size(b)
   do i=1,n
      x(i)=b(i)
      do k=1,i-1
         x(i)=x(i)-l(i,k)*x(k)
      end do
      x(i)=x(i)/l(i,i)
   end do
end subroutine solve_lower

subroutine solve_upper(u, b, x)
   real(dp), intent(in) :: u(:,:), b(:)
   real(dp), intent(out) :: x(:)
   integer :: i, k, n
   n=size(b)
   do i=n,1,-1
      x(i)=b(i)
      do k=i+1,n
         x(i)=x(i)-u(i,k)*x(k)
      end do
      x(i)=x(i)/u(i,i)
   end do
end subroutine solve_upper

subroutine solve_spd(a,b,x,info)
   real(dp), intent(in) :: a(:,:), b(:)
   real(dp), intent(out) :: x(:)
   integer, intent(out) :: info
   real(dp), allocatable :: l(:,:), y(:)
   integer :: n
   n=size(b); allocate(l(n,n),y(n))
   call chol_lower(a,l,info)
   if(info/=0) then; x=0.0_dp; return; end if
   call solve_lower(l,b,y)
   call solve_upper(transpose(l),y,x)
end subroutine solve_spd

subroutine inverse_spd(a, ainv, info)
   real(dp), intent(in) :: a(:,:)
   real(dp), intent(out) :: ainv(:,:)
   integer, intent(out) :: info
   integer :: i,n
   real(dp), allocatable :: e(:),x(:)
   n=size(a,1); allocate(e(n),x(n)); ainv=0.0_dp
   do i=1,n
      e=0.0_dp; e(i)=1.0_dp
      call solve_spd(a,e,x,info)
      if(info/=0) return
      ainv(:,i)=x
   end do
   ainv=0.5_dp*(ainv+transpose(ainv))
end subroutine inverse_spd

function logdet_spd(a, info) result(v)
   real(dp), intent(in) :: a(:,:)
   integer, intent(out), optional :: info
   real(dp) :: v
   real(dp), allocatable :: l(:,:)
   integer :: i,ifail,n
   n=size(a,1); allocate(l(n,n)); call chol_lower(a,l,ifail)
   if(present(info)) info=ifail
   if(ifail/=0) then; v=-huge(1.0_dp); return; end if
   v=0.0_dp
   do i=1,n; v=v+2.0_dp*log(l(i,i)); end do
end function logdet_spd

subroutine cov_to_cor(s, r, sd)
   real(dp), intent(in) :: s(:,:)
   real(dp), intent(out) :: r(:,:)
   real(dp), intent(out), optional :: sd(:)
   integer :: i,j,n
   real(dp), allocatable :: d(:)
   n=size(s,1); allocate(d(n))
   do i=1,n; d(i)=sqrt(max(s(i,i),0.0_dp)); end do
   do i=1,n
      do j=1,n
         if(d(i)>0.0_dp .and. d(j)>0.0_dp) then
            r(i,j)=s(i,j)/(d(i)*d(j))
         else
            r(i,j)=0.0_dp
         end if
      end do
   end do
   if(present(sd)) sd=d
end subroutine cov_to_cor

subroutine symmetrize(a)
   real(dp), intent(inout) :: a(:,:)
   a=0.5_dp*(a+transpose(a))
end subroutine symmetrize

subroutine make_positive_definite(a, jitter_used, info)
   real(dp), intent(inout) :: a(:,:)
   real(dp), intent(out), optional :: jitter_used
   integer, intent(out), optional :: info
   real(dp), allocatable :: l(:,:)
   real(dp) :: jitter, scale
   integer :: k,ifail,n,i
   n=size(a,1); allocate(l(n,n)); call symmetrize(a)
   scale=max(1.0_dp,maxval(abs([(a(i,i),i=1,n)])))
   jitter=0.0_dp
   do k=0,12
      call chol_lower(a,l,ifail)
      if(ifail==0) exit
      if(k==0) then; jitter=1.0e-12_dp*scale; else; jitter=jitter*10.0_dp; end if
      do i=1,n; a(i,i)=a(i,i)+jitter; end do
   end do
   if(present(jitter_used)) jitter_used=jitter
   if(present(info)) info=ifail
end subroutine make_positive_definite

pure function outer_product(x,y) result(a)
   real(dp), intent(in) :: x(:),y(:)
   real(dp) :: a(size(x),size(y))
   integer :: i,j
   do i=1,size(x); do j=1,size(y); a(i,j)=x(i)*y(j); end do; end do
end function outer_product

subroutine sample_covariance(x,cov,mean)
   real(dp), intent(in) :: x(:,:)
   real(dp), intent(out) :: cov(:,:)
   real(dp), intent(out), optional :: mean(:)
   integer :: n,p,i
   real(dp), allocatable :: m(:),d(:)
   n=size(x,1); p=size(x,2); allocate(m(p),d(p))
   m=sum(x,dim=1)/real(n,dp); cov=0.0_dp
   do i=1,n
      d=x(i,:)-m; cov=cov+outer_product(d,d)
   end do
   if(n>1) cov=cov/real(n-1,dp)
   if(present(mean)) mean=m
end subroutine sample_covariance

function weighted_mean(x,w) result(m)
   real(dp), intent(in) :: x(:,:),w(:)
   real(dp) :: m(size(x,2)), sw
   integer :: i
   sw=sum(w); m=0.0_dp
   if(sw<=0.0_dp) return
   do i=1,size(x,1); m=m+w(i)*x(i,:); end do
   m=m/sw
end function weighted_mean

function weighted_covariance(x,w,m) result(c)
   real(dp), intent(in) :: x(:,:),w(:),m(:)
   real(dp) :: c(size(x,2),size(x,2)),sw
   real(dp) :: d(size(x,2))
   integer :: i
   c=0.0_dp; sw=sum(w)
   if(sw<=0.0_dp) return
   do i=1,size(x,1); d=x(i,:)-m; c=c+w(i)*outer_product(d,d); end do
   c=c/sw
end function weighted_covariance

pure function trace_matrix(a) result(v)
   real(dp), intent(in) :: a(:,:)
   real(dp) :: v
   integer :: i
   v=0.0_dp; do i=1,min(size(a,1),size(a,2)); v=v+a(i,i); end do
end function trace_matrix
end module ld_linalg
