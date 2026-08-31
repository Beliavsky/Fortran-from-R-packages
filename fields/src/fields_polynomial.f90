! GPL-2.0-or-later. Polynomial/null-space utilities translated from fields.
module fields_polynomial
use fields_kinds, only: dp
use fields_native, only: dmaket
implicit none
private
public :: polynomial_basis, polynomial_power_table, evaluate_polynomial, &
          evaluate_multivariate_polynomial, polynomial_gradient, make_multi_index, choose_int

contains

pure integer function choose_int(n,k) result(v)
integer, intent(in) :: n,k
integer :: j,kk
if(k<0 .or. k>n) then; v=0; return; end if
kk=min(k,n-k); v=1
do j=1,kk; v=v*(n-kk+j)/j; end do
end function choose_int

subroutine polynomial_basis(x,m,t,ptab,info)
real(dp), intent(in) :: x(:,:)
integer, intent(in) :: m
real(dp), allocatable, intent(out) :: t(:,:)
integer, allocatable, intent(out), optional :: ptab(:,:)
integer, intent(out), optional :: info
integer :: n,d,npoly,ierr
integer, allocatable :: wp(:),pt(:,:)
if(m<0) error stop 'polynomial_basis: m must be nonnegative'
n=size(x,1); d=size(x,2)
if(m==0) then
   allocate(t(n,0)); if(present(ptab)) allocate(ptab(0,d)); if(present(info)) info=0; return
end if
npoly=choose_int(m+d-1,d)
allocate(t(n,npoly),pt(npoly,d),wp(d)); t=0.0_dp; pt=0; wp=0; ierr=0
call dmaket(m,n,d,x,n,npoly,t,n,wp,ierr,pt,npoly)
if(present(ptab)) then; allocate(ptab(npoly,d)); ptab=pt; end if
if(present(info)) info=ierr
end subroutine polynomial_basis

function polynomial_power_table(d,m) result(ptab)
integer, intent(in) :: d,m
integer, allocatable :: ptab(:,:)
real(dp), allocatable :: x(:,:),t(:,:)
integer :: info
allocate(x(1,d)); x=1.0_dp
call polynomial_basis(x,m,t,ptab,info)
if(info/=0) error stop 'polynomial_power_table: dmaket failed'
end function polynomial_power_table

function evaluate_polynomial(x,coef) result(y)
real(dp), intent(in) :: x(:),coef(:)
real(dp), allocatable :: y(:)
integer :: i,j
allocate(y(size(x))); y=0.0_dp
do j=size(coef),1,-1
   do i=1,size(x); y(i)=y(i)*x(i)+coef(j); end do
end do
end function evaluate_polynomial

function evaluate_multivariate_polynomial(x,coef,ptab) result(y)
real(dp), intent(in) :: x(:,:),coef(:)
integer, intent(in) :: ptab(:,:)
real(dp), allocatable :: y(:)
integer :: i,j,k
real(dp) :: term
if(size(ptab,1)/=size(coef) .or. size(ptab,2)/=size(x,2)) error stop 'evaluate_multivariate_polynomial: dimension mismatch'
allocate(y(size(x,1))); y=0.0_dp
do i=1,size(x,1)
   do j=1,size(coef)
      term=coef(j)
      do k=1,size(x,2)
         if(ptab(j,k)>0) term=term*x(i,k)**ptab(j,k)
      end do
      y(i)=y(i)+term
   end do
end do
end function evaluate_multivariate_polynomial

function polynomial_gradient(x,coef,ptab) result(g)
real(dp), intent(in) :: x(:,:),coef(:)
integer, intent(in) :: ptab(:,:)
real(dp), allocatable :: g(:,:)
integer :: i,j,k,l,pow
real(dp) :: term
if(size(ptab,1)/=size(coef) .or. size(ptab,2)/=size(x,2)) error stop 'polynomial_gradient: dimension mismatch'
allocate(g(size(x,1),size(x,2))); g=0.0_dp
do k=1,size(x,2)
   do i=1,size(x,1)
      do j=1,size(coef)
         if(ptab(j,k)==0) cycle
         term=coef(j)*real(ptab(j,k),dp)
         do l=1,size(x,2)
            pow=ptab(j,l); if(l==k) pow=pow-1
            if(pow>0) term=term*x(i,l)**pow
         end do
         g(i,k)=g(i,k)+term
      end do
   end do
end do
end function polynomial_gradient

function make_multi_index(m) result(idx)
integer, intent(in) :: m(:)
integer, allocatable :: idx(:,:)
integer :: l,n,i,k,rep1,rep2,block,pos,val
if(any(m<=0)) error stop 'make_multi_index: all dimensions must be positive'
l=size(m); n=product(m); allocate(idx(n,l))
do k=1,l
   if(k==1) then; rep1=1; else; rep1=product(m(:k-1)); end if
   if(k==l) then; rep2=1; else; rep2=product(m(k+1:)); end if
   block=rep1*m(k)
   do i=1,n
      pos=mod(i-1,block)
      val=pos/rep1+1
      idx(i,k)=val
   end do
end do
end function make_multi_index

end module fields_polynomial
