! GPL-2.0-or-later. Grid/interpolation computational utilities from fields.
module fields_grid
use fields_kinds, only: dp
use fields_native, only: multwendlandg
use fields_linalg, only: inverse_spd
use fields_covariance, only: exponential, matern, powered_exponential, gaussian_covariance, double_exponential
use spam_types, only: csr_matrix
use spam_csr, only: csr_from_triplet
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
implicit none
private
public :: bilinear_interpolate, find_grid_box, make_surface_grid, image_smooth_direct, mult_wendland_grid
public :: offgrid_weights_1d, offgrid_weights_2d, offgrid_weight_result

type :: offgrid_weight_result
   type(csr_matrix) :: b
   real(dp), allocatable :: prediction_variance(:)
end type

contains

function make_surface_grid(xgrid,ygrid) result(x)
real(dp), intent(in) :: xgrid(:),ygrid(:)
real(dp), allocatable :: x(:,:)
integer :: i,j,k
allocate(x(size(xgrid)*size(ygrid),2)); k=0
do j=1,size(ygrid)
   do i=1,size(xgrid)
      k=k+1; x(k,1)=xgrid(i); x(k,2)=ygrid(j)
   end do
end do
end function make_surface_grid

function bilinear_interpolate(xgrid,ygrid,z,loc) result(v)
real(dp), intent(in) :: xgrid(:),ygrid(:),z(:,:),loc(:,:)
real(dp), allocatable :: v(:)
real(dp) :: ex,ey,nanv
integer :: q,i,j
if(size(z,1)/=size(xgrid) .or. size(z,2)/=size(ygrid) .or. size(loc,2)/=2) error stop 'bilinear_interpolate: dimension mismatch'
if(size(xgrid)<2 .or. size(ygrid)<2) error stop 'bilinear_interpolate: grid too small'
nanv=ieee_value(0.0_dp,ieee_quiet_nan); allocate(v(size(loc,1)))
do q=1,size(loc,1)
   i=lower_cell(xgrid,loc(q,1)); j=lower_cell(ygrid,loc(q,2))
   if(i==0 .or. j==0) then; v(q)=nanv; cycle; end if
   ex=(loc(q,1)-xgrid(i))/(xgrid(i+1)-xgrid(i)); ey=(loc(q,2)-ygrid(j))/(ygrid(j+1)-ygrid(j))
   v(q)=z(i,j)*(1.0_dp-ex)*(1.0_dp-ey)+z(i+1,j)*ex*(1.0_dp-ey)+ &
        z(i,j+1)*(1.0_dp-ex)*ey+z(i+1,j+1)*ex*ey
end do
end function bilinear_interpolate

subroutine find_grid_box(xgrid,ygrid,s,index,ij)
real(dp), intent(in) :: xgrid(:),ygrid(:),s(:,:)
integer, allocatable, intent(out) :: index(:),ij(:,:)
real(dp) :: dx,dy,x0,y0
integer :: q,i,j,m,n
if(size(s,2)/=2 .or. size(xgrid)<2 .or. size(ygrid)<2) error stop 'find_grid_box: invalid input'
m=size(xgrid); n=size(ygrid); dx=xgrid(2)-xgrid(1); dy=ygrid(2)-ygrid(1)
x0=xgrid(1)-0.5_dp*dx; y0=ygrid(1)-0.5_dp*dy
allocate(index(size(s,1)),ij(size(s,1),2)); index=0; ij=0
do q=1,size(s,1)
   i=int((s(q,1)-x0)/dx)+1; j=int((s(q,2)-y0)/dy)+1
   if(i>=1 .and. i<=m .and. j>=1 .and. j<=n) then
      ij(q,:)=[i,j]; index(q)=i+(j-1)*m
   end if
end do
end subroutine find_grid_box

function image_smooth_direct(y,dx,dy,a_range,kernel,weights,tol) result(out)
real(dp), intent(in) :: y(:,:)
real(dp), intent(in), optional :: dx,dy,a_range,weights(:,:),tol
character(len=*), intent(in), optional :: kernel
real(dp), allocatable :: out(:,:)
real(dp) :: ddx,ddy,ar,eps,dist2,w,sumw,sumy,obsw,nanv
character(len=32) :: ker
integer :: i,j,ii,jj
if(present(weights)) then
   if(any(shape(weights)/=shape(y))) error stop 'image_smooth_direct: weights shape mismatch'
end if
ddx=1.0_dp; if(present(dx)) ddx=dx
ddy=1.0_dp; if(present(dy)) ddy=dy
ar=1.0_dp; if(present(a_range)) ar=a_range
eps=1.0e-8_dp; if(present(tol)) eps=tol
ker='double_exponential'; if(present(kernel)) ker=lower(trim(kernel))
nanv=ieee_value(0.0_dp,ieee_quiet_nan); allocate(out(size(y,1),size(y,2)))
do j=1,size(y,2)
   do i=1,size(y,1)
      sumw=0.0_dp; sumy=0.0_dp
      do jj=1,size(y,2)
         do ii=1,size(y,1)
            if(ieee_is_nan(y(ii,jj))) cycle
            dist2=((real(ii-i,dp)*ddx)/ar)**2+((real(jj-j,dp)*ddy)/ar)**2
            select case(trim(ker))
            case('gaussian','gauss'); w=exp(-dist2)
            case('exponential','exp'); w=exp(-sqrt(dist2))
            case default; w=double_exponential(dist2)
            end select
            obsw=1.0_dp; if(present(weights)) obsw=weights(ii,jj)
            sumw=sumw+w*obsw; sumy=sumy+w*obsw*y(ii,jj)
         end do
      end do
      if(sumw>eps) then; out(i,j)=sumy/sumw; else; out(i,j)=nanv; end if
   end do
end do
end function image_smooth_direct

function offgrid_weights_1d(s,xgrid,a_range,sigma2,nn_size,model,smoothness,power) result(out)
real(dp), intent(in) :: s(:),xgrid(:),a_range,sigma2
integer, intent(in), optional :: nn_size
character(len=*), intent(in), optional :: model
real(dp), intent(in), optional :: smoothness,power
type(offgrid_weight_result) :: out
integer :: np,nnb,m,q,j,k,start,idx,nnz
integer, allocatable :: ir(:),jc(:)
real(dp), allocatable :: val(:),coords(:),sigma11(:,:),sigma11i(:,:),cross(:)
real(dp) :: dx,nu,pp
character(len=32) :: mod
np=2; if(present(nn_size)) np=nn_size
nnb=2*np; m=size(xgrid); dx=xgrid(2)-xgrid(1)
nu=0.5_dp; if(present(smoothness)) nu=smoothness
pp=1.0_dp; if(present(power)) pp=power
mod='matern'; if(present(model)) mod=lower(trim(model))
allocate(coords(nnb)); do j=1,nnb; coords(j)=real(j-np,dp)*dx; end do
allocate(sigma11(nnb,nnb))
do j=1,nnb; do k=1,nnb; sigma11(j,k)=sigma2*cov_value(abs(coords(j)-coords(k)),a_range,mod,nu,pp); end do; end do
sigma11i=inverse_spd(sigma11,idx); if(idx/=0) error stop 'offgrid_weights_1d: local covariance singular'
nnz=size(s)*nnb; allocate(ir(nnz),jc(nnz),val(nnz),out%prediction_variance(size(s)),cross(nnb)); k=0
do q=1,size(s)
   start=int((s(q)-xgrid(1))/dx)+1-(np-1)
   if(start<1 .or. start+nnb-1>m) error stop 'offgrid_weights_1d: grid does not extend far enough'
   do j=1,nnb; cross(j)=sigma2*cov_value(abs(xgrid(start+j-1)-s(q)),a_range,mod,nu,pp); end do
   val(k+1:k+nnb)=matmul(cross,sigma11i)
   do j=1,nnb; k=k+1; ir(k)=q; jc(k)=start+j-1; end do
   out%prediction_variance(q)=max(0.0_dp,sigma2-dot_product(cross,matmul(sigma11i,cross)))
end do
out%b=csr_from_triplet(size(s),m,ir,jc,val)
end function offgrid_weights_1d

function offgrid_weights_2d(s,xgrid,ygrid,a_range,sigma2,nn_size,model,smoothness,power) result(out)
real(dp), intent(in) :: s(:,:),xgrid(:),ygrid(:),a_range,sigma2
integer, intent(in), optional :: nn_size
character(len=*), intent(in), optional :: model
real(dp), intent(in), optional :: smoothness,power
type(offgrid_weight_result) :: out
integer :: np,nside,nnb,m,n,q,a,b,c,d,k,ix0,iy0,ii,jj,info,nnz
integer, allocatable :: ir(:),jc(:),gx(:),gy(:)
real(dp), allocatable :: val(:),sigma11(:,:),sigma11i(:,:),cross(:),cx(:),cy(:)
real(dp) :: dx,dy,nu,pp,dist
character(len=32) :: mod
if(size(s,2)/=2) error stop 'offgrid_weights_2d: s must have two columns'
np=2; if(present(nn_size)) np=nn_size
nside=2*np; nnb=nside*nside; m=size(xgrid); n=size(ygrid); dx=xgrid(2)-xgrid(1); dy=ygrid(2)-ygrid(1)
nu=0.5_dp; if(present(smoothness)) nu=smoothness
pp=1.0_dp; if(present(power)) pp=power
mod='matern'; if(present(model)) mod=lower(trim(model))
allocate(gx(nnb),gy(nnb),cx(nnb),cy(nnb)); k=0
do b=0,nside-1; do a=0,nside-1; k=k+1; gx(k)=a; gy(k)=b; cx(k)=real(a,dp)*dx; cy(k)=real(b,dp)*dy; end do; end do
allocate(sigma11(nnb,nnb))
do b=1,nnb; do a=1,nnb
   dist=sqrt((cx(a)-cx(b))**2+(cy(a)-cy(b))**2); sigma11(a,b)=sigma2*cov_value(dist,a_range,mod,nu,pp)
end do; end do
sigma11i=inverse_spd(sigma11,info); if(info/=0) error stop 'offgrid_weights_2d: local covariance singular'
nnz=size(s,1)*nnb; allocate(ir(nnz),jc(nnz),val(nnz),cross(nnb),out%prediction_variance(size(s,1))); k=0
do q=1,size(s,1)
   ix0=int((s(q,1)-xgrid(1))/dx)+1-(np-1); iy0=int((s(q,2)-ygrid(1))/dy)+1-(np-1)
   if(ix0<1 .or. iy0<1 .or. ix0+nside-1>m .or. iy0+nside-1>n) error stop 'offgrid_weights_2d: grid too small around point'
   do a=1,nnb
      ii=ix0+gx(a); jj=iy0+gy(a); dist=sqrt((xgrid(ii)-s(q,1))**2+(ygrid(jj)-s(q,2))**2)
      cross(a)=sigma2*cov_value(dist,a_range,mod,nu,pp)
   end do
   val(k+1:k+nnb)=matmul(cross,sigma11i)
   do a=1,nnb
      k=k+1; ir(k)=q; ii=ix0+gx(a); jj=iy0+gy(a); jc(k)=ii+(jj-1)*m
   end do
   out%prediction_variance(q)=max(0.0_dp,sigma2-dot_product(cross,matmul(sigma11i,cross)))
end do
out%b=csr_from_triplet(size(s,1),m*n,ir,jc,val)
end function offgrid_weights_2d

real(dp) function cov_value(d,ar,model,nu,p) result(v)
real(dp), intent(in) :: d,ar,nu,p
character(len=*), intent(in) :: model
select case(trim(model))
case('exponential','exp'); v=exponential(d,ar)
case('gaussian','gauss'); v=gaussian_covariance(d,ar)
case('powered_exponential','power'); v=powered_exponential(d,ar,p)
case default; v=matern(d,ar,nu)
end select
end function cov_value

integer function lower_cell(g,x) result(i)
real(dp), intent(in) :: g(:),x
integer :: lo,hi,mid,n
n=size(g); if(x<g(1) .or. x>g(n)) then; i=0; return; end if
if(x==g(n)) then; i=n-1; return; end if
lo=1; hi=n-1
do while(lo<=hi)
   mid=(lo+hi)/2
   if(x<g(mid)) then; hi=mid-1
   else if(x>=g(mid+1)) then; lo=mid+1
   else; i=mid; return; end if
end do
i=0
end function lower_cell

function mult_wendland_grid(xgrid,ygrid,center,delta,coef) result(h)
real(dp),intent(in)::xgrid(:),ygrid(:),center(:,:),delta,coef(:)
real(dp),allocatable::h(:,:),cs(:,:)
real(dp)::dx,dy,delta_x,delta_y
integer::mx,my,nc,flag
mx=size(xgrid);my=size(ygrid);nc=size(center,1)
if(mx<2 .or. my<2 .or. size(center,2)/=2 .or. size(coef)/=nc .or. delta<=0.0_dp) &
 error stop 'mult_wendland_grid: invalid input'
dx=(xgrid(mx)-xgrid(1))/real(mx-1,dp);dy=(ygrid(my)-ygrid(1))/real(my-1,dp)
if(dx<=0.0_dp .or. dy<=0.0_dp)error stop 'mult_wendland_grid: grids must increase'
allocate(cs(nc,2),h(mx,my));h=0.0_dp
cs(:,1)=(center(:,1)-xgrid(1))/dx+1.0_dp
cs(:,2)=(center(:,2)-ygrid(1))/dy+1.0_dp
delta_x=delta/dx;delta_y=delta/dy;flag=1
call multwendlandg(mx,my,delta_x,delta_y,nc,cs,coef,h,flag)
if(flag/=0)error stop 'mult_wendland_grid: native kernel failed'
end function mult_wendland_grid

pure function lower(s) result(t)
character(len=*),intent(in)::s
character(len=len(s))::t
integer::i,c
t=s; do i=1,len(s); c=iachar(t(i:i)); if(c>=65.and.c<=90)t(i:i)=achar(c+32); end do
end function lower

end module fields_grid
