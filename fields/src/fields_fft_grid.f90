! GPL-2.0-or-later. FFT/circulant gridded computations corresponding to fields.
module fields_fft_grid
use fields_kinds, only: dp
use fields_covariance, only: exponential, matern, powered_exponential, gaussian_covariance, double_exponential
use r_mod, only: rnorm1
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
implicit none
private
public :: circulant_embedding_2d, circulant_setup_2d, circulant_sample_2d, image_smooth_fft
public :: fft_interp_surface, fft_interp_result

type :: circulant_embedding_2d
   integer :: nx=0,ny=0,mx=0,my=0
   real(dp) :: dx=1.0_dp,dy=1.0_dp
   complex(dp), allocatable :: sqrt_spectrum(:,:)
end type

type :: fft_interp_result
   real(dp),allocatable::x(:),y(:),z(:,:)
end type fft_interp_result

contains

subroutine circulant_setup_2d(nx,ny,dx,dy,embedding,model,a_range,smoothness,power,phi,info)
integer,intent(in)::nx,ny
real(dp),intent(in)::dx,dy
type(circulant_embedding_2d),intent(out)::embedding
character(len=*),intent(in),optional::model
real(dp),intent(in),optional::a_range,smoothness,power,phi
integer,intent(out),optional::info
complex(dp),allocatable::spec(:,:)
integer::mx,my,i,j,ii,jj,istat,tries
real(dp)::ar,nu,pp,ph,h,tol
character(len=32)::mod
if(nx<1 .or. ny<1 .or. dx<=0.0_dp .or. dy<=0.0_dp) error stop 'circulant_setup_2d: invalid grid'
ar=1.0_dp;if(present(a_range))ar=a_range
nu=0.5_dp;if(present(smoothness))nu=smoothness
pp=1.0_dp;if(present(power))pp=power
ph=1.0_dp;if(present(phi))ph=phi
mod='exponential';if(present(model))mod=lower(trim(model))
mx=next_power_two(max(2,2*(nx-1))); my=next_power_two(max(2,2*(ny-1))); istat=1
do tries=1,8
   allocate(spec(mx,my))
   do j=1,my
      jj=j-1;if(jj>my/2)jj=jj-my
      do i=1,mx
         ii=i-1;if(ii>mx/2)ii=ii-mx
         h=hypot(dx*real(ii,dp),dy*real(jj,dp))
         spec(i,j)=cmplx(cov_scalar(h,mod,ar,nu,pp,ph),0.0_dp,dp)
      end do
   end do
   call fft2(spec,.false.)
   tol=1.0e-10_dp*max(1.0_dp,maxval(abs(real(spec,dp))))
   if(maxval(abs(aimag(spec)))<=100.0_dp*tol .and. minval(real(spec,dp))>=-tol)then
      istat=0; exit
   end if
   deallocate(spec); mx=2*mx; my=2*my
end do
if(istat==0)then
   embedding%nx=nx;embedding%ny=ny;embedding%mx=mx;embedding%my=my;embedding%dx=dx;embedding%dy=dy
   allocate(embedding%sqrt_spectrum(mx,my));embedding%sqrt_spectrum=cmplx(sqrt(max(real(spec,dp),0.0_dp)),0.0_dp,dp)
   deallocate(spec)
end if
if(present(info))info=istat
end subroutine circulant_setup_2d

function circulant_sample_2d(embedding,nsim,info) result(draws)
type(circulant_embedding_2d),intent(in)::embedding
integer,intent(in),optional::nsim
integer,intent(out),optional::info
real(dp),allocatable::draws(:,:,:)
complex(dp),allocatable::work(:,:)
integer::n,k,i,j
n=1;if(present(nsim))n=nsim
if(.not.allocated(embedding%sqrt_spectrum))then
   allocate(draws(0,0,0));if(present(info))info=1;return
end if
allocate(draws(embedding%nx,embedding%ny,n),work(embedding%mx,embedding%my))
do k=1,n
   do j=1,embedding%my;do i=1,embedding%mx;work(i,j)=cmplx(rnorm1(),0.0_dp,dp);end do;end do
   call fft2(work,.false.);work=work*embedding%sqrt_spectrum;call fft2(work,.true.)
   draws(:,:,k)=real(work(:embedding%nx,:embedding%ny),dp)
end do
if(present(info))info=0
end function circulant_sample_2d

function image_smooth_fft(y,dx,dy,a_range,kernel,weights,tol) result(out)
real(dp),intent(in)::y(:,:)
real(dp),intent(in),optional::dx,dy,a_range,weights(:,:),tol
character(len=*),intent(in),optional::kernel
real(dp),allocatable::out(:,:)
complex(dp),allocatable::num(:,:),den(:,:),kerf(:,:)
real(dp)::ddx,ddy,ar,eps,h2,wobs,nanv
integer::nx,ny,mx,my,i,j,ii,jj
character(len=32)::ker
nx=size(y,1);ny=size(y,2);ddx=1.0_dp;if(present(dx))ddx=dx;ddy=1.0_dp;if(present(dy))ddy=dy
ar=1.0_dp;if(present(a_range))ar=a_range;eps=1.0e-8_dp;if(present(tol))eps=tol
ker='double_exponential';if(present(kernel))ker=lower(trim(kernel))
if(present(weights))then;if(any(shape(weights)/=shape(y)))error stop 'image_smooth_fft: weights shape mismatch';end if
mx=next_power_two(2*nx);my=next_power_two(2*ny);allocate(num(mx,my),den(mx,my),kerf(mx,my));num=0;den=0
! Data/missingness numerator and denominator.
do j=1,ny;do i=1,nx
   if(.not.ieee_is_nan(y(i,j)))then
      wobs=1.0_dp;if(present(weights))wobs=weights(i,j)
      num(i,j)=cmplx(y(i,j)*wobs,0.0_dp,dp);den(i,j)=cmplx(wobs,0.0_dp,dp)
   end if
end do;end do
! Periodic kernel with enough zero padding to make the retained convolution linear.
do j=1,my
   jj=j-1;if(jj>my/2)jj=jj-my
   do i=1,mx
      ii=i-1;if(ii>mx/2)ii=ii-mx
      h2=((ddx*real(ii,dp))/ar)**2+((ddy*real(jj,dp))/ar)**2
      select case(trim(ker))
      case('gaussian','gauss');kerf(i,j)=cmplx(exp(-h2),0.0_dp,dp)
      case('exponential','exp');kerf(i,j)=cmplx(exp(-sqrt(h2)),0.0_dp,dp)
      case default;kerf(i,j)=cmplx(double_exponential(h2),0.0_dp,dp)
      end select
   end do
end do
call fft2(num,.false.);call fft2(den,.false.);call fft2(kerf,.false.)
num=num*kerf;den=den*kerf;call fft2(num,.true.);call fft2(den,.true.)
nanv=ieee_value(0.0_dp,ieee_quiet_nan);allocate(out(nx,ny))
do j=1,ny;do i=1,nx
   if(real(den(i,j),dp)>eps)then;out(i,j)=real(num(i,j),dp)/real(den(i,j),dp);else;out(i,j)=nanv;end if
end do;end do
end function image_smooth_fft

real(dp) function cov_scalar(d,model,ar,nu,p,ph) result(v)
real(dp),intent(in)::d,ar,nu,p,ph
character(len=*),intent(in)::model
select case(trim(model))
case('matern');v=matern(d,ar,nu,ph)
case('gaussian','gauss');v=gaussian_covariance(d,ar,ph)
case('powered_exponential','power');v=powered_exponential(d,ar,p,ph)
case default;v=exponential(d,ar,ph)
end select
end function cov_scalar

pure integer function next_power_two(n) result(m)
integer,intent(in)::n
m=1;do while(m<n);m=2*m;end do
end function next_power_two

subroutine fft2(a,inverse)
complex(dp),intent(inout)::a(:,:)
logical,intent(in)::inverse
complex(dp),allocatable::v(:)
integer::i,j,n1,n2
n1=size(a,1);n2=size(a,2);allocate(v(max(n1,n2)))
do j=1,n2;v(1:n1)=a(:,j);call fft1(v(1:n1),inverse);a(:,j)=v(1:n1);end do
do i=1,n1;v(1:n2)=a(i,:);call fft1(v(1:n2),inverse);a(i,:)=v(1:n2);end do
end subroutine fft2

subroutine fft1(a,inverse)
complex(dp),intent(inout)::a(:)
logical,intent(in)::inverse
integer::n,i,j,k,m,half
complex(dp)::tmp,w,wm,u,t
real(dp)::sgn,pi
n=size(a);if(iand(n,n-1)/=0)error stop 'fft1: length must be power of two';j=1
do i=1,n
   if(i<j)then;tmp=a(i);a(i)=a(j);a(j)=tmp;end if
   k=n/2;do while(k>=1 .and. j>k);j=j-k;k=k/2;end do;j=j+k
end do
pi=acos(-1.0_dp);m=2;sgn=-1.0_dp;if(inverse)sgn=1.0_dp
do while(m<=n)
   half=m/2;wm=exp(cmplx(0.0_dp,sgn*2.0_dp*pi/real(m,dp),dp))
   do k=1,n,m
      w=cmplx(1.0_dp,0.0_dp,dp)
      do j=0,half-1
         u=a(k+j);t=w*a(k+j+half);a(k+j)=u+t;a(k+j+half)=u-t;w=w*wm
      end do
   end do
   m=2*m
end do
if(inverse)a=a/real(n,dp)
end subroutine fft1

pure function lower(s) result(t)
character(len=*),intent(in)::s
character(len=len(s))::t
integer::i,c
t=s;do i=1,len(s);c=iachar(t(i:i));if(c>=65.and.c<=90)t(i:i)=achar(c+32);end do
end function lower

function fft_interp_surface(x,y,z,factor) result(out)
real(dp),intent(in)::x(:),y(:),z(:,:)
integer,intent(in)::factor
type(fft_interp_result)::out
complex(dp),allocatable::spec(:,:)
complex(dp)::acc,phase
real(dp)::twopi,angle
integer::m1,m2,n1,n2,o1,o2,i,j,k,l,pk,pl
m1=size(z,1);m2=size(z,2)
if(size(x)/=m1 .or. size(y)/=m2 .or. m1<1 .or. m2<1)error stop 'fft_interp_surface: dimensions'
if(mod(m1,2)/=1 .or. mod(m2,2)/=1)error stop 'fft_interp_surface: source grid dimensions must be odd'
if(factor<1)error stop 'fft_interp_surface: factor must be positive'
n1=factor*m1;n2=factor*m2;o1=n1-factor+1;o2=n2-factor+1
twopi=2.0_dp*acos(-1.0_dp);allocate(spec(m1,m2));spec=(0.0_dp,0.0_dp)
! R's fft is unnormalized in both directions. Build the source spectrum exactly.
do k=0,m1-1
 do l=0,m2-1
  acc=(0.0_dp,0.0_dp)
  do i=0,m1-1
   do j=0,m2-1
    angle=-twopi*(real(k*i,dp)/real(m1,dp)+real(l*j,dp)/real(m2,dp))
    phase=cmplx(cos(angle),sin(angle),kind=dp);acc=acc+z(i+1,j+1)*phase
   end do
  end do
  spec(k+1,l+1)=acc
 end do
end do
allocate(out%z(o1,o2),out%x(o1),out%y(o2))
do i=0,o1-1
 do j=0,o2-1
  acc=(0.0_dp,0.0_dp)
  do k=0,m1-1
   pk=k;if(k>m1/2)pk=k-m1
   do l=0,m2-1
    pl=l;if(l>m2/2)pl=l-m2
    angle=twopi*(real(pk*i,dp)/real(n1,dp)+real(pl*j,dp)/real(n2,dp))
    acc=acc+spec(k+1,l+1)*cmplx(cos(angle),sin(angle),kind=dp)
   end do
  end do
  out%z(i+1,j+1)=real(acc,dp)/real(m1*m2,dp)
 end do
end do
if(o1==1)then;out%x=x(1);else;do i=1,o1;out%x(i)=x(1)+real(i-1,dp)*(x(m1)-x(1))/real(o1-1,dp);end do;end if
if(o2==1)then;out%y=y(1);else;do j=1,o2;out%y(j)=y(1)+real(j-1,dp)*(y(m2)-y(1))/real(o2-1,dp);end do;end if
end function fft_interp_surface

end module fields_fft_grid
