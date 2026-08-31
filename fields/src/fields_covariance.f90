! GPL-2.0-or-later. Translation of covariance kernels from fields.
module fields_covariance
use fields_kinds, only: dp
use fields_distance, only: fields_rdist, fields_rdist_earth
use fields_linalg, only: solve_general, identity_matrix
use spam_types, only: csr_matrix
use spam_distance, only: nearest_dist
use r_mod, only: besselK
implicit none
private
public :: exponential, matern, powered_exponential, gaussian_covariance, matern_cor_to_range, &
          stationary_covariance, stationary_earth_covariance, cubic_covariance, &
          radial_basis, radial_basis_constant, radial_covariance, &
          wendland, wendland_beta, wendland_covariance, wendland_covariance_sparse, &
          stationary_taper_sparse, paciorek_covariance, double_exponential

contains

pure elemental real(dp) function double_exponential(x) result(v)
real(dp), intent(in) :: x
v=0.5_dp*exp(-abs(x))
end function double_exponential

pure elemental real(dp) function exponential(d,a_range,phi) result(v)
real(dp), intent(in) :: d
real(dp), intent(in), optional :: a_range,phi
real(dp) :: ar,ph
ar=1.0_dp; if(present(a_range)) ar=a_range
ph=1.0_dp; if(present(phi)) ph=phi
if(d<0.0_dp .or. ar<=0.0_dp) then
   v=huge(1.0_dp)
else
   v=ph*exp(-d/ar)
end if
end function exponential

pure elemental real(dp) function powered_exponential(d,a_range,p,phi) result(v)
real(dp), intent(in) :: d
real(dp), intent(in), optional :: a_range,p,phi
real(dp) :: ar,pp,ph
ar=1.0_dp; if(present(a_range)) ar=a_range
pp=1.0_dp; if(present(p)) pp=p
ph=1.0_dp; if(present(phi)) ph=phi
if(d<0.0_dp .or. ar<=0.0_dp .or. pp<=0.0_dp) then
   v=huge(1.0_dp)
else
   v=ph*exp(-(d/ar)**pp)
end if
end function powered_exponential

pure elemental real(dp) function gaussian_covariance(d,a_range,phi) result(v)
real(dp), intent(in) :: d
real(dp), intent(in), optional :: a_range,phi
real(dp) :: ar,ph
ar=1.0_dp; if(present(a_range)) ar=a_range
ph=1.0_dp; if(present(phi)) ph=phi
v=ph*exp(-(d/ar)**2)
end function gaussian_covariance

pure elemental real(dp) function matern(d,a_range,smoothness,phi) result(v)
real(dp), intent(in) :: d
real(dp), intent(in), optional :: a_range,smoothness,phi
real(dp) :: ar,nu,ph,x,con
ar=1.0_dp; if(present(a_range)) ar=a_range
nu=0.5_dp; if(present(smoothness)) nu=smoothness
ph=1.0_dp; if(present(phi)) ph=phi
if(d<0.0_dp .or. ar<=0.0_dp .or. nu<=0.0_dp) then
   v=huge(1.0_dp); return
end if
x=d/ar
if(x<=1.0e-12_dp) then
   v=ph
else if(abs(nu-0.5_dp)<=1.0e-14_dp) then
   v=ph*exp(-x)
else if(abs(nu-1.5_dp)<=1.0e-14_dp) then
   v=ph*(1.0_dp+x)*exp(-x)
else if(abs(nu-2.5_dp)<=1.0e-14_dp) then
   v=ph*(1.0_dp+x+x*x/3.0_dp)*exp(-x)
else
   con=1.0_dp/(2.0_dp**(nu-1.0_dp)*gamma(nu))
   v=ph*con*x**nu*besselK(x,nu)
end if
end function matern

function stationary_covariance(x1,x2,model,a_range,smoothness,p,phi,vmat) result(k)
real(dp), intent(in) :: x1(:,:),x2(:,:)
character(len=*), intent(in), optional :: model
real(dp), intent(in), optional :: a_range,smoothness,p,phi,vmat(:,:)
real(dp), allocatable :: k(:,:)
real(dp), allocatable :: z1(:,:),z2(:,:),d(:,:),vinv(:,:)
real(dp) :: ar,nu,pp,ph
character(len=32) :: mod
integer :: i,j,info,nd
ar=1.0_dp; if(present(a_range)) ar=a_range
nu=0.5_dp; if(present(smoothness)) nu=smoothness
pp=1.0_dp; if(present(p)) pp=p
ph=1.0_dp; if(present(phi)) ph=phi
mod='exponential'; if(present(model)) mod=lower(trim(model))
if(size(x1,2)/=size(x2,2)) error stop 'stationary_covariance: coordinate dimensions differ'
nd=size(x1,2)
allocate(z1(size(x1,1),nd),z2(size(x2,1),nd)); z1=x1; z2=x2
if(present(vmat)) then
   if(size(vmat,1)/=nd .or. size(vmat,2)/=nd) error stop 'stationary_covariance: V dimension mismatch'
   if(abs(ar-1.0_dp)>epsilon(1.0_dp)) error stop 'stationary_covariance: specify V or a_range, not both'
   vinv=solve_general(vmat,identity_matrix(nd),info)
   if(info/=0) error stop 'stationary_covariance: singular V'
   z1=matmul(z1,transpose(vinv)); z2=matmul(z2,transpose(vinv))
end if
d=fields_rdist(z1,z2); allocate(k(size(d,1),size(d,2)))
select case(trim(mod))
case('exponential','exp')
   do j=1,size(d,2); do i=1,size(d,1); k(i,j)=exponential(d(i,j),ar,ph); end do; end do
case('gaussian','gauss')
   do j=1,size(d,2); do i=1,size(d,1); k(i,j)=gaussian_covariance(d(i,j),ar,ph); end do; end do
case('powered_exponential','powered-exponential','power','exp.cov')
   do j=1,size(d,2); do i=1,size(d,1); k(i,j)=powered_exponential(d(i,j),ar,pp,ph); end do; end do
case('matern')
   do j=1,size(d,2); do i=1,size(d,1); k(i,j)=matern(d(i,j),ar,nu,ph); end do; end do
case default
   error stop 'stationary_covariance: unsupported model'
end select
end function stationary_covariance

function stationary_earth_covariance(x1,x2,a_range,miles,phi) result(k)
real(dp), intent(in) :: x1(:,:),x2(:,:)
real(dp), intent(in), optional :: a_range,phi
logical, intent(in), optional :: miles
real(dp), allocatable :: k(:,:),d(:,:)
real(dp) :: ar,ph
integer :: i,j
ar=1.0_dp; if(present(a_range)) ar=a_range
ph=1.0_dp; if(present(phi)) ph=phi
d=fields_rdist_earth(x1,x2,miles); allocate(k(size(d,1),size(d,2)))
do j=1,size(d,2); do i=1,size(d,1); k(i,j)=ph*exp(-d(i,j)/ar); end do; end do
end function stationary_earth_covariance

pure elemental real(dp) function cubic_covariance(u,v) result(c)
real(dp), intent(in) :: u,v
if(u<v) then
   c=1.0_dp+v*u*u/2.0_dp-u**3/6.0_dp
else
   c=1.0_dp+u*v*v/2.0_dp-v**3/6.0_dp
end if
end function cubic_covariance

pure real(dp) function gamma_local(xin) result(g)
real(dp), intent(in) :: xin
real(dp) :: x,tmp
x=xin; tmp=1.0_dp
if(x<0.0_dp) then
   do while(x<0.0_dp)
      tmp=tmp*x; x=x+1.0_dp
   end do
   g=gamma(x)/tmp
else
   g=gamma(x)
end if
end function gamma_local

pure real(dp) function radial_basis_constant(m,d) result(a)
integer, intent(in) :: m,d
real(dp) :: pi
pi=acos(-1.0_dp)
if(mod(d,2)==0) then
   a=(-1.0_dp)**(1+m+d/2)*2.0_dp**(1-2*m)*pi**(-0.5_dp*d) / &
     (gamma(real(m,dp))*gamma_local(real(m-d/2+1,dp)))
else
   a=gamma_local(0.5_dp*d-m)*2.0_dp**(-2*m)*pi**(-0.5_dp*d)/gamma(real(m,dp))
end if
end function radial_basis_constant

pure elemental real(dp) function radial_basis(d,m,dimension,derivative,with_constant,with_log) result(v)
real(dp), intent(in) :: d
integer, intent(in) :: m,dimension
integer, intent(in), optional :: derivative
logical, intent(in), optional :: with_constant,with_log
integer :: der,powr
logical :: wc,wl
real(dp) :: con
powr=2*m-dimension
if(powr<=0) then; v=huge(1.0_dp); return; end if
der=0; if(present(derivative)) der=derivative
wc=.true.; if(present(with_constant)) wc=with_constant
wl=.true.; if(present(with_log)) wl=with_log
con=1.0_dp; if(wc) con=radial_basis_constant(m,dimension)
if(d<=1.0e-14_dp) then; v=0.0_dp; return; end if
if(der==0) then
   if(mod(dimension,2)==0 .and. wl) then; v=con*d**powr*log(d); else; v=con*d**powr; end if
else if(der==1) then
   if(mod(dimension,2)==0 .and. wl) then
      v=con*d**(powr-1)*(powr*log(d)+1.0_dp)
   else
      v=con*powr*d**(powr-1)
   end if
else
   ! Higher radial derivatives are not exposed by upstream RadialBasis.
   v=huge(1.0_dp)
end if
end function radial_basis

function radial_covariance(x1,x2,m,with_constant,with_log) result(k)
real(dp), intent(in) :: x1(:,:),x2(:,:)
integer, intent(in) :: m
logical, intent(in), optional :: with_constant,with_log
real(dp), allocatable :: k(:,:),d(:,:)
integer :: i,j,nd
logical :: wc,wl
nd=size(x1,2); if(size(x2,2)/=nd) error stop 'radial_covariance: coordinate dimensions differ'
wc=.true.; if(present(with_constant)) wc=with_constant
wl=.true.; if(present(with_log)) wl=with_log
d=fields_rdist(x1,x2); allocate(k(size(d,1),size(d,2)))
do j=1,size(d,2); do i=1,size(d,1)
   k(i,j)=radial_basis(d(i,j),m,nd,with_constant=wc,with_log=wl)
end do; end do
end function radial_covariance

pure real(dp) function falling(q,k) result(v)
integer, intent(in) :: q,k
integer :: j
if(k==0) then; v=1.0_dp; return; end if
if(k>q .and. q>=0) then; v=0.0_dp; return; end if
v=1.0_dp
do j=0,k-1; v=v*real(q-j,dp); end do
end function falling

pure real(dp) function rising_real(q,k) result(v)
real(dp), intent(in) :: q
integer, intent(in) :: k
integer :: j
v=1.0_dp
do j=0,k-1; v=v*(q+real(j,dp)); end do
end function rising_real

function wendland_beta(n,k) result(beta)
integer, intent(in) :: n,k
real(dp), allocatable :: beta(:,:)
integer :: l,col,row,m
real(dp) :: s
if(n<1 .or. k<0) error stop 'wendland_beta: invalid n or k'
l=n/2+k+1; allocate(beta(k+1,k+1)); beta=0.0_dp; beta(1,1)=1.0_dp
if(k==0) return
do col=0,k-1
   row=0; s=0.0_dp
   do m=0,col
      s=s+beta(m+1,col+1)*falling(m+1,m-row+1)/rising_real(real(l+2*col-m+1,dp),m-row+2)
   end do
   beta(row+1,col+2)=s
   do row=1,col+1
      s=0.0_dp
      do m=row-1,col
         s=s+beta(m+1,col+1)*falling(m+1,m-row+1)/rising_real(real(l+2*col-m+1,dp),m-row+2)
      end do
      beta(row+1,col+2)=s
   end do
end do
end function wendland_beta

pure real(dp) function binomial_int(n,k) result(v)
integer, intent(in) :: n,k
integer :: j,kk
if(k<0 .or. k>n) then; v=0.0_dp; return; end if
kk=min(k,n-k); v=1.0_dp
do j=1,kk; v=v*real(n-kk+j,dp)/real(j,dp); end do
end function binomial_int

pure real(dp) function monomial_product_derivative(r,m,p,der) result(v)
real(dp), intent(in) :: r
integer, intent(in) :: m,p,der
integer :: j
real(dp) :: term
v=0.0_dp
do j=0,der
   if(j>m .or. der-j>p) cycle
   term=binomial_int(der,j)*falling(m,j)*(-1.0_dp)**(der-j)*falling(p,der-j)
   if(m-j>0) term=term*r**(m-j)
   if(p-(der-j)>0) term=term*(1.0_dp-r)**(p-(der-j))
   v=v+term
end do
end function monomial_product_derivative

real(dp) function wendland_eval(r,n,k,derivative) result(phi)
real(dp), intent(in) :: r
integer, intent(in) :: n,k,derivative
real(dp), allocatable :: beta(:,:)
integer :: l,m,p
if(r<0.0_dp .or. r>1.0_dp) then; phi=0.0_dp; return; end if
beta=wendland_beta(n,k); l=n/2+k+1; phi=0.0_dp
do m=0,k
   p=l+2*k-m
   phi=phi+beta(m+1,k+1)*monomial_product_derivative(r,m,p,derivative)
end do
end function wendland_eval

real(dp) function wendland(d,a_range,dimension,k,derivative) result(v)
real(dp), intent(in) :: d
real(dp), intent(in), optional :: a_range
integer, intent(in) :: dimension,k
integer, intent(in), optional :: derivative
real(dp) :: ar,r,scale
integer :: der
ar=1.0_dp; if(present(a_range)) ar=a_range
der=0; if(present(derivative)) der=derivative
if(d<0.0_dp .or. ar<=0.0_dp) then; v=huge(1.0_dp); return; end if
r=d/ar
if(r>=1.0_dp) then; v=0.0_dp; return; end if
if(k==2 .and. dimension==2 .and. der==0) then
   v=(1.0_dp-r)**6*(35.0_dp*r*r+18.0_dp*r+3.0_dp)/3.0_dp
   return
end if
scale=wendland_eval(0.0_dp,dimension,k,0)*ar**der
v=wendland_eval(r,dimension,k,der)/scale
end function wendland

function wendland_covariance(x1,x2,a_range,korder) result(cov)
real(dp), intent(in) :: x1(:,:),x2(:,:),a_range
integer, intent(in), optional :: korder
real(dp), allocatable :: cov(:,:),d(:,:)
integer :: i,j,k,nd
k=2; if(present(korder)) k=korder
nd=size(x1,2); if(size(x2,2)/=nd) error stop 'wendland_covariance: coordinate dimensions differ'
d=fields_rdist(x1,x2); allocate(cov(size(d,1),size(d,2)))
do j=1,size(d,2); do i=1,size(d,1); cov(i,j)=wendland(d(i,j),a_range,nd,k); end do; end do
end function wendland_covariance

function wendland_covariance_sparse(x1,x2,a_range,korder) result(cov)
real(dp), intent(in) :: x1(:,:),x2(:,:),a_range
integer, intent(in), optional :: korder
type(csr_matrix) :: cov
integer :: k,idx,nd
k=2; if(present(korder)) k=korder
nd=size(x1,2)
cov=nearest_dist(x1,x2,method='euclidean',delta=a_range,full=.true.)
do idx=1,cov%nnz(); cov%entries(idx)=wendland(cov%entries(idx),a_range,nd,k); end do
end function wendland_covariance_sparse

function stationary_taper_sparse(x,a_range,taper_range,model,smoothness,p,korder) result(cov)
real(dp), intent(in) :: x(:,:),a_range,taper_range
character(len=*), intent(in), optional :: model
real(dp), intent(in), optional :: smoothness,p
integer, intent(in), optional :: korder
type(csr_matrix) :: cov
character(len=32) :: mod
real(dp) :: nu,pp,dv
integer :: k,idx,nd
mod='exponential'; if(present(model)) mod=lower(trim(model))
nu=0.5_dp; if(present(smoothness)) nu=smoothness
pp=1.0_dp; if(present(p)) pp=p
k=2; if(present(korder)) k=korder
nd=size(x,2); cov=nearest_dist(x,method='euclidean',delta=taper_range,full=.true.)
do idx=1,cov%nnz()
   dv=cov%entries(idx)
   select case(trim(mod))
   case('exponential','exp'); cov%entries(idx)=exponential(dv,a_range)
   case('gaussian','gauss'); cov%entries(idx)=gaussian_covariance(dv,a_range)
   case('matern'); cov%entries(idx)=matern(dv,a_range,nu)
   case('powered_exponential','power'); cov%entries(idx)=powered_exponential(dv,a_range,pp)
   case default; error stop 'stationary_taper_sparse: unsupported model'
   end select
   cov%entries(idx)=cov%entries(idx)*wendland(dv,taper_range,nd,k)
end do
end function stationary_taper_sparse

function paciorek_covariance(x1,x2,range1,range2,smoothness,rho1,rho2) result(cov)
real(dp), intent(in) :: x1(:,:),x2(:,:),range1(:),range2(:)
real(dp), intent(in), optional :: smoothness,rho1(:),rho2(:)
real(dp), allocatable :: cov(:,:),d(:,:)
real(dp) :: nu,a2,normc,s1,s2
integer :: i,j,dim
if(size(x1,1)/=size(range1) .or. size(x2,1)/=size(range2)) error stop 'paciorek_covariance: range length mismatch'
if(size(x1,2)/=size(x2,2)) error stop 'paciorek_covariance: coordinate dimensions differ'
nu=0.5_dp; if(present(smoothness)) nu=smoothness
d=fields_rdist(x1,x2); dim=size(x1,2); allocate(cov(size(x1,1),size(x2,1)))
do j=1,size(x2,1)
   do i=1,size(x1,1)
      a2=0.5_dp*(range1(i)**2+range2(j)**2)
      normc=(range1(i)*range2(j))**(0.5_dp*dim)/a2**(0.5_dp*dim)
      cov(i,j)=normc*matern(d(i,j)/sqrt(a2),smoothness=nu)
      s1=1.0_dp; if(present(rho1)) s1=sqrt(max(0.0_dp,rho1(i)))
      s2=1.0_dp; if(present(rho2)) s2=sqrt(max(0.0_dp,rho2(j)))
      cov(i,j)=cov(i,j)*s1*s2
   end do
end do
end function paciorek_covariance

pure function lower(s) result(t)
character(len=*), intent(in) :: s
character(len=len(s)) :: t
integer :: i,c
t=s
do i=1,len(s); c=iachar(t(i:i)); if(c>=65 .and. c<=90) t(i:i)=achar(c+32); end do
end function lower

real(dp) function matern_cor_to_range(d,nu,cor_target,tol,maxiter) result(a_range)
real(dp),intent(in)::d,nu
real(dp),intent(in),optional::cor_target,tol
integer,intent(in),optional::maxiter
real(dp)::target,eps,lo,hi,mid,fm
integer::it,nit
if(d<=0.0_dp .or. nu<=0.0_dp)error stop 'matern_cor_to_range: d and nu must be positive'
target=0.5_dp;if(present(cor_target))target=cor_target
if(target<=0.0_dp .or. target>=1.0_dp)error stop 'matern_cor_to_range: cor_target must be in (0,1)'
eps=1.0e-8_dp;if(present(tol))eps=tol
nit=100;if(present(maxiter))nit=maxiter
lo=-d/log(target);hi=lo
do it=1,100
   if(matern(d,a_range=hi,smoothness=nu)>=target)exit
   hi=2.0_dp*hi
end do
do it=1,100
   if(matern(d,a_range=lo,smoothness=nu)<=target)exit
   lo=0.5_dp*lo
end do
do it=1,nit
   mid=0.5_dp*(lo+hi);fm=matern(d,a_range=mid,smoothness=nu)-target
   if(abs(fm)<=eps)exit
   if(fm<0.0_dp)then;lo=mid;else;hi=mid;end if
end do
a_range=0.5_dp*(lo+hi)
end function matern_cor_to_range

end module fields_covariance
