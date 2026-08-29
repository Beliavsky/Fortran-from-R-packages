module spam_covariance
use spam_kinds, only: dp
use spam_types, only: csr_matrix
use r_compat, only: besselK
implicit none
private
public :: cov_exp, cov_sph, cor_sph, cov_nug, cov_wu1, cov_wu2, cov_wu3
public :: cov_wend1, cov_wend2, cov_mat, cov_finnmat, cov_mat12, cov_mat32, cov_mat52

interface cov_exp
   module procedure cov_exp_vec, cov_exp_csr
end interface
interface cov_sph
   module procedure cov_sph_vec, cov_sph_csr
end interface
interface cor_sph
   module procedure cor_sph_vec, cor_sph_csr
end interface
interface cov_nug
   module procedure cov_nug_vec, cov_nug_csr
end interface
interface cov_wu1
   module procedure cov_wu1_vec, cov_wu1_csr
end interface
interface cov_wu2
   module procedure cov_wu2_vec, cov_wu2_csr
end interface
interface cov_wu3
   module procedure cov_wu3_vec, cov_wu3_csr
end interface
interface cov_wend1
   module procedure cov_wend1_vec, cov_wend1_csr
end interface
interface cov_wend2
   module procedure cov_wend2_vec, cov_wend2_csr
end interface
interface cov_mat
   module procedure cov_mat_vec, cov_mat_csr
end interface
interface cov_finnmat
   module procedure cov_finnmat_vec, cov_finnmat_csr
end interface
interface cov_mat12
   module procedure cov_mat12_vec, cov_mat12_csr
end interface
interface cov_mat32
   module procedure cov_mat32_vec, cov_mat32_csr
end interface
interface cov_mat52
   module procedure cov_mat52_vec, cov_mat52_csr
end interface

contains

function theta_at(theta,i,default) result(x)
real(dp),intent(in)::theta(:),default
integer,intent(in)::i
real(dp)::x
if(i<=size(theta))then
x=abs(theta(i))
else
x=default
end if
end function theta_at

function cov_apply_vec(h,theta,kind,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
integer,intent(in)::kind
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
integer::i
allocate(v(size(h)))
do i=1,size(h)
v(i)=cov_value(h(i),theta,kind,eps)
end do
end function cov_apply_vec

function cov_apply_csr(h,theta,kind,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
integer,intent(in)::kind
real(dp),intent(in),optional::eps
type(csr_matrix)::v
integer::k
v=h
if(.not.allocated(v%entries))return
do k=1,size(v%entries)
v%entries(k)=cov_value(v%entries(k),theta,kind,eps)
end do
end function cov_apply_csr

real(dp) function cov_value(h,theta,kind,eps) result(v)
real(dp),intent(in)::h,theta(:)
integer,intent(in)::kind
real(dp),intent(in),optional::eps
real(dp)::tol,range,sill,nug,nu,x
real(dp),parameter::one=1.0_dp
tol=epsilon(one)
if(present(eps))tol=eps
select case(kind)
case(4) ! nugget has a one-parameter convention
   sill=theta_at(theta,1,1.0_dp)
   if(h<tol)then
   v=sill
   else
   v=0.0_dp
   end if
   return
case default
   range=theta_at(theta,1,1.0_dp)
   sill=theta_at(theta,2,1.0_dp)
   if(kind==10 .or. kind==11) then
      nug=theta_at(theta,4,0.0_dp)
   else
      nug=theta_at(theta,3,0.0_dp)
   end if
end select
if(range<tol)then
   if(h<tol)then
   v=sill+nug
   else
   v=0.0_dp
   end if
   return
end if
x=h/range
if(h<tol)then
   v=sill+nug
   return
end if
select case(kind)
case(1) ! exponential / Matern 1/2
   v=sill*exp(-x)
case(2) ! spherical
   if(x>=1.0_dp)then
   v=0.0_dp
   else
   v=sill*(1.0_dp-1.5_dp*x+0.5_dp*x*x*x)
   end if
case(3) ! spherical correlation
   if(x>=1.0_dp)then
   v=0.0_dp
   else
   v=1.0_dp-1.5_dp*x+0.5_dp*x*x*x
   end if
case(5) ! Wu1
   if(x>=1.0_dp)then
   v=0.0_dp
   else
   v=sill*(1.0_dp-x)**3*(1.0_dp+3.0_dp*x+x*x)
   end if
case(6) ! Wu2
   if(x>=1.0_dp)then
   v=0.0_dp
   else
   v=sill*(1.0_dp-x)**4*(4.0_dp+16.0_dp*x+12.0_dp*x*x+3.0_dp*x**3)/4.0_dp
   end if
case(7) ! Wu3
   if(x>=1.0_dp)then
      v=0.0_dp
   else
      v=sill*(1.0_dp-x)**6*(1.0_dp+6.0_dp*x+(41.0_dp/3.0_dp)*x*x+12.0_dp*x**3+5.0_dp*x**4+(5.0_dp/6.0_dp)*x**5)
   end if
case(8) ! Wendland 1
   if(x>=1.0_dp)then
   v=0.0_dp
   else
   v=sill*(1.0_dp-x)**4*(4.0_dp*x+1.0_dp)
   end if
case(9) ! Wendland 2
   if(x>=1.0_dp)then
   v=0.0_dp
   else
   v=sill*(1.0_dp-x)**6*(35.0_dp*x*x+18.0_dp*x+3.0_dp)/3.0_dp
   end if
case(10) ! general Matern: theta=(range,sill,nu,nugget)
   nu=theta_at(theta,3,1.0_dp)
   nug=theta_at(theta,4,0.0_dp)
   if(h<tol)then
      v=sill+nug
   else
      v=sill*(2.0_dp**(1.0_dp-nu)/gamma(nu))*x**nu*besselK(x,nu)
   end if
case(11) ! Furrer/INLA Matern range parametrization
   nu=theta_at(theta,3,1.0_dp)
   nug=theta_at(theta,4,0.0_dp)
   x=sqrt(8.0_dp*nu)*h/range
   if(h<tol)then
      v=sill+nug
   else
      v=sill*(2.0_dp**(1.0_dp-nu)/gamma(nu))*x**nu*besselK(x,nu)
   end if
case(12) ! Matern 3/2 under spam's range convention
   v=sill*exp(-x)*(1.0_dp+x)
case(13) ! Matern 5/2 under spam's range convention
   v=sill*exp(-x)*(1.0_dp+x+x*x/3.0_dp)
case default
   error stop 'cov_value: unknown covariance kind'
end select
end function cov_value

function cov_exp_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,1,eps)
end function
function cov_exp_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,1,eps)
end function
function cov_sph_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,2,eps)
end function
function cov_sph_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,2,eps)
end function
function cor_sph_vec(h,range,eps) result(v)
real(dp),intent(in)::h(:)
real(dp),intent(in)::range
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
real(dp)::th(3)
th=[range,1.0_dp,0.0_dp]
v=cov_apply_vec(h,th,3,eps)
end function
function cor_sph_csr(h,range,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::range
real(dp),intent(in),optional::eps
type(csr_matrix)::v
real(dp)::th(3)
th=[range,1.0_dp,0.0_dp]
v=cov_apply_csr(h,th,3,eps)
end function
function cov_nug_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,4,eps)
end function
function cov_nug_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,4,eps)
end function
function cov_wu1_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,5,eps)
end function
function cov_wu1_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,5,eps)
end function
function cov_wu2_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,6,eps)
end function
function cov_wu2_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,6,eps)
end function
function cov_wu3_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,7,eps)
end function
function cov_wu3_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,7,eps)
end function
function cov_wend1_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,8,eps)
end function
function cov_wend1_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,8,eps)
end function
function cov_wend2_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,9,eps)
end function
function cov_wend2_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,9,eps)
end function
function cov_mat_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,10,eps)
end function
function cov_mat_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,10,eps)
end function
function cov_finnmat_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,11,eps)
end function
function cov_finnmat_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,11,eps)
end function
function cov_mat12_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,1,eps)
end function
function cov_mat12_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,1,eps)
end function
function cov_mat32_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,12,eps)
end function
function cov_mat32_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,12,eps)
end function
function cov_mat52_vec(h,theta,eps) result(v)
real(dp),intent(in)::h(:),theta(:)
real(dp),intent(in),optional::eps
real(dp),allocatable::v(:)
v=cov_apply_vec(h,theta,13,eps)
end function
function cov_mat52_csr(h,theta,eps) result(v)
type(csr_matrix),intent(in)::h
real(dp),intent(in)::theta(:)
real(dp),intent(in),optional::eps
type(csr_matrix)::v
v=cov_apply_csr(h,theta,13,eps)
end function

end module spam_covariance
