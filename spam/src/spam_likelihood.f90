module spam_likelihood
use spam_kinds,only:dp
use spam_types,only:csr_matrix,mle_result,spam_chol
use spam_covariance,only:cov_exp,cov_sph,cov_nug,cov_wu1,cov_wu2,cov_wu3,cov_wend1,cov_wend2, &
                         cov_mat,cov_finnmat,cov_mat12,cov_mat32,cov_mat52
use spam_cholesky,only:spam_chol_factor,spam_solve,spam_logdet
implicit none
private
public::neg2loglikelihood_spam,neg2loglikelihood_nomean_spam,mle_spam,mle_nomean_spam,make_covariance
public::neg2loglikelihood_custom,mle_spam_custom,mle_nomean_spam_custom
abstract interface
 function covariance_builder(dist,theta) result(s)
  import dp,csr_matrix
  type(csr_matrix),intent(in)::dist
  real(dp),intent(in)::theta(:)
  type(csr_matrix)::s
 end function covariance_builder
end interface
contains

function make_covariance(dist,theta,model) result(s)
type(csr_matrix),intent(in)::dist
real(dp),intent(in)::theta(:)
character(len=*),intent(in)::model
type(csr_matrix)::s
select case(lower(trim(model)))
case('exp','exponential','mat12');s=cov_exp(dist,theta)
case('sph','spherical');s=cov_sph(dist,theta)
case('nug','nugget');s=cov_nug(dist,theta)
case('wu1')
s=cov_wu1(dist,theta)
case('wu2')
s=cov_wu2(dist,theta)
case('wu3')
s=cov_wu3(dist,theta)
case('wend1','wendland1')
s=cov_wend1(dist,theta)
case('wend2','wendland2')
s=cov_wend2(dist,theta)
case('mat','matern')
s=cov_mat(dist,theta)
case('finnmat')
s=cov_finnmat(dist,theta)
case('mat32')
s=cov_mat32(dist,theta)
case('mat52')
s=cov_mat52(dist,theta)
case default;error stop 'make_covariance: unknown model'
end select
end function

real(dp) function neg2loglikelihood_spam(y,x,dist,beta,theta,model) result(v)
real(dp),intent(in)::y(:),x(:,:),beta(:),theta(:)
type(csr_matrix),intent(in)::dist
character(len=*),intent(in)::model
type(csr_matrix)::s
type(spam_chol)::f
real(dp),allocatable::r(:),z(:)
if(size(x,1)/=size(y).or.size(x,2)/=size(beta))error stop 'neg2loglikelihood_spam: dimension mismatch'
s=make_covariance(dist,theta,model)
f=spam_chol_factor(s)
if(f%info/=0)then
v=huge(1.0_dp)/16.0_dp
return
end if
r=y-matmul(x,beta)
z=spam_solve(f,r)
v=real(size(y),dp)*log(2.0_dp*acos(-1.0_dp))+spam_logdet(f)+dot_product(r,z)
end function
real(dp) function neg2loglikelihood_nomean_spam(y,dist,theta,model) result(v)
real(dp),intent(in)::y(:),theta(:)
type(csr_matrix),intent(in)::dist
character(len=*),intent(in)::model
type(csr_matrix)::s
type(spam_chol)::f
real(dp),allocatable::z(:)
s=make_covariance(dist,theta,model)
f=spam_chol_factor(s)
if(f%info/=0)then
v=huge(1.0_dp)/16.0_dp
return
end if
z=spam_solve(f,y)
v=real(size(y),dp)*log(2.0_dp*acos(-1.0_dp))+spam_logdet(f)+dot_product(y,z)
end function

real(dp) function neg2loglikelihood_custom(y,x,dist,beta,theta,covfun) result(v)
real(dp),intent(in)::y(:),x(:,:),beta(:),theta(:)
type(csr_matrix),intent(in)::dist
procedure(covariance_builder)::covfun
type(csr_matrix)::s
type(spam_chol)::f
real(dp),allocatable::r(:),z(:)
if(size(x,1)/=size(y).or.size(x,2)/=size(beta))error stop 'neg2loglikelihood_custom: dimension mismatch'
s=covfun(dist,theta)
f=spam_chol_factor(s)
if(f%info/=0)then
v=huge(1.0_dp)/16.0_dp
return
end if
r=y-matmul(x,beta)
z=spam_solve(f,r)
v=real(size(y),dp)*log(2.0_dp*acos(-1.0_dp))+spam_logdet(f)+dot_product(r,z)
end function

function mle_spam(y,x,dist,beta0,theta0,thetalower,thetaupper,model,maxit,reltol) result(out)
real(dp),intent(in)::y(:),x(:,:),beta0(:),theta0(:),thetalower(:),thetaupper(:)
type(csr_matrix),intent(in)::dist
character(len=*),intent(in)::model
integer,intent(in),optional::maxit
real(dp),intent(in),optional::reltol
type(mle_result)::out
real(dp),allocatable::p(:),lo(:),hi(:)
integer::nb
nb=size(beta0)
allocate(p(nb+size(theta0)),lo(nb+size(theta0)),hi(nb+size(theta0)))
p=[beta0,theta0]
lo(:nb)=-huge(1.0_dp)/100
hi(:nb)=huge(1.0_dp)/100
lo(nb+1:)=thetalower
hi(nb+1:)=thetaupper
call bounded_bfgs(obj,p,lo,hi,maxit,reltol,out)
contains
 real(dp) function obj(z) result(f)
 real(dp),intent(in)::z(:)
 f=neg2loglikelihood_spam(y,x,dist,z(:nb),z(nb+1:),model)
 end function
end function

function mle_nomean_spam(y,dist,theta0,thetalower,thetaupper,model,maxit,reltol) result(out)
real(dp),intent(in)::y(:),theta0(:),thetalower(:),thetaupper(:)
type(csr_matrix),intent(in)::dist
character(len=*),intent(in)::model
integer,intent(in),optional::maxit
real(dp),intent(in),optional::reltol
type(mle_result)::out
real(dp),allocatable::p(:)
p=theta0
call bounded_bfgs(obj,p,thetalower,thetaupper,maxit,reltol,out)
contains
 real(dp) function obj(z) result(f)
 real(dp),intent(in)::z(:)
 f=neg2loglikelihood_nomean_spam(y,dist,z,model)
 end function
end function

function mle_spam_custom(y,x,dist,beta0,theta0,thetalower,thetaupper,covfun,maxit,reltol) result(out)
real(dp),intent(in)::y(:),x(:,:),beta0(:),theta0(:),thetalower(:),thetaupper(:)
type(csr_matrix),intent(in)::dist
procedure(covariance_builder)::covfun
integer,intent(in),optional::maxit
real(dp),intent(in),optional::reltol
type(mle_result)::out
real(dp),allocatable::p(:),lo(:),hi(:)
integer::nb
nb=size(beta0)
allocate(p(nb+size(theta0)),lo(nb+size(theta0)),hi(nb+size(theta0)))
p=[beta0,theta0]
lo(:nb)=-huge(1.0_dp)/100
hi(:nb)=huge(1.0_dp)/100
lo(nb+1:)=thetalower
hi(nb+1:)=thetaupper
call bounded_bfgs(obj,p,lo,hi,maxit,reltol,out)
contains
 real(dp) function obj(z) result(f)
 real(dp),intent(in)::z(:)
 f=neg2loglikelihood_custom(y,x,dist,z(:nb),z(nb+1:),covfun)
 end function
end function mle_spam_custom

function mle_nomean_spam_custom(y,dist,theta0,thetalower,thetaupper,covfun,maxit,reltol) result(out)
real(dp),intent(in)::y(:),theta0(:),thetalower(:),thetaupper(:)
type(csr_matrix),intent(in)::dist
procedure(covariance_builder)::covfun
integer,intent(in),optional::maxit
real(dp),intent(in),optional::reltol
type(mle_result)::out
real(dp),allocatable::p(:),z0(:,:)
integer::n
n=size(y)
allocate(z0(n,0))
p=theta0
call bounded_bfgs(obj,p,thetalower,thetaupper,maxit,reltol,out)
contains
 real(dp) function obj(z) result(f)
 real(dp),intent(in)::z(:)
 type(csr_matrix)::s
 type(spam_chol)::cf
 real(dp),allocatable::q(:)
 s=covfun(dist,z)
 cf=spam_chol_factor(s)
 if(cf%info/=0)then
 f=huge(1.0_dp)/16.0_dp
 return
 end if
 q=spam_solve(cf,y)
 f=real(n,dp)*log(2.0_dp*acos(-1.0_dp))+spam_logdet(cf)+dot_product(y,q)
 end function
end function mle_nomean_spam_custom

subroutine bounded_bfgs(fn,p,lo,hi,maxit,reltol,out)
interface
 function fn(x) result(f)
  import dp
  real(dp),intent(in)::x(:)
  real(dp)::f
 end function
end interface
real(dp),intent(inout)::p(:)
real(dp),intent(in)::lo(:),hi(:)
integer,intent(in),optional::maxit
real(dp),intent(in),optional::reltol
type(mle_result),intent(out)::out
real(dp),allocatable::g(:),gn(:),h(:,:),d(:),pn(:),s(:),y(:),hy(:)
real(dp)::f,fnv,step,tol,ys,yhy,epsi
integer::n,it,mi,j
logical::ok
n=size(p)
mi=100
if(present(maxit))mi=maxit
tol=1e-7_dp
if(present(reltol))tol=reltol
p=max(lo,min(hi,p))
allocate(g(n),gn(n),h(n,n),d(n),pn(n),s(n),y(n),hy(n))
h=0
do j=1,n
h(j,j)=1.0_dp
end do
it=0
f=fn(p)
call fdgrad(fn,p,lo,hi,g)
out%convergence=1
if(.not.isfinite(f))goto 90
do it=1,mi
 if(maxval(abs(g)*max(1.0_dp,abs(p)))<=tol)then
 out%convergence=0
 exit
 end if
 d=-matmul(h,g)
 if(dot_product(d,g)>=0.0_dp)d=-g
 step=1.0_dp
 ok=.false.
 do j=1,25
  pn=max(lo,min(hi,p+step*d))
  fnv=fn(pn)
  if(isfinite(fnv).and.fnv<=f+1e-4_dp*dot_product(g,pn-p))then
  ok=.true.
  exit
  end if
  step=0.5_dp*step
 end do
 if(.not.ok)exit
 call fdgrad(fn,pn,lo,hi,gn)
 s=pn-p
 y=gn-g
 ys=dot_product(y,s)
 if(ys>sqrt(epsilon(1.0_dp))*sqrt(dot_product(s,s)*dot_product(y,y)))then
  hy=matmul(h,y)
  yhy=dot_product(y,hy)
  h=h+((ys+yhy)/(ys*ys))*outer(s,s)-(outer(hy,s)+outer(s,hy))/ys
 else
  h=0
  do j=1,n
  h(j,j)=1.0_dp
  end do
 end if
 if(abs(f-fnv)<=tol*(1.0_dp+abs(f)))then
 p=pn
 f=fnv
 g=gn
 out%convergence=0
 exit
 end if
 p=pn
 f=fnv
 g=gn
end do
90 continue
allocate(out%par(n))
out%par=p
out%objective=f
out%iterations=min(it,mi)
end subroutine

subroutine fdgrad(fn,p,lo,hi,g)
interface
 function fn(x) result(f)
 import dp
 real(dp),intent(in)::x(:)
 real(dp)::f
 end function
end interface
real(dp),intent(in)::p(:),lo(:),hi(:)
real(dp),intent(out)::g(:)
real(dp),allocatable::q(:)
real(dp)::hh,fp,fm
integer::j
q=p
do j=1,size(p)
 hh=1e-5_dp*max(1.0_dp,abs(p(j)))
 if(p(j)+hh<=hi(j).and.p(j)-hh>=lo(j))then
  q(j)=p(j)+hh
  fp=fn(q)
  q(j)=p(j)-hh
  fm=fn(q)
  g(j)=(fp-fm)/(2*hh)
 else if(p(j)+hh<=hi(j))then
  q(j)=p(j)+hh
  fp=fn(q)
  q(j)=p(j)
  fm=fn(q)
  g(j)=(fp-fm)/hh
 else
  q(j)=p(j)-hh
  fm=fn(q)
  q(j)=p(j)
  fp=fn(q)
  g(j)=(fp-fm)/hh
 end if
 q(j)=p(j)
end do
end subroutine
pure function outer(a,b) result(c)
real(dp),intent(in)::a(:),b(:)
real(dp)::c(size(a),size(b))
integer::i
do i=1,size(a)
c(i,:)=a(i)*b
end do
end function
pure logical function isfinite(x) result(ok)
use,intrinsic::ieee_arithmetic,only:ieee_is_finite
real(dp),intent(in)::x
ok=ieee_is_finite(x)
end function
pure function lower(s) result(t)
character(len=*),intent(in)::s
character(len=len(s))::t
integer::i,c
t=s
do i=1,len(s)
c=iachar(t(i:i))
if(c>=65.and.c<=90)t(i:i)=achar(c+32)
end do
end function
end module spam_likelihood
