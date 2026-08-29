! Minimax-tilting probability and simulation algorithms from TruncatedNormal 2.3.
! SPDX-License-Identifier: GPL-3.0-only
module truncated_normal_core
 use, intrinsic :: iso_fortran_env, only: int32
 use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
 use r_compat, only: dp, runif1, sd, pt, qt
 use nleqslv_fortran, only: nleq_options,nleq_result,solve_nleqslv,NLEQ_BROYDEN,NLEQ_NEWTON,NLEQ_PWLDOG
 use alabama, only: auglag3, alabama_result_t, alabama_outer_control_t
 use spacefillr_sobol, only: generate_sobol_owen_set
 use qrng_sobol_mod, only: sobol
 use truncated_normal_math, only: lnNpr,norminvp,trandn
 use truncated_normal_linear, only: cholperm,cholperm_result
 use truncated_normal_tilting, only: gradpsi_eval,jacpsi_eval,psy_eval,gradpsi_t_eval,psy_t_eval
 implicit none
 private
 real(dp),parameter::pi=acos(-1.0_dp)
 type, public :: prob_result
   real(dp) :: prob = 0.0_dp, err = 0.0_dp, relerr = 0.0_dp, upbnd = 0.0_dp
   integer :: status = 0
 end type
 type, public :: tregress_result
   real(dp), allocatable :: r(:)
   real(dp), allocatable :: z(:,:)
   integer :: status = 0
 end type
 public :: mvncdf, mvnqmc, mvtcdf, mvtqmc, mvrandn, mvrandt, tregress
 public::mvnpr,mvnprqmc,mvtpr,mvtprqmc
contains

subroutine prep(sigma,lo,hi,lmat,l,u,lfull,perm,status)
 real(dp),intent(in)::sigma(:,:),lo(:),hi(:)
 real(dp),allocatable,intent(out)::lmat(:,:),l(:),u(:),lfull(:,:)
 integer,allocatable,intent(out)::perm(:)
 integer,intent(out)::status
 type(cholperm_result)::cp
 real(dp),allocatable::d(:)
 integer::i,n
 n=size(lo)
 call cholperm(sigma,lo,hi,cp,'GGE')
 status=cp%status
 if(status/=0)return
 allocate(lfull(n,n),lmat(n,n),l(n),u(n),perm(n),d(n))
 lfull=cp%lmat
 perm=cp%perm
 d=[(lfull(i,i),i=1,n)]
 l=cp%lower/d
 u=cp%upper/d
 lmat=lfull
 do i=1,n
 lmat(i,:)=lmat(i,:)/d(i)
 lmat(i,i)=0.0_dp
 end do
end subroutine

subroutine solve_normal_tilt(lmat,l,u,x,mu,status,use_newton)
 real(dp),intent(in)::lmat(:,:),l(:),u(:)
 real(dp),allocatable,intent(out)::x(:),mu(:)
 integer,intent(out)::status
 logical,intent(in),optional::use_newton
 integer::d
 real(dp),allocatable::x0(:)
 type(nleq_options)::opt
 type(nleq_result)::res
 logical::nw
 d=size(l)
 allocate(x0(2*d-2))
 x0=0
 nw=.false.
 if(present(use_newton))nw=use_newton
 opt=nleq_options()
 opt%global=NLEQ_PWLDOG
 opt%maxit=500
 if(nw)then
 opt%method=NLEQ_NEWTON
 call solve_nleqslv(x0,fn,res,opt,jac)
 else
 opt%method=NLEQ_BROYDEN
 call solve_nleqslv(x0,fn,res,opt,jac)
 end if
 allocate(x(d-1),mu(d-1))
 x=res%x(1:d-1)
 mu=res%x(d:2*d-2)
 status=merge(0,res%termcd,res%termcd==1.or.res%termcd==2)
 if(maxval(abs(res%fvec))>1e-6_dp)status=100+res%termcd
 if(status/=0) call fallback_alabama(res%x,x,mu,status)
contains
 subroutine fallback_alabama(y0,xo,mo,istat)
  real(dp),intent(in)::y0(:)
  real(dp),intent(out)::xo(:),mo(:)
  integer,intent(out)::istat
  real(dp)::par(size(y0))
  type(alabama_result_t)::ar
  type(alabama_outer_control_t)::oc
  par=y0
  oc=alabama_outer_control_t()
  oc%trace=.false.
  oc%method='BFGS'
  oc%itmax=40
  call auglag3(par,obj,hin,heq,ar,gr=ogr,control_outer=oc)
  if(ar%convergence==0)then
  xo=par(1:d-1)
  mo=par(d:2*d-2)
  istat=0
  else
  istat=200+ar%convergence
  end if
 end subroutine
 function obj(y) result(v)
  real(dp),intent(in)::y(:)
  real(dp)::v
  v=-psy_eval(y(1:d-1),lmat,l,u,y(d:2*d-2))
 end function
 subroutine ogr(y,g)
  real(dp),intent(in)::y(:)
  real(dp),intent(out)::g(:)
  call gradpsi_eval(y,lmat,l,u,g)
  g=-g
 end subroutine
 function heq(y) result(h)
  real(dp),intent(in)::y(:)
  real(dp),allocatable::h(:)
  real(dp)::gg(2*d-2)
  call gradpsi_eval(y,lmat,l,u,gg)
  allocate(h(d-1))
  h=gg(d:2*d-2)
 end function
 function hin(y) result(h)
  real(dp),intent(in)::y(:)
  real(dp),allocatable::h(:)
  real(dp)::xx(d),cv(d)
  allocate(h(2*d-2))
  xx=0
  xx(1:d-1)=y(1:d-1)
  cv=xx+matmul(lmat,xx)
  h(1:d-1)=u(1:d-1)-cv(1:d-1)
  h(d:)=cv(1:d-1)-l(1:d-1)
 end function
 subroutine fn(y,f)
 real(dp),intent(in)::y(:)
 real(dp),intent(out)::f(:)
 call gradpsi_eval(y,lmat,l,u,f)
 end subroutine
 subroutine jac(y,jm)
 real(dp),intent(in)::y(:)
 real(dp),intent(out)::jm(:,:)
 call jacpsi_eval(y,lmat,l,u,jm)
 end subroutine
end subroutine

subroutine solve_student_tilt(lmat,l,u,nu,x,mu,status)
 real(dp),intent(in)::lmat(:,:),l(:),u(:),nu
 real(dp),allocatable,intent(out)::x(:),mu(:)
 integer,intent(out)::status
 integer::d
 real(dp),allocatable::x0(:)
 type(nleq_options)::opt
 type(nleq_result)::res
 d=size(l)
 allocate(x0(2*d))
 x0=0
 x0(d)=log(sqrt(nu))
 x0(2*d)=sqrt(nu)
 opt=nleq_options()
 opt%global=NLEQ_PWLDOG
 opt%method=NLEQ_BROYDEN
 opt%maxit=500
 call solve_nleqslv(x0,fn,res,opt)
 allocate(x(d),mu(d))
 x=res%x(1:d)
 x(d)=exp(x(d))
 mu=res%x(d+1:2*d)
 status=merge(0,res%termcd,res%termcd==1.or.res%termcd==2)
 if(maxval(abs(res%fvec))>1e-6_dp)status=100+res%termcd
contains
 subroutine fn(y,f)
 real(dp),intent(in)::y(:)
 real(dp),intent(out)::f(:)
 call gradpsi_t_eval(y,lmat,l,u,nu,f)
 end subroutine
end subroutine

function mvnpr(n,lmat,l,u,muin) result(res)
 integer,intent(in)::n
 real(dp),intent(in)::lmat(:,:),l(:),u(:),muin(:)
 type(prob_result)::res
 integer::d,k,j
 real(dp),allocatable::z(:,:),p(:),col(:),tl(:),tu(:),mu(:),w(:)
 real(dp)::m,ss
 d=size(l)
 allocate(z(d,n),p(n),col(n),tl(n),tu(n),mu(d))
 z=0
 p=0
 mu=0
 mu(1:min(d-1,size(muin)))=muin(1:min(d-1,size(muin)))
 do k=1,d-1
  col=0
  do j=1,k-1
  col=col+lmat(k,j)*z(j,:)
  end do
  tl=l(k)-mu(k)-col
  tu=u(k)-mu(k)-col
  z(k,:)=mu(k)+trandn(tl,tu)
  w=lnNpr(tl,tu)
  p=p+w+0.5_dp*mu(k)**2-mu(k)*z(k,:)
 end do
 col=0
 do j=1,d-1
 col=col+lmat(d,j)*z(j,:)
 end do
 p=p+lnNpr(l(d)-col,u(d)-col)
 p=exp(p)
 m=sum(p)/n
 ss=sqrt(sum((p-m)**2)/max(1,n-1))
 res%prob=m
 res%err=ss
 res%relerr=ss/sqrt(real(n,dp))/max(m,tiny(m))
end function

function mvnprqmc(n,lmat,l,u,muin,seed) result(prob)
 integer,intent(in)::n
 real(dp),intent(in)::lmat(:,:),l(:),u(:),muin(:)
 integer,intent(in),optional::seed
 real(dp)::prob
 integer::d,k,j,s
 real(dp),allocatable::z(:,:),p(:),col(:),tl(:),tu(:),mu(:),w(:),q(:,:)
 integer(int32)::s32
 d=size(l)
 allocate(z(d,n),p(n),col(n),tl(n),tu(n),mu(d))
 z=0
 p=0
 mu=0
 mu(1:min(d-1,size(muin)))=muin(1:min(d-1,size(muin)))
 s=ceiling(1e6_dp*runif1())
 if(present(seed))s=seed
 s32=int(s,int32)
 allocate(q(n,d-1))
 call generate_sobol_owen_set(n,d-1,q,s32)
 do k=1,d-1
  col=0
  do j=1,k-1
  col=col+lmat(k,j)*z(j,:)
  end do
  tl=l(k)-mu(k)-col
  tu=u(k)-mu(k)-col
  z(k,:)=mu(k)+norminvp(q(:,k),tl,tu)
  w=lnNpr(tl,tu)
  p=p+w+0.5_dp*mu(k)**2-mu(k)*z(k,:)
 end do
 col=0
 do j=1,d-1
 col=col+lmat(d,j)*z(j,:)
 end do
 p=p+lnNpr(l(d)-col,u(d)-col)
 prob=sum(exp(p))/n
end function

function mvncdf(lo,hi,sigma,n) result(res)
 real(dp),intent(in)::lo(:),hi(:),sigma(:,:)
 integer,intent(in),optional::n
 type(prob_result)::res
 real(dp),allocatable::lmat(:,:),l(:),u(:),lfull(:,:),x(:),mu(:)
 integer,allocatable::perm(:)
 integer::nn,st,d
 real(dp)::tmp(1)
 nn=100000
 if(present(n))nn=n
 d=size(lo)
 if(d==1)then
 tmp=lnNpr([lo(1)/sqrt(sigma(1,1))],[hi(1)/sqrt(sigma(1,1))])
 res%prob=exp(tmp(1))
 return
 end if
 call prep(sigma,lo,hi,lmat,l,u,lfull,perm,st)
 if(st/=0)then
 res%status=st
 return
 end if
 call solve_normal_tilt(lmat,l,u,x,mu,st,.false.)
 res=mvnpr(nn,lmat,l,u,mu)
 res%status=st
 res%upbnd=exp(psy_eval(x,lmat,l,u,mu))
end function

function mvnqmc(lo,hi,sigma,n) result(res)
 real(dp),intent(in)::lo(:),hi(:),sigma(:,:)
 integer,intent(in),optional::n
 type(prob_result)::res
 real(dp),allocatable::lmat(:,:),l(:),u(:),lfull(:,:),x(:),mu(:)
 integer,allocatable::perm(:)
 integer::nn,st,i,m,d
 real(dp)::p(12),av,ss,tmp(1)
 nn=100000
 if(present(n))nn=n
 d=size(lo)
 if(d==1)then
 tmp=lnNpr([lo(1)/sqrt(sigma(1,1))],[hi(1)/sqrt(sigma(1,1))])
 res%prob=exp(tmp(1))
 return
 end if
 call prep(sigma,lo,hi,lmat,l,u,lfull,perm,st)
 if(st/=0)then
 res%status=st
 return
 end if
 call solve_normal_tilt(lmat,l,u,x,mu,st,.true.)
 m=(nn+11)/12
 do i=1,12
 p(i)=mvnprqmc(m,lmat,l,u,mu)
 end do
 av=sum(p)/12
 ss=sqrt(sum((p-av)**2)/11)
 res%prob=av
 res%relerr=ss/sqrt(12.0_dp)/max(av,tiny(av))
 res%upbnd=exp(psy_eval(x,lmat,l,u,mu))
 res%status=st
end function

function mvtpr(n,lmat,l,u,nu,muin) result(res)
 integer,intent(in)::n
 real(dp),intent(in)::lmat(:,:),l(:),u(:),nu,muin(:)
 type(prob_result)::res
 integer::d,k,j
 real(dp),allocatable::z(:,:),p(:),r(:),col(:),tl(:),tu(:),mu(:),w(:)
 real(dp)::eta,c,m,ss,inf
 real(dp)::t1(1)
 d=size(l)
 eta=muin(d)
 allocate(z(d,n),p(n),r(n),col(n),tl(n),tu(n),mu(d))
 z=0
 mu=muin
 mu(d)=0
 inf=ieee_value(1.0_dp,ieee_positive_inf)
 t1=lnNpr([-eta],[inf])
 c=0.5_dp*log(2*pi)-log_gamma(nu/2)-(nu/2-1)*log(2.0_dp)+t1(1)+0.5_dp*eta*eta
 r=eta+trandn(spread(-eta,1,n),spread(inf,1,n))
 p=(nu-1)*log(r)-eta*r
 r=r/sqrt(nu)
 do k=1,d-1
 col=0
 do j=1,k-1
 col=col+lmat(k,j)*z(j,:)
 end do
 tl=r*l(k)-mu(k)-col
 tu=r*u(k)-mu(k)-col
 z(k,:)=mu(k)+trandn(tl,tu)
 w=lnNpr(tl,tu)
 p=p+w+0.5_dp*mu(k)**2-mu(k)*z(k,:)
 end do
 col=0
 do j=1,d-1
 col=col+lmat(d,j)*z(j,:)
 end do
 p=exp(p+lnNpr(r*l(d)-col,r*u(d)-col))
 m=exp(c)*sum(p)/n
 ss=exp(c)*sqrt(sum((p-sum(p)/n)**2)/max(1,n-1))
 res%prob=m
 res%err=ss
 res%relerr=ss/sqrt(real(n,dp))/max(m,tiny(m))
end function

function mvtprqmc(n,lmat,l,u,nu,muin,seed) result(prob)
 integer,intent(in)::n
 real(dp),intent(in)::lmat(:,:),l(:),u(:),nu,muin(:)
 integer,intent(in),optional::seed
 real(dp)::prob
 integer::d,k,j,s
 real(dp),allocatable::z(:,:),p(:),r(:),col(:),tl(:),tu(:),mu(:),q(:,:),w(:)
 real(dp)::eta,c,inf,t1(1)
 d=size(l)
 eta=muin(d)
 allocate(z(d,n),p(n),r(n),col(n),tl(n),tu(n),mu(d))
 z=0
 mu=muin
 mu(d)=0
 inf=ieee_value(1.0_dp,ieee_positive_inf)
 t1=lnNpr([-eta],[inf])
 c=0.5_dp*log(2*pi)-log_gamma(nu/2)-(nu/2-1)*log(2.0_dp)+t1(1)+0.5_dp*eta*eta
 r=eta+trandn(spread(-eta,1,n),spread(inf,1,n))
 p=(nu-1)*log(r)-eta*r
 r=r/sqrt(nu)
 s=ceiling(1e6_dp*runif1())
 if(present(seed))s=seed
 q=sobol(n,d-1,randomize=.true.,seed=s)
 do k=1,d-1
 col=0
 do j=1,k-1
 col=col+lmat(k,j)*z(j,:)
 end do
 tl=r*l(k)-mu(k)-col
 tu=r*u(k)-mu(k)-col
 z(k,:)=mu(k)+norminvp(q(:,k),tl,tu)
 w=lnNpr(tl,tu)
 p=p+w+0.5_dp*mu(k)**2-mu(k)*z(k,:)
 end do
 col=0
 do j=1,d-1
 col=col+lmat(d,j)*z(j,:)
 end do
 p=p+lnNpr(r*l(d)-col,r*u(d)-col)
 prob=exp(c)*sum(exp(p))/n
end function

function mvtcdf(lo,hi,sigma,nu,n) result(res)
 real(dp),intent(in)::lo(:),hi(:),sigma(:,:),nu
 integer,intent(in),optional::n
 type(prob_result)::res
 real(dp),allocatable::lmat(:,:),l(:),u(:),lfull(:,:),x(:),mu(:)
 integer,allocatable::perm(:)
 integer::nn,st,d
 nn=100000
 if(present(n))nn=n
 d=size(lo)
 if(d==1)then
 res%prob=pt(hi(1)/sqrt(sigma(1,1)),nu)-pt(lo(1)/sqrt(sigma(1,1)),nu)
 return
 end if
 call prep(sigma,lo,hi,lmat,l,u,lfull,perm,st)
 if(st/=0)then
 res%status=st
 return
 end if
 call solve_student_tilt(lmat,l,u,nu,x,mu,st)
 res=mvtpr(nn,lmat,l,u,nu,mu)
 res%upbnd=exp(psy_t_eval(x,lmat,l,u,nu,mu))
 res%status=st
end function
function mvtqmc(lo,hi,sigma,nu,n) result(res)
 real(dp),intent(in)::lo(:),hi(:),sigma(:,:),nu
 integer,intent(in),optional::n
 type(prob_result)::res
 real(dp),allocatable::lmat(:,:),l(:),u(:),lfull(:,:),x(:),mu(:)
 integer,allocatable::perm(:)
 integer::nn,st,d,i,m
 real(dp)::p(12),av,ss
 nn=100000
 if(present(n))nn=n
 d=size(lo)
 if(d==1)then
 res%prob=pt(hi(1)/sqrt(sigma(1,1)),nu)-pt(lo(1)/sqrt(sigma(1,1)),nu)
 return
 end if
 call prep(sigma,lo,hi,lmat,l,u,lfull,perm,st)
 if(st/=0)then
 res%status=st
 return
 end if
 call solve_student_tilt(lmat,l,u,nu,x,mu,st)
 m=(nn+11)/12
 do i=1,12
 p(i)=mvtprqmc(m,lmat,l,u,nu,mu)
 end do
 av=sum(p)/12
 ss=sqrt(sum((p-av)**2)/11)
 res%prob=av
 res%relerr=ss/sqrt(12.0_dp)/max(av,tiny(av))
 res%upbnd=exp(psy_t_eval(x,lmat,l,u,nu,mu))
 res%status=st
end function

subroutine proposal_normal(n,lmat,l,u,muin,z,logpr)
 integer,intent(in)::n
 real(dp),intent(in)::lmat(:,:),l(:),u(:),muin(:)
 real(dp),intent(out)::z(:,:),logpr(:)
 integer::d,k,j
 real(dp)::mu(size(l)),col(n),tl(n),tu(n)
 d=size(l)
 mu=0
 mu(1:min(d-1,size(muin)))=muin(1:min(d-1,size(muin)))
 z=0
 logpr=0
 do k=1,d
 col=0
 do j=1,k-1
 col=col+lmat(k,j)*z(j,:)
 end do
 tl=l(k)-mu(k)-col
 tu=u(k)-mu(k)-col
 z(k,:)=mu(k)+trandn(tl,tu)
 logpr=logpr+lnNpr(tl,tu)+0.5_dp*mu(k)**2-mu(k)*z(k,:)
 end do
end subroutine

function mvrandn(lo,hi,sigma,n,mean) result(out)
 real(dp),intent(in)::lo(:),hi(:),sigma(:,:)
 integer,intent(in)::n
 real(dp),intent(in),optional::mean(:)
 real(dp),allocatable::out(:,:)
 integer::d,st,acc,m,i,j,k,tries
 real(dp),allocatable::l0(:),u0(:),lmat(:,:),l(:),u(:),lfull(:,:),x(:),mu(:),z(:,:),lp(:),cand(:,:)
 integer,allocatable::perm(:),invp(:)
 real(dp)::ps,ex
 d=size(lo)
 allocate(l0(d),u0(d))
 l0=lo
 u0=hi
 if(present(mean))then
 l0=l0-mean
 u0=u0-mean
 end if
 if(d==1)then
 allocate(out(n,1))
 out(:,1)=sqrt(sigma(1,1))*trandn(spread(l0(1)/sqrt(sigma(1,1)),1,n),spread(u0(1)/sqrt(sigma(1,1)),1,n))
 if(present(mean))out(:,1)=out(:,1)+mean(1)
 return
 end if
 call prep(sigma,l0,u0,lmat,l,u,lfull,perm,st)
 if(st/=0)error stop 'mvrandn: Cholesky failure'
 call solve_normal_tilt(lmat,l,u,x,mu,st,.false.)
 ps=psy_eval(x,lmat,l,u,mu)
 allocate(out(n,d),z(d,n),lp(n),invp(d))
 acc=0
 tries=0
 do while(acc<n)
 call proposal_normal(n,lmat,l,u,mu,z,lp)
 do j=1,n
 if(-log(max(runif1(),tiny(1.0_dp)))>ps-lp(j))then
 acc=acc+1
 out(acc,:)=matmul(lfull,z(:,j))
 if(acc==n)exit
 end if
 end do
 tries=tries+1
 if(tries>100000)error stop 'mvrandn: acceptance loop failed'
 end do
 do i=1,d
 invp(perm(i))=i
 end do
 out=out(:,invp)
 if(present(mean)) then
  do i=1,n
  out(i,:)=out(i,:)+mean
  end do
 end if
end function

subroutine proposal_student(n,lmat,l,u,nu,muin,z,r,logpr)
 integer,intent(in)::n
 real(dp),intent(in)::lmat(:,:),l(:),u(:),nu,muin(:)
 real(dp),intent(out)::z(:,:),r(:),logpr(:)
 integer::d,k,j
 real(dp)::mu(size(l)),eta,inf,c,t1(1),col(n),tl(n),tu(n)
 d=size(l)
 mu=muin
 eta=mu(d)
 mu(d)=0
 inf=ieee_value(1.0_dp,ieee_positive_inf)
 t1=lnNpr([-eta],[inf])
 c=.5_dp*log(2*pi)-log_gamma(nu/2)-(nu/2-1)*log(2.0_dp)+t1(1)+.5_dp*eta*eta
 r=eta+trandn(spread(-eta,1,n),spread(inf,1,n))
 logpr=(nu-1)*log(r)-eta*r+c
 z=0
 do k=1,d
 col=0
 do j=1,k-1
 col=col+lmat(k,j)*z(j,:)
 end do
 tl=r*l(k)/sqrt(nu)-mu(k)-col
 tu=r*u(k)/sqrt(nu)-mu(k)-col
 z(k,:)=mu(k)+trandn(tl,tu)
 logpr=logpr+lnNpr(tl,tu)+.5_dp*mu(k)**2-mu(k)*z(k,:)
 end do
end subroutine

function tregress(n, lo, hi, sigma, nu) result(res)
 integer, intent(in) :: n
 real(dp), intent(in) :: lo(:), hi(:), sigma(:,:), nu
 type(tregress_result) :: res
 integer :: d, st, acc, i, j, tries
 real(dp), allocatable :: lmat(:,:), l(:), u(:), lfull(:,:), x(:), mu(:)
 real(dp), allocatable :: z(:,:), r(:), lp(:)
 integer, allocatable :: perm(:), invp(:)
 real(dp) :: ps

 d = size(lo)
 call prep(sigma, lo, hi, lmat, l, u, lfull, perm, st)
 if (st /= 0) then
   res%status = st
   return
 end if
 call solve_student_tilt(lmat, l, u, nu, x, mu, st)
 res%status = st
 ps = psy_t_eval(x, lmat, l, u, nu, mu)
 allocate(res%r(n), res%z(n,d), z(d,n), r(n), lp(n), invp(d))
 acc = 0
 tries = 0
 do while (acc < n)
   call proposal_student(n, lmat, l, u, nu, mu, z, r, lp)
   do j = 1, n
     if (-log(max(runif1(), tiny(1.0_dp))) > ps - lp(j)) then
       acc = acc + 1
       res%r(acc) = r(j)
       res%z(acc,:) = matmul(lfull, z(:,j))
       if (acc == n) exit
     end if
   end do
   tries = tries + 1
   if (tries > 100000) then
     res%status = max(res%status, 300)
     exit
   end if
 end do
 if (acc < n) then
   if (acc == 0) then
     res%r = 0.0_dp
     res%z = 0.0_dp
     return
   end if
   res%r(acc+1:) = 0.0_dp
   res%z(acc+1:,:) = 0.0_dp
 end if
 do i = 1, d
   invp(perm(i)) = i
 end do
 res%z = res%z(:,invp)
end function tregress

function mvrandt(lo,hi,sigma,nu,n,mean) result(out)
 real(dp),intent(in)::lo(:),hi(:),sigma(:,:),nu
 integer,intent(in)::n
 real(dp),intent(in),optional::mean(:)
 real(dp),allocatable::out(:,:)
 integer::d,st,acc,i,j,tries
 real(dp),allocatable::l0(:),u0(:),lmat(:,:),l(:),u(:),lfull(:,:),x(:),mu(:),z(:,:),r(:),lp(:)
 integer,allocatable::perm(:),invp(:)
 real(dp)::ps,p0,p1,uu,sdv
 d=size(lo)
 allocate(l0(d),u0(d))
 l0=lo
 u0=hi
 if(present(mean))then
 l0=l0-mean
 u0=u0-mean
 end if
 if(d==1)then
 allocate(out(n,1))
 sdv=sqrt(sigma(1,1))
 do i=1,n
 p0=pt(l0(1)/sdv,nu)
 p1=pt(u0(1)/sdv,nu)
 uu=p0+runif1()*(p1-p0)
 out(i,1)=sdv*qt(uu,nu)
 end do
 if(present(mean))out(:,1)=out(:,1)+mean(1)
 return
 end if
 call prep(sigma,l0,u0,lmat,l,u,lfull,perm,st)
 if(st/=0)error stop 'mvrandt: Cholesky failure'
 call solve_student_tilt(lmat,l,u,nu,x,mu,st)
 ps=psy_t_eval(x,lmat,l,u,nu,mu)
 allocate(out(n,d),z(d,n),r(n),lp(n),invp(d))
 acc=0
 tries=0
 do while(acc<n)
 call proposal_student(n,lmat,l,u,nu,mu,z,r,lp)
 do j=1,n
 if(-log(max(runif1(),tiny(1.0_dp)))>ps-lp(j))then
 acc=acc+1
 out(acc,:)=sqrt(nu)*matmul(lfull,z(:,j))/r(j)
 if(acc==n)exit
 end if
 end do
 tries=tries+1
 if(tries>100000)error stop 'mvrandt: acceptance loop failed'
 end do
 do i=1,d
 invp(perm(i))=i
 end do
 out=out(:,invp)
 if(present(mean)) then
  do i=1,n
  out(i,:)=out(i,:)+mean
  end do
 end if
end function
end module
