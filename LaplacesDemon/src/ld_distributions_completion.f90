module ld_distributions_completion
use ld_kinds, only: dp, pi
use ld_random, only: rand_uniform, rand_normal, rand_chisq, rand_dirichlet, rand_mvn, rand_mv_t, rand_categorical
use ld_linalg, only: inverse_spd, chol_lower, logdet_spd
use ld_distributions, only: normal_logpdf, normal_cdf, normal_quantile, dbern, dcat, dlaplace, plaplace, qlaplace, rlaplace
use ld_distributions, only: dhalfnorm, dhalfcauchy, ddirichlet, dgpd, dmvn, dmvt, dst, dstp, dwishart, dinvwishart
use ld_distributions, only: dmatrixnorm
use ld_distributions_extra, only: dhalft, dpe, dsdlaplace, psdlaplace, qsdlaplace, dslaplace, pslaplace, qslaplace, rslaplace
use ld_distributions_catalog, only: dmvc, rmvc, dmvl, rmvl, dmvpe, rmvpe, dmvpolya, dhorseshoe, rhorseshoe
use ld_distributions_catalog, only: dhuangwand, rhuangwand, dlasso, dlaplace_mixture, rlaplace_mixture
use ld_distributions_catalog, only: dnormal_mixture, rnormal_mixture, dnorminvwishart, dnormwishart, dyangberger, dzellner
use ld_distributions_catalog, only: dmvnormal_precision, dmvstudent_precision, dmvpower_precision, stick_density, rstick
implicit none
private
public :: pbern, qbern, rbern, qcat, rcat, rdirichlet, rgpd
public :: rhalfcauchy, rhalfnorm, phalft, qhalft, dhs, rhs
public :: rinvwishart, dinvwishartc, rinvwishartc, rwishart, dwishartc, rwishartc
public :: dlaplacep, plaplacep, qlaplacep, rlaplacep, dlaplacem, plaplacem, rlaplacem, rlasso
public :: dlnormp, plnormp, qlnormp, rlnormp, rmatrixnorm
public :: dmvcc, rmvcc, dmvcp, rmvcp, dmvcpc, rmvcpc
public :: dmvlc, rmvlc, rmvn, dmvnc, rmvnc, dmvnp, rmvnp, dmvnpc, rmvnpc
public :: rmvpolya, dmvpec, rmvpec, rmvt, dmvtc, rmvtc, dmvtp, rmvtp, dmvtpc, rmvtpc
public :: dnormm, pnormm, rnormm, dnormp, pnormp, qnormp, rnormp, dnormv, pnormv, qnormv, rnormv
public :: rnorminvwishart, rnormwishart, ppe, qpe, dsiw, rsiw, rsdlaplace
public :: dstick, pst, qst, rst, pstp, qstp, rstp
public :: dtrunc, ptrunc, qtrunc, rtrunc, extrunc, vartrunc, dyangbergerc, rzellner
public :: dcrmrf, rcrmrf, dhuangwandc, rhuangwandc

abstract interface
   function scalar_pdf_iface_c(x) result(v)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: v
   end function scalar_pdf_iface_c
   function scalar_cdf_iface_c(x) result(v)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: v
   end function scalar_cdf_iface_c
   function scalar_quantile_iface_c(p) result(v)
      import dp
      real(dp), intent(in) :: p
      real(dp) :: v
   end function scalar_quantile_iface_c
end interface

contains

function dcrmrf(x,alpha,omega,log_density) result(v)
   real(dp),intent(in)::x(:),alpha(:),omega(:,:)
   logical,intent(in),optional::log_density
   real(dp)::v,oinv(size(x),size(x)),z
   integer::i,info
   logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(size(alpha)/=size(x) .or. size(omega,1)/=size(x) .or. size(omega,2)/=size(x)) then
      v=merge(-huge(1.0_dp),0.0_dp,lg); return
   end if
   call inverse_spd(omega,oinv,info)
   if(info/=0) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=-0.5_dp*dot_product(x,matmul(oinv,x))
   do i=1,size(x)
      z=x(i)+alpha(i)
      if(z>0.0_dp) then
         v=v+z+log(1.0_dp+exp(-z))
      else
         v=v+log(1.0_dp+exp(z))
      end if
   end do
   if(.not.lg) v=exp(v)
end function dcrmrf

integer function rcrmrf(alpha,omega) result(outcome)
   real(dp),intent(in)::alpha(:),omega(:,:)
   real(dp)::mean(size(alpha)),z(size(alpha)),prob(size(alpha)),tot
   integer::j,k
   k=size(alpha)
   do j=1,k
      mean=omega(:,j)
      call rand_mvn(mean,omega,z)
      prob(j)=dcrmrf(z,alpha,omega)
   end do
   tot=sum(prob)
   if(tot<=0.0_dp .or. .not.(tot<huge(1.0_dp))) then
      prob=1.0_dp/real(k,dp)
   else
      prob=prob/tot
   end if
   outcome=rand_categorical(prob)
end function rcrmrf

function dhuangwandc(u,nu,alatent,ascale,log_density) result(v)
   real(dp),intent(in)::u(:,:),nu,alatent(:),ascale(:)
   logical,intent(in),optional::log_density
   real(dp)::v,x(size(u,1),size(u,2))
   x=matmul(transpose(u),u)
   if(present(log_density)) then
      v=dhuangwand(x,nu,alatent,ascale,log_density)
   else
      v=dhuangwand(x,nu,alatent,ascale)
   end if
end function dhuangwandc

subroutine rhuangwandc(nu,alatent,ascale,u)
   real(dp),intent(in)::nu,alatent(:),ascale(:)
   real(dp),intent(out)::u(:,:)
   real(dp)::x(size(u,1),size(u,2)),l(size(u,1),size(u,2))
   integer::info
   call rhuangwand(nu,alatent,ascale,x)
   call chol_lower(x,l,info)
   if(info/=0) then; u=0.0_dp; else; u=transpose(l); end if
end subroutine rhuangwandc

function dlaplacem(x,p,location,scale,log_density) result(v)
   real(dp),intent(in)::x,p(:),location(:),scale(:)
   logical,intent(in),optional::log_density
   real(dp)::v
   if(present(log_density)) then
      v=dlaplace_mixture(x,p,location,scale,log_density)
   else
      v=dlaplace_mixture(x,p,location,scale)
   end if
end function dlaplacem

pure function pbern(q,prob) result(v)
   integer,intent(in)::q; real(dp),intent(in)::prob; real(dp)::v
   if(q<0) then; v=0.0_dp; else if(q==0) then; v=1.0_dp-prob; else; v=1.0_dp; end if
end function pbern
pure function qbern(p,prob) result(v)
   real(dp),intent(in)::p,prob; integer::v
   v=merge(0,1,p<=1.0_dp-prob)
end function qbern
integer function rbern(prob) result(v)
   real(dp),intent(in)::prob; v=merge(1,0,rand_uniform()<prob)
end function rbern
integer function qcat(pr,p) result(v)
   real(dp),intent(in)::pr,p(:); real(dp)::c; integer::i
   c=0.0_dp; v=size(p)
   do i=1,size(p); c=c+p(i); if(pr<=c) then; v=i; return; end if; end do
end function qcat
integer function rcat(p) result(v)
   real(dp),intent(in)::p(:); v=rand_categorical(p)
end function rcat
subroutine rdirichlet(alpha,x)
   real(dp),intent(in)::alpha(:); real(dp),intent(out)::x(:); call rand_dirichlet(alpha,x)
end subroutine rdirichlet
function rgpd(mu,sigma,xi) result(x)
   real(dp),intent(in)::mu,sigma,xi; real(dp)::x,u
   u=max(rand_uniform(),0.001_dp)
   if(abs(xi)<1.0e-14_dp) then; x=mu-sigma*log(u); else; x=mu+sigma*(u**(-xi)-1.0_dp)/xi; end if
end function rgpd
function rhalfcauchy(scale) result(x)
   real(dp),intent(in)::scale; real(dp)::x; x=scale*tan(0.5_dp*pi*rand_uniform())
end function rhalfcauchy
function rhalfnorm(scale) result(x)
   real(dp),intent(in)::scale; real(dp)::x; x=abs(sqrt(pi/2.0_dp)/scale*rand_normal())
end function rhalfnorm
function phalft(q,scale,nu) result(v)
   real(dp),intent(in)::q,scale,nu; real(dp)::v
   if(q<=0.0_dp) then; v=0.0_dp; else; v=2.0_dp*student_t_cdf(q/scale,nu)-1.0_dp; end if
end function phalft
function qhalft(p,scale,nu) result(v)
   real(dp),intent(in)::p,scale,nu; real(dp)::v
   v=scale*student_t_quantile(0.5_dp*(p+1.0_dp),nu)
end function qhalft
function dhs(x,lambda,tau,log_density) result(v)
   real(dp),intent(in)::x,lambda,tau; logical,intent(in),optional::log_density; real(dp)::v
   if(present(log_density)) then; v=dhorseshoe(x,lambda,tau,log_density); else; v=dhorseshoe(x,lambda,tau); end if
end function dhs
function rhs(lambda,tau) result(x)
   real(dp),intent(in)::lambda,tau; real(dp)::x; x=rhorseshoe(lambda,tau)
end function rhs

subroutine rwishart(nu,s,x)
   real(dp),intent(in)::nu,s(:,:); real(dp),intent(out)::x(:,:)
   real(dp)::l(size(s,1),size(s,2)),a(size(s,1),size(s,2)),b(size(s,1),size(s,2)); integer::i,j,info,p
   p=size(s,1); call chol_lower(s,l,info); if(info/=0) then; x=0.0_dp; return; end if
   a=0.0_dp
   do i=1,p
      a(i,i)=sqrt(rand_chisq(nu-real(i-1,dp)))
      do j=1,i-1; a(i,j)=rand_normal(); end do
   end do
   b=matmul(l,a); x=matmul(b,transpose(b))
end subroutine rwishart
subroutine rinvwishart(nu,s,x)
   real(dp),intent(in)::nu,s(:,:); real(dp),intent(out)::x(:,:)
   real(dp)::sinv(size(s,1),size(s,2)),w(size(s,1),size(s,2)); integer::info
   call inverse_spd(s,sinv,info); if(info/=0) then; x=0.0_dp; return; end if
   call rwishart(nu,sinv,w); call inverse_spd(w,x,info)
end subroutine rinvwishart
function dinvwishartc(u,nu,s,log_density) result(v)
   real(dp),intent(in)::u(:,:),nu,s(:,:); logical,intent(in),optional::log_density; real(dp)::v,x(size(u,1),size(u,2))
   x=matmul(transpose(u),u)
   if(present(log_density)) then; v=dinvwishart(x,nu,s,log_density); else; v=dinvwishart(x,nu,s); end if
end function dinvwishartc
subroutine rinvwishartc(nu,s,u)
   real(dp),intent(in)::nu,s(:,:)
   real(dp),intent(out)::u(:,:)
   real(dp)::x(size(s,1),size(s,2)),l(size(s,1),size(s,2))
   integer::info
   call rinvwishart(nu,s,x); call chol_lower(x,l,info); u=transpose(l)
end subroutine rinvwishartc
function dwishartc(u,nu,s,log_density) result(v)
   real(dp),intent(in)::u(:,:),nu,s(:,:); logical,intent(in),optional::log_density; real(dp)::v,x(size(u,1),size(u,2))
   x=matmul(transpose(u),u)
   if(present(log_density)) then; v=dwishart(x,nu,s,log_density); else; v=dwishart(x,nu,s); end if
end function dwishartc
subroutine rwishartc(nu,s,u)
   real(dp),intent(in)::nu,s(:,:)
   real(dp),intent(out)::u(:,:)
   real(dp)::x(size(s,1),size(s,2)),l(size(s,1),size(s,2))
   integer::info
   call rwishart(nu,s,x); call chol_lower(x,l,info); u=transpose(l)
end subroutine rwishartc

function dlaplacep(x,mu,tau,log_density) result(v)
   real(dp),intent(in)::x,mu,tau; logical,intent(in),optional::log_density; real(dp)::v
   if(present(log_density)) then; v=dlaplace(x,mu,1.0_dp/tau,log_density); else; v=dlaplace(x,mu,1.0_dp/tau); end if
end function dlaplacep
pure function plaplacep(q,mu,tau) result(v)
   real(dp),intent(in)::q,mu,tau; real(dp)::v; v=plaplace(q,mu,1.0_dp/tau)
end function plaplacep
pure function qlaplacep(p,mu,tau) result(v)
   real(dp),intent(in)::p,mu,tau; real(dp)::v; v=qlaplace(p,mu,1.0_dp/tau)
end function qlaplacep
function rlaplacep(mu,tau) result(v)
   real(dp),intent(in)::mu,tau; real(dp)::v; v=rlaplace(mu,1.0_dp/tau)
end function rlaplacep
function plaplacem(q,p,location,scale) result(v)
   real(dp),intent(in)::q,p(:),location(:),scale(:); real(dp)::v; integer::i
   v=0.0_dp; do i=1,size(p); v=v+p(i)*plaplace(q,location(i),scale(i)); end do
end function plaplacem
function rlaplacem(p,location,scale) result(v)
   real(dp),intent(in)::p(:),location(:),scale(:); real(dp)::v; v=rlaplace_mixture(p,location,scale)
end function rlaplacem
subroutine rlasso(sigma,tau,x)
   real(dp),intent(in)::sigma,tau(:); real(dp),intent(out)::x(:); integer::i
   do i=1,size(x); x(i)=sigma*tau(i)*rand_normal(); end do
end subroutine rlasso

function dlnormp(x,mu,tau,log_density) result(v)
   real(dp),intent(in)::x,mu,tau; logical,intent(in),optional::log_density; real(dp)::v
   logical::lg; lg=.false.; if(present(log_density)) lg=log_density
   if(x<=0.0_dp .or. tau<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=0.5_dp*(log(tau)-log(2.0_dp*pi))-log(x)-0.5_dp*tau*(log(x)-mu)**2; if(.not.lg) v=exp(v)
end function dlnormp
function plnormp(q,mu,tau) result(v)
   real(dp),intent(in)::q,mu,tau; real(dp)::v
   if(q<=0.0_dp) then; v=0.0_dp; else; v=normal_cdf(log(q),mu,1.0_dp/sqrt(tau)); end if
end function plnormp
function qlnormp(p,mu,tau) result(v)
   real(dp),intent(in)::p,mu,tau; real(dp)::v; v=exp(normal_quantile(p,mu,1.0_dp/sqrt(tau)))
end function qlnormp
function rlnormp(mu,tau) result(v)
   real(dp),intent(in)::mu,tau; real(dp)::v; v=exp(mu+rand_normal()/sqrt(tau))
end function rlnormp
subroutine rmatrixnorm(m,u,v,x)
   real(dp),intent(in)::m(:,:),u(:,:),v(:,:); real(dp),intent(out)::x(:,:)
   real(dp)::lu(size(u,1),size(u,2)),lv(size(v,1),size(v,2)),z(size(m,1),size(m,2)); integer::i,j,info
   call chol_lower(u,lu,info)
   if(info/=0) then
   x=m
   return
   end if
   call chol_lower(v,lv,info)
   if(info/=0) then
   x=m
   return
   end if
   do i=1,size(m,1); do j=1,size(m,2); z(i,j)=rand_normal(); end do; end do
   x=m+matmul(lu,matmul(z,transpose(lv)))
end subroutine rmatrixnorm

function dmvcc(x,mu,u,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),u(:,:); logical,intent(in),optional::log_density; real(dp)::v,s(size(u,1),size(u,2))
   s=matmul(transpose(u),u); if(present(log_density)) then; v=dmvc(x,mu,s,log_density); else; v=dmvc(x,mu,s); end if
end function dmvcc
subroutine rmvcc(mu,u,x)
   real(dp),intent(in)::mu(:),u(:,:)
   real(dp),intent(out)::x(:)
   real(dp)::s(size(u,1),size(u,2))
   s=matmul(transpose(u),u)
   call rmvc(mu,s,x)
end subroutine rmvcc
function dmvcp(x,mu,omega,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),omega(:,:)
   logical,intent(in),optional::log_density
   real(dp)::v,s(size(omega,1),size(omega,2))
   integer::info
   call inverse_spd(omega,s,info); if(info/=0) then; v=0.0_dp; return; end if
   if(present(log_density)) then; v=dmvc(x,mu,s,log_density); else; v=dmvc(x,mu,s); end if
end function dmvcp
subroutine rmvcp(mu,omega,x)
   real(dp),intent(in)::mu(:),omega(:,:); real(dp),intent(out)::x(:); real(dp)::s(size(omega,1),size(omega,2)); integer::info
   call inverse_spd(omega,s,info); if(info/=0) then; x=mu; else; call rmvc(mu,s,x); end if
end subroutine rmvcp
function dmvcpc(x,mu,u,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),u(:,:); logical,intent(in),optional::log_density; real(dp)::v,omega(size(u,1),size(u,2))
   omega=matmul(transpose(u),u); if(present(log_density)) then; v=dmvcp(x,mu,omega,log_density); else; v=dmvcp(x,mu,omega); end if
end function dmvcpc
subroutine rmvcpc(mu,u,x)
   real(dp),intent(in)::mu(:),u(:,:)
   real(dp),intent(out)::x(:)
   real(dp)::omega(size(u,1),size(u,2))
   omega=matmul(transpose(u),u)
   call rmvcp(mu,omega,x)
end subroutine rmvcpc
function dmvlc(x,mu,u,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),u(:,:); logical,intent(in),optional::log_density; real(dp)::v,s(size(u,1),size(u,2))
   s=matmul(transpose(u),u); if(present(log_density)) then; v=dmvl(x,mu,s,log_density); else; v=dmvl(x,mu,s); end if
end function dmvlc
subroutine rmvlc(mu,u,x)
   real(dp),intent(in)::mu(:),u(:,:)
   real(dp),intent(out)::x(:)
   real(dp)::s(size(u,1),size(u,2))
   s=matmul(transpose(u),u)
   call rmvl(mu,s,x)
end subroutine rmvlc
subroutine rmvn(mu,sigma,x)
   real(dp),intent(in)::mu(:),sigma(:,:); real(dp),intent(out)::x(:); integer::info; call rand_mvn(mu,sigma,x,info)
end subroutine rmvn
function dmvnc(x,mu,u,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),u(:,:); logical,intent(in),optional::log_density; real(dp)::v,s(size(u,1),size(u,2))
   s=matmul(transpose(u),u); if(present(log_density)) then; v=dmvn(x,mu,s,log_density); else; v=dmvn(x,mu,s); end if
end function dmvnc
subroutine rmvnc(mu,u,x)
   real(dp),intent(in)::mu(:),u(:,:)
   real(dp),intent(out)::x(:)
   real(dp)::s(size(u,1),size(u,2))
   s=matmul(transpose(u),u)
   call rmvn(mu,s,x)
end subroutine rmvnc
function dmvnp(x,mu,omega,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),omega(:,:); logical,intent(in),optional::log_density; real(dp)::v
   if(present(log_density)) then; v=dmvnormal_precision(x,mu,omega,log_density); else; v=dmvnormal_precision(x,mu,omega); end if
end function dmvnp
subroutine rmvnp(mu,omega,x)
   real(dp),intent(in)::mu(:),omega(:,:); real(dp),intent(out)::x(:); real(dp)::s(size(omega,1),size(omega,2)); integer::info
   call inverse_spd(omega,s,info); if(info/=0) then; x=mu; else; call rmvn(mu,s,x); end if
end subroutine rmvnp
function dmvnpc(x,mu,u,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),u(:,:); logical,intent(in),optional::log_density; real(dp)::v,omega(size(u,1),size(u,2))
   omega=matmul(transpose(u),u); if(present(log_density)) then; v=dmvnp(x,mu,omega,log_density); else; v=dmvnp(x,mu,omega); end if
end function dmvnpc
subroutine rmvnpc(mu,u,x)
   real(dp),intent(in)::mu(:),u(:,:)
   real(dp),intent(out)::x(:)
   real(dp)::omega(size(u,1),size(u,2))
   omega=matmul(transpose(u),u)
   call rmvnp(mu,omega,x)
end subroutine rmvnpc
integer function rmvpolya(alpha) result(x)
   real(dp),intent(in)::alpha(:); real(dp)::p(size(alpha)); call rand_dirichlet(alpha,p); x=rand_categorical(p)
end function rmvpolya
function dmvpec(x,mu,u,kappa,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),u(:,:),kappa; logical,intent(in),optional::log_density; real(dp)::v,s(size(u,1),size(u,2))
   s=matmul(transpose(u),u); if(present(log_density)) then; v=dmvpe(x,mu,s,kappa,log_density); else; v=dmvpe(x,mu,s,kappa); end if
end function dmvpec
subroutine rmvpec(mu,u,kappa,x)
   real(dp),intent(in)::mu(:),u(:,:),kappa
   real(dp),intent(out)::x(:)
   real(dp)::s(size(u,1),size(u,2))
   s=matmul(transpose(u),u)
   call rmvpe(mu,s,kappa,x)
end subroutine rmvpec
subroutine rmvt(mu,s,df,x)
   real(dp),intent(in)::mu(:),s(:,:),df; real(dp),intent(out)::x(:); integer::info; call rand_mv_t(mu,s,df,x,info)
end subroutine rmvt
function dmvtc(x,mu,u,df,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),u(:,:),df; logical,intent(in),optional::log_density; real(dp)::v,s(size(u,1),size(u,2))
   s=matmul(transpose(u),u); if(present(log_density)) then; v=dmvt(x,mu,s,df,log_density); else; v=dmvt(x,mu,s,df); end if
end function dmvtc
subroutine rmvtc(mu,u,df,x)
   real(dp),intent(in)::mu(:),u(:,:),df
   real(dp),intent(out)::x(:)
   real(dp)::s(size(u,1),size(u,2))
   s=matmul(transpose(u),u)
   call rmvt(mu,s,df,x)
end subroutine rmvtc
function dmvtp(x,mu,omega,nu,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),omega(:,:),nu; logical,intent(in),optional::log_density; real(dp)::v
   if(present(log_density)) then
   v=dmvstudent_precision(x,mu,omega,nu,log_density)
   else
   v=dmvstudent_precision(x,mu,omega,nu)
   end if
end function dmvtp
subroutine rmvtp(mu,omega,nu,x)
   real(dp),intent(in)::mu(:),omega(:,:),nu; real(dp),intent(out)::x(:); real(dp)::s(size(omega,1),size(omega,2)); integer::info
   call inverse_spd(omega,s,info); if(info/=0) then; x=mu; else; call rmvt(mu,s,nu,x); end if
end subroutine rmvtp
function dmvtpc(x,mu,u,nu,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),u(:,:),nu; logical,intent(in),optional::log_density; real(dp)::v,omega(size(u,1),size(u,2))
   omega=matmul(transpose(u),u)
   if(present(log_density)) then
   v=dmvtp(x,mu,omega,nu,log_density)
   else
   v=dmvtp(x,mu,omega,nu)
   end if
end function dmvtpc
subroutine rmvtpc(mu,u,nu,x)
   real(dp),intent(in)::mu(:),u(:,:),nu
   real(dp),intent(out)::x(:)
   real(dp)::omega(size(u,1),size(u,2))
   omega=matmul(transpose(u),u)
   call rmvtp(mu,omega,nu,x)
end subroutine rmvtpc

function dnormm(x,p,mu,sigma,log_density) result(v)
   real(dp),intent(in)::x,p(:),mu(:),sigma(:); logical,intent(in),optional::log_density; real(dp)::v
   if(present(log_density)) then; v=dnormal_mixture(x,p,mu,sigma,log_density); else; v=dnormal_mixture(x,p,mu,sigma); end if
end function dnormm
function pnormm(q,p,mu,sigma) result(v)
   real(dp),intent(in)::q,p(:),mu(:),sigma(:); real(dp)::v; integer::i
   v=0.0_dp; do i=1,size(p); v=v+p(i)*normal_cdf(q,mu(i),sigma(i)); end do
end function pnormm
function rnormm(p,mu,sigma) result(v)
   real(dp),intent(in)::p(:),mu(:),sigma(:); real(dp)::v; v=rnormal_mixture(p,mu,sigma)
end function rnormm
function dnormp(x,mean,prec,log_density) result(v)
   real(dp),intent(in)::x,mean,prec; logical,intent(in),optional::log_density; real(dp)::v
   if(present(log_density)) then
   v=merge(normal_logpdf(x,mean,1.0_dp/sqrt(prec)),exp(normal_logpdf(x,mean,1.0_dp/sqrt(prec))),log_density)
   else; v=exp(normal_logpdf(x,mean,1.0_dp/sqrt(prec))); end if
end function dnormp
pure function pnormp(q,mean,prec) result(v)
   real(dp),intent(in)::q,mean,prec; real(dp)::v; v=normal_cdf(q,mean,1.0_dp/sqrt(prec))
end function pnormp
pure function qnormp(p,mean,prec) result(v)
   real(dp),intent(in)::p,mean,prec; real(dp)::v; v=normal_quantile(p,mean,1.0_dp/sqrt(prec))
end function qnormp
function rnormp(mean,prec) result(v)
   real(dp),intent(in)::mean,prec; real(dp)::v; v=mean+rand_normal()/sqrt(prec)
end function rnormp
function dnormv(x,mean,var,log_density) result(v)
   real(dp),intent(in)::x,mean,var; logical,intent(in),optional::log_density; real(dp)::v
   if(present(log_density)) then; v=merge(normal_logpdf(x,mean,sqrt(var)),exp(normal_logpdf(x,mean,sqrt(var))),log_density)
   else; v=exp(normal_logpdf(x,mean,sqrt(var))); end if
end function dnormv
pure function pnormv(q,mean,var) result(v)
   real(dp),intent(in)::q,mean,var; real(dp)::v; v=normal_cdf(q,mean,sqrt(var))
end function pnormv
pure function qnormv(p,mean,var) result(v)
   real(dp),intent(in)::p,mean,var; real(dp)::v; v=normal_quantile(p,mean,sqrt(var))
end function qnormv
function rnormv(mean,var) result(v)
   real(dp),intent(in)::mean,var; real(dp)::v; v=mean+sqrt(var)*rand_normal()
end function rnormv
subroutine rnorminvwishart(mu0,lambda,s,nu,mu,sigma)
   real(dp),intent(in)::mu0(:),lambda,s(:,:),nu; real(dp),intent(out)::mu(:),sigma(:,:); integer::info
   call rinvwishart(nu,s,sigma); call rand_mvn(mu0,sigma/lambda,mu,info)
end subroutine rnorminvwishart
subroutine rnormwishart(mu0,lambda,s,nu,mu,omega)
   real(dp),intent(in)::mu0(:),lambda,s(:,:),nu; real(dp),intent(out)::mu(:),omega(:,:)
   real(dp)::cov(size(s,1),size(s,2)); integer::info
   call rwishart(nu,s,omega); call inverse_spd(lambda*omega,cov,info); call rand_mvn(mu0,cov,mu,info)
end subroutine rnormwishart

function ppe(q,mu,sigma,kappa) result(v)
   real(dp),intent(in)::q,mu,sigma,kappa; real(dp)::v
   v=numeric_cdf_pe(q,mu,sigma,kappa)
end function ppe
function qpe(p,mu,sigma,kappa) result(v)
   real(dp),intent(in)::p,mu,sigma,kappa; real(dp)::v,lo,hi,mid; integer::i
   lo=mu-20.0_dp*sigma; hi=mu+20.0_dp*sigma
   do i=1,100; mid=0.5_dp*(lo+hi); if(ppe(mid,mu,sigma,kappa)<p) then; lo=mid; else; hi=mid; end if; end do
   v=0.5_dp*(lo+hi)
end function qpe
function dsiw(q,nu,s,zeta,mu,delta,log_density) result(v)
   real(dp),intent(in)::q(:,:),nu,s(:,:),zeta(:),mu(:),delta(:); logical,intent(in),optional::log_density
   real(dp)::v,lp; integer::i
   lp=dinvwishart(q,nu,s,.true.)
   do i=1,size(zeta); lp=lp+normal_logpdf(log(zeta(i)),mu(i),sqrt(delta(i))); end do
   if(present(log_density)) then; v=merge(lp,exp(lp),log_density); else; v=exp(lp); end if
end function dsiw
subroutine rsiw(nu,s,mu,delta,x)
   real(dp),intent(in)::nu,s(:,:),mu(:),delta(:); real(dp),intent(out)::x(:,:)
   real(dp)::q(size(s,1),size(s,2)),z(size(mu)),d(size(s,1),size(s,2)); integer::i
   call rinvwishart(nu,s,q); d=0.0_dp; do i=1,size(mu); z(i)=exp(mu(i)+sqrt(delta(i))*rand_normal()); d(i,i)=z(i); end do
   x=matmul(d,matmul(q,d))
end subroutine rsiw
function rsdlaplace(p,q) result(x)
   real(dp),intent(in)::p,q; integer::x
   x=q_sdlaplace(rand_uniform(),p,q)
end function rsdlaplace
function dstick(theta,gamma0,log_density) result(v)
   real(dp),intent(in)::theta(:),gamma0; logical,intent(in),optional::log_density; real(dp)::v
   if(present(log_density)) then; v=stick_density(theta,gamma0,log_density); else; v=stick_density(theta,gamma0); end if
end function dstick
function pst(q,mu,sigma,nu) result(v)
   real(dp),intent(in)::q,mu,sigma,nu; real(dp)::v; v=student_t_cdf((q-mu)/sigma,nu)
end function pst
function qst(p,mu,sigma,nu) result(v)
   real(dp),intent(in)::p,mu,sigma,nu; real(dp)::v; v=mu+sigma*student_t_quantile(p,nu)
end function qst
function rst(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu; real(dp)::v; v=mu+sigma*student_t_random(nu)
end function rst
function pstp(q,mu,tau,nu) result(v)
   real(dp),intent(in)::q,mu,tau,nu; real(dp)::v; v=pst(q,mu,1.0_dp/sqrt(tau),nu)
end function pstp
function qstp(p,mu,tau,nu) result(v)
   real(dp),intent(in)::p,mu,tau,nu; real(dp)::v; v=qst(p,mu,1.0_dp/sqrt(tau),nu)
end function qstp
function rstp(mu,tau,nu) result(v)
   real(dp),intent(in)::mu,tau,nu; real(dp)::v; v=rst(mu,1.0_dp/sqrt(tau),nu)
end function rstp

function dtrunc(x,pdf,cdf,a,b,log_density) result(v)
   real(dp),intent(in)::x,a,b; procedure(scalar_pdf_iface_c)::pdf; procedure(scalar_cdf_iface_c)::cdf
   logical,intent(in),optional::log_density; real(dp)::v,z; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density; z=cdf(b)-cdf(a)
   if(x<a .or. x>b .or. z<=0.0_dp) then
   v=merge(-huge(1.0_dp),0.0_dp,lg)
   else if(lg) then
   v=log(max(pdf(x),tiny(1.0_dp)))-log(z)
   else
   v=pdf(x)/z
   end if
end function dtrunc
function ptrunc(x,cdf,a,b) result(v)
   real(dp),intent(in)::x,a,b; procedure(scalar_cdf_iface_c)::cdf; real(dp)::v,z,xx
   xx=min(max(x,a),b); z=cdf(b)-cdf(a); v=(cdf(xx)-cdf(a))/z
end function ptrunc
function qtrunc(p,cdf,quantile,a,b) result(v)
   real(dp),intent(in)::p,a,b; procedure(scalar_cdf_iface_c)::cdf; procedure(scalar_quantile_iface_c)::quantile; real(dp)::v
   v=quantile(cdf(a)+p*(cdf(b)-cdf(a)))
end function qtrunc
function rtrunc(cdf,quantile,a,b) result(v)
   real(dp),intent(in)::a,b; procedure(scalar_cdf_iface_c)::cdf; procedure(scalar_quantile_iface_c)::quantile; real(dp)::v
   v=qtrunc(rand_uniform(),cdf,quantile,a,b)
end function rtrunc
function extrunc(pdf,cdf,a,b) result(v)
   real(dp),intent(in)::a,b; procedure(scalar_pdf_iface_c)::pdf; procedure(scalar_cdf_iface_c)::cdf; real(dp)::v
   v=trunc_moment(pdf,cdf,a,b,1)
end function extrunc
function vartrunc(pdf,cdf,a,b) result(v)
   real(dp),intent(in)::a,b; procedure(scalar_pdf_iface_c)::pdf; procedure(scalar_cdf_iface_c)::cdf; real(dp)::v,m1,m2
   m1=trunc_moment(pdf,cdf,a,b,1); m2=trunc_moment(pdf,cdf,a,b,2); v=max(0.0_dp,m2-m1*m1)
end function vartrunc
function dyangbergerc(u,log_density) result(v)
   real(dp),intent(in)::u(:,:); logical,intent(in),optional::log_density; real(dp)::v,x(size(u,1),size(u,2))
   x=matmul(transpose(u),u); if(present(log_density)) then; v=dyangberger(x,log_density); else; v=dyangberger(x); end if
end function dyangbergerc
subroutine rzellner(g,sigma,xmat,beta)
   real(dp),intent(in)::g,sigma,xmat(:,:); real(dp),intent(out)::beta(:)
   real(dp)::xtx(size(beta),size(beta)),inv(size(beta),size(beta)); integer::info
   xtx=matmul(transpose(xmat),xmat)
   call inverse_spd(xtx,inv,info)
   if(info/=0) then
   beta=0.0_dp
   else
   call rand_mvn(0.0_dp*beta,g*sigma*sigma*inv,beta,info)
   end if
end subroutine rzellner

function student_t_random(nu) result(x)
   real(dp),intent(in)::nu; real(dp)::x; x=rand_normal()/sqrt(rand_chisq(nu)/nu)
end function student_t_random
function student_t_pdf(x,nu) result(v)
   real(dp),intent(in)::x,nu; real(dp)::v
   v=exp(log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)-0.5_dp*log(nu*pi)-0.5_dp*(nu+1.0_dp)*log(1.0_dp+x*x/nu))
end function student_t_pdf
function student_t_cdf(x,nu) result(v)
   real(dp),intent(in)::x,nu; real(dp)::v,a,h,s,t; integer,parameter::n=800; integer::i
   if(nu>1.0e6_dp) then; v=normal_cdf(x,0.0_dp,1.0_dp); return; end if
   if(x==0.0_dp) then; v=0.5_dp; return; end if
   a=abs(x); h=a/real(n,dp); s=student_t_pdf(0.0_dp,nu)+student_t_pdf(a,nu)
   do i=1,n-1; t=h*real(i,dp); s=s+merge(4.0_dp,2.0_dp,mod(i,2)==1)*student_t_pdf(t,nu); end do
   v=0.5_dp+sign(1.0_dp,x)*h*s/3.0_dp; v=min(max(v,0.0_dp),1.0_dp)
end function student_t_cdf
function student_t_quantile(p,nu) result(v)
   real(dp),intent(in)::p,nu; real(dp)::v,lo,hi,mid; integer::i
   if(nu>1.0e6_dp) then; v=normal_quantile(p,0.0_dp,1.0_dp); return; end if
   lo=-100.0_dp; hi=100.0_dp
   do i=1,100; mid=0.5_dp*(lo+hi); if(student_t_cdf(mid,nu)<p) then; lo=mid; else; hi=mid; end if; end do
   v=0.5_dp*(lo+hi)
end function student_t_quantile
function numeric_cdf_pe(q,mu,sigma,kappa) result(v)
   real(dp),intent(in)::q,mu,sigma,kappa; real(dp)::v,lo,hi,h,s,x; integer,parameter::n=1000; integer::i
   lo=mu-20.0_dp*sigma
   hi=min(q,mu+20.0_dp*sigma)
   if(q<=lo) then
   v=0.0_dp
   return
   end if
   if(q>=mu+20.0_dp*sigma) then
   v=1.0_dp
   return
   end if
   h=(hi-lo)/real(n,dp); s=dpe(lo,mu,sigma,kappa)+dpe(hi,mu,sigma,kappa)
   do i=1,n-1; x=lo+h*real(i,dp); s=s+merge(4.0_dp,2.0_dp,mod(i,2)==1)*dpe(x,mu,sigma,kappa); end do
   v=min(max(h*s/3.0_dp,0.0_dp),1.0_dp)
end function numeric_cdf_pe
integer function q_sdlaplace(u,p,q) result(x)
   real(dp),intent(in)::u,p,q
   x=qsdlaplace(u,p,q)
end function q_sdlaplace

function trunc_moment(pdf,cdf,a,b,power) result(v)
   procedure(scalar_pdf_iface_c) :: pdf
   procedure(scalar_cdf_iface_c) :: cdf
   real(dp),intent(in)::a,b
   integer,intent(in)::power
   real(dp)::v,lo,hi,h,s,x,z
   integer,parameter::n=1600
   integer::i
   lo=a; hi=b
   if(lo < -1.0e100_dp) lo=-20.0_dp
   if(hi >  1.0e100_dp) hi= 20.0_dp
   z=cdf(b)-cdf(a)
   if(z<=0.0_dp .or. hi<=lo) then; v=0.0_dp; return; end if
   h=(hi-lo)/real(n,dp)
   s=(lo**power)*pdf(lo)+(hi**power)*pdf(hi)
   do i=1,n-1
      x=lo+h*real(i,dp)
      s=s+merge(4.0_dp,2.0_dp,mod(i,2)==1)*(x**power)*pdf(x)
   end do
   v=h*s/(3.0_dp*z)
end function trunc_moment

end module ld_distributions_completion
