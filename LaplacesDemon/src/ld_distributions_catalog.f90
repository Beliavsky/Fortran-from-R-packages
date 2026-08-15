module ld_distributions_catalog
use ld_kinds, only: dp, pi, log2pi
use ld_random, only: rand_uniform, rand_normal, rand_gamma, rand_chisq, rand_mvn, rand_categorical
use ld_linalg, only: inverse_spd, logdet_spd, trace_matrix, chol_lower, make_positive_definite
use ld_distributions, only: normal_logpdf, normal_cdf, normal_quantile, dmvn, dmvt, dlaplace, rlaplace
use ld_distributions, only: dinvgamma, rinvgamma, dinvwishart, dwishart, log_multivariate_gamma
use ld_distributions_extra, only: dslaplace, rslaplace
implicit none
private
public :: daml, raml, dmvc, rmvc, dmvl, rmvl, dmvpe, rmvpe
public :: dmvpolya, dmatrixgamma, rinvmatrixgamma, dinvmatrixgamma, rmatrixgamma
public :: dinvbeta, rinvbeta, dhorseshoe, rhorseshoe, dhuangwand, rhuangwand
public :: dlasso, dlaplace_mixture, rlaplace_mixture, dnormal_mixture, rnormal_mixture
public :: dnormlaplace, rnormlaplace, dnorminvwishart, dnormwishart
public :: dyangberger, dhyperg, dzellner, dmvnormal_precision
public :: dmvstudent_precision, dmvpower_precision, stick_density, rstick
public :: dtrunc_generic, ptrunc_generic, qtrunc_generic, rtrunc_generic

abstract interface
   function scalar_pdf_iface(x) result(v)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: v
   end function scalar_pdf_iface
   function scalar_cdf_iface(x) result(v)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: v
   end function scalar_cdf_iface
   function scalar_quantile_iface(p) result(v)
      import dp
      real(dp), intent(in) :: p
      real(dp) :: v
   end function scalar_quantile_iface
end interface

contains

function daml(x,mu,sigma,log_density) result(v)
   real(dp), intent(in) :: x(:),mu(:),sigma(:,:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,omega(size(x),size(x)),xox,xom,mom,arg,nu
   integer :: info,k
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x); call inverse_spd(sigma,omega,info)
   if(info/=0) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   xox=dot_product(x,matmul(omega,x)); xom=dot_product(x,matmul(omega,mu)); mom=dot_product(mu,matmul(omega,mu))
   xox=max(xox,1.0e-300_dp); arg=sqrt((2.0_dp+mom)*xox); nu=0.5_dp*real(2-k,dp)
   v=log(2.0_dp)+xom-0.5_dp*(real(k,dp)*log(2.0_dp*pi)+logdet_spd(sigma))
   v=v+0.25_dp*real(2-k,dp)*(log(xox)-log(2.0_dp+mom))+log(max(bessel_k_numeric(arg,nu),tiny(1.0_dp)))
   if(.not.lg) v=exp(v)
end function daml

subroutine raml(mu,sigma,x)
   real(dp), intent(in) :: mu(:),sigma(:,:)
   real(dp), intent(out) :: x(:)
   real(dp) :: z(size(mu)),e(size(mu))
   integer :: i,info
   call rand_mvn(0.0_dp*mu,sigma,z,info)
   do i=1,size(mu); e(i)=-log(rand_uniform()); end do
   x=mu*e+sqrt(e)*z
end subroutine raml

function dmvc(x,mu,s,log_density) result(v)
   real(dp), intent(in) :: x(:),mu(:),s(:,:)
   logical, intent(in), optional :: log_density
   real(dp) :: v
   if(present(log_density)) then; v=dmvt(x,mu,s,1.0_dp,log_density); else; v=dmvt(x,mu,s,1.0_dp); end if
end function dmvc

subroutine rmvc(mu,s,x)
   real(dp), intent(in) :: mu(:),s(:,:)
   real(dp), intent(out) :: x(:)
   real(dp) :: z(size(mu))
   integer :: info
   call rand_mvn(0.0_dp*mu,s,z,info); x=mu+z/abs(rand_normal())
end subroutine rmvc

function dmvl(x,mu,sigma,log_density) result(v)
   real(dp), intent(in) :: x(:),mu(:),sigma(:,:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,omega(size(x),size(x)),d(size(x)),z,arg,nu
   integer :: info,k
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x); call inverse_spd(sigma,omega,info)
   if(info/=0) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   d=x-mu; z=max(dot_product(d,matmul(omega,d)),1.0e-300_dp); arg=sqrt(2.0_dp*z); nu=1.0_dp-0.5_dp*real(k,dp)
   v=log(2.0_dp)-0.5_dp*real(k,dp)*log(2.0_dp*pi)-0.5_dp*logdet_spd(sigma)
   v=v+0.5_dp*nu*(log(z)-log(2.0_dp))+log(max(bessel_k_numeric(arg,nu),tiny(1.0_dp)))
   if(.not.lg) v=exp(v)
end function dmvl

subroutine rmvl(mu,sigma,x)
   real(dp), intent(in) :: mu(:),sigma(:,:)
   real(dp), intent(out) :: x(:)
   real(dp) :: z(size(mu)),e
   integer :: info
   e=-log(rand_uniform()); call rand_mvn(0.0_dp*mu,sigma,z,info); x=mu+sqrt(e)*z
end subroutine rmvl

function dmvpe(x,mu,sigma,kappa,log_density) result(v)
   real(dp), intent(in) :: x(:),mu(:),sigma(:,:),kappa
   logical, intent(in), optional :: log_density
   real(dp) :: v,omega(size(x),size(x)),d(size(x)),temp
   integer :: info,k
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x); call inverse_spd(sigma,omega,info)
   if(info/=0 .or. kappa<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   d=x-mu; temp=dot_product(d,matmul(omega,d))
   v=log(real(k,dp))+log_gamma(0.5_dp*real(k,dp))-0.5_dp*real(k,dp)*log(pi)
   v=v-0.5_dp*logdet_spd(sigma)-log_gamma(1.0_dp+real(k,dp)/(2.0_dp*kappa))
   v=v-(1.0_dp+real(k,dp)/(2.0_dp*kappa))*log(2.0_dp)-0.5_dp*kappa*temp
   if(.not.lg) v=exp(v)
end function dmvpe

subroutine rmvpe(mu,sigma,kappa,x)
   real(dp), intent(in) :: mu(:),sigma(:,:),kappa
   real(dp), intent(out) :: x(:)
   real(dp) :: z(size(mu)),u(size(mu)),normu,radius,l(size(mu),size(mu))
   integer :: i,info,k
   k=size(mu); do i=1,k; u(i)=rand_normal(); end do; normu=sqrt(dot_product(u,u)); u=u/max(normu,tiny(1.0_dp))
   radius=rand_gamma(real(k,dp)/(2.0_dp*kappa),2.0_dp)**(1.0_dp/(2.0_dp*kappa))
   call chol_lower(sigma,l,info); if(info/=0) then; x=mu; else; z=matmul(l,u); x=mu+radius*z; end if
end subroutine rmvpe

function dmvpolya(x,alpha,log_density) result(v)
   integer, intent(in) :: x(:),alpha(:)
   logical, intent(in), optional :: log_density
   real(dp) :: v
   integer :: i
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(size(x)/=size(alpha) .or. any(x<0) .or. any(alpha<1)) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log_gamma(real(sum(x)+1,dp))-sum([(log_gamma(real(x(i)+1,dp)),i=1,size(x))])
   v=v+log_gamma(real(sum(alpha),dp))-log_gamma(real(sum(x)+sum(alpha),dp))
   do i=1,size(x); v=v+log_gamma(real(x(i)+alpha(i),dp))-log_gamma(real(alpha(i),dp)); end do
   if(.not.lg) v=exp(v)
end function dmvpolya

function dmatrixgamma(x,alpha,beta,sigma,log_density) result(v)
   real(dp), intent(in) :: x(:,:),alpha,beta,sigma(:,:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,omega(size(x,1),size(x,2))
   integer :: info,k
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x,1); call inverse_spd(sigma,omega,info)
   if(info/=0 .or. alpha<=2.0_dp .or. beta<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=alpha*logdet_spd(omega)+(alpha-0.5_dp*real(k+1,dp))*logdet_spd(x)
   v=v-real(k,dp)*alpha*log(beta)-log_multivariate_gamma(alpha,k)-trace_matrix(matmul(omega,x))/beta
   if(.not.lg) v=exp(v)
end function dmatrixgamma

function dinvmatrixgamma(x,alpha,beta,psi,log_density) result(v)
   real(dp), intent(in) :: x(:,:),alpha,beta,psi(:,:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,xinv(size(x,1),size(x,2))
   integer :: info,k
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x,1); call inverse_spd(x,xinv,info)
   if(info/=0 .or. alpha<=2.0_dp .or. beta<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=alpha*logdet_spd(psi)-real(k,dp)*alpha*log(beta)-log_multivariate_gamma(alpha,k)
   v=v+(-alpha-0.5_dp*real(k+1,dp))*logdet_spd(x)-trace_matrix(matmul(psi,xinv))/beta
   if(.not.lg) v=exp(v)
end function dinvmatrixgamma

subroutine rmatrixgamma(alpha,beta,sigma,x)
   real(dp), intent(in) :: alpha,beta,sigma(:,:)
   real(dp), intent(out) :: x(:,:)
   call rwishart_local(2.0_dp*alpha,0.5_dp*beta*sigma,x)
end subroutine rmatrixgamma

subroutine rinvmatrixgamma(alpha,beta,psi,x)
   real(dp), intent(in) :: alpha,beta,psi(:,:)
   real(dp), intent(out) :: x(:,:)
   real(dp) :: scale(size(psi,1),size(psi,2)),scale_inv(size(psi,1),size(psi,2))
   real(dp) :: w(size(psi,1),size(psi,2))
   integer :: info
   scale=2.0_dp*psi/beta; call inverse_spd(scale,scale_inv,info)
   call rwishart_local(2.0_dp*alpha,scale_inv,w); call inverse_spd(w,x,info)
end subroutine rinvmatrixgamma

function dinvbeta(x,a,b,log_density) result(v)
   real(dp), intent(in) :: x,a,b
   logical, intent(in), optional :: log_density
   real(dp) :: v
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<=0.0_dp .or. a<=0.0_dp .or. b<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log_gamma(a+b)-log_gamma(a)-log_gamma(b)+(a-1.0_dp)*log(x)-(a+b)*log(1.0_dp+x)
   if(.not.lg) v=exp(v)
end function dinvbeta

function rinvbeta(a,b) result(x)
   use ld_random, only: rand_beta
   real(dp), intent(in) :: a,b
   real(dp) :: x,u
   u=rand_beta(a,b); x=u/max(1.0_dp-u,tiny(1.0_dp))
end function rinvbeta

function dhorseshoe(x,lambda,tau,log_density) result(v)
   real(dp), intent(in) :: x,lambda,tau
   logical, intent(in), optional :: log_density
   real(dp) :: v
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(lambda<=0.0_dp .or. tau<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=normal_logpdf(x,0.0_dp,lambda*tau); if(.not.lg) v=exp(v)
end function dhorseshoe

function rhorseshoe(lambda,tau) result(x)
   real(dp), intent(in) :: lambda,tau
   real(dp) :: x
   x=lambda*tau*rand_normal()
end function rhorseshoe

function dhuangwand(x,nu,alatent,ascale,log_density) result(v)
   real(dp), intent(in) :: x(:,:),nu,alatent(:),ascale(:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,s(size(x,1),size(x,2))
   integer :: i,k
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x,1)
   if(size(alatent)/=k .or. size(ascale)/=k .or. any(alatent<=0.0_dp) .or. any(ascale<=0.0_dp)) then
      v=merge(-huge(1.0_dp),0.0_dp,lg); return
   end if
   v=0.0_dp; s=0.0_dp
   do i=1,k
      v=v+dinvgamma(alatent(i),0.5_dp,1.0_dp/(ascale(i)*ascale(i)),.true.)
      s(i,i)=2.0_dp*nu/alatent(i)
   end do
   v=v+dinvwishart(x,nu+real(k-1,dp),s,.true.); if(.not.lg) v=exp(v)
end function dhuangwand

subroutine rhuangwand(nu,alatent,ascale,x)
   real(dp), intent(in) :: nu,alatent(:),ascale(:)
   real(dp), intent(out) :: x(:,:)
   real(dp) :: aa(size(alatent)),s(size(x,1),size(x,2)),sinv(size(x,1),size(x,2))
   real(dp) :: w(size(x,1),size(x,2))
   integer :: i,info,k
   k=size(alatent); s=0.0_dp
   do i=1,k
      aa(i)=rinvgamma(0.5_dp,1.0_dp/(ascale(i)*ascale(i)))
      s(i,i)=2.0_dp*nu/aa(i)
   end do
   call inverse_spd(s,sinv,info); call rwishart_local(nu+real(k-1,dp),sinv,w); call inverse_spd(w,x,info)
end subroutine rhuangwand

function dlasso(x,sigma,tau,lambda,a,b,log_density) result(v)
   real(dp), intent(in) :: x(:),sigma,tau(:),lambda,a,b
   logical, intent(in), optional :: log_density
   real(dp) :: v,cov(size(x),size(x)),lam2
   integer :: i
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(sigma<=0.0_dp .or. lambda<=0.0_dp .or. a<=0.0_dp .or. b<=0.0_dp .or. any(tau<=0.0_dp)) then
      v=merge(-huge(1.0_dp),0.0_dp,lg); return
   end if
   cov=0.0_dp; do i=1,size(x); cov(i,i)=sigma*sigma*tau(i)*tau(i); end do
   lam2=lambda*lambda; v=dmvn(x,0.0_dp*x,cov,.true.)-2.0_dp*log(sigma)
   do i=1,size(tau); v=v+log(0.5_dp*lam2)-0.5_dp*lam2*tau(i)*tau(i); end do
   v=v+a*log(b)-log_gamma(a)+(a-1.0_dp)*log(lam2)-b*lam2
   if(.not.lg) v=exp(v)
end function dlasso

function dlaplace_mixture(x,p,location,scale,log_density) result(v)
   real(dp), intent(in) :: x,p(:),location(:),scale(:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,lw(size(p))
   integer :: i
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   do i=1,size(p); lw(i)=log(max(p(i),tiny(1.0_dp)))+dlaplace(x,location(i),scale(i),.true.); end do
   v=logsumexp_local(lw); if(.not.lg) v=exp(v)
end function dlaplace_mixture

function rlaplace_mixture(p,location,scale) result(x)
   real(dp), intent(in) :: p(:),location(:),scale(:)
   real(dp) :: x
   integer :: k
   k=rand_categorical(p); x=rlaplace(location(k),scale(k))
end function rlaplace_mixture

function dnormal_mixture(x,p,mu,sigma,log_density) result(v)
   real(dp), intent(in) :: x,p(:),mu(:),sigma(:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,lw(size(p))
   integer :: i
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   do i=1,size(p); lw(i)=log(max(p(i),tiny(1.0_dp)))+normal_logpdf(x,mu(i),sigma(i)); end do
   v=logsumexp_local(lw); if(.not.lg) v=exp(v)
end function dnormal_mixture

function rnormal_mixture(p,mu,sigma) result(x)
   real(dp), intent(in) :: p(:),mu(:),sigma(:)
   real(dp) :: x
   integer :: k
   k=rand_categorical(p); x=mu(k)+sigma(k)*rand_normal()
end function rnormal_mixture

function dnormlaplace(x,mu,sigma,alpha,beta,log_density) result(v)
   real(dp), intent(in) :: x,mu,sigma,alpha,beta
   logical, intent(in), optional :: log_density
   real(dp) :: v,z1,z2,t1,t2
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(sigma<=0.0_dp .or. alpha<=0.0_dp .or. beta<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   z1=(x-mu)/sigma-alpha*sigma; z2=-(x-mu)/sigma-beta*sigma
   t1=log(alpha*beta/(alpha+beta))+alpha*(mu-x)+0.5_dp*(alpha*sigma)**2+log(max(normal_cdf(z1,0.0_dp,1.0_dp),tiny(1.0_dp)))
   t2=log(alpha*beta/(alpha+beta))+beta*(x-mu)+0.5_dp*(beta*sigma)**2+log(max(normal_cdf(z2,0.0_dp,1.0_dp),tiny(1.0_dp)))
   v=logsumexp_local([t1,t2]); if(.not.lg) v=exp(v)
end function dnormlaplace

function rnormlaplace(mu,sigma,alpha,beta) result(x)
   real(dp), intent(in) :: mu,sigma,alpha,beta
   real(dp) :: x
   x=sigma*rand_normal()+rslaplace(mu,1.0_dp/beta,1.0_dp/alpha)
end function rnormlaplace

function dnorminvwishart(mu,mu0,lambda,sigma,s,nu,log_density) result(v)
   real(dp), intent(in) :: mu(:),mu0(:),lambda,sigma(:,:),s(:,:),nu
   logical, intent(in), optional :: log_density
   real(dp) :: v,cov(size(mu),size(mu))
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; cov=sigma/lambda
   v=dinvwishart(sigma,nu,s,.true.)+dmvn(mu,mu0,cov,.true.); if(.not.lg) v=exp(v)
end function dnorminvwishart

function dnormwishart(mu,mu0,lambda,omega,s,nu,log_density) result(v)
   real(dp), intent(in) :: mu(:),mu0(:),lambda,omega(:,:),s(:,:),nu
   logical, intent(in), optional :: log_density
   real(dp) :: v,cov(size(mu),size(mu))
   integer :: info
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; call inverse_spd(lambda*omega,cov,info)
   if(info/=0) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=dwishart(omega,nu,s,.true.)+dmvn(mu,mu0,cov,.true.); if(.not.lg) v=exp(v)
end function dnormwishart

function dyangberger(x,log_density) result(v)
   real(dp), intent(in) :: x(:,:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,eig(size(x,1)),prodgap
   integer :: i,j
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; call jacobi_eigenvalues(x,eig); call sort_vector(eig); prodgap=1.0_dp
   do i=1,size(eig)-1; do j=i+1,size(eig); prodgap=prodgap*max(abs(eig(j)-eig(i)),tiny(1.0_dp)); end do; end do
   v=-logdet_spd(x)*prodgap; if(.not.lg) v=exp(v)
end function dyangberger

function dhyperg(g,alpha,log_density) result(v)
   real(dp), intent(in) :: g,alpha
   logical, intent(in), optional :: log_density
   real(dp) :: v
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(g<=0.0_dp .or. alpha<=2.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log((alpha-2.0_dp)/2.0_dp)-0.5_dp*alpha*log(1.0_dp+g); if(.not.lg) v=exp(v)
end function dhyperg

function dzellner(beta,g,sigma,xmat,log_density) result(v)
   real(dp), intent(in) :: beta(:),g,sigma,xmat(:,:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,xtx(size(beta),size(beta)),inv(size(beta),size(beta)),cov(size(beta),size(beta))
   integer :: info
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; xtx=matmul(transpose(xmat),xmat); call inverse_spd(xtx,inv,info)
   if(info/=0 .or. g<=0.0_dp .or. sigma<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   cov=g*sigma*sigma*inv; v=dmvn(beta,0.0_dp*beta,cov,.true.); if(.not.lg) v=exp(v)
end function dzellner

function dmvnormal_precision(x,mu,omega,log_density) result(v)
   real(dp), intent(in) :: x(:),mu(:),omega(:,:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,d(size(x))
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; d=x-mu
   v=-0.5_dp*real(size(x),dp)*log(2.0_dp*pi)+0.5_dp*logdet_spd(omega)-0.5_dp*dot_product(d,matmul(omega,d))
   if(.not.lg) v=exp(v)
end function dmvnormal_precision

function dmvstudent_precision(x,mu,omega,nu,log_density) result(v)
   real(dp), intent(in) :: x(:),mu(:),omega(:,:),nu
   logical, intent(in), optional :: log_density
   real(dp) :: v,d(size(x)),z
   integer :: k
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x); d=x-mu; z=dot_product(d,matmul(omega,d))
   v=log_gamma(0.5_dp*(nu+real(k,dp)))-log_gamma(0.5_dp*nu)+0.5_dp*logdet_spd(omega)
   v=v-0.5_dp*real(k,dp)*log(nu*pi)-0.5_dp*(nu+real(k,dp))*log(1.0_dp+z/nu); if(.not.lg) v=exp(v)
end function dmvstudent_precision

function dmvpower_precision(x,mu,omega,kappa,log_density) result(v)
   real(dp), intent(in) :: x(:),mu(:),omega(:,:),kappa
   logical, intent(in), optional :: log_density
   real(dp) :: v,d(size(x)),z
   integer :: k
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x); d=x-mu; z=dot_product(d,matmul(omega,d))
   v=log(real(k,dp))+log_gamma(0.5_dp*real(k,dp))-0.5_dp*real(k,dp)*log(pi)+0.5_dp*logdet_spd(omega)
   v=v-log_gamma(1.0_dp+real(k,dp)/(2.0_dp*kappa))-(1.0_dp+real(k,dp)/(2.0_dp*kappa))*log(2.0_dp)-0.5_dp*kappa*z
   if(.not.lg) v=exp(v)
end function dmvpower_precision

function stick_density(theta,gamma0,log_density) result(v)
   real(dp), intent(in) :: theta(:),gamma0
   logical, intent(in), optional :: log_density
   real(dp) :: v
   integer :: i
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(gamma0<=0.0_dp .or. any(theta<=0.0_dp) .or. any(theta>=1.0_dp)) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=0.0_dp; do i=1,size(theta); v=v+log(gamma0)+(gamma0-1.0_dp)*log(1.0_dp-theta(i)); end do
   if(.not.lg) v=exp(v)
end function stick_density

subroutine rstick(m,gamma0,w)
   integer, intent(in) :: m
   real(dp), intent(in) :: gamma0
   real(dp), intent(out) :: w(m+1)
   real(dp) :: rem,beta
   integer :: i
   rem=1.0_dp
   do i=1,m+1; beta=1.0_dp-rand_uniform()**(1.0_dp/gamma0); w(i)=rem*beta; rem=rem*(1.0_dp-beta); end do
   w=w/sum(w)
end subroutine rstick

function dtrunc_generic(x,pdf,cdf,a,b,log_density) result(v)
   real(dp), intent(in) :: x,a,b
   procedure(scalar_pdf_iface) :: pdf
   procedure(scalar_cdf_iface) :: cdf
   logical, intent(in), optional :: log_density
   real(dp) :: v,z
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; z=cdf(b)-cdf(a)
   if(a>=b .or. x<a .or. x>b .or. z<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   if(lg) then; v=log(max(pdf(x),tiny(1.0_dp)))-log(z); else; v=pdf(x)/z; end if
end function dtrunc_generic

function ptrunc_generic(x,cdf,a,b) result(v)
   real(dp), intent(in) :: x,a,b
   procedure(scalar_cdf_iface) :: cdf
   real(dp) :: v,z,xx
   xx=min(max(x,a),b); z=cdf(b)-cdf(a); v=(cdf(xx)-cdf(a))/z
end function ptrunc_generic

function qtrunc_generic(p,cdf,quantile,a,b) result(v)
   real(dp), intent(in) :: p,a,b
   procedure(scalar_cdf_iface) :: cdf
   procedure(scalar_quantile_iface) :: quantile
   real(dp) :: v
   v=quantile(cdf(a)+p*(cdf(b)-cdf(a)))
end function qtrunc_generic

function rtrunc_generic(cdf,quantile,a,b) result(v)
   real(dp), intent(in) :: a,b
   procedure(scalar_cdf_iface) :: cdf
   procedure(scalar_quantile_iface) :: quantile
   real(dp) :: v
   v=qtrunc_generic(rand_uniform(),cdf,quantile,a,b)
end function rtrunc_generic

function bessel_k_numeric(x,nu) result(v)
   real(dp), intent(in) :: x,nu
   real(dp) :: v,h,t,s,term
   integer, parameter :: n=600
   integer :: i
   if(x<=0.0_dp) then; v=huge(1.0_dp); return; end if; h=14.0_dp/real(n,dp); s=0.0_dp
   do i=0,n
      t=real(i,dp)*h; term=exp(-x*cosh(t))*cosh(abs(nu)*t)
      if(i==0 .or. i==n) then; s=s+term; else if(mod(i,2)==0) then; s=s+2.0_dp*term; else; s=s+4.0_dp*term; end if
   end do
   v=h*s/3.0_dp
end function bessel_k_numeric

subroutine rwishart_local(nu,s,x)
   real(dp), intent(in) :: nu,s(:,:)
   real(dp), intent(out) :: x(:,:)
   integer :: p,i,j,info
   real(dp) :: a(size(s,1),size(s,2)),l(size(s,1),size(s,2)),b(size(s,1),size(s,2))
   p=size(s,1); a=0.0_dp; call chol_lower(s,l,info); if(info/=0) then; x=0.0_dp; return; end if
   do i=1,p
      a(i,i)=sqrt(rand_chisq(nu-real(i-1,dp))); do j=1,i-1; a(i,j)=rand_normal(); end do
   end do
   b=matmul(l,a); x=matmul(b,transpose(b))
end subroutine rwishart_local

subroutine jacobi_eigenvalues(a,eig)
   real(dp), intent(in) :: a(:,:)
   real(dp), intent(out) :: eig(:)
   real(dp) :: b(size(a,1),size(a,2)),app,aqq,apq,phi,c,s,bip,biq
   integer :: n,it,p,q,i
   n=size(a,1); b=0.5_dp*(a+transpose(a))
   do it=1,100*n*n
      call max_offdiag(b,p,q,apq); if(abs(apq)<1.0e-12_dp) exit; app=b(p,p); aqq=b(q,q)
      phi=0.5_dp*atan2(2.0_dp*apq,aqq-app); c=cos(phi); s=sin(phi)
      do i=1,n
         if(i/=p .and. i/=q) then
            bip=b(i,p); biq=b(i,q); b(i,p)=c*bip-s*biq; b(p,i)=b(i,p); b(i,q)=s*bip+c*biq; b(q,i)=b(i,q)
         end if
      end do
      b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq; b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq; b(p,q)=0.0_dp; b(q,p)=0.0_dp
   end do
   eig=[(b(i,i),i=1,n)]
end subroutine jacobi_eigenvalues

subroutine max_offdiag(a,p,q,val)
   real(dp), intent(in) :: a(:,:)
   integer, intent(out) :: p,q
   real(dp), intent(out) :: val
   integer :: i,j
   p=1; q=min(2,size(a,1)); val=0.0_dp
   do i=1,size(a,1)-1; do j=i+1,size(a,1); if(abs(a(i,j))>abs(val)) then; val=a(i,j); p=i; q=j; end if; end do; end do
end subroutine max_offdiag

subroutine sort_vector(x)
   real(dp), intent(inout) :: x(:)
   integer :: i,j
   real(dp) :: t
   do i=2,size(x)
      t=x(i); j=i-1
      do while(j>=1)
         if(x(j)<=t) exit
         x(j+1)=x(j); j=j-1
      end do
      x(j+1)=t
   end do
end subroutine sort_vector

pure function logsumexp_local(x) result(v)
   real(dp), intent(in) :: x(:)
   real(dp) :: v,m
   m=maxval(x); v=m+log(sum(exp(x-m)))
end function logsumexp_local

end module ld_distributions_catalog
