module ld_distributions
use ld_kinds, only: dp, pi, log2pi
use ld_linalg, only: inverse_spd, logdet_spd, chol_lower, trace_matrix
use ld_random, only: rand_uniform, rand_normal, rand_gamma, rand_chisq, rand_mvn, rand_mv_t, rand_dirichlet
implicit none
private
public :: normal_pdf, normal_logpdf, normal_cdf, normal_quantile
public :: dbern, dcat, dlaplace, plaplace, qlaplace, rlaplace
public :: dhalfnorm, phalfnorm, qhalfnorm, dhalfcauchy, phalfcauchy, qhalfcauchy
public :: dinvgamma, rinvgamma, dpareto, ppareto, qpareto, rpareto
public :: ddirichlet, dmvn, dmvt, dgpd, dgpois, dst, dstp
public :: dwishart, dinvwishart, dmatrixnorm, log_multivariate_gamma
contains
pure function normal_logpdf(x,mu,sigma) result(v)
   real(dp), intent(in) :: x,mu,sigma
   real(dp) :: v,z
   if(sigma<=0.0_dp) then; v=-huge(1.0_dp); return; end if
   z=(x-mu)/sigma; v=-0.5_dp*log2pi-log(sigma)-0.5_dp*z*z
end function normal_logpdf
pure function normal_pdf(x,mu,sigma) result(v)
   real(dp), intent(in) :: x,mu,sigma
   real(dp) :: v
   v=exp(normal_logpdf(x,mu,sigma))
end function normal_pdf
pure function normal_cdf(x,mu,sigma) result(v)
   real(dp), intent(in) :: x,mu,sigma
   real(dp) :: v
   if(sigma<=0.0_dp) then; v=merge(1.0_dp,0.0_dp,x>=mu); return; end if
   v=0.5_dp*erfc(-(x-mu)/(sigma*sqrt(2.0_dp)))
end function normal_cdf
pure function normal_quantile(p,mu,sigma) result(x)
   real(dp), intent(in) :: p,mu,sigma
   real(dp) :: x,q,r,z
   real(dp), parameter :: a1=-3.969683028665376e+01_dp,a2=2.209460984245205e+02_dp
   real(dp), parameter :: a3=-2.759285104469687e+02_dp,a4=1.383577518672690e+02_dp
   real(dp), parameter :: a5=-3.066479806614716e+01_dp,a6=2.506628277459239e+00_dp
   real(dp), parameter :: b1=-5.447609879822406e+01_dp,b2=1.615858368580409e+02_dp
   real(dp), parameter :: b3=-1.556989798598866e+02_dp,b4=6.680131188771972e+01_dp
   real(dp), parameter :: b5=-1.328068155288572e+01_dp,c1=-7.784894002430293e-03_dp
   real(dp), parameter :: c2=-3.223964580411365e-01_dp,c3=-2.400758277161838e+00_dp
   real(dp), parameter :: c4=-2.549732539343734e+00_dp,c5=4.374664141464968e+00_dp
   real(dp), parameter :: c6=2.938163982698783e+00_dp,d1=7.784695709041462e-03_dp
   real(dp), parameter :: d2=3.224671290700398e-01_dp,d3=2.445134137142996e+00_dp
   real(dp), parameter :: d4=3.754408661907416e+00_dp,pl=0.02425_dp,ph=1.0_dp-pl
   if(p<=0.0_dp) then; x=-huge(1.0_dp); return; end if
   if(p>=1.0_dp) then; x=huge(1.0_dp); return; end if
   if(p<pl) then
      q=sqrt(-2.0_dp*log(p)); z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
   else if(p<=ph) then
      q=p-0.5_dp; r=q*q; z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
   else
      q=sqrt(-2.0_dp*log(1.0_dp-p)); z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
   end if
   x=mu+sigma*z
end function normal_quantile

pure function dbern(x,prob,log_density) result(v)
   integer, intent(in) :: x
   real(dp), intent(in) :: prob
   logical, intent(in), optional :: log_density
   real(dp) :: v
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(prob<0.0_dp .or. prob>1.0_dp .or. (x/=0 .and. x/=1)) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   if(x==1) then; v=log(max(prob,tiny(1.0_dp))); else; v=log(max(1.0_dp-prob,tiny(1.0_dp))); end if
   if(.not.lg) v=exp(v)
end function dbern

pure function dcat(x,p,log_density) result(v)
   integer, intent(in) :: x
   real(dp), intent(in) :: p(:)
   logical, intent(in), optional :: log_density
   real(dp) :: v,s
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density; s=sum(p)
   if(x<1 .or. x>size(p) .or. s<=0.0_dp .or. p(x)<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log(p(x)/s); if(.not.lg) v=exp(v)
end function dcat

pure function dlaplace(x,location,scale,log_density) result(v)
   real(dp), intent(in) :: x,location,scale
   logical, intent(in), optional :: log_density
   real(dp) :: v
   logical :: lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(scale<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=-abs(x-location)/scale-log(2.0_dp*scale); if(.not.lg) v=exp(v)
end function dlaplace
pure function plaplace(q,location,scale) result(p)
   real(dp),intent(in)::q,location,scale; real(dp)::p,z
   z=(q-location)/scale; if(q<location) then; p=0.5_dp*exp(z); else; p=1.0_dp-0.5_dp*exp(-z); end if
end function plaplace
pure function qlaplace(p,location,scale) result(q)
   real(dp),intent(in)::p,location,scale; real(dp)::q
   if(p<=0.0_dp) then; q=-huge(1.0_dp); else if(p>=1.0_dp) then; q=huge(1.0_dp)
   else if(p<=0.5_dp) then; q=location+scale*log(2.0_dp*p); else; q=location-scale*log(2.0_dp*(1.0_dp-p)); end if
end function qlaplace
function rlaplace(location,scale) result(x)
   real(dp),intent(in)::location,scale; real(dp)::x
   x=qlaplace(rand_uniform(),location,scale)
end function rlaplace

pure function dhalfnorm(x,scale,log_density) result(v)
   real(dp),intent(in)::x,scale; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<0.0_dp .or. scale<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log(2.0_dp)+normal_logpdf(x,0.0_dp,scale); if(.not.lg) v=exp(v)
end function dhalfnorm
pure function phalfnorm(q,scale) result(p)
   real(dp),intent(in)::q,scale; real(dp)::p
   if(q<=0.0_dp) then; p=0.0_dp; else; p=2.0_dp*normal_cdf(q,0.0_dp,scale)-1.0_dp; end if
end function phalfnorm
pure function qhalfnorm(p,scale) result(q)
   real(dp),intent(in)::p,scale; real(dp)::q
   q=normal_quantile(0.5_dp*(p+1.0_dp),0.0_dp,scale)
end function qhalfnorm
pure function dhalfcauchy(x,scale,log_density) result(v)
   real(dp),intent(in)::x,scale; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<0.0_dp .or. scale<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log(2.0_dp/(pi*scale))-log(1.0_dp+(x/scale)**2); if(.not.lg) v=exp(v)
end function dhalfcauchy
pure function phalfcauchy(q,scale) result(p)
   real(dp),intent(in)::q,scale; real(dp)::p
   if(q<=0.0_dp) then; p=0.0_dp; else; p=2.0_dp*atan(q/scale)/pi; end if
end function phalfcauchy
pure function qhalfcauchy(p,scale) result(q)
   real(dp),intent(in)::p,scale; real(dp)::q
   q=scale*tan(0.5_dp*pi*p)
end function qhalfcauchy

pure function dinvgamma(x,shape,scale,log_density) result(v)
   real(dp),intent(in)::x,shape,scale; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<=0.0_dp .or. shape<=0.0_dp .or. scale<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=shape*log(scale)-log_gamma(shape)-(shape+1.0_dp)*log(x)-scale/x; if(.not.lg) v=exp(v)
end function dinvgamma
function rinvgamma(shape,scale) result(x)
   real(dp),intent(in)::shape,scale; real(dp)::x
   x=1.0_dp/rand_gamma(shape,1.0_dp/scale)
end function rinvgamma

pure function dpareto(x,alpha,log_density) result(v)
   real(dp),intent(in)::x,alpha; logical,intent(in),optional::log_density; real(dp)::v; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<1.0_dp .or. alpha<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log(alpha)-(alpha+1.0_dp)*log(x); if(.not.lg) v=exp(v)
end function dpareto
pure function ppareto(q,alpha) result(p)
   real(dp),intent(in)::q,alpha; real(dp)::p
   if(q<1.0_dp) then; p=0.0_dp; else; p=1.0_dp-q**(-alpha); end if
end function ppareto
pure function qpareto(p,alpha) result(q)
   real(dp),intent(in)::p,alpha; real(dp)::q
   q=(1.0_dp-p)**(-1.0_dp/alpha)
end function qpareto
function rpareto(alpha) result(x)
   real(dp),intent(in)::alpha; real(dp)::x
   x=rand_uniform()**(-1.0_dp/alpha)
end function rpareto

pure function ddirichlet(x,alpha,log_density) result(v)
   real(dp),intent(in)::x(:),alpha(:); logical,intent(in),optional::log_density; real(dp)::v; logical::lg; integer::i
   lg=.false.; if(present(log_density)) lg=log_density
   if(size(x)/=size(alpha) .or. any(x<=0.0_dp) .or. any(alpha<=0.0_dp) .or. abs(sum(x)-1.0_dp)>1e-8_dp) then
      v=merge(-huge(1.0_dp),0.0_dp,lg); return
   end if
   v=log_gamma(sum(alpha)); do i=1,size(x); v=v-log_gamma(alpha(i))+(alpha(i)-1.0_dp)*log(x(i)); end do
   if(.not.lg) v=exp(v)
end function ddirichlet

function dmvn(x,mu,sigma,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),sigma(:,:); logical,intent(in),optional::log_density; real(dp)::v
   real(dp)::omega(size(x),size(x)),d(size(x)); integer::info; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density; call inverse_spd(sigma,omega,info)
   if(info/=0) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   d=x-mu; v=-0.5_dp*(real(size(x),dp)*log2pi+logdet_spd(sigma)+dot_product(d,matmul(omega,d)))
   if(.not.lg) v=exp(v)
end function dmvn

function dmvt(x,mu,s,df,log_density) result(v)
   real(dp),intent(in)::x(:),mu(:),s(:,:),df; logical,intent(in),optional::log_density; real(dp)::v
   real(dp)::omega(size(x),size(x)),d(size(x)),z; integer::info,k; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density; k=size(x); call inverse_spd(s,omega,info)
   if(info/=0 .or. df<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   d=x-mu; z=dot_product(d,matmul(omega,d))
   v=log_gamma(0.5_dp*(df+k))-log_gamma(0.5_dp*df)-0.5_dp*real(k,dp)*log(df*pi)&
      -0.5_dp*logdet_spd(s)-0.5_dp*(df+k)*log(1.0_dp+z/df)
   if(.not.lg) v=exp(v)
end function dmvt

pure function dgpd(x,mu,sigma,xi,log_density) result(v)
   real(dp),intent(in)::x,mu,sigma,xi; logical,intent(in),optional::log_density; real(dp)::v,z; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(sigma<=0.0_dp .or. x<mu) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   if(abs(xi)<1e-12_dp) then
      v=-log(sigma)-(x-mu)/sigma
   else
      z=1.0_dp+xi*(x-mu)/sigma
      if(z<=0.0_dp) then
         v=merge(-huge(1.0_dp),0.0_dp,lg)
         return
      end if
      v=-log(sigma)-(1.0_dp/xi+1.0_dp)*log(z)
   end if
   if(.not.lg) v=exp(v)
end function dgpd

pure function dgpois(x,lambda,omega,log_density) result(v)
   integer,intent(in)::x; real(dp),intent(in)::lambda,omega; logical,intent(in),optional::log_density; real(dp)::v,t; logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(x<0 .or. lambda<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   if(x==0) then
      v=-lambda
   else
      t=lambda+omega*real(x,dp)
      if(t<=0.0_dp) then
         v=merge(-huge(1.0_dp),0.0_dp,lg)
         return
      end if
      v=log(lambda)+real(x-1,dp)*log(t)-t-log_gamma(real(x+1,dp))
   end if
   if(.not.lg) v=exp(v)
end function dgpois

pure function dst(x,mu,sigma,nu,log_density) result(v)
   real(dp),intent(in)::x,mu,sigma,nu
   logical,intent(in),optional::log_density
   real(dp)::v
   logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(sigma<=0.0_dp .or. nu<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)-0.5_dp*log(pi*nu)-log(sigma)&
      -0.5_dp*(nu+1.0_dp)*log(1.0_dp+((x-mu)/sigma)**2/nu)
   if(.not.lg) v=exp(v)
end function dst
pure function dstp(x,mu,tau,nu,log_density) result(v)
   real(dp),intent(in)::x,mu,tau,nu
   logical,intent(in),optional::log_density
   real(dp)::v
   logical::lg
   lg=.false.; if(present(log_density)) lg=log_density
   if(tau<=0.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=dst(x,mu,1.0_dp/sqrt(tau),nu,.true.); if(.not.lg) v=exp(v)
end function dstp

pure function log_multivariate_gamma(a,p) result(v)
   real(dp),intent(in)::a
   integer,intent(in)::p
   real(dp)::v
   integer::j
   v=0.25_dp*real(p*(p-1),dp)*log(pi); do j=1,p; v=v+log_gamma(a+0.5_dp*real(1-j,dp)); end do
end function log_multivariate_gamma

function dwishart(omega,nu,s,log_density) result(v)
   real(dp),intent(in)::omega(:,:),nu,s(:,:)
   logical,intent(in),optional::log_density
   real(dp)::v
   real(dp)::sinv(size(s,1),size(s,2))
   integer::info,p
   logical::lg
   lg=.false.; if(present(log_density)) lg=log_density; p=size(s,1); call inverse_spd(s,sinv,info)
   if(info/=0 .or. nu<=real(p-1,dp)) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=0.5_dp*(nu-real(p+1,dp))*logdet_spd(omega)-0.5_dp*sum(sinv*transpose(omega))&
      -0.5_dp*nu*real(p,dp)*log(2.0_dp)-0.5_dp*nu*logdet_spd(s)-log_multivariate_gamma(0.5_dp*nu,p)
   if(.not.lg) v=exp(v)
end function dwishart

function dinvwishart(sigma,nu,s,log_density) result(v)
   real(dp),intent(in)::sigma(:,:),nu,s(:,:)
   logical,intent(in),optional::log_density
   real(dp)::v
   real(dp)::sig_inv(size(s,1),size(s,2))
   integer::info,p
   logical::lg
   lg=.false.; if(present(log_density)) lg=log_density; p=size(s,1); call inverse_spd(sigma,sig_inv,info)
   if(info/=0 .or. nu<=real(p-1,dp)) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   v=0.5_dp*nu*logdet_spd(s)-0.5_dp*nu*real(p,dp)*log(2.0_dp)-log_multivariate_gamma(0.5_dp*nu,p)&
      -0.5_dp*(nu+real(p+1,dp))*logdet_spd(sigma)-0.5_dp*sum(s*transpose(sig_inv))
   if(.not.lg) v=exp(v)
end function dinvwishart

function dmatrixnorm(x,m,u,vmat,log_density) result(v)
   real(dp),intent(in)::x(:,:),m(:,:),u(:,:),vmat(:,:)
   logical,intent(in),optional::log_density
   real(dp)::v
   real(dp)::ui(size(u,1),size(u,2)),vi(size(vmat,1),size(vmat,2))
   real(dp)::d(size(x,1),size(x,2))
   integer::info,n,p
   logical::lg
   lg=.false.; if(present(log_density)) lg=log_density; n=size(x,1); p=size(x,2); call inverse_spd(u,ui,info)
   if(info/=0) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if; call inverse_spd(vmat,vi,info)
   if(info/=0) then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
   d=x-m; v=-0.5_dp*real(n*p,dp)*log2pi-0.5_dp*real(p,dp)*logdet_spd(u)-0.5_dp*real(n,dp)*logdet_spd(vmat)&
      -0.5_dp*trace_matrix(matmul(ui,matmul(d,matmul(vi,transpose(d)))))
   if(.not.lg) v=exp(v)
end function dmatrixnorm
end module ld_distributions
