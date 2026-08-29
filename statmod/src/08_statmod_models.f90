module statmod_models
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf, ieee_quiet_nan
use r_compat, only: dp, r_digamma, r_trigamma
use statmod_linalg, only: least_squares, weighted_least_squares, solve_spd, symmetric_inverse, &
   column_space_basis, weighted_hat_basis, svd_full_u, logdet_xtwx, variance_vec
use statmod_special, only: canonic_digamma
implicit none
private
public :: glm_fit_result_t, mixed_fit_result_t, reml_fit_result_t
public :: glmgam_fit, glmnb_fit, fit_nbp, mixed_model2_fit, remlscore, remlscoregamma

type :: glm_fit_result_t
   real(dp),allocatable::coefficients(:),fitted(:)
   real(dp)::deviance=0
   integer::iter=0
   character(len=16)::convergence='converged'
end type

type :: mixed_fit_result_t
   real(dp)::varcomp(2)=0,se_varcomp(2)=0
   real(dp),allocatable::coefficients(:),se_coefficients(:)
   integer::info=0
end type

type :: reml_fit_result_t
   real(dp),allocatable::beta(:),se_beta(:),gamma(:),se_gamma(:),mu(:),phi(:),h(:),cov_beta(:,:),cov_gamma(:,:)
   real(dp)::deviance=0
   integer::iter=0,info=0
end type
contains

pure function gamma_deviance(y,mu) result(dev)
real(dp),intent(in)::y(:),mu(:)
real(dp)::dev
integer::i
dev=0
if(any(mu<0))then
dev=ieee_value(dev,ieee_positive_inf)
return
end if
do i=1,size(y)
   if(y(i)<1e-15_dp.and.mu(i)<1e-15_dp)cycle
   if(y(i)<=0.or.mu(i)<=0)then
   dev=ieee_value(dev,ieee_positive_inf)
   return
   end if
   dev=dev+2*((y(i)-mu(i))/mu(i)-log(y(i)/mu(i)))
end do
end function

function glmgam_fit(x,y,coef_start,tol,maxit) result(out)
real(dp),intent(in)::x(:,:),y(:)
real(dp),intent(in),optional::coef_start(:),tol
integer,intent(in),optional::maxit
type(glm_fit_result_t)::out
real(dp),allocatable::beta(:),mu(:),fit(:),res(:),w(:),xvx(:,:),dl(:),db(:),a(:,:),b0(:),ainv(:,:)
real(dp)::dev,devold,lambda,maxinfo,tolerance,maxy
integer::n,p,rank,info,iter,lev,mx,i,j
n=size(y)
p=size(x,2)
tolerance=1e-6_dp
if(present(tol))tolerance=tol
mx=50
if(present(maxit))mx=maxit
if(n==0)then
allocate(out%coefficients(0),out%fitted(0))
return
end if
maxy=maxval(y)
if(maxy==0)then
   allocate(out%coefficients(p),out%fitted(n))
   out%coefficients=0
   out%fitted=0
   out%deviance=ieee_value(out%deviance,ieee_quiet_nan)
   return
end if
if(present(coef_start))then
beta=coef_start
mu=matmul(x,beta)
else
   call least_squares(x,y,beta,fit,res,rank,info)
   mu=fit
   if(any(mu<0))then
      w=1/max(y,maxy*1e-3_dp)**2
      call weighted_least_squares(x,y,w,beta,fit,res,rank,info)
      mu=fit
   end if
   if(any(mu<0))then
      b0=[sum(y)/real(n,dp)]
      call least_squares(x,spread(sum(y)/real(n,dp),1,n),beta,fit,res,rank,info)
      mu=fit
   end if
   if(any(mu<0))then
      ! Final upstream fallback: use the first column whose entries all have
      ! the same sign, otherwise return an infinite-deviance failure fit.
      block
         integer :: jj
         real(dp) :: den
         logical :: found
         found=.false.
         do jj=1,p
            if(all(x(:,jj)>0.0_dp).or.all(x(:,jj)<0.0_dp))then
               beta=0.0_dp
               w=1.0_dp/max(y,maxy*1.0e-3_dp)**2
               den=sum(w*x(:,jj)*x(:,jj))
               if(den>0.0_dp)beta(jj)=sum(w*x(:,jj)*y)/den
               mu=x(:,jj)*beta(jj)
               found=.true.
               exit
            end if
         end do
         if(.not.found)then
            allocate(out%coefficients(p),out%fitted(n))
            out%coefficients=0.0_dp
            out%fitted=0.0_dp
            out%deviance=ieee_value(out%deviance,ieee_positive_inf)
            return
         end if
      end block
   end if
end if
dev=gamma_deviance(y,mu)
allocate(xvx(p,p),dl(p),db(p),a(p,p))
do iter=1,mx
   w=1/max(mu*mu,maxval(mu*mu)/1e3_dp)
   xvx=0
   do i=1,n
   do j=1,p
   xvx(j,:)=xvx(j,:)+w(i)*x(i,j)*x(i,:)
   end do
   end do
   maxinfo=maxval([(xvx(i,i),i=1,p)])
   if(iter==1)lambda=abs(sum([(xvx(i,i),i=1,p)])/p)/p
   dl=matmul(transpose(x),(y-mu)*w)
   devold=dev
   do lev=1,100
      a=xvx
      do i=1,p
      a(i,i)=a(i,i)+lambda
      end do
      call solve_spd(a,dl,db,info)
      if(info/=0)then
      lambda=lambda*2
      cycle
      end if
      fit=beta+db
      res=matmul(x,fit)
      dev=gamma_deviance(y,res)
      if(dev<=devold.or.dev/maxval(max(res,1e-300_dp))<1e-15_dp)then
      beta=fit
      mu=res
      exit
      end if
      lambda=lambda*2
      if(lambda/max(maxinfo,tiny(1.0_dp))>1e15_dp)exit
   end do
   if(lambda/max(maxinfo,tiny(1.0_dp))>1e15_dp)then
   out%convergence='lambdabig'
   exit
   end if
   if(lev==1)lambda=lambda/10
   if(dot_product(dl,db)<tolerance.or.dev/maxval(max(mu,1e-300_dp))<1e-15_dp)exit
end do
out%coefficients=beta
out%fitted=mu
out%deviance=dev
out%iter=iter
if(iter>mx)out%convergence='maxit'
end function

pure function nb_total_deviance(y,mu,phi,w) result(dev)
real(dp),intent(in)::y(:),mu(:),phi(:),w(:)
real(dp)::dev,yy,mm,ph,b,b2,ud,alpha
integer::i
dev=0
if(any(.not.ieee_is_finite(mu)))then
dev=huge(1.0_dp)
return
end if
do i=1,size(y)
 yy=y(i)+1e-8_dp
 mm=mu(i)+1e-8_dp
 ph=phi(i)
 if(ph<1e-4_dp)then
 b=yy-mm
 b2=0.5_dp*b*b*ph*(1+ph*(2.0_dp/3*b-yy))
 ud=2*(yy*log(yy/mm)-b-b2)
 else if(ph*mm>1e6_dp)then
 alpha=mm/(1+ph*mm)
 ud=2*((yy-mm)/mm-log(yy/mm))*alpha
 else
 ud=2*(yy*log(yy/mm)-(yy+1/ph)*log((1+yy*ph)/(1+mm*ph)))
 end if
 dev=dev+w(i)*ud
end do
end function

function glmnb_fit(x,y,dispersion,weights,offset,coef_start,tol,maxit) result(out)
real(dp),intent(in)::x(:,:),y(:),dispersion(:)
real(dp),intent(in),optional::weights(:),offset(:),coef_start(:),tol
integer,intent(in),optional::maxit
type(glm_fit_result_t)::out
real(dp),allocatable::phi(:),wt(:),off(:),beta(:),mu(:),fit(:),res(:),iv(:),xvx(:,:),dl(:),db(:),a(:,:),rate(:),nw(:)
real(dp)::dev,devnew,lambda,maxinfo,ceiling,tolerance,ymax,bmean
integer::n,p,i,j,rank,info,iter,lev,mx
n=size(y)
p=size(x,2)
ymax=maxval(y)
tolerance=1e-6_dp
if(present(tol))tolerance=tol
mx=50
if(present(maxit))mx=maxit
allocate(phi(n),wt(n),off(n))
do i=1,n
phi(i)=dispersion(1+mod(i-1,size(dispersion)))
end do
wt=1
if(present(weights))wt=weights
off=0
if(present(offset))off=offset
if(ymax==0.and.all(abs(off)<1e-14_dp))then
   ! Upstream all-zero special case: if X spans the intercept, drive its
   ! linear predictor to a very large negative value and return mu=0.
   call least_squares(x,spread(1.0_dp,1,n),beta,fit,res,rank,info)
   allocate(out%coefficients(p),out%fitted(n))
   if(info==0.and.maxval(abs(res))<1.0_dp)then
      out%coefficients=-1.0e10_dp*beta
   else
      out%coefficients=0.0_dp
   end if
   out%fitted=0
   out%deviance=0
   out%iter=0
   return
end if
if(present(coef_start))then
beta=coef_start
mu=exp(matmul(x,beta)+off)
else
   rate=y/exp(off)
   nw=wt*exp(off)/(1+phi*exp(off))
   bmean=log(max(sum(nw*rate)/sum(nw),tiny(1.0_dp)))
   call least_squares(x,spread(bmean,1,n),beta,fit,res,rank,info)
   mu=exp(matmul(x,beta)+off)
end if
dev=nb_total_deviance(y,mu,phi,wt)
allocate(xvx(p,p),dl(p),db(p),a(p,p))
do iter=1,mx
   iv=1+phi*mu
   xvx=0
   do i=1,n
   do j=1,p
   xvx(j,:)=xvx(j,:)+wt(i)*mu(i)/iv(i)*x(i,j)*x(i,:)
   end do
   end do
   maxinfo=maxval([(xvx(i,i),i=1,p)])
   if(iter==1)then
   lambda=max(maxinfo*1e-6_dp,1e-13_dp)
   ceiling=maxinfo*1e13_dp
   end if
   dl=matmul(transpose(x),wt*(y-mu)/iv)
   do lev=1,100
      a=xvx
      do i=1,p
      a(i,i)=a(i,i)+lambda
      end do
      call solve_spd(a,dl,db,info)
      if(info/=0)then
      lambda=lambda*10
      cycle
      end if
      fit=beta+db
      res=exp(matmul(x,fit)+off)
      devnew=nb_total_deviance(y,res,phi,wt)
      if(devnew<=dev.or.devnew/max(ymax,tiny(1.0_dp))<1e-13_dp)then
      beta=fit
      mu=res
      dev=devnew
      exit
      end if
      lambda=lambda*2
      if(lambda>ceiling)exit
   end do
   if(lambda>ceiling)then
   out%convergence='lambdabig'
   exit
   end if
   if(dot_product(dl,db)<tolerance.or.dev/max(ymax,tiny(1.0_dp))<1e-12_dp)exit
   if(lev==1)lambda=lambda/10
end do
out%coefficients=beta
out%fitted=mu
out%deviance=dev
out%iter=iter
if(iter>=mx)out%convergence='maxit'
end function

subroutine fit_nbp(y,group,lib_size,coefficients,fitted,dispersion,tol,maxit)
real(dp),intent(in)::y(:,:),lib_size(:)
integer,intent(in)::group(:)
real(dp),allocatable,intent(out)::coefficients(:,:),fitted(:,:)
real(dp),intent(out)::dispersion
real(dp),intent(in),optional::tol
integer,intent(in),optional::maxit
integer,allocatable::levels(:),idx(:)
real(dp),allocatable::mu(:,:),off(:,:),w(:,:),z(:,:),eta(:,:),e2(:,:),v(:,:),dv(:,:)
real(dp)::phi,x2,dx2,step,conv,tolerance,low
integer::ng,nlib,ngene,g,i,j,iter,inner,mx
call unique_int_local(group,levels)
ng=size(levels)
nlib=size(y,2)
ngene=size(y,1)
tolerance=1e-5_dp
if(present(tol))tolerance=tol
mx=40
if(present(maxit))mx=maxit
allocate(coefficients(ngene,ng),fitted(ngene,nlib),off(ngene,nlib),mu(ngene,nlib),w(ngene,nlib),z(ngene,nlib),eta(ngene,nlib))
do j=1,nlib
off(:,j)=log(lib_size(j))
end do
mu=max(y,0.5_dp)
phi=0
w=mu
z=w*(log(mu)-off)
coefficients=0
eta=off
do g=1,ng
   idx=pack([(j,j=1,nlib)],group==levels(g))
   do i=1,ngene
   coefficients(i,g)=sum(z(i,idx))/sum(w(i,idx))
   eta(i,idx)=eta(i,idx)+coefficients(i,g)
   end do
end do
mu=exp(eta)
do iter=1,mx
   e2=(y-mu)**2
   dv=mu*mu
   do inner=1,10
      v=mu*(1+phi*mu)
      x2=sum(e2/v)/real(nlib-ng,dp)-ngene
      if(x2>=0)then
      low=phi
      exit
      end if
      if(phi==0)exit
      if(inner>4)then
      phi=0.9_dp*phi
      else
      phi=(low+phi)/2
      end if
   end do
   if(x2<0)exit
   dx2=sum(e2/(v*v)*dv)/real(nlib-ng,dp)
   step=x2/max(dx2,1e-6_dp)
   phi=phi+step
   conv=step/(phi+1)
   if(conv<tolerance)exit
   w=mu/(1+phi*mu)
   z=(y-mu)/v*mu
   eta=off
   do g=1,ng
      idx=pack([(j,j=1,nlib)],group==levels(g))
      do i=1,ngene
      coefficients(i,g)=coefficients(i,g)+sum(z(i,idx))/sum(w(i,idx))
      eta(i,idx)=eta(i,idx)+coefficients(i,g)
      end do
   end do
   mu=exp(eta)
end do
fitted=mu
dispersion=phi
end subroutine

function mixed_model2_fit(y,x,z,w,only_varcomp,tol,maxit) result(out)
real(dp),intent(in)::y(:),x(:,:),z(:,:)
real(dp),intent(in),optional::w(:),tol
logical,intent(in),optional::only_varcomp
integer,intent(in),optional::maxit
type(mixed_fit_result_t)::out
real(dp),allocatable::ys(:),xs(:,:),zs(:,:),u(:,:),s(:),qres(:,:),qtz(:,:),qty(:),us(:,:),sv(:),uqy(:),d(:),dx(:,:),dy(:)
real(dp),allocatable::beta(:),fit(:),res(:),vc(:),dfit(:),dres(:),dinv(:,:),vv(:),uz(:,:),xt(:,:),ainv(:,:)
real(dp)::tolerance
integer::n,p,nz,rank,info,mq,i,j,mx
logical::ov
n=size(y)
p=size(x,2)
nz=size(z,2)
ov=.false.
if(present(only_varcomp))ov=only_varcomp
mx=50
if(present(maxit))mx=maxit
tolerance=1.0e-6_dp
if(present(tol))tolerance=tol
ys=y
xs=x
zs=z
if(present(w))then
   do i=1,n
      ys(i)=sqrt(w(i))*ys(i)
      xs(i,:)=sqrt(w(i))*xs(i,:)
      zs(i,:)=sqrt(w(i))*zs(i,:)
   end do
end if
call svd_full_u(xs,u,s,rank,info)
if(info/=0)then
out%info=info
return
end if
mq=n-rank
if(mq<=0)then
out%varcomp=[0.0_dp,0.0_dp]
out%info=1
return
end if
qres=u(:,rank+1:n)
qtz=matmul(transpose(qres),zs)
qty=matmul(transpose(qres),ys)
call svd_full_u(qtz,us,sv,rank,info)
uqy=matmul(transpose(us),qty)
allocate(d(mq))
d=0
if(size(sv)>0)d(1:size(sv))=sv*sv
allocate(dx(mq,2))
dx(:,1)=1
dx(:,2)=d
dy=uqy*uqy
call least_squares(dx,dy,vc,fit,dres,rank,info)
if(mq>2.and.count(abs(d)>1e-15_dp)>1.and.variance_vec(d)>1e-15_dp)then
   if(any(fit<0))vc=[sum(dy)/mq,0.0_dp]
   block
      type(glm_fit_result_t)::gf
      gf=glmgam_fit(dx,dy,coef_start=vc,tol=tolerance,maxit=mx)
      vc=gf%coefficients
      fit=gf%fitted
   end block
end if
out%varcomp=vc
if(ov)return
call info_inverse(dx,1/max(fit*fit,tiny(1.0_dp)),ainv,info)
if(info==0)out%se_varcomp=sqrt(2*[(ainv(i,i),i=1,2)])
! V = residual I + block Z Z'. Use eigen/SVD of Z through full U.
call svd_full_u(zs,u,sv,rank,info)
allocate(vv(n))
vv=vc(1)
if(size(sv)>0)vv(1:size(sv))=vv(1:size(sv))+vc(2)*sv*sv
xt=matmul(transpose(u),xs)
qty=matmul(transpose(u),ys)
call weighted_least_squares(xt,qty,1/max(vv,tiny(1.0_dp)),beta,fit,res,rank,info)
out%coefficients=beta
call info_inverse(xt,1/max(vv,tiny(1.0_dp)),ainv,info)
if(info==0)out%se_coefficients=sqrt([(ainv(i,i),i=1,size(beta))])
end function

subroutine info_inverse(x,w,ainv,info)
real(dp),intent(in)::x(:,:),w(:)
real(dp),allocatable,intent(out)::ainv(:,:)
integer,intent(out)::info
real(dp),allocatable::a(:,:)
integer::i,j,p
p=size(x,2)
allocate(a(p,p))
a=0
do i=1,size(x,1)
do j=1,p
a(j,:)=a(j,:)+w(i)*x(i,j)*x(i,:)
end do
end do
call symmetric_inverse(a,ainv,info)
end subroutine

function remlscore(y,x,z,tol,maxit) result(out)
real(dp),intent(in)::y(:),x(:,:),z(:,:)
real(dp),intent(in),optional::tol
integer,intent(in),optional::maxit
type(reml_fit_result_t)::out
real(dp),allocatable::beta(:),fit(:),res(:),q(:,:),h(:),d(:),wd(:),zd(:),gam(:),gfit(:),gres(:),phi(:),wm(:)
real(dp),allocatable::q2(:,:),q2z(:,:),zvz(:,:),dl(:),dgam(:),a(:,:),covg(:,:),covb(:,:)
real(dp)::dev,devold,lambda,maxinfo,tolerance,ldet
integer::n,p,nq,rank,info,iter,lev,mx,i,j,k,j0
n=size(y)
p=size(x,2)
nq=size(z,2)
tolerance=1e-5_dp
if(present(tol))tolerance=tol
mx=40
if(present(maxit))mx=maxit
call least_squares(x,y,beta,fit,res,rank,info)
if(rank<p)then
out%info=1
return
end if
call column_space_basis(x,q,rank,info)
h=sum(q*q,dim=2)
d=res*res
wd=1-h
zd=log(max(d/max(wd,tiny(1.0_dp)),tiny(1.0_dp)))+1.27_dp
call weighted_least_squares(z,zd,max(wd,0.0_dp),gam,gfit,gres,rank,info)
phi=exp(matmul(z,gam))
wm=1/phi
call weighted_least_squares(x,y,wm,beta,fit,res,rank,info)
d=res*res
ldet=logdet_xtwx(x,wm,info)
dev=sum(d/phi)+sum(log(phi))+n*log(2*acos(-1.0_dp))+ldet
allocate(q2(n,p*(p+1)/2),zvz(nq,nq),dl(nq),dgam(nq),a(nq,nq))
do iter=1,mx
   call weighted_hat_basis(x,wm,q,rank,info)
   h=sum(q*q,dim=2)
   j0=0
   do k=0,p-1
      do j=1,p-k
      j0=j0+1
      q2(:,j0)=q(:,j)*q(:,j+k)
      if(k>0)q2(:,j0)=sqrt(2.0_dp)*q2(:,j0)
      end do
   end do
   q2z=matmul(transpose(q2),z)
   zvz=0
   do i=1,n
   do j=1,nq
   zvz(j,:)=zvz(j,:)+(1-2*h(i))*z(i,j)*z(i,:)
   end do
   end do
   zvz=(zvz+matmul(transpose(q2z),q2z))/2
   maxinfo=maxval([(zvz(i,i),i=1,nq)])
   if(iter==1)lambda=abs(sum([(zvz(i,i),i=1,nq)])/nq)/nq
   zd=(d-(1-h)*phi)/phi
   dl=matmul(transpose(z),zd)/2
   devold=dev
   do lev=1,100
      a=zvz
      do i=1,nq
      a(i,i)=a(i,i)+lambda
      end do
      call solve_spd(a,dl,dgam,info)
      if(info/=0)then
      lambda=lambda*2
      cycle
      end if
      gfit=gam+dgam
      phi=exp(matmul(z,gfit))
      wm=1/phi
      call weighted_least_squares(x,y,wm,beta,fit,res,rank,info)
      d=res*res
      ldet=logdet_xtwx(x,wm,info)
      dev=sum(d/phi)+sum(log(phi))+n*log(2*acos(-1.0_dp))+ldet
      if(dev<devold-1e-15_dp)then
      gam=gfit
      exit
      end if
      lambda=lambda*2
      if(lambda/max(maxinfo,tiny(1.0_dp))>1e15_dp)exit
   end do
   if(lambda/max(maxinfo,tiny(1.0_dp))>1e15_dp)exit
   if(lev==1)lambda=lambda/10
   if(dot_product(dl,dgam)<tolerance)exit
end do
call symmetric_inverse(zvz,covg,info)
call info_inverse(x,wm,covb,info)
out%beta=beta
out%gamma=gam
out%mu=fit
out%phi=phi
out%h=h
out%deviance=dev
out%iter=iter
out%cov_gamma=covg
out%cov_beta=covb
out%se_gamma=sqrt([(covg(i,i),i=1,nq)])
out%se_beta=sqrt([(covb(i,i),i=1,p)])
end function

function remlscoregamma(y,x,z,tol,maxit) result(out)
! Default upstream link choices (log mean, log dispersion) are the numerical API here.
real(dp),intent(in)::y(:),x(:,:),z(:,:)
real(dp),intent(in),optional::tol
integer,intent(in),optional::maxit
type(reml_fit_result_t)::out
real(dp),allocatable::beta(:),mu(:),phi(:),gam(:),fit(:),res(:),d(:),q(:,:),h(:),z2(:,:)
real(dp),allocatable::q2(:,:),q2z(:,:),info_m(:,:),dl(:),dgam(:),covg(:,:),covb(:,:),w(:)
real(dp)::dev,tolerance,ldet,mean_d,phi0,extradisp,deltah
integer::n,p,nq,rank,info,iter,mx,i,j,k,j0
n=size(y)
p=size(x,2)
nq=size(z,2)
tolerance=1e-5_dp
if(present(tol))tolerance=tol
mx=40
if(present(maxit))mx=maxit
! Gamma GLM log link by IRLS.
call gamma_log_glm(x,y,spread(1.0_dp,1,n),beta,mu,rank,info)
d=2*((y-mu)/mu-log(y/mu))
mean_d=sum(d)/n
phi0=-1/canonic_digamma(mean_d)*real(n,dp)/real(n-p,dp)
call least_squares(z,spread(log(phi0),1,n),gam,fit,res,rank,info)
phi=exp(matmul(z,gam))
w=1/phi
call gamma_log_glm(x,y,w,beta,mu,rank,info)
d=2*((y-mu)/mu-log(y/mu))
ldet=logdet_xtwx(x,w,info)
dev=4*sum(log(y))+sum(2*(log_gamma(1/phi)+(1+log(phi))/phi)+d/phi)+ldet
allocate(info_m(nq,nq),dl(nq),dgam(nq),q2(n,p*(p+1)/2))
do iter=1,mx
   call weighted_hat_basis(x,w,q,rank,info)
   h=sum(q*q,dim=2)
   z2=z/sqrt(2.0_dp) ! log link: phidot=phi Z, so Z2=Z/sqrt(2)
   j0=0
   do k=0,p-1
   do j=1,p-k
   j0=j0+1
   q2(:,j0)=q(:,j)*q(:,j+k)
   if(k>0)q2(:,j0)=sqrt(2.0_dp)*q2(:,j0)
   end do
   end do
   q2z=matmul(transpose(q2),z2)
   info_m=0
   do i=1,n
      if(h(i)>tiny(1.0_dp))then
      extradisp=2*(r_trigamma(1/phi(i))-r_trigamma(1/(h(i)*phi(i)))/h(i))/phi(i)**2-(1-h(i))
      else
      extradisp=0
      end if
      do j=1,nq
      info_m(j,:)=info_m(j,:)+(extradisp+1-2*h(i))*z2(i,j)*z2(i,:)
      end do
   end do
   info_m=info_m+matmul(transpose(q2z),q2z)
   dl=0
   do i=1,n
      if(h(i)>tiny(1.0_dp))then
      deltah=2*(r_digamma(1/(h(i)*phi(i)))+log(h(i))-r_digamma(1/phi(i)))
      else
      deltah=0
      end if
      dl=dl+z(i,:)*(d(i)-deltah)/(2*phi(i)) ! phidot=phi*z, divide phi^2
   end do
   call solve_spd(info_m,dl,dgam,info)
   if(info/=0)exit
   gam=gam+dgam
   phi=exp(matmul(z,gam))
   w=1/phi
   call gamma_log_glm(x,y,w,beta,mu,rank,info)
   d=2*((y-mu)/mu-log(y/mu))
   ldet=logdet_xtwx(x,w,info)
   dev=4*sum(log(y))+sum(2*(log_gamma(1/phi)+(1+log(phi))/phi)+d/phi)+ldet
   if(dot_product(dl,dgam)<tolerance)exit
end do
call symmetric_inverse(info_m,covg,info)
call info_inverse(x,w,covb,info)
out%beta=beta
out%gamma=gam
out%mu=mu
out%phi=phi
out%h=h
out%deviance=dev
out%iter=iter
out%cov_gamma=covg
out%cov_beta=covb
out%se_gamma=sqrt([(covg(i,i),i=1,nq)])
out%se_beta=sqrt([(covb(i,i),i=1,p)])
end function

subroutine gamma_log_glm(x,y,prior,beta,mu,rank,info)
real(dp),intent(in)::x(:,:),y(:),prior(:)
real(dp),allocatable,intent(out)::beta(:),mu(:)
integer,intent(out)::rank,info
real(dp),allocatable::eta(:),z(:),fit(:),res(:),w(:),old(:)
integer::iter
call least_squares(x,log(max(y,1e-6_dp)),beta,fit,res,rank,info)
do iter=1,50
   old=beta
   eta=matmul(x,beta)
   mu=exp(eta)
   z=eta+(y-mu)/mu
   w=prior
   call weighted_least_squares(x,z,w,beta,fit,res,rank,info)
   if(maxval(abs(beta-old))<1e-10_dp)exit
end do
mu=exp(matmul(x,beta))
end subroutine

subroutine unique_int_local(x,u)
integer,intent(in)::x(:)
integer,allocatable,intent(out)::u(:)
integer::i,n
allocate(u(size(x)))
n=0
do i=1,size(x)
if(n==0.or..not.any(u(1:n)==x(i)))then
n=n+1
u(n)=x(i)
end if
end do
u=u(:n)
end subroutine

end module statmod_models
