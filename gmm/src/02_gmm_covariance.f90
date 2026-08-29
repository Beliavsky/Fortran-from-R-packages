! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_covariance
use r_compat, only: dp, ar_fit_t, ar_fit, arima_fit_t, arima_fit
use gmm_linalg, only: center_columns, colmeans_mat, invert_matrix, symmetrize
implicit none
private
public :: moment_covariance, hac_covariance, kernel_weight, smooth_moments
public :: bw_andrews, bw_wilhelm

contains

pure function kernel_weight(u,kernel) result(w)
real(dp),intent(in)::u
character(len=*),intent(in)::kernel
real(dp)::w,z,pi
pi=acos(-1.0_dp)
z=abs(u)
select case(trim(kernel))
case('Bartlett','bartlett')
   if(z<=1.0_dp) then
   w=1.0_dp-z
   else
   w=0.0_dp
   end if
case('Parzen','parzen')
   if(z<=0.5_dp) then
      w=1.0_dp-6.0_dp*z*z+6.0_dp*z*z*z
   else if(z<=1.0_dp) then
      w=2.0_dp*(1.0_dp-z)**3
   else
      w=0.0_dp
   end if
case('Tukey-Hanning','tukey-hanning','Tukey')
   if(z<=1.0_dp) then
   w=0.5_dp*(1.0_dp+cos(pi*z))
   else
   w=0.0_dp
   end if
case('Truncated','truncated')
   if(z<=1.0_dp) then
   w=1.0_dp
   else
   w=0.0_dp
   end if
case('Quadratic Spectral','QS','qs')
   if(z<=sqrt(epsilon(1.0_dp))) then
      w=1.0_dp
   else
      z=6.0_dp*pi*z/5.0_dp
      w=25.0_dp/(12.0_dp*pi*pi*(abs(u)**2))*(sin(z)/z-cos(z))
   end if
case default
   w=0.0_dp
end select
end function kernel_weight

pure function moment_covariance(gt,centered) result(s)
real(dp),intent(in)::gt(:,:)
logical,intent(in),optional::centered
real(dp)::s(size(gt,2),size(gt,2)),x(size(gt,1),size(gt,2))
logical::ctr
ctr=.true.
if(present(centered)) ctr=centered
if(ctr) then
x=center_columns(gt)
else
x=gt
end if
s=matmul(transpose(x),x)/real(size(x,1),dp)
s=symmetrize(s)
end function moment_covariance

pure function hac_covariance(gt,bw,kernel,centered) result(s)
real(dp),intent(in)::gt(:,:),bw
character(len=*),intent(in)::kernel
logical,intent(in),optional::centered
real(dp)::s(size(gt,2),size(gt,2)),x(size(gt,1),size(gt,2))
real(dp)::gam(size(gt,2),size(gt,2)),w
integer::n,l,lmax
logical::ctr
ctr=.true.
if(present(centered)) ctr=centered
if(ctr) then
x=center_columns(gt)
else
x=gt
end if
n=size(x,1)
s=matmul(transpose(x),x)/real(n,dp)
if(trim(kernel)=='Quadratic Spectral' .or. trim(kernel)=='QS' .or. trim(kernel)=='qs') then
   lmax=n-1
else
   lmax=min(n-1,max(0,ceiling(max(bw,0.0_dp))))
end if
if(bw<=0.0_dp) return
do l=1,lmax
   w=kernel_weight(real(l,dp)/bw,kernel)
   if(abs(w)<=epsilon(1.0_dp)) cycle
   gam=matmul(transpose(x(l+1:n,:)),x(1:n-l,:))/real(n,dp)
   s=s+w*(gam+transpose(gam))
end do
s=symmetrize(s)
end function hac_covariance

pure function smooth_moments(x,bw,kernel) result(y)
real(dp),intent(in)::x(:,:),bw
character(len=*),intent(in)::kernel
real(dp)::y(size(x,1),size(x,2)),den,w
integer::i,j,l
! Symmetric kernel smoother corresponding to the computational part of smoothG.
y=0.0_dp
do i=1,size(x,1)
   den=0.0_dp
   do j=1,size(x,1)
      l=abs(i-j)
      if(bw<=0.0_dp) then
         w=merge(1.0_dp,0.0_dp,l==0)
      else
         w=kernel_weight(real(l,dp)/bw,kernel)
      end if
      if(w/=0.0_dp) then
      y(i,:)=y(i,:)+w*x(j,:)
      den=den+w
      end if
   end do
   if(den>0.0_dp) y(i,:)=y(i,:)/den
end do
end function smooth_moments

function bw_andrews(x,kernel,approx,weights) result(bw)
real(dp),intent(in)::x(:,:)
character(len=*),intent(in)::kernel
character(len=*),intent(in),optional::approx
real(dp),intent(in),optional::weights(:)
real(dp)::bw
real(dp),allocatable::wt(:)
real(dp)::rho(size(x,2)),sig(size(x,2)),psi(size(x,2)),denom,a1,a2,c,powr
integer::j,n,k
character(len=16)::ap
n=size(x,1)
k=size(x,2)
ap='AR(1)'
if(present(approx)) ap=trim(approx)
allocate(wt(k))
wt=1.0_dp
if(present(weights)) wt=weights(1:k)
rho=0.0_dp
sig=0.0_dp
psi=0.0_dp
if(trim(ap)=='ARMA(1,1)') then
   do j=1,k
      block
         type(arima_fit_t)::fit
         fit=arima_fit(x(:,j),[1,0,1],.false.)
         if(size(fit%coef)>=2) then
         rho(j)=fit%coef(1)
         psi(j)=fit%coef(2)
         end if
         sig(j)=sqrt(max(fit%sigma2,tiny(1.0_dp)))
      end block
   end do
   ! Andrews plug-in using the implied ARMA(1,1) zero-frequency spectrum.
   denom=sum(wt*((1.0_dp+psi)**4)*sig**4/max((1.0_dp-rho)**4,tiny(1.0_dp)))
   a1=sum(wt*4.0_dp*(1.0_dp+psi*rho)**2*(psi+rho)**2*sig**4 / &
      max((1.0_dp-rho)**6*(1.0_dp+rho)**2,tiny(1.0_dp)))/max(denom,tiny(1.0_dp))
   a2=sum(wt*4.0_dp*(1.0_dp+psi*rho)**2*(psi+rho)**2*sig**4 / &
      max((1.0_dp-rho)**8,tiny(1.0_dp)))/max(denom,tiny(1.0_dp))
else
   do j=1,k
      block
         type(ar_fit_t)::fit
         fit=ar_fit(x(:,j),order_max=1,aic=.false.,method='ols')
         if(size(fit%ar)>=1) rho(j)=fit%ar(1)
         sig(j)=sqrt(max(fit%var_pred,tiny(1.0_dp)))
      end block
   end do
   denom=sum(wt*sig**4/max((1.0_dp-rho)**4,tiny(1.0_dp)))
   a1=sum(wt*4.0_dp*rho**2*sig**4/max((1.0_dp-rho)**6*(1.0_dp+rho)**2,tiny(1.0_dp)))/max(denom,tiny(1.0_dp))
   a2=sum(wt*4.0_dp*rho**2*sig**4/max((1.0_dp-rho)**8,tiny(1.0_dp)))/max(denom,tiny(1.0_dp))
end if
select case(trim(kernel))
case('Bartlett','bartlett')
c=1.1447_dp
powr=1.0_dp/3.0_dp
bw=c*(real(n,dp)*max(a1,0.0_dp))**powr
case('Parzen','parzen')
c=2.6614_dp
powr=1.0_dp/5.0_dp
bw=c*(real(n,dp)*max(a2,0.0_dp))**powr
case('Tukey-Hanning','Tukey')
c=1.7462_dp
powr=1.0_dp/5.0_dp
bw=c*(real(n,dp)*max(a2,0.0_dp))**powr
case default
c=1.3221_dp
powr=1.0_dp/5.0_dp
bw=c*(real(n,dp)*max(a2,0.0_dp))**powr
end select
end function bw_andrews

function bw_wilhelm(umat,G,kernel,approx,weights) result(bw)
real(dp),intent(in)::umat(:,:),G(:,:)
character(len=*),intent(in)::kernel
character(len=*),intent(in),optional::approx
real(dp),intent(in),optional::weights(:)
real(dp)::bw
real(dp),allocatable::wt(:),omega0(:),omegaq(:,:),oinv(:,:),sigma0(:,:),h0(:,:),p0(:,:),W(:,:),tmp(:,:)
real(dp)::rho(size(umat,2)),sig(size(umat,2)),psi(size(umat,2)),qv,gq,mu1,mu2,nu2,nu3,c0
integer::j,k,p,n,info
character(len=16)::ap
n=size(umat,1)
k=size(umat,2)
p=size(G,2)
ap='AR(1)'
if(present(approx)) ap=trim(approx)
if(p==k) then
bw=bw_andrews(umat,kernel,ap,weights)
return
end if
allocate(wt(p))
wt=1.0_dp
if(present(weights)) wt=weights(1:p)
rho=0.0_dp
sig=0.0_dp
psi=0.0_dp
if(trim(ap)=='ARMA(1,1)') then
   do j=1,k
      block
         type(arima_fit_t)::fit
         fit=arima_fit(umat(:,j),[1,0,1],.false.)
         if(size(fit%coef)>=2) then
         rho(j)=fit%coef(1)
         psi(j)=fit%coef(2)
         end if
         sig(j)=sqrt(max(fit%sigma2,tiny(1.0_dp)))
      end block
   end do
else
   do j=1,k
      block
         type(ar_fit_t)::fit
         fit=ar_fit(umat(:,j),order_max=1,aic=.false.,method='ols')
         if(size(fit%ar)>=1) rho(j)=fit%ar(1)
         sig(j)=sqrt(max(fit%var_pred,tiny(1.0_dp)))
      end block
   end do
end if
select case(trim(kernel))
case('Bartlett')
qv=1
gq=1
mu1=1
mu2=2.0_dp/3.0_dp
case('Parzen')
qv=2
gq=6
mu1=0.75_dp
mu2=0.539286_dp
case('Tukey-Hanning')
qv=2
gq=acos(-1.0_dp)**2/4
mu1=1
mu2=0.75_dp
case default
qv=2
gq=1.421223_dp
mu1=1.25003_dp
mu2=0.999985_dp
end select
allocate(omega0(k),omegaq(k,k))
omegaq=0.0_dp
if(trim(ap)=='ARMA(1,1)') then
   omega0=(1.0_dp+psi)**2*sig**2/max((1.0_dp-rho)**2,tiny(1.0_dp))
   do j=1,k
      if(qv==1) then
         omegaq(j,j)=2*(1+psi(j)*rho(j))*(psi(j)+rho(j))*sig(j)**2 / &
            max((1-rho(j))**3*(1+rho(j)),tiny(1.0_dp))
      else
         omegaq(j,j)=2*(1+psi(j)*rho(j))*(psi(j)+rho(j))*sig(j)**2/max((1-rho(j))**4,tiny(1.0_dp))
      end if
   end do
else
   omega0=sig**2/max((1.0_dp-rho)**2,tiny(1.0_dp))
   do j=1,k
      if(qv==1) then
         omegaq(j,j)=2*sig(j)**2*rho(j)/max((1-rho(j))**3*(1+rho(j)),tiny(1.0_dp))
      else
         omegaq(j,j)=2*sig(j)**2*rho(j)/max((1-rho(j))**4,tiny(1.0_dp))
      end if
   end do
end if
allocate(oinv(k,k))
oinv=0.0_dp
do j=1,k
oinv(j,j)=1.0_dp/max(omega0(j),tiny(1.0_dp))
end do
allocate(tmp(p,p))
tmp=matmul(transpose(G),matmul(oinv,G))
call invert_matrix(tmp,sigma0,info)
if(info/=0) then
bw=bw_andrews(umat,kernel,ap)
return
end if
h0=matmul(sigma0,matmul(transpose(G),oinv))
p0=oinv-matmul(oinv,matmul(G,h0))
allocate(W(p,p))
W=0.0_dp
do j=1,p
W(j,j)=wt(j)
end do
nu2=(2*mu1+mu2)*real(k-p,dp)*sum([(sum(sigma0(j,:)*W(:,j)),j=1,p)])
tmp=matmul(transpose(omegaq),matmul(transpose(h0),matmul(W,matmul(h0,matmul(omegaq,p0)))))
nu3=gq**2*sum([(tmp(j,j),j=1,k)])
if(nu2*nu3>0) then
c0=2*qv
else
c0=-1
end if
if(c0*nu3/nu2>0) then
   bw=(c0*nu3/nu2*real(n,dp))**(1.0_dp/(1.0_dp+2*qv))
else
   bw=0.0_dp
end if
end function bw_wilhelm

end module gmm_covariance
