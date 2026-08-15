module ld_diagnostics
use ld_kinds, only: dp
implicit none
private
public :: acceptance_rate, iat, ess, mcse_imps, mcse_batch_means, geweke_z
public :: gelman_rhat, waic_result_t, waic, kld

type :: waic_result_t
   real(dp) :: waic=0.0_dp, lppd=0.0_dp, p_waic=0.0_dp, p_waic1=0.0_dp
end type
contains
function mean1(x) result(m)
   real(dp),intent(in)::x(:); real(dp)::m
   m=sum(x)/real(size(x),dp)
end function mean1
function var1(x) result(v)
   real(dp),intent(in)::x(:); real(dp)::v,m
   if(size(x)<2) then; v=0.0_dp; return; end if
   m=mean1(x); v=sum((x-m)**2)/real(size(x)-1,dp)
end function var1
function autocov(x,lag) result(v)
   real(dp),intent(in)::x(:); integer,intent(in)::lag; real(dp)::v,m
   integer::n
   n=size(x); if(lag>=n) then; v=0.0_dp; return; end if
   m=mean1(x); v=sum((x(1:n-lag)-m)*(x(1+lag:n)-m))/real(n,dp)
end function autocov

function acceptance_rate(chain) result(r)
   real(dp),intent(in)::chain(:,:); real(dp)::r(size(chain,2)); integer::i,j,n
   n=size(chain,1); r=0.0_dp; if(n<2) return
   do j=1,size(chain,2); do i=1,n-1; if(chain(i,j)/=chain(i+1,j)) r(j)=r(j)+1.0_dp; end do; r(j)=r(j)/real(n-1,dp); end do
end function acceptance_rate

function iat(x) result(tau)
   real(dp),intent(in)::x(:); real(dp)::tau,s2,gprev,gcur,s
   integer::n,m,maxlag
   n=size(x); if(n<4) then; tau=1.0_dp; return; end if
   s2=var1(x); if(s2<=0.0_dp) then; tau=1.0_dp; return; end if
   maxlag=max(3,n/2)
   gprev=s2+autocov(x,1); s=gprev/s2; m=1
   do
      if(2*m+1>maxlag) exit
      gcur=autocov(x,2*m)+autocov(x,2*m+1)
      if(gcur<=0.0_dp .or. gcur>=gprev) exit
      s=s+gcur/s2; gprev=gcur; m=m+1
   end do
   tau=max(1.0_dp,-1.0_dp+2.0_dp*s)
end function iat

function ess(x) result(v)
   real(dp),intent(in)::x(:); real(dp)::v
   v=min(real(size(x),dp),max(epsilon(1.0_dp),real(size(x),dp)/iat(x)))
end function ess

function mcse_imps(x) result(se)
   real(dp),intent(in)::x(:); real(dp)::se
   se=sqrt(max(var1(x),0.0_dp)*iat(x)/real(size(x),dp))
end function mcse_imps

subroutine mcse_batch_means(x,estimate,se,batch_size)
   real(dp),intent(in)::x(:); real(dp),intent(out)::estimate,se; integer,intent(in),optional::batch_size
   integer::n,b,a,k; real(dp),allocatable::y(:)
   n=size(x); b=max(2,int(sqrt(real(n,dp)))); if(present(batch_size)) b=max(2,batch_size); a=n/b
   if(a<2) then; estimate=mean1(x); se=sqrt(var1(x)/real(max(n,1),dp)); return; end if
   allocate(y(a)); do k=1,a; y(k)=mean1(x((k-1)*b+1:k*b)); end do
   estimate=mean1(y); se=sqrt(real(b,dp)*sum((y-estimate)**2)/real(a-1,dp)/real(n,dp))
end subroutine mcse_batch_means

function geweke_z(x) result(z)
   real(dp),intent(in)::x(:); real(dp)::z
   integer::n,n1,n2,s2
   real(dp)::m1,m2,v1,v2
   n=size(x); if(n<20) then; z=0.0_dp; return; end if
   n1=max(2,int(0.1_dp*n)); n2=max(2,int(0.5_dp*n)); s2=n-n2+1
   m1=mean1(x(1:n1)); m2=mean1(x(s2:n)); v1=var1(x(1:n1))*iat(x(1:n1))/real(n1,dp); v2=var1(x(s2:n))*iat(x(s2:n))/real(n2,dp)
   if(v1+v2>0.0_dp) then; z=(m1-m2)/sqrt(v1+v2); else; z=0.0_dp; end if
end function geweke_z

subroutine gelman_rhat(chains,rhat)
   real(dp),intent(in)::chains(:,:,:)
   real(dp),intent(out)::rhat(:)
   integer::n,m,p,i,j
   real(dp),allocatable::means(:),vars(:)
   real(dp)::w,b,varhat,grand
   n=size(chains,1); m=size(chains,2); p=size(chains,3); allocate(means(m),vars(m))
   do j=1,p
      do i=1,m; means(i)=mean1(chains(:,i,j)); vars(i)=var1(chains(:,i,j)); end do
      w=mean1(vars); grand=mean1(means); b=real(n,dp)*sum((means-grand)**2)/real(max(1,m-1),dp)
      varhat=(real(n-1,dp)/real(n,dp))*w+b/real(n,dp)
      if(w>0.0_dp) then; rhat(j)=sqrt(varhat/w); else; rhat(j)=1.0_dp; end if
   end do
end subroutine gelman_rhat

subroutine waic(log_lik,res)
   real(dp),intent(in)::log_lik(:,:)
   type(waic_result_t),intent(out)::res
   integer::i,n,s
   real(dp)::m,lv,lr
   n=size(log_lik,1); s=size(log_lik,2); res%lppd=0.0_dp; res%p_waic=0.0_dp; res%p_waic1=0.0_dp
   do i=1,n
      m=maxval(log_lik(i,:)); lr=m+log(sum(exp(log_lik(i,:)-m))/real(s,dp)); lv=var1(log_lik(i,:))
      res%lppd=res%lppd+lr; res%p_waic=res%p_waic+lv; res%p_waic1=res%p_waic1+2.0_dp*(lr-mean1(log_lik(i,:)))
   end do
   res%waic=-2.0_dp*res%lppd+2.0_dp*res%p_waic
end subroutine waic

function kld(px,py,base) result(v)
   real(dp),intent(in)::px(:),py(:); real(dp),intent(in),optional::base; real(dp)::v,b
   integer::i
   b=exp(1.0_dp); if(present(base)) b=base; v=0.0_dp
   do i=1,size(px); if(px(i)>0.0_dp .and. py(i)>0.0_dp) v=v+px(i)*log(px(i)/py(i))/log(b); end do
end function kld
end module ld_diagnostics
