module ld_diagnostics_catalog
use ld_kinds, only: dp, pi
implicit none
private
public :: heidelberger_result_t, bmk_diagnostic, heidelberger_diagnostic, ks_diagnostic

type :: heidelberger_result_t
   logical, allocatable :: stationarity_pass(:), halfwidth_pass(:)
   integer, allocatable :: start(:)
   real(dp), allocatable :: pvalue(:), mean(:), halfwidth(:)
end type heidelberger_result_t

contains

subroutine bmk_diagnostic(chain,batches,hellinger)
   real(dp), intent(in) :: chain(:,:)
   integer, intent(in) :: batches
   real(dp), allocatable, intent(out) :: hellinger(:,:)
   integer :: n,p,bs,j,b
   n=size(chain,1); p=size(chain,2); if(batches<2) error stop 'bmk_diagnostic: batches must be >=2'
   if(mod(n,batches)/=0) error stop 'bmk_diagnostic: equal batch sizes are required'
   bs=n/batches; allocate(hellinger(p,batches-1))
   do j=1,p
      do b=1,batches-1
         hellinger(j,b)=hellinger_distance(chain((b-1)*bs+1:b*bs,j),chain(b*bs+1:(b+1)*bs,j))
      end do
   end do
   hellinger=min(1.0_dp,max(0.0_dp,hellinger))
end subroutine bmk_diagnostic

subroutine heidelberger_diagnostic(chain,res,eps,pvalue_cut)
   real(dp), intent(in) :: chain(:,:)
   type(heidelberger_result_t), intent(out) :: res
   real(dp), intent(in), optional :: eps,pvalue_cut
   integer :: n,p,j,step,s,ncur,i
   real(dp) :: ee,pcut,s0,ybar,istat,pv,hw
   real(dp), allocatable :: y(:),bbridge(:)
   n=size(chain,1); p=size(chain,2); ee=0.1_dp; if(present(eps)) ee=eps
   pcut=0.05_dp; if(present(pvalue_cut)) pcut=pvalue_cut
   allocate(res%stationarity_pass(p),res%halfwidth_pass(p),res%start(p),res%pvalue(p),res%mean(p),res%halfwidth(p))
   res%stationarity_pass=.false.; res%halfwidth_pass=.false.; res%start=0
   res%pvalue=0.0_dp; res%mean=0.0_dp; res%halfwidth=huge(1.0_dp)
   step=max(1,n/10)
   do j=1,p
      s0=long_run_variance(chain(max(1,n/2):n,j)); if(s0<=0.0_dp) cycle
      do s=1,max(1,n/2),step
         y=chain(s:n,j); ncur=size(y); ybar=sum(y)/real(ncur,dp); allocate(bbridge(ncur))
         bbridge=cumsum_local(y-ybar); istat=sum((bbridge*bbridge)/(real(ncur,dp)*s0))/real(ncur,dp)
         pv=1.0_dp-cramer_cdf(istat); deallocate(bbridge)
         if(pv>pcut) then
            res%stationarity_pass(j)=.true.; res%start(j)=s; res%pvalue(j)=pv
            s0=long_run_variance(y); hw=1.96_dp*sqrt(max(s0,0.0_dp)/real(ncur,dp))
            res%mean(j)=ybar; res%halfwidth(j)=hw
            if(abs(ybar)>tiny(1.0_dp)) then; res%halfwidth_pass(j)=abs(hw/ybar)<=ee
            else; res%halfwidth_pass(j)=hw<=ee; end if
            exit
         end if
      end do
   end do
end subroutine heidelberger_diagnostic

subroutine ks_diagnostic(chain,statistic,pvalue)
   real(dp), intent(in) :: chain(:,:)
   real(dp), intent(out) :: statistic(:),pvalue(:)
   integer :: n,p,j,n1,n2
   real(dp), allocatable :: a(:),b(:)
   n=size(chain,1); p=size(chain,2); n1=n/2; n2=n-n1
   if(size(statistic)/=p .or. size(pvalue)/=p) error stop 'ks_diagnostic: output size mismatch'
   do j=1,p
      a=chain(1:n1,j); b=chain(n1+1:n,j); call sort_real(a); call sort_real(b)
      statistic(j)=ks_two_sample_stat(a,b)
      pvalue(j)=ks_asymptotic_p(statistic(j)*sqrt(real(n1*n2,dp)/real(n1+n2,dp)))
   end do
end subroutine ks_diagnostic

function hellinger_distance(a,b) result(h)
   real(dp), intent(in) :: a(:),b(:)
   real(dp) :: h,lo,hi,dx,bw1,bw2,x,pd,qd,ss
   integer :: m,i
   lo=min(minval(a),minval(b)); hi=max(maxval(a),maxval(b)); if(hi<=lo) then; h=0.0_dp; return; end if
   m=max(64,min(512,max(size(a),size(b)))); dx=(hi-lo)/real(m-1,dp); bw1=kde_bandwidth(a); bw2=kde_bandwidth(b); ss=0.0_dp
   do i=1,m
      x=lo+real(i-1,dp)*dx; pd=kde_value(a,x,bw1); qd=kde_value(b,x,bw2)
      ss=ss+(sqrt(max(pd,0.0_dp))-sqrt(max(qd,0.0_dp)))**2*dx
   end do
   h=sqrt(0.5_dp*ss)
end function hellinger_distance

function kde_bandwidth(x) result(bw)
   real(dp), intent(in) :: x(:)
   real(dp) :: bw,m,sd
   m=sum(x)/real(size(x),dp)
   if(size(x)>1) then; sd=sqrt(sum((x-m)**2)/real(size(x)-1,dp)); else; sd=1.0_dp; end if
   bw=1.06_dp*max(sd,1.0e-8_dp)*real(size(x),dp)**(-0.2_dp); bw=max(bw,1.0e-8_dp)
end function kde_bandwidth

function kde_value(x,z,bw) result(v)
   real(dp), intent(in) :: x(:),z,bw
   real(dp) :: v
   v=sum(exp(-0.5_dp*((z-x)/bw)**2))/(real(size(x),dp)*bw*sqrt(2.0_dp*pi))
end function kde_value

function long_run_variance(x) result(v)
   real(dp), intent(in) :: x(:)
   real(dp) :: v,m,g0,g1,g2,pair
   integer :: n,k
   n=size(x); if(n<3) then; v=0.0_dp; return; end if; m=sum(x)/real(n,dp); g0=sum((x-m)**2)/real(n,dp); v=g0
   k=1
   do while(k+1<n)
      g1=sum((x(1:n-k)-m)*(x(1+k:n)-m))/real(n,dp)
      g2=sum((x(1:n-k-1)-m)*(x(2+k:n)-m))/real(n,dp); pair=g1+g2
      if(pair<=0.0_dp) exit; v=v+2.0_dp*pair; k=k+2
   end do
   v=max(v,0.0_dp)
end function long_run_variance

function cramer_cdf(q) result(p)
   real(dp), intent(in) :: q
   real(dp) :: p,z,u,term
   integer :: k
   if(q<=0.0_dp) then; p=0.0_dp; return; end if
   p=0.0_dp
   do k=0,3
      z=gamma(real(k,dp)+0.5_dp)*sqrt(real(4*k+1,dp))/(gamma(real(k+1,dp))*pi**1.5_dp*sqrt(q))
      u=real((4*k+1)**2,dp)/(16.0_dp*q)
      if(u<80.0_dp) then; term=z*exp(-u)*bessel_k_numeric(u,0.25_dp); p=p+term; end if
   end do
   p=min(1.0_dp,max(0.0_dp,p))
end function cramer_cdf

function bessel_k_numeric(x,nu) result(v)
   real(dp), intent(in) :: x,nu
   real(dp) :: v,h,t,sumv,term
   integer, parameter :: n=400
   integer :: i
   if(x<=0.0_dp) then; v=huge(1.0_dp); return; end if
   h=12.0_dp/real(n,dp); sumv=0.0_dp
   do i=0,n
      t=real(i,dp)*h; term=exp(-x*cosh(t))*cosh(nu*t)
      if(i==0 .or. i==n) then; sumv=sumv+term; else if(mod(i,2)==0) then; sumv=sumv+2.0_dp*term; else; sumv=sumv+4.0_dp*term; end if
   end do
   v=h*sumv/3.0_dp
end function bessel_k_numeric

pure function cumsum_local(x) result(y)
   real(dp), intent(in) :: x(:)
   real(dp) :: y(size(x))
   integer :: i
   if(size(x)==0) return; y(1)=x(1); do i=2,size(x); y(i)=y(i-1)+x(i); end do
end function cumsum_local

subroutine sort_real(x)
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
end subroutine sort_real

function ks_two_sample_stat(a,b) result(d)
   real(dp), intent(in) :: a(:),b(:)
   real(dp) :: d,fa,fb
   integer :: i,j,n1,n2
   n1=size(a); n2=size(b); i=1; j=1; fa=0.0_dp; fb=0.0_dp; d=0.0_dp
   do while(i<=n1 .or. j<=n2)
      if(j>n2) then
         i=i+1; fa=real(i-1,dp)/real(n1,dp)
      else if(i>n1) then
         j=j+1; fb=real(j-1,dp)/real(n2,dp)
      else if(a(i)<=b(j)) then
         i=i+1; fa=real(i-1,dp)/real(n1,dp)
      else
         j=j+1; fb=real(j-1,dp)/real(n2,dp)
      end if
      d=max(d,abs(fa-fb))
   end do
end function ks_two_sample_stat

function ks_asymptotic_p(lambda) result(p)
   real(dp), intent(in) :: lambda
   real(dp) :: p,term
   integer :: k
   if(lambda<=0.0_dp) then; p=1.0_dp; return; end if; p=0.0_dp
   do k=1,100
      term=2.0_dp*(-1.0_dp)**(k-1)*exp(-2.0_dp*real(k*k,dp)*lambda*lambda); p=p+term
      if(abs(term)<1.0e-14_dp) exit
   end do
   p=min(1.0_dp,max(0.0_dp,p))
end function ks_asymptotic_p

end module ld_diagnostics_catalog
