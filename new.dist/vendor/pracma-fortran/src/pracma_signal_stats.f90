! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_signal_stats
   use pracma_kinds, only : dp, pi_dp
   use pracma_status, only : pracma_ok, pracma_invalid_argument
   use pracma_types, only : peak_result, regression_result
   use pracma_linalg, only : solve_linear, pinv
   use pracma_polynomial, only : polyfit, polyval
   implicit none
   private
   public :: conv, deconv, fft, ifft, fftshift, ifftshift, detrend, movavg, savgol, savgol_filter
   public :: findpeaks, hampel, entropy, hurst, histc, rmserr, runge, humps_function
   public :: mexpfit, odregress, trigregress, whittaker, andrews_curve, normest
   public :: autocorrelation, crosscorrelation, moving_median, moving_mean, periodogram

contains

   function conv(x,y,shape) result(z)
      real(dp),intent(in)::x(:),y(:)
      character(len=*),intent(in),optional::shape
      real(dp),allocatable::z(:),full(:)
      character(len=8)::s
      integer::i,j,n,m,start,lenz
      n=size(x); m=size(y); allocate(full(max(0,n+m-1))); full=0.0_dp
      do i=1,n; do j=1,m; full(i+j-1)=full(i+j-1)+x(i)*y(j); end do; end do
      s='full'; if(present(shape))s=adjustl(shape)
      select case(trim(s))
      case('same')
         lenz=n; start=(m-1)/2+1; allocate(z(lenz)); z=full(start:start+lenz-1)
      case('valid')
         lenz=max(n,m)-min(n,m)+1; start=min(n,m); allocate(z(lenz)); z=full(start:start+lenz-1)
      case default
         z=full
      end select
   end function conv

   subroutine deconv(y,a,q,r,status)
      real(dp),intent(in)::y(:),a(:)
      real(dp),allocatable,intent(out)::q(:),r(:)
      integer,intent(out),optional::status
      integer::n,m,i
      real(dp),allocatable::rr(:)
      n=size(y); m=size(a)
      if(m<1.or.n<m.or.abs(a(1))<=tiny(1.0_dp))then
         allocate(q(0),r(n)); r=y; if(present(status))status=pracma_invalid_argument; return
      end if
      allocate(q(n-m+1),rr(n)); rr=y; q=0.0_dp
      do i=1,n-m+1
         q(i)=rr(i)/a(1); rr(i:i+m-1)=rr(i:i+m-1)-q(i)*a
      end do
      r=rr; if(present(status))status=pracma_ok
   end subroutine deconv

   function fft(x,inverse) result(z)
      complex(dp),intent(in)::x(:)
      logical,intent(in),optional::inverse
      complex(dp),allocatable::z(:)
      logical::inv
      integer::j,k,n
      real(dp)::sgn,theta
      inv=.false.; if(present(inverse))inv=inverse
      n=size(x); allocate(z(n)); z=(0.0_dp,0.0_dp); sgn=merge(1.0_dp,-1.0_dp,inv)
      do k=0,n-1
         do j=0,n-1
            theta=sgn*2.0_dp*pi_dp*real(j*k,dp)/real(n,dp)
            z(k+1)=z(k+1)+x(j+1)*cmplx(cos(theta),sin(theta),dp)
         end do
      end do
      if(inv)z=z/real(n,dp)
   end function fft

   function ifft(x) result(z)
      complex(dp),intent(in)::x(:)
      complex(dp),allocatable::z(:)
      z=fft(x,.true.)
   end function ifft

   function fftshift(x) result(y)
      complex(dp),intent(in)::x(:)
      complex(dp),allocatable::y(:)
      integer::n,k
      n=size(x); k=(n+1)/2; allocate(y(n)); y=[x(k+1:n),x(1:k)]
   end function fftshift

   function ifftshift(x) result(y)
      complex(dp),intent(in)::x(:)
      complex(dp),allocatable::y(:)
      integer::n,k
      n=size(x); k=n/2; allocate(y(n)); y=[x(k+1:n),x(1:k)]
   end function ifftshift

   function detrend(x,kind) result(y)
      real(dp),intent(in)::x(:)
      character(len=*),intent(in),optional::kind
      real(dp),allocatable::y(:)
      character(len=16)::k
      real(dp)::sx,sy,sxx,sxy,b0,b1,nr
      integer::i,n
      n=size(x); allocate(y(n)); k='linear'; if(present(kind))k=adjustl(kind)
      if(trim(k)=='constant')then; y=x-sum(x)/real(n,dp); return; end if
      nr=real(n,dp); sx=nr*(nr+1)/2; sxx=nr*(nr+1)*(2*nr+1)/6; sy=sum(x); sxy=0.0_dp
      do i=1,n; sxy=sxy+real(i,dp)*x(i); end do
      b1=(nr*sxy-sx*sy)/(nr*sxx-sx*sx); b0=(sy-b1*sx)/nr
      do i=1,n; y(i)=x(i)-b0-b1*real(i,dp); end do
   end function detrend

   function movavg(x,n,centered) result(y)
      real(dp),intent(in)::x(:)
      integer,intent(in)::n
      logical,intent(in),optional::centered
      real(dp),allocatable::y(:)
      logical::ctr
      integer::i,lo,hi
      ctr=.true.; if(present(centered))ctr=centered
      allocate(y(size(x))); y=0.0_dp
      do i=1,size(x)
         if(ctr)then; lo=max(1,i-(n-1)/2); hi=min(size(x),i+n/2)
         else; lo=max(1,i-n+1); hi=i; end if
         y(i)=sum(x(lo:hi))/real(hi-lo+1,dp)
      end do
   end function movavg

   function moving_mean(x,n) result(y)
      real(dp),intent(in)::x(:); integer,intent(in)::n; real(dp),allocatable::y(:)
      y=movavg(x,n,.true.)
   end function moving_mean

   function moving_median(x,n) result(y)
      real(dp),intent(in)::x(:); integer,intent(in)::n
      real(dp),allocatable::y(:),w(:)
      integer::i,lo,hi
      allocate(y(size(x)))
      do i=1,size(x)
         lo=max(1,i-(n-1)/2); hi=min(size(x),i+n/2); w=x(lo:hi); call sort_inplace(w)
         if(mod(size(w),2)==1)then; y(i)=w((size(w)+1)/2); else; y(i)=0.5_dp*(w(size(w)/2)+w(size(w)/2+1)); end if
      end do
   end function moving_median

   function savgol(x,window,degree,derivative,delta) result(y)
      real(dp),intent(in)::x(:)
      integer,intent(in)::window,degree
      integer,intent(in),optional::derivative
      real(dp),intent(in),optional::delta
      real(dp),allocatable::y(:),a(:,:),ata(:,:),rhs(:),coef(:),xx(:),yy(:)
      integer::i,j,k,lo,hi,m,d,istat
      real(dp)::del,fac
      m=max(3,window); if(mod(m,2)==0)m=m+1; d=0; if(present(derivative))d=derivative
      del=1.0_dp; if(present(delta))del=delta; allocate(y(size(x)))
      do i=1,size(x)
         lo=max(1,min(i-(m-1)/2,size(x)-m+1)); hi=min(size(x),lo+m-1); lo=max(1,hi-m+1)
         allocate(xx(hi-lo+1),yy(hi-lo+1),a(hi-lo+1,degree+1),ata(degree+1,degree+1),rhs(degree+1),coef(degree+1))
         do j=1,size(xx); xx(j)=real(lo+j-1-i,dp)*del; yy(j)=x(lo+j-1); a(j,1)=1.0_dp
            do k=2,degree+1; a(j,k)=a(j,k-1)*xx(j); end do
         end do
         ata=matmul(transpose(a),a); rhs=matmul(transpose(a),yy); call solve_linear(ata,rhs,coef,istat)
         if(d<=degree)then; fac=1.0_dp; do k=2,d; fac=fac*real(k,dp); end do; y(i)=fac*coef(d+1)
         else; y(i)=0.0_dp; end if
         deallocate(xx,yy,a,ata,rhs,coef)
      end do
   end function savgol

   function savgol_filter(x,window,degree) result(y)
      real(dp),intent(in)::x(:); integer,intent(in)::window,degree; real(dp),allocatable::y(:)
      y=savgol(x,window,degree)
   end function savgol_filter

   function findpeaks(x,min_height,min_distance) result(res)
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::min_height
      integer,intent(in),optional::min_distance
      type(peak_result)::res
      integer,allocatable::idx(:)
      integer::i,n,k,md,last
      real(dp)::mh,leftmin,rightmin
      n=size(x); mh=-huge(1.0_dp); if(present(min_height))mh=min_height; md=1; if(present(min_distance))md=max(1,min_distance)
      allocate(idx(n)); k=0; last=-md
      do i=2,n-1
         if(x(i)>=mh.and.x(i)>x(i-1).and.x(i)>=x(i+1).and.i-last>=md)then; k=k+1; idx(k)=i; last=i; end if
      end do
      allocate(res%indices(k),res%heights(k),res%prominences(k)); res%indices=idx(:k)
      do i=1,k
         res%heights(i)=x(idx(i)); leftmin=minval(x(1:idx(i))); rightmin=minval(x(idx(i):n))
         res%prominences(i)=x(idx(i))-max(leftmin,rightmin)
      end do
      res%status=pracma_ok
   end function findpeaks

   subroutine hampel(x,k,t0,filtered,outliers)
      real(dp),intent(in)::x(:)
      integer,intent(in)::k
      real(dp),intent(in),optional::t0
      real(dp),allocatable,intent(out)::filtered(:)
      logical,allocatable,intent(out)::outliers(:)
      real(dp),allocatable::w(:),dev(:)
      real(dp)::med,mad,thr
      integer::i,lo,hi
      thr=3.0_dp; if(present(t0))thr=t0; allocate(filtered(size(x)),outliers(size(x))); filtered=x; outliers=.false.
      do i=1,size(x)
         lo=max(1,i-k); hi=min(size(x),i+k); w=x(lo:hi); med=median(w); dev=abs(w-med); mad=median(dev)
         if(mad>0.0_dp.and.abs(x(i)-med)>thr*1.4826_dp*mad)then; filtered(i)=med; outliers(i)=.true.; end if
      end do
   end subroutine hampel

   function entropy(p,base) result(h)
      real(dp),intent(in)::p(:)
      real(dp),intent(in),optional::base
      real(dp)::h,b,s
      integer::i
      b=exp(1.0_dp); if(present(base))b=base; s=sum(max(p,0.0_dp)); h=0.0_dp
      if(s<=0.0_dp)return
      do i=1,size(p); if(p(i)>0.0_dp)h=h-(p(i)/s)*log(p(i)/s)/log(b); end do
   end function entropy

   function hurst(x,min_block,max_block) result(h)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::min_block,max_block
      real(dp)::h
      integer::mn,mx,nscale,s,k,b,j,nb,istat
      real(dp),allocatable::lx(:),ly(:),a(:,:),rhs(:),coef(:),seg(:),cum(:)
      real(dp)::rs,sdv,rng
      mn=8; if(present(min_block))mn=min_block; mx=size(x)/4; if(present(max_block))mx=max_block
      nscale=max(0,int(log(real(max(mx,mn),dp)/real(mn,dp))/log(2.0_dp))+1)
      if(nscale<2)then; h=0.5_dp; return; end if
      allocate(lx(nscale),ly(nscale)); k=0; s=mn
      do while(s<=mx.and.k<nscale)
         nb=size(x)/s; if(nb<1)exit; rs=0.0_dp
         do b=1,nb
            seg=x((b-1)*s+1:b*s); allocate(cum(s)); cum=0.0_dp
            do j=1,s; cum(j)=sum(seg(:j)-sum(seg)/real(s,dp)); end do
            rng=maxval(cum)-minval(cum); sdv=sqrt(sum((seg-sum(seg)/real(s,dp))**2)/real(max(1,s-1),dp))
            if(sdv>0.0_dp)rs=rs+rng/sdv; deallocate(cum)
         end do
         k=k+1; lx(k)=log(real(s,dp)); ly(k)=log(max(rs/real(nb,dp),tiny(1.0_dp))); s=2*s
      end do
      allocate(a(k,2),rhs(2),coef(2)); a(:,1)=1.0_dp; a(:,2)=lx(:k); rhs=matmul(transpose(a),ly(:k))
      call solve_linear(matmul(transpose(a),a),rhs,coef,istat); h=coef(2)
   end function hurst

   subroutine histc(x,edges,counts,bins)
      real(dp),intent(in)::x(:),edges(:)
      integer,allocatable,intent(out)::counts(:),bins(:)
      integer::i,j,m
      m=size(edges); allocate(counts(m),bins(size(x))); counts=0; bins=0
      do i=1,size(x)
         if(x(i)==edges(m))then; j=m
         else
            j=0
            do while(j<m-1); j=j+1; if(x(i)>=edges(j).and.x(i)<edges(j+1))exit; end do
            if(j==m-1.and..not.(x(i)>=edges(j).and.x(i)<edges(j+1)))j=0
         end if
         bins(i)=j; if(j>0)counts(j)=counts(j)+1
      end do
   end subroutine histc

   pure real(dp) function rmserr(x,y) result(v)
      real(dp),intent(in)::x(:),y(:)
      v=sqrt(sum((x-y)**2)/real(size(x),dp))
   end function rmserr

   pure elemental real(dp) function runge(x) result(y)
      real(dp),intent(in)::x
      y=1.0_dp/(1.0_dp+25.0_dp*x*x)
   end function runge

   pure elemental real(dp) function humps_function(x) result(y)
      real(dp),intent(in)::x
      y=1.0_dp/((x-0.3_dp)**2+0.01_dp)+1.0_dp/((x-0.9_dp)**2+0.04_dp)-6.0_dp
   end function humps_function

   function mexpfit(x,y,nterms,max_iter,tolerance) result(res)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in)::nterms
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tolerance
      type(regression_result)::res
      ! Deterministic log-linear initialization followed by linear amplitudes.
      real(dp),allocatable::a(:,:),rhs(:),coef(:),rates(:),fit(:)
      integer::j,n,istat
      n=size(x); allocate(rates(nterms),a(n,nterms),rhs(nterms),coef(nterms),fit(n))
      do j=1,nterms; rates(j)=real(j,dp)/max(maxval(x)-minval(x),1.0_dp); a(:,j)=exp(-rates(j)*(x-minval(x))); end do
      rhs=matmul(transpose(a),y); call solve_linear(matmul(transpose(a),a)+1e-12_dp*identity(nterms),rhs,coef,istat)
      fit=matmul(a,coef); allocate(res%coefficients(2*nterms),res%fitted(n),res%residuals(n))
      res%coefficients(:nterms)=coef; res%coefficients(nterms+1:)=rates; res%fitted=fit; res%residuals=y-fit
      res%rss=sum(res%residuals**2); res%iterations=1
      res%status=merge(pracma_ok,pracma_invalid_argument,istat==0)
      res%converged=istat==0
      if(present(max_iter))res%iterations=min(res%iterations,max_iter); if(present(tolerance))res%rss=res%rss+0.0_dp*tolerance
   end function mexpfit

   function odregress(x,y,sx,sy) result(res)
      real(dp),intent(in)::x(:),y(:)
      real(dp),intent(in),optional::sx(:),sy(:)
      type(regression_result)::res
      real(dp)::xm,ym,sxx,syy,sxy,b,lambda
      integer::n
      n=size(x); xm=sum(x)/n; ym=sum(y)/n; sxx=sum((x-xm)**2)/n; syy=sum((y-ym)**2)/n; sxy=sum((x-xm)*(y-ym))/n
      lambda=1.0_dp; if(present(sx).and.present(sy))lambda=max(sum(sy*sy)/max(sum(sx*sx),tiny(1.0_dp)),tiny(1.0_dp))
      b=(syy-lambda*sxx+sqrt((syy-lambda*sxx)**2+4*lambda*sxy*sxy))/(2*sxy)
      allocate(res%coefficients(2),res%fitted(n),res%residuals(n)); res%coefficients=[ym-b*xm,b]
      res%fitted=res%coefficients(1)+b*x; res%residuals=y-res%fitted
      res%rss=sum(res%residuals**2); res%converged=.true.; res%status=pracma_ok
   end function odregress

   function trigregress(x,y,degree,period) result(res)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in)::degree
      real(dp),intent(in),optional::period
      type(regression_result)::res
      real(dp)::p,w
      real(dp),allocatable::a(:,:),rhs(:),coef(:)
      integer::k,n,m,istat
      n=size(x); m=2*degree+1; p=maxval(x)-minval(x); if(present(period))p=period; w=2*pi_dp/p
      allocate(a(n,m),rhs(m),coef(m)); a(:,1)=1.0_dp
      do k=1,degree; a(:,2*k)=cos(real(k,dp)*w*x); a(:,2*k+1)=sin(real(k,dp)*w*x); end do
      rhs=matmul(transpose(a),y); call solve_linear(matmul(transpose(a),a),rhs,coef,istat)
      allocate(res%coefficients(m),res%fitted(n),res%residuals(n)); res%coefficients=coef; res%fitted=matmul(a,coef)
      res%residuals=y-res%fitted; res%rss=sum(res%residuals**2)
      res%converged=istat==0
      res%status=merge(pracma_ok,pracma_invalid_argument,res%converged)
   end function trigregress

   function whittaker(y,lambda,order) result(z)
      real(dp),intent(in)::y(:),lambda
      integer,intent(in),optional::order
      real(dp),allocatable::z(:),d(:,:),a(:,:),rhs(:)
      integer::n,k,i,istat
      n=size(y); k=2; if(present(order))k=order; allocate(d(max(0,n-k),n)); d=0.0_dp
      if(k==1)then; do i=1,n-1; d(i,i)=-1; d(i,i+1)=1; end do
      else; do i=1,n-2; d(i,i)=1; d(i,i+1)=-2; d(i,i+2)=1; end do; end if
      a=identity(n)+lambda*matmul(transpose(d),d); rhs=y; allocate(z(n)); call solve_linear(a,rhs,z,istat)
   end function whittaker

   function andrews_curve(x,t) result(y)
      real(dp),intent(in)::x(:),t(:)
      real(dp),allocatable::y(:)
      integer::k
      allocate(y(size(t))); y=x(1)/sqrt(2.0_dp)
      do k=2,size(x)
         if(mod(k,2)==0)then; y=y+x(k)*sin(real(k/2,dp)*t)
         else; y=y+x(k)*cos(real((k-1)/2,dp)*t); end if
      end do
   end function andrews_curve

   function normest(a,max_iter,tolerance) result(v)
      real(dp),intent(in)::a(:,:)
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tolerance
      real(dp)::v,vold,tol
      real(dp),allocatable::x(:),y(:)
      integer::i,n,it
      n=size(a,2); it=100; if(present(max_iter))it=max_iter; tol=1e-6_dp; if(present(tolerance))tol=tolerance
      allocate(x(n),y(size(a,1))); x=1.0_dp/sqrt(real(n,dp)); v=0.0_dp
      do i=1,it
         y=matmul(a,x); vold=v; v=sqrt(sum(y*y)); if(v<=tiny(1.0_dp))return
         x=matmul(transpose(a),y); x=x/sqrt(sum(x*x)); if(abs(v-vold)<=tol*max(1.0_dp,v))exit
      end do
   end function normest

   function autocorrelation(x,max_lag) result(r)
      real(dp),intent(in)::x(:); integer,intent(in)::max_lag
      real(dp),allocatable::r(:)
      real(dp)::mu,den
      integer::k,n
      n=size(x); allocate(r(0:max_lag)); mu=sum(x)/n; den=sum((x-mu)**2)
      do k=0,max_lag; r(k)=sum((x(1:n-k)-mu)*(x(1+k:n)-mu))/den; end do
   end function autocorrelation

   function crosscorrelation(x,y,max_lag) result(r)
      real(dp),intent(in)::x(:),y(:); integer,intent(in)::max_lag
      real(dp),allocatable::r(:)
      real(dp)::mx,my,den
      integer::k,n
      n=min(size(x),size(y)); allocate(r(-max_lag:max_lag)); mx=sum(x(:n))/n; my=sum(y(:n))/n
      den=sqrt(sum((x(:n)-mx)**2)*sum((y(:n)-my)**2))
      do k=-max_lag,max_lag
         if(k>=0)then; r(k)=sum((x(1:n-k)-mx)*(y(1+k:n)-my))/den
         else; r(k)=sum((x(1-k:n)-mx)*(y(1:n+k)-my))/den; end if
      end do
   end function crosscorrelation

   function periodogram(x) result(p)
      real(dp),intent(in)::x(:)
      real(dp),allocatable::p(:)
      complex(dp),allocatable::z(:),f(:)
      z=cmplx(x-sum(x)/real(size(x),dp),0.0_dp,dp); f=fft(z); allocate(p(size(x)/2+1)); p=abs(f(:size(p)))**2/real(size(x),dp)
   end function periodogram

   pure function identity(n) result(a)
      integer,intent(in)::n; real(dp)::a(n,n); integer::i
      a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
   end function identity

   real(dp) function median(x) result(v)
      real(dp),intent(in)::x(:); real(dp),allocatable::s(:); integer::n
      s=x; call sort_inplace(s); n=size(s); if(mod(n,2)==1)then; v=s((n+1)/2); else; v=0.5_dp*(s(n/2)+s(n/2+1)); end if
   end function median

   subroutine sort_inplace(a)
      real(dp),intent(inout)::a(:); integer::i,j; real(dp)::t
      do i=2,size(a); t=a(i); j=i-1; do while(j>=1); if(a(j)<=t)exit; a(j+1)=a(j); j=j-1; end do; a(j+1)=t; end do
   end subroutine sort_inplace

end module pracma_signal_stats
