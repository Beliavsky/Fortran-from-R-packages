module forecast_arima
   use forecast_kinds, only : dp, pi
   use forecast_types, only : arima_model, forecast_result
   use forecast_features, only : difference_series, seasonal_difference
   use forecast_stats, only : normal_quantile
   use forecast_optimize, only : pattern_search
   use forecast_unitroot, only : ndiffs, nsdiffs
   use forecast_linalg, only : least_squares
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private
   public :: arima_fit, arima_refit, arima_forecast, auto_arima, arima_errors, arima_order, arima_simulate
   public :: arima_innovations_loglik, arima_diffuse_loglik, arima_impulse_weights

   real(dp), allocatable, save :: ctx_z(:), ctx_y(:), ctx_xreg(:,:)
   integer, save :: cp=0,cq=0,csp=0,csq=0,cm=1,cdiff=0,csdiff=0,cxreg=0
   logical, save :: cmean=.false., cdiffuse=.false.
contains
   function polynomial_multiply(a,b) result(c)
      ! Coefficients are stored in ascending powers with the constant at index 1.
      real(dp),intent(in)::a(:),b(:)
      real(dp),allocatable::c(:)
      integer::i,j
      allocate(c(size(a)+size(b)-1))
      c=0.0_dp
      do i=1,size(a)
         do j=1,size(b)
            c(i+j-1)=c(i+j-1)+a(i)*b(j)
         end do
      end do
   end function polynomial_multiply

   subroutine effective_coefficients(ar,ma,sar,sma,m,phi,theta)
      real(dp),intent(in)::ar(:),ma(:),sar(:),sma(:)
      integer,intent(in)::m
      real(dp),allocatable,intent(out)::phi(:),theta(:)
      real(dp),allocatable::a(:),sa(:),b(:),sb(:),ap(:),bp(:)
      integer::i
      allocate(a(size(ar)+1),sa(m*size(sar)+1),b(size(ma)+1),sb(m*size(sma)+1))
      a=0.0_dp
      sa=0.0_dp
      b=0.0_dp
      sb=0.0_dp
      a(1)=1.0_dp
      sa(1)=1.0_dp
      b(1)=1.0_dp
      sb(1)=1.0_dp
      do i=1,size(ar)
         a(i+1)=-ar(i)
      end do
      do i=1,size(sar)
         sa(i*m+1)=-sar(i)
      end do
      do i=1,size(ma)
         b(i+1)=ma(i)
      end do
      do i=1,size(sma)
         sb(i*m+1)=sma(i)
      end do
      ap=polynomial_multiply(a,sa)
      bp=polynomial_multiply(b,sb)
      allocate(phi(max(0,size(ap)-1)),theta(max(0,size(bp)-1)))
      if(size(phi)>0)phi=-ap(2:)
      if(size(theta)>0)theta=bp(2:)
   end subroutine effective_coefficients

   subroutine unpack_params(x,ar,ma,sar,sma,mu,beta)
      real(dp),intent(in)::x(:)
      real(dp),allocatable,intent(out)::ar(:),ma(:),sar(:),sma(:),beta(:)
      real(dp),intent(out)::mu
      integer::k
      allocate(ar(cp),ma(cq),sar(csp),sma(csq),beta(cxreg))
      k=1
      if(cp>0)then
         ar=x(k:k+cp-1)
         k=k+cp
      end if
      if(cq>0)then
         ma=x(k:k+cq-1)
         k=k+cq
      end if
      if(csp>0)then
         sar=x(k:k+csp-1)
         k=k+csp
      end if
      if(csq>0)then
         sma=x(k:k+csq-1)
         k=k+csq
      end if
      mu=0.0_dp
      if(cmean)then
         mu=x(k)
         k=k+1
      end if
      if(cxreg>0)beta=x(k:k+cxreg-1)
   end subroutine unpack_params


   function objective_series(beta) result(z)
      real(dp),intent(in)::beta(:)
      real(dp),allocatable::z(:),raw(:)
      if(cxreg>0)then
         raw=ctx_y-matmul(ctx_xreg,beta)
         z=transformed_series(raw,cdiff,cm,csdiff)
      else
         z=ctx_z
      end if
   end function objective_series

   logical function impulse_stable(phi) result(ok)
      real(dp),intent(in)::phi(:)
      real(dp),allocatable::psi(:)
      integer::k,j,ncheck
      real(dp)::a
      ok=.true.
      if(size(phi)==0)return
      ncheck=max(256,32*size(phi))
      allocate(psi(0:ncheck))
      psi=0.0_dp
      psi(0)=1.0_dp
      do k=1,ncheck
         a=0.0_dp
         do j=1,min(size(phi),k)
            a=a+phi(j)*psi(k-j)
         end do
         psi(k)=a
         if(abs(a)>1.0e6_dp)then
            ok=.false.
            return
         end if
      end do
      if(sum(psi**2)>1.0e8_dp)ok=.false.
   end function impulse_stable

   logical function ma_invertible(theta) result(ok)
      real(dp),intent(in)::theta(:)
      real(dp),allocatable::a(:)
      integer::k,j,ncheck
      ok=.true.
      if(size(theta)==0)return
      ncheck=max(256,32*size(theta))
      allocate(a(0:ncheck))
      a=0.0_dp
      a(0)=1.0_dp
      do k=1,ncheck
         do j=1,min(size(theta),k)
            a(k)=a(k)-theta(j)*a(k-j)
         end do
         if(abs(a(k))>1.0e6_dp)then
            ok=.false.
            return
         end if
      end do
      if(sum(a**2)>1.0e8_dp)ok=.false.
   end function ma_invertible

   subroutine arma_psi(phi,theta,nmax,psi,info)
      real(dp),intent(in)::phi(:),theta(:)
      integer,intent(in)::nmax
      real(dp),allocatable,intent(out)::psi(:)
      integer,intent(out)::info
      integer::k,j
      allocate(psi(0:nmax))
      psi=0.0_dp
      psi(0)=1.0_dp
      info=0
      do k=1,nmax
         if(k<=size(theta))psi(k)=theta(k)
         do j=1,min(size(phi),k)
            psi(k)=psi(k)+phi(j)*psi(k-j)
         end do
         if(abs(psi(k))>1.0e100_dp)then
            info=1
            return
         end if
      end do
   end subroutine arma_psi

   subroutine arma_missing_loglik(z,phi,theta,mu,loglik,sigma2,residuals,vfactor,info)
      ! State-space Gaussian likelihood for stationary ARMA data containing NaNs.
      ! Missing observations skip the measurement update but retain time propagation.
      real(dp),intent(in)::z(:),phi(:),theta(:),mu
      real(dp),intent(out)::loglik,sigma2
      real(dp),allocatable,intent(out)::residuals(:),vfactor(:)
      integer,intent(out)::info
      real(dp),allocatable::tt(:,:),rr(:),pw(:,:),pnew(:,:),a(:),p0(:,:),kvec(:),prow(:)
      real(dp)::f,v,ss,ldiff,diff,scale
      integer::pord,qord,rwlen,nw,i,it,t,nuse

      info=0
      loglik=-huge(1.0_dp)
      sigma2=huge(1.0_dp)
      if(size(z)<1 .or. .not.impulse_stable(phi) .or. .not.ma_invertible(theta))then
         info=1
         return
      end if
      pord=size(phi)
      qord=size(theta)
      rwlen=max(1,pord)
      nw=rwlen+qord
      allocate(tt(nw,nw),rr(nw))
      tt=0.0_dp
      rr=0.0_dp
      if(pord>0)tt(1,1:pord)=phi
      do i=2,rwlen
         tt(i,i-1)=1.0_dp
      end do
      if(qord>0)then
         tt(1,rwlen+1:rwlen+qord)=theta
         do i=2,qord
            tt(rwlen+i,rwlen+i-1)=1.0_dp
         end do
         rr(rwlen+1)=1.0_dp
      end if
      rr(1)=1.0_dp
      allocate(pw(nw,nw))
      pw=0.0_dp
      do it=1,20000
         pnew=matmul(matmul(tt,pw),transpose(tt))+matmul(reshape(rr,[nw,1]),reshape(rr,[1,nw]))
         diff=maxval(abs(pnew-pw))
         scale=max(1.0_dp,maxval(abs(pnew)))
         pw=pnew
         if(diff<1.0e-13_dp*scale)exit
      end do
      if(it>20000)then
      info=2
      return
      end if
      allocate(a(nw),p0(nw,nw),residuals(size(z)),vfactor(size(z)))
      a=0.0_dp
      p0=pw
      residuals=0.0_dp
      vfactor=0.0_dp
      ss=0.0_dp
      ldiff=0.0_dp
      nuse=0
      do t=1,size(z)
         f=max(p0(1,1),tiny(1.0_dp))
         if(.not.ieee_is_nan(z(t)))then
            v=z(t)-mu-a(1)
            residuals(t)=v
            vfactor(t)=f
            kvec=p0(:,1)/f
            prow=p0(1,:)
            a=a+kvec*v
            p0=p0-matmul(reshape(kvec,[nw,1]),reshape(prow,[1,nw]))
            p0=0.5_dp*(p0+transpose(p0))
            ss=ss+v*v/f
            ldiff=ldiff+log(f)
            nuse=nuse+1
         end if
         if(t<size(z))then
            a=matmul(tt,a)
            p0=matmul(matmul(tt,p0),transpose(tt))+matmul(reshape(rr,[nw,1]),reshape(rr,[1,nw]))
            p0=0.5_dp*(p0+transpose(p0))
         end if
      end do
      if(nuse<1 .or. ss<=0.0_dp)then
      info=3
      return
      end if
      sigma2=max(ss/real(nuse,dp),tiny(1.0_dp))
      loglik=-0.5_dp*(real(nuse,dp)*(log(2.0_dp*pi*sigma2)+1.0_dp)+ldiff)
   end subroutine arma_missing_loglik

   subroutine arima_innovations_loglik(z,phi,theta,mu,loglik,sigma2,residuals,vfactor,info)
      ! Exact Gaussian likelihood for a stationary ARMA process using the
      ! innovations/Durbin-Levinson algorithm. sigma2 is profiled analytically.
      real(dp),intent(in)::z(:),phi(:),theta(:),mu
      real(dp),intent(out)::loglik,sigma2
      real(dp),allocatable,intent(out)::residuals(:),vfactor(:)
      integer,intent(out)::info
      real(dp),allocatable::psi(:),gamma(:),prev(:),curr(:)
      real(dp)::num,kappa,vprev,vcur,pred,q
      integer::n,j,k,jmax,istat

      n=size(z)
      info=0
      loglik=-huge(1.0_dp)
      sigma2=huge(1.0_dp)
      if(any(ieee_is_nan(z)))then
         call arma_missing_loglik(z,phi,theta,mu,loglik,sigma2,residuals,vfactor,info)
         return
      end if
      allocate(residuals(n),vfactor(n))
      residuals=0.0_dp
      vfactor=0.0_dp
      if(n<1 .or. .not.impulse_stable(phi) .or. .not.ma_invertible(theta))then
         info=1
         return
      end if
      jmax=max(n-1,max(512,min(32768,20*n+20*(size(phi)+size(theta))+256)))
      call arma_psi(phi,theta,jmax,psi,istat)
      if(istat/=0)then
         info=2
         return
      end if
      allocate(gamma(0:n-1))
      do k=0,n-1
         gamma(k)=dot_product(psi(0:jmax-k),psi(k:jmax))
      end do
      if(gamma(0)<=tiny(1.0_dp))then
         info=3
         return
      end if
      allocate(prev(max(1,n-1)),curr(max(1,n-1)))
      prev=0.0_dp
      curr=0.0_dp
      residuals(1)=z(1)-mu
      vfactor(1)=gamma(0)
      vprev=gamma(0)
      do k=1,n-1
         num=gamma(k)
         do j=1,k-1
            num=num-prev(j)*gamma(k-j)
         end do
         kappa=num/max(vprev,tiny(1.0_dp))
         if(abs(kappa)>=1.0_dp-1.0e-10_dp)then
            info=4
            return
         end if
         curr=0.0_dp
         do j=1,k-1
            curr(j)=prev(j)-kappa*prev(k-j)
         end do
         curr(k)=kappa
         vcur=vprev*(1.0_dp-kappa*kappa)
         if(vcur<=tiny(1.0_dp))then
            info=5
            return
         end if
         pred=mu
         do j=1,k
            pred=pred+curr(j)*(z(k+1-j)-mu)
         end do
         residuals(k+1)=z(k+1)-pred
         vfactor(k+1)=vcur
         prev(1:k)=curr(1:k)
         vprev=vcur
      end do
      q=sum(residuals**2/max(vfactor,tiny(1.0_dp)))
      sigma2=max(q/real(n,dp),tiny(1.0_dp))
      loglik=-0.5_dp*(real(n,dp)*(log(2.0_dp*pi*sigma2)+1.0_dp)+sum(log(vfactor)))
   end subroutine arima_innovations_loglik

   subroutine arima_diffuse_loglik(y,phi,theta,d,m,sd,drift,loglik,sigma2,residuals,vfactor,info)
      ! Gaussian likelihood for integrated ARIMA using an innovations state-space model
      ! with explicit diffuse states for the inverse differencing polynomial.
      real(dp),intent(in)::y(:),phi(:),theta(:),drift
      integer,intent(in)::d,m,sd
      real(dp),intent(out)::loglik,sigma2
      real(dp),allocatable,intent(out)::residuals(:),vfactor(:)
      integer,intent(out)::info
      real(dp),allocatable::delta(:),kernel(:),tmp(:),tw(:,:),rw(:),pw(:,:),pnew(:,:)
      real(dp),allocatable::tt(:,:),rr(:),a(:),p0(:,:),kvec(:),prow(:),yc(:)
      integer::rep,l,rwlen,pord,qord,nw,nstate,i,j,t,it,nuse
      real(dp)::f,v,ss,ldiff,kappa,scale,diff

      info=0
      loglik=-huge(1.0_dp)
      sigma2=huge(1.0_dp)
      if(size(y)<2)then
      info=1
      return
      end if
      if(.not.impulse_stable(phi).or..not.ma_invertible(theta))then
      info=2
      return
      end if
      allocate(delta(1))
      delta=1.0_dp
      do rep=1,d
         kernel=[1.0_dp,-1.0_dp]
         tmp=polynomial_multiply(delta,kernel)
         delta=tmp
      end do
      do rep=1,sd
         allocate(kernel(m+1))
         kernel=0.0_dp
         kernel(1)=1.0_dp
         kernel(m+1)=-1.0_dp
         tmp=polynomial_multiply(delta,kernel)
         delta=tmp
         deallocate(kernel)
      end do
      l=size(delta)-1
      if(l<=0)then
         call arima_innovations_loglik(y,phi,theta,0.0_dp,loglik,sigma2,residuals,vfactor,info)
         return
      end if
      if(size(y)<=l)then
      info=3
      return
      end if

      pord=size(phi)
      qord=size(theta)
      rwlen=max(1,pord)
      nw=rwlen+qord
      allocate(tw(nw,nw),rw(nw))
      tw=0.0_dp
      rw=0.0_dp
      if(pord>0)tw(1,1:pord)=phi
      do i=2,rwlen
         tw(i,i-1)=1.0_dp
      end do
      if(qord>0)then
         tw(1,rwlen+1:rwlen+qord)=theta
         do i=2,qord
            tw(rwlen+i,rwlen+i-1)=1.0_dp
         end do
         rw(rwlen+1)=1.0_dp
      end if
      rw(1)=1.0_dp
      allocate(pw(nw,nw))
      pw=0.0_dp
      do it=1,20000
         pnew=matmul(matmul(tw,pw),transpose(tw))+matmul(reshape(rw,[nw,1]),reshape(rw,[1,nw]))
         diff=maxval(abs(pnew-pw))
         scale=max(1.0_dp,maxval(abs(pnew)))
         pw=pnew
         if(diff<1.0e-13_dp*scale)exit
      end do
      if(it>20000)then
      info=4
      return
      end if

      nstate=l+nw
      allocate(tt(nstate,nstate),rr(nstate),a(nstate),p0(nstate,nstate))
      tt=0.0_dp
      rr=0.0_dp
      a=0.0_dp
      p0=0.0_dp
      ! y_{t+1} = w_{t+1} - sum(delta_k*y_{t+1-k}).
      tt(1,1:l)=-delta(2:)
      tt(1,l+1:l+nw)=tw(1,:)
      rr(1)=rw(1)
      do i=2,l
         tt(i,i-1)=1.0_dp
      end do
      tt(l+1:l+nw,l+1:l+nw)=tw
      rr(l+1:l+nw)=rw
      kappa=1.0e7_dp
      do i=1,l
         p0(i,i)=kappa
      end do
      p0(l+1:l+nw,l+1:l+nw)=pw
      allocate(yc(size(y)))
      yc=y
      if(abs(drift)>0.0_dp)then
         do i=1,size(y)
            yc(i)=yc(i)-drift*real(i,dp)
         end do
      end if
      allocate(residuals(size(y)),vfactor(size(y)))
      residuals=0.0_dp
      vfactor=0.0_dp
      ss=0.0_dp
      ldiff=0.0_dp
      nuse=0
      do t=1,size(y)
         f=max(p0(1,1),tiny(1.0_dp))
         if(.not.ieee_is_nan(yc(t)))then
            v=yc(t)-a(1)
            residuals(t)=v
            vfactor(t)=f
            kvec=p0(:,1)/f
            prow=p0(1,:)
            a=a+kvec*v
            p0=p0-matmul(reshape(kvec,[nstate,1]),reshape(prow,[1,nstate]))
            p0=0.5_dp*(p0+transpose(p0))
            if(t>l)then
               ss=ss+v*v/f
               ldiff=ldiff+log(f)
               nuse=nuse+1
            end if
         end if
         if(t<size(y))then
            a=matmul(tt,a)
            p0=matmul(matmul(tt,p0),transpose(tt))+matmul(reshape(rr,[nstate,1]),reshape(rr,[1,nstate]))
            p0=0.5_dp*(p0+transpose(p0))
         end if
      end do
      if(nuse<1.or.ss<=0.0_dp)then
      info=5
      return
      end if
      sigma2=max(ss/real(nuse,dp),tiny(1.0_dp))
      loglik=-0.5_dp*(real(nuse,dp)*(log(2.0_dp*pi*sigma2)+1.0_dp)+ldiff)
   end subroutine arima_diffuse_loglik

   function arima_css_objective(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v,mu
      real(dp),allocatable::ar(:),ma(:),sar(:),sma(:),beta(:),phi(:),theta(:),e(:),f(:),z(:)
      integer::t,j,start,n
      call unpack_params(x,ar,ma,sar,sma,mu,beta)
      call effective_coefficients(ar,ma,sar,sma,cm,phi,theta)
      if(.not.impulse_stable(phi) .or. .not.ma_invertible(theta))then
         v=1.0e40_dp
         return
      end if
      z=objective_series(beta)
      n=size(z)
      allocate(e(n),f(n))
      e=0.0_dp
      f=0.0_dp
      start=max(size(phi),size(theta))+1
      do t=1,n
         f(t)=mu
         do j=1,min(size(phi),t-1)
            f(t)=f(t)+phi(j)*(z(t-j)-mu)
         end do
         do j=1,min(size(theta),t-1)
            f(t)=f(t)+theta(j)*e(t-j)
         end do
         e(t)=z(t)-f(t)
      end do
      if(start>n)start=1
      v=sum(e(start:)**2)
   end function arima_css_objective

   function arima_ml_objective(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v,mu,ll,sig
      real(dp),allocatable::ar(:),ma(:),sar(:),sma(:),beta(:),phi(:),theta(:),e(:),vf(:),z(:)
      integer::info
      call unpack_params(x,ar,ma,sar,sma,mu,beta)
      call effective_coefficients(ar,ma,sar,sma,cm,phi,theta)
      z=objective_series(beta)
      if(cdiffuse .and. cdiff+csdiff>0)then
         if(cxreg>0)then
            z=ctx_y-matmul(ctx_xreg,beta)
         else
            z=ctx_y
         end if
         call arima_diffuse_loglik(z,phi,theta,cdiff,cm,csdiff,mu,ll,sig,e,vf,info)
      else
         call arima_innovations_loglik(z,phi,theta,mu,ll,sig,e,vf,info)
      end if
      if(info/=0)then
         v=1.0e40_dp
      else
         v=-ll
      end if
   end function arima_ml_objective

   function transformed_series(y,d,m,sd) result(z)
      real(dp),intent(in)::y(:)
      integer,intent(in)::d,m,sd
      real(dp),allocatable::z(:)
      integer::k
      z=y
      do k=1,sd
         z=seasonal_difference(z,m,1)
      end do
      do k=1,d
         z=difference_series(z,1)
      end do
   end function transformed_series

   function arima_fit(y,p,d,q,m,sp,sd,sq,include_mean,optimize,method,include_drift,xreg) result(model)
      real(dp),intent(in)::y(:)
      integer,intent(in)::p,d,q
      integer,intent(in),optional::m,sp,sd,sq
      logical,intent(in),optional::include_mean,optimize,include_drift
      character(len=*),intent(in),optional::method
      real(dp),intent(in),optional::xreg(:,:)
      type(arima_model)::model
      integer::mm,spp,sdd,sqq,npar,k,start,n,nobs,info,info_store,nx,rank
      logical::inc,opt,drift
      character(len=8)::meth
      real(dp),allocatable::x(:),lo(:),hi(:),ar(:),ma(:),sar(:),sma(:),beta(:),phi(:),theta(:)
      real(dp),allocatable::e(:),fit(:),vf(:),z(:),b0(:),olsres(:),zfinite(:)
      real(dp)::mu,rss,ll,sig,obj,sy,sx,span,ll_store,sig_store

      mm=1
      if(present(m))mm=max(1,m)
      spp=0
      if(present(sp))spp=sp
      sdd=0
      if(present(sd))sdd=sd
      sqq=0
      if(present(sq))sqq=sq
      inc=(d+sdd==0)
      if(present(include_mean))inc=include_mean
      drift=.false.
      if(present(include_drift))drift=include_drift
      if(drift .and. .not.(d==1 .and. sdd==0))drift=.false.
      if(d+sdd>0 .and. .not.drift)inc=.false.
      opt=.true.
      if(present(optimize))opt=optimize
      meth='ml'
      if(present(method))meth=adjustl(method)
      cdiffuse=(trim(meth)=='ml'.or.trim(meth)=='ML'.or.trim(meth)=='diffuse'.or.trim(meth)=='DIFFUSE')
      if((trim(meth)=='css'.or.trim(meth)=='CSS') .and. any(ieee_is_nan(y))) &
         error stop 'arima_fit: CSS does not support missing observations; use ML'

      nx=0
      if(present(xreg))then
         if(size(xreg,1)/=size(y))error stop 'arima_fit: xreg rows must equal length(y)'
         nx=size(xreg,2)
      end if
      ctx_y=y
      cdiff=d
      csdiff=sdd
      if(allocated(ctx_xreg))deallocate(ctx_xreg)
      if(present(xreg))then
         ctx_xreg=xreg
      else
         allocate(ctx_xreg(size(y),0))
      end if
      ctx_z=transformed_series(y,d,mm,sdd)
      cp=p
      cq=q
      csp=spp
      csq=sqq
      cm=mm
      cmean=inc.or.drift
      cxreg=nx
      npar=p+q+spp+sqq+merge(1,0,cmean)+nx
      allocate(x(npar),lo(npar),hi(npar))
      x=0.0_dp
      lo=-0.98_dp
      hi=0.98_dp
      k=p+q+spp+sqq+1
      if(cmean)then
         zfinite=pack(ctx_z,.not.ieee_is_nan(ctx_z))
         if(size(zfinite)==0)error stop 'arima_fit: all transformed observations are missing'
         x(k)=sum(zfinite)/real(size(zfinite),dp)
         lo(k)=minval(zfinite)-2.0_dp*max(1.0_dp,abs(x(k)))
         hi(k)=maxval(zfinite)+2.0_dp*max(1.0_dp,abs(x(k)))
         k=k+1
      end if
      if(nx>0)then
         call least_squares(xreg,y,b0,olsres,rank,info)
         if(info/=0)error stop 'arima_fit: xreg least-squares initialization failed'
         sy=sqrt(max(sum((y-sum(y)/real(size(y),dp))**2)/real(max(1,size(y)-1),dp),tiny(1.0_dp)))
         do start=1,nx
            sx=sqrt(max(sum((xreg(:,start)-sum(xreg(:,start))/real(size(y),dp))**2)/ &
               real(max(1,size(y)-1),dp),tiny(1.0_dp)))
            if(sx<=sqrt(tiny(1.0_dp)))error stop 'arima_fit: constant xreg column'
            x(k)=b0(start)
            span=max(2.0_dp,abs(b0(start))+10.0_dp*sy/sx)
            lo(k)=b0(start)-span
            hi(k)=b0(start)+span
            k=k+1
         end do
      end if

      if(opt .and. npar>0)then
         ! CSS gives a good inexpensive starting point when the series is complete.
         if(.not.any(ieee_is_nan(ctx_y)))then
            call pattern_search(arima_css_objective,x,lo,hi,maxit=300,tol=5.0e-5_dp)
         end if
         if(trim(meth)/='css' .and. trim(meth)/='CSS')then
            call pattern_search(arima_ml_objective,x,lo,hi,maxit=380,tol=2.0e-5_dp,fval=obj)
         end if
      end if
      call unpack_params(x,ar,ma,sar,sma,mu,beta)
      call effective_coefficients(ar,ma,sar,sma,mm,phi,theta)
      if(nx>0)then
         z=transformed_series(y-matmul(xreg,beta),d,mm,sdd)
      else
         z=ctx_z
      end if
      n=size(z)
      if(trim(meth)=='css' .or. trim(meth)=='CSS')then
         allocate(e(n),fit(n))
         e=0.0_dp
         fit=0.0_dp
         do k=1,n
            fit(k)=mu
            if(k>1)then
               do start=1,min(size(phi),k-1)
                  fit(k)=fit(k)+phi(start)*(z(k-start)-mu)
               end do
               do start=1,min(size(theta),k-1)
                  fit(k)=fit(k)+theta(start)*e(k-start)
               end do
            end if
            e(k)=z(k)-fit(k)
         end do
         start=max(size(phi),size(theta))+1
         if(start>n)start=1
         rss=sum(e(start:)**2)
         sig=rss/real(max(1,n-start+1),dp)
         ll=-0.5_dp*real(n-start+1,dp)*(log(2.0_dp*pi*max(sig,tiny(1.0_dp)))+1.0_dp)
      else
         if(cdiffuse .and. d+sdd>0)then
            if(nx>0)then
               call arima_diffuse_loglik(y-matmul(xreg,beta),phi,theta,d,mm,sdd,mu,ll,sig,e,vf,info)
            else
               call arima_diffuse_loglik(y,phi,theta,d,mm,sdd,mu,ll,sig,e,vf,info)
            end if
            ! Store ARMA innovations on the differenced scale for forecasting/diagnostics.
            if(info==0)then
               call arima_innovations_loglik(z,phi,theta,mu,ll_store,sig_store,e,vf,info_store)
               if(info_store/=0)info=info_store
            end if
         else
            call arima_innovations_loglik(z,phi,theta,mu,ll,sig,e,vf,info)
         end if
         if(info/=0)then
            ! A pathological exact-likelihood endpoint falls back to CSS rather than returning an invalid model.
            allocate(fit(n))
            fit=0.0_dp
            e=0.0_dp
            do k=1,n
               fit(k)=mu
               do start=1,min(size(phi),k-1)
                  fit(k)=fit(k)+phi(start)*(z(k-start)-mu)
               end do
               do start=1,min(size(theta),k-1)
                  fit(k)=fit(k)+theta(start)*e(k-start)
               end do
               e(k)=z(k)-fit(k)
            end do
            sig=sum(e**2)/real(n,dp)
            ll=-0.5_dp*real(n,dp)*(log(2.0_dp*pi*max(sig,tiny(1.0_dp)))+1.0_dp)
         else
            allocate(fit(n))
            fit=z-e
         end if
      end if

      model%sigma2=sig
      model%loglik=ll
      model%p=p
      model%d=d
      model%q=q
      model%sp=spp
      model%sd=sdd
      model%sq=sqq
      model%m=mm
      model%include_mean=inc
      model%include_drift=drift
      model%intercept=merge(mu,0.0_dp,inc)
      model%drift=merge(mu,0.0_dp,drift)
      model%ar=ar
      model%ma=ma
      model%sar=sar
      model%sma=sma
      model%xreg_coef=beta
      if(present(xreg))then
         model%xreg=xreg
      else
         allocate(model%xreg(size(y),0))
      end if
      model%residuals=e
      model%fitted=fit
      npar=npar+1
      nobs=count(.not.ieee_is_nan(z))
      model%aic=-2.0_dp*model%loglik+2.0_dp*real(npar,dp)
      if(nobs>0)then
         model%bic=-2.0_dp*model%loglik+log(real(nobs,dp))*real(npar,dp)
      else
         model%bic=huge(1.0_dp)
      end if
      if(nobs>npar+1)then
         model%aicc=model%aic+2.0_dp*real(npar*(npar+1),dp)/real(nobs-npar-1,dp)
      else
         model%aicc=huge(1.0_dp)
      end if
   end function arima_fit

   function arima_refit(y,model,xreg) result(refit)
      ! Re-evaluate a fitted ARIMA structure on new data while holding all
      ! coefficients fixed, analogous to forecast::Arima(..., model=old).
      real(dp),intent(in)::y(:)
      type(arima_model),intent(in)::model
      real(dp),intent(in),optional::xreg(:,:)
      type(arima_model)::refit
      real(dp),allocatable::raw(:),z(:),phi(:),theta(:),e(:),vf(:)
      real(dp)::ll,sig,ll_store,sig_store,mu
      integer::info,info_store,npar,nobs

      refit=model
      raw=y
      if(allocated(model%xreg_coef) .and. size(model%xreg_coef)>0)then
         if(.not.present(xreg))error stop 'arima_refit: regressors required by supplied model'
         if(size(xreg,1)/=size(y) .or. size(xreg,2)/=size(model%xreg_coef)) &
            error stop 'arima_refit: regressor dimensions do not match model'
         raw=y-matmul(xreg,model%xreg_coef)
         refit%xreg=xreg
      else
         if(present(xreg))then
            if(size(xreg,1)/=size(y) .or. size(xreg,2)/=0) &
               error stop 'arima_refit: supplied model has no regressors'
         end if
         allocate(refit%xreg(size(y),0))
      end if
      call effective_coefficients(model%ar,model%ma,model%sar,model%sma,model%m,phi,theta)
      z=transformed_series(raw,model%d,model%m,model%sd)
      mu=0.0_dp
      if(model%include_mean)mu=model%intercept
      if(model%include_drift)mu=model%drift
      if(model%d+model%sd>0)then
         call arima_diffuse_loglik(raw,phi,theta,model%d,model%m,model%sd,mu,ll,sig,e,vf,info)
         if(info/=0)error stop 'arima_refit: diffuse likelihood failed'
         call arima_innovations_loglik(z,phi,theta,mu,ll_store,sig_store,e,vf,info_store)
         if(info_store/=0)error stop 'arima_refit: innovations reconstruction failed'
      else
         call arima_innovations_loglik(z,phi,theta,mu,ll,sig,e,vf,info)
         if(info/=0)error stop 'arima_refit: likelihood failed'
      end if
      refit%loglik=ll
      refit%sigma2=sig
      refit%residuals=e
      refit%fitted=z-e
      nobs=count(.not.ieee_is_nan(z))
      npar=model%p+model%q+model%sp+model%sq+merge(1,0,model%include_mean.or.model%include_drift)+1
      if(allocated(model%xreg_coef))npar=npar+size(model%xreg_coef)
      refit%aic=-2.0_dp*ll+2.0_dp*real(npar,dp)
      if(nobs>0)then
         refit%bic=-2.0_dp*ll+log(real(nobs,dp))*real(npar,dp)
      else
         refit%bic=huge(1.0_dp)
      end if
      if(nobs>npar+1)then
         refit%aicc=refit%aic+2.0_dp*real(npar*(npar+1),dp)/real(nobs-npar-1,dp)
      else
         refit%aicc=huge(1.0_dp)
      end if
   end function arima_refit

   function auto_arima(y,m,max_p,max_q,max_sp,max_sq,seasonal,stepwise,approximation,max_order, &
      allowmean,allowdrift,seasonal_test,unitroot_test,unitroot_type,alpha,xreg,ic,start_p,start_q,start_sp,start_sq, &
      nmodels,d_fixed,sd_fixed,truncate) result(best)
      real(dp),intent(in)::y(:)
      integer,intent(in),optional::m,max_p,max_q,max_sp,max_sq,max_order,start_p,start_q,start_sp,start_sq,nmodels
      integer,intent(in),optional::d_fixed,sd_fixed,truncate
      logical,intent(in),optional::seasonal,stepwise,approximation,allowmean,allowdrift
      character(len=*),intent(in),optional::seasonal_test,unitroot_test,unitroot_type,ic
      real(dp),intent(in),optional::alpha,xreg(:,:)
      type(arima_model)::best,cand
      integer::mm,mp,mq,msp,msq,mo,d,sd,p,q,sp,sq,ps,qs,Ps0,Qs0,nmod,nfit
      integer::curp,curq,cursp,cursq,dpn,dqn,dspn,dsqn,i,j,k,rank,info,nvisited,nsel
      integer,allocatable::visited(:,:)
      logical::seas,sw,approx,amean,adrift,inc,dr,curconst,tryconst,improved
      character(len=8)::smethod,fitmethod,urmethod,urtype,icname
      real(dp)::ualpha,bestscore,score
      real(dp),allocatable::xx(:),dx(:),breg(:),rreg(:),ysel(:),xsel(:,:)

      mm=1
      if(present(m))mm=max(1,m)
      mp=5
      if(present(max_p))mp=max(0,max_p)
      mq=5
      if(present(max_q))mq=max(0,max_q)
      msp=2
      if(present(max_sp))msp=max(0,max_sp)
      msq=2
      if(present(max_sq))msq=max(0,max_sq)
      mo=5
      if(present(max_order))mo=max(0,max_order)
      seas=mm>1
      if(present(seasonal))seas=seasonal
      if(.not.seas)then
         msp=0
         msq=0
      end if
      sw=.true.
      if(present(stepwise))sw=stepwise
      approx=.false.
      if(present(approximation))approx=approximation
      amean=.true.
      if(present(allowmean))amean=allowmean
      adrift=.true.
      if(present(allowdrift))adrift=allowdrift
      smethod='seas'
      if(present(seasonal_test))smethod=adjustl(seasonal_test)
      urmethod='kpss'
      if(present(unitroot_test))urmethod=adjustl(unitroot_test)
      urtype='level'
      if(present(unitroot_type))urtype=adjustl(unitroot_type)
      icname='aicc'
      if(present(ic))icname=adjustl(ic)
      ualpha=0.05_dp
      if(present(alpha))ualpha=alpha
      nmod=94
      if(present(nmodels))nmod=max(1,nmodels)
      if(present(xreg))then
         if(size(xreg,1)/=size(y))error stop 'auto_arima: xreg rows must equal length(y)'
         if(size(xreg,2)>0)then
            call least_squares(xreg,y,breg,rreg,rank,info)
            if(info/=0)error stop 'auto_arima: xreg regression failed'
            xx=rreg
         else
            xx=y
         end if
      else
         xx=y
      end if

      if(present(sd_fixed))then
         sd=max(0,sd_fixed)
      else
         sd=merge(nsdiffs(xx,mm,test=smethod),0,seas)
      end if
      dx=xx
      if(sd>0)dx=transformed_series(xx,0,mm,sd)
      if(present(d_fixed))then
         d=max(0,d_fixed)
      else
         d=ndiffs(dx,test=urmethod,type=urtype,alpha=ualpha)
      end if
      fitmethod=merge('css     ','ml      ',approx)
      ysel=y
      if(present(xreg))then
         xsel=xreg
      else
         allocate(xsel(size(y),0))
      end if
      if(approx .and. present(truncate))then
         nsel=max(1,truncate)
         if(size(y)>nsel)then
            ysel=y(size(y)-nsel+1:)
            if(present(xreg))xsel=xreg(size(y)-nsel+1:,:)
         end if
      end if
      best%aicc=huge(1.0_dp)
      best%aic=huge(1.0_dp)
      best%bic=huge(1.0_dp)
      bestscore=huge(1.0_dp)
      nfit=0

      if(.not.sw)then
         do p=0,mp
            do q=0,mq
               do sp=0,msp
                  do sq=0,msq
                     if(p+q+sp+sq>mo)cycle
                     inc=(d+sd==0 .and. amean)
                     dr=(d==1 .and. sd==0 .and. adrift)
                     call consider_model(p,q,sp,sq,inc.or.dr)
                     if(inc .or. dr)call consider_model(p,q,sp,sq,.false.)
                  end do
               end do
            end do
         end do
      else
         ps=min(2,mp)
         if(present(start_p))ps=min(max(0,start_p),mp)
         qs=min(2,mq)
         if(present(start_q))qs=min(max(0,start_q),mq)
         Ps0=min(1,msp)
         if(present(start_sp))Ps0=min(max(0,start_sp),msp)
         Qs0=min(1,msq)
         if(present(start_sq))Qs0=min(max(0,start_sq),msq)
         curp=ps
         curq=qs
         cursp=Ps0
         cursq=Qs0
         curconst=(d+sd==0 .and. amean) .or. (d==1 .and. sd==0 .and. adrift)
         allocate(visited(5,nmod))
         visited=-999999
         nvisited=0

         call consider_model(curp,curq,cursp,cursq,curconst)
         curp=best%p
         curq=best%q
         cursp=best%sp
         cursq=best%sq
         curconst=best%include_mean.or.best%include_drift
         call consider_model(0,0,0,0,curconst)
         call consider_model(0,0,0,0,.false.)
         if(mp>0 .or. msp>0)call consider_model(merge(1,0,mp>0),0,merge(1,0,msp>0),0,curconst)
         if(mq>0 .or. msq>0)call consider_model(0,merge(1,0,mq>0),0,merge(1,0,msq>0),curconst)
         curp=best%p
         curq=best%q
         cursp=best%sp
         cursq=best%sq
         curconst=best%include_mean.or.best%include_drift

         improved=.true.
         do while(improved .and. nfit<nmod)
            improved=.false.
            ! Match forecast::auto.arima stepwise candidate order. An
            ! improvement restarts the sequence from the new incumbent.
            call try_step(curp,curq,cursp-1,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq,cursp,cursq-1,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq,cursp+1,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq,cursp,cursq+1,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq,cursp-1,cursq-1,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq,cursp-1,cursq+1,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq,cursp+1,cursq-1,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq,cursp+1,cursq+1,curconst,improved)
            if(improved)cycle
            call try_step(curp-1,curq,cursp,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq-1,cursp,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp+1,curq,cursp,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp,curq+1,cursp,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp-1,curq-1,cursp,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp-1,curq+1,cursp,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp+1,curq-1,cursp,cursq,curconst,improved)
            if(improved)cycle
            call try_step(curp+1,curq+1,cursp,cursq,curconst,improved)
            if(improved)cycle
            if((d+sd==0 .and. amean).or.(d==1.and.sd==0.and.adrift))then
               call try_step(curp,curq,cursp,cursq,.not.curconst,improved)
               if(improved)cycle
            end if
         end do
      end if

      ! Upstream approximation mode uses CSS for selection but refits the selected order by ML.
      if(approx)then
         if(present(xreg))then
            best=arima_fit(y,best%p,best%d,best%q,best%m,best%sp,best%sd,best%sq, &
               best%include_mean,.true.,'ml',best%include_drift,xreg)
         else
            best=arima_fit(y,best%p,best%d,best%q,best%m,best%sp,best%sd,best%sq, &
               best%include_mean,.true.,'ml',best%include_drift)
         end if
      end if
   contains
      real(dp) function model_score(a) result(v)
         type(arima_model),intent(in)::a
         select case(trim(icname))
         case('aic','AIC'); v=a%aic
         case('bic','BIC'); v=a%bic
         case default; v=a%aicc
         end select
      end function model_score

      logical function already_seen(pp,qq,PPs,QQs,cflag) result(seen)
         integer,intent(in)::pp,qq,PPs,QQs
         logical,intent(in)::cflag
         integer::ii,cf
         cf=merge(1,0,cflag)
         seen=.false.
         do ii=1,nvisited
            if(all(visited(:,ii)==[pp,qq,PPs,QQs,cf]))then
               seen=.true.
               return
            end if
         end do
         if(nvisited<nmod)then
            nvisited=nvisited+1
            visited(:,nvisited)=[pp,qq,PPs,QQs,cf]
         end if
      end function already_seen

      subroutine try_step(pp,qq,PPs,QQs,cflag,changed)
         integer,intent(in)::pp,qq,PPs,QQs
         logical,intent(in)::cflag
         logical,intent(out)::changed
         real(dp)::oldscore
         changed=.false.
         if(nfit>=nmod)return
         if(pp<0.or.pp>mp.or.qq<0.or.qq>mq.or.PPs<0.or.PPs>msp.or.QQs<0.or.QQs>msq)return
         oldscore=bestscore
         call consider_model(pp,qq,PPs,QQs,cflag)
         if(bestscore<oldscore-1.0e-10_dp)then
            curp=best%p
            curq=best%q
            cursp=best%sp
            cursq=best%sq
            curconst=best%include_mean.or.best%include_drift
            changed=.true.
         end if
      end subroutine try_step

      subroutine consider_model(pp,qq,PPs,QQs,cflag)
         integer,intent(in)::pp,qq,PPs,QQs
         logical,intent(in)::cflag
         logical::imean,idrift
         if(pp<0.or.pp>mp.or.qq<0.or.qq>mq.or.PPs<0.or.PPs>msp.or.QQs<0.or.QQs>msq)return
         if(sw)then
            if(nfit>=nmod)return
            if(already_seen(pp,qq,PPs,QQs,cflag))return
         else if(pp+qq+PPs+QQs>mo)then
            return
         end if
         imean=cflag.and.(d+sd==0).and.amean
         idrift=cflag.and.(d==1).and.(sd==0).and.adrift
         if(present(xreg))then
            cand=arima_fit(ysel,pp,d,qq,mm,PPs,sd,QQs,imean,.true.,fitmethod,idrift,xsel)
         else
            cand=arima_fit(ysel,pp,d,qq,mm,PPs,sd,QQs,imean,.true.,fitmethod,idrift)
         end if
         nfit=nfit+1
         score=model_score(cand)
         if(score<bestscore)then
            best=cand
            bestscore=score
         end if
      end subroutine consider_model
   end function auto_arima

   function arima_errors(model) result(e)
      type(arima_model),intent(in)::model
      real(dp),allocatable::e(:)
      e=model%residuals
   end function arima_errors

   function arima_order(model) result(ord)
      type(arima_model),intent(in)::model
      integer::ord(7)
      ord=[model%p,model%d,model%q,model%sp,model%sd,model%sq,model%m]
   end function arima_order

   subroutine arima_impulse_weights(phi,theta,d,m,sd,h,weights)
      real(dp),intent(in)::phi(:),theta(:)
      integer,intent(in)::d,m,sd,h
      real(dp),allocatable,intent(out)::weights(:)
      real(dp),allocatable::psi(:),delta(:),tmp(:),kernel(:)
      integer::info,k,j,rep
      if(h<=0)then
         allocate(weights(0))
         return
      end if
      call arma_psi(phi,theta,h-1,psi,info)
      allocate(delta(0:h-1))
      delta=0.0_dp
      delta(0)=1.0_dp
      do rep=1,d
         allocate(kernel(0:h-1))
         kernel=1.0_dp
         allocate(tmp(0:h-1))
         tmp=0.0_dp
         do k=0,h-1
            do j=0,k
               tmp(k)=tmp(k)+delta(j)*kernel(k-j)
            end do
         end do
         call move_alloc(tmp,delta)
         deallocate(kernel)
      end do
      do rep=1,sd
         allocate(kernel(0:h-1))
         kernel=0.0_dp
         do k=0,h-1,m
            kernel(k)=1.0_dp
         end do
         allocate(tmp(0:h-1))
         tmp=0.0_dp
         do k=0,h-1
            do j=0,k
               tmp(k)=tmp(k)+delta(j)*kernel(k-j)
            end do
         end do
         call move_alloc(tmp,delta)
         deallocate(kernel)
      end do
      allocate(weights(0:h-1))
      weights=0.0_dp
      do k=0,h-1
         do j=0,k
            weights(k)=weights(k)+psi(j)*delta(k-j)
         end do
      end do
   end subroutine arima_impulse_weights

   function arima_forecast(model,y,h,levels,future_xreg) result(fc)
      type(arima_model),intent(in)::model
      real(dp),intent(in)::y(:)
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:),future_xreg(:,:)
      type(forecast_result)::fc
      real(dp),allocatable::z(:),work(:),e(:),phi(:),theta(:),lev(:),base(:),seasbase(:),weights(:),yadj(:)
      real(dp)::f,zq,mu
      integer::i,j,n,nz,k,nx
      nx=0
      if(allocated(model%xreg_coef))nx=size(model%xreg_coef)
      if(nx>0)then
         if(.not.allocated(model%xreg))error stop 'arima_forecast: fitted model lacks xreg history'
         if(size(model%xreg,1)/=size(y).or.size(model%xreg,2)/=nx)error stop 'arima_forecast: xreg history mismatch'
         if(.not.present(future_xreg))error stop 'arima_forecast: future_xreg is required for an ARIMAX model'
         if(size(future_xreg,1)/=h.or.size(future_xreg,2)/=nx)error stop 'arima_forecast: future_xreg shape mismatch'
         yadj=y-matmul(model%xreg,model%xreg_coef)
      else
         yadj=y
      end if
      call effective_coefficients(model%ar,model%ma,model%sar,model%sma,model%m,phi,theta)
      z=transformed_series(yadj,model%d,model%m,model%sd)
      nz=size(z)
      allocate(work(nz+h),e(nz+h))
      work(1:nz)=z
      e=0.0_dp
      e(1:min(nz,size(model%residuals)))=model%residuals(1:min(nz,size(model%residuals)))
      mu=0.0_dp
      if(model%include_mean)mu=model%intercept
      if(model%include_drift)mu=model%drift
      do i=1,h
         n=nz+i
         f=mu
         do j=1,min(size(phi),n-1)
            f=f+phi(j)*(work(n-j)-mu)
         end do
         do j=1,min(size(theta),n-1)
            f=f+theta(j)*e(n-j)
         end do
         work(n)=f
         e(n)=0.0_dp
      end do
      base=work(nz+1:)
      if(model%sd>0)then
         seasbase=seasonal_difference(yadj,model%m,model%sd)
      else
         seasbase=yadj
      end if
      do k=1,model%d
         call integrate_regular(base,seasbase)
         if(k<model%d)seasbase=[seasbase,base]
      end do
      do k=1,model%sd
         call integrate_seasonal(base,yadj,model%m)
      end do
      if(nx>0)base=base+matmul(future_xreg,model%xreg_coef)
      if(present(levels))then
         lev=levels
      else
         lev=[80.0_dp,95.0_dp]
      end if
      allocate(fc%mean(h),fc%se(h),fc%level(size(lev)),fc%lower(h,size(lev)),fc%upper(h,size(lev)))
      fc%mean=base
      fc%level=lev
      call arima_impulse_weights(phi,theta,model%d,model%m,model%sd,h,weights)
      do i=1,h
         fc%se(i)=sqrt(max(model%sigma2*sum(weights(0:i-1)**2),0.0_dp))
      end do
      do j=1,size(lev)
         zq=normal_quantile(0.5_dp+lev(j)/200.0_dp)
         fc%lower(:,j)=fc%mean-zq*fc%se
         fc%upper(:,j)=fc%mean+zq*fc%se
      end do
   contains
      subroutine integrate_regular(pred,hist)
         real(dp),intent(inout)::pred(:)
         real(dp),intent(in)::hist(:)
         integer::ii
         real(dp)::last
         last=hist(size(hist))
         do ii=1,size(pred)
            pred(ii)=last+pred(ii)
            last=pred(ii)
         end do
      end subroutine integrate_regular
      subroutine integrate_seasonal(pred,hist,m)
         real(dp),intent(inout)::pred(:)
         real(dp),intent(in)::hist(:)
         integer,intent(in)::m
         real(dp),allocatable::all(:)
         integer::ii,n0
         n0=size(hist)
         allocate(all(n0+size(pred)))
         all(1:n0)=hist
         do ii=1,size(pred)
            all(n0+ii)=pred(ii)+all(n0+ii-m)
            pred(ii)=all(n0+ii)
         end do
      end subroutine integrate_seasonal
   end function arima_forecast

   function arima_simulate(model,n,burnin,initial,future_xreg) result(y)
      type(arima_model),intent(in)::model
      integer,intent(in)::n
      integer,intent(in),optional::burnin
      real(dp),intent(in),optional::initial(:),future_xreg(:,:)
      real(dp),allocatable::y(:),z(:),e(:),phi(:),theta(:),u(:),delta(:),hist(:),all(:),tmp(:),kernel(:)
      integer::b,t,j,tot,rep,k,degree,nx
      real(dp)::r1,r2,mu
      b=200
      if(present(burnin))b=max(0,burnin)
      tot=n+b
      allocate(z(tot),e(tot),u(tot))
      z=0.0_dp
      e=0.0_dp
      call effective_coefficients(model%ar,model%ma,model%sar,model%sma,model%m,phi,theta)
      mu=0.0_dp
      if(model%include_mean)mu=model%intercept
      if(model%include_drift)mu=model%drift
      t=1
      do while(t<=tot)
         call random_number(r1)
         call random_number(r2)
         r1=max(r1,1.0e-12_dp)
         u(t)=sqrt(-2.0_dp*log(r1))*cos(2*pi*r2)*sqrt(model%sigma2)
         t=t+1
      end do
      do t=1,tot
         z(t)=mu+u(t)
         do j=1,min(size(phi),t-1)
            z(t)=z(t)+phi(j)*(z(t-j)-mu)
         end do
         do j=1,min(size(theta),t-1)
            z(t)=z(t)+theta(j)*e(t-j)
         end do
         e(t)=u(t)
      end do
      allocate(y(n))
      y=z(b+1:)

      if(model%d+model%sd>0)then
         allocate(delta(1))
         delta=1.0_dp
         do rep=1,model%d
            allocate(kernel(2))
            kernel=[1.0_dp,-1.0_dp]
            tmp=polynomial_multiply(delta,kernel)
            call move_alloc(tmp,delta)
            deallocate(kernel)
         end do
         do rep=1,model%sd
            allocate(kernel(model%m+1))
            kernel=0.0_dp
            kernel(1)=1.0_dp
            kernel(model%m+1)=-1.0_dp
            tmp=polynomial_multiply(delta,kernel)
            call move_alloc(tmp,delta)
            deallocate(kernel)
         end do
         degree=size(delta)-1
         if(.not.present(initial))error stop 'arima_simulate: initial history required for integrated model'
         if(size(initial)<degree)error stop 'arima_simulate: initial history is too short for differencing order'
         hist=initial(size(initial)-degree+1:)
         nx=0
         if(allocated(model%xreg_coef))nx=size(model%xreg_coef)
         if(nx>0 .and. allocated(model%xreg))then
            if(size(model%xreg,1)>=degree)then
               hist=hist-matmul(model%xreg(size(model%xreg,1)-degree+1:,:),model%xreg_coef)
            end if
         end if
         allocate(all(degree+n))
         all(1:degree)=hist
         do t=1,n
            all(degree+t)=y(t)
            do k=1,degree
               all(degree+t)=all(degree+t)-delta(k+1)*all(degree+t-k)
            end do
         end do
         y=all(degree+1:)
      end if
      nx=0
      if(allocated(model%xreg_coef))nx=size(model%xreg_coef)
      if(nx>0)then
         if(.not.present(future_xreg))error stop 'arima_simulate: future_xreg required for ARIMAX model'
         if(size(future_xreg,1)/=n.or.size(future_xreg,2)/=nx)error stop 'arima_simulate: future_xreg shape mismatch'
         y=y+matmul(future_xreg,model%xreg_coef)
      end if
   end function arima_simulate
end module forecast_arima
