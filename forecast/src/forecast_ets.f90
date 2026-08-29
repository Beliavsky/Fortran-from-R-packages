module forecast_ets
   use forecast_kinds, only : dp
   use forecast_types, only : ets_model, forecast_result
   use forecast_stats, only : mean_value, variance_value, normal_quantile, quantile_type8
   use forecast_optimize, only : pattern_search
   use forecast_linalg, only : polynomial_max_root_mod
   implicit none
   private
   integer, parameter, public :: ETS_NONE=0, ETS_ADD=1, ETS_MULT=2
   public :: ets_calc, ets_fit, ets_auto, ets_forecast, ets_forecast_simulated, ets_simulate
   public :: ses_fit, holt_fit, hw_fit
   real(dp), allocatable, save :: ctx_y(:),ctx_init(:)
   integer, save :: ctx_m=1,ctx_err=1,ctx_trend=0,ctx_season=0
   logical, save :: ctx_damped=.false.,ctx_opt_init=.true.
   integer, save :: ctx_nsmooth=0
contains
   subroutine ets_state_forecast(l,b,s,m,trend,season,phi,f)
      real(dp),intent(in)::l,b,s(:),phi
      integer,intent(in)::m,trend,season
      real(dp),intent(out)::f(:)
      real(dp)::phistar
      integer::i,j
      phistar=phi
      do i=1,size(f)
         select case(trend)
         case(ETS_NONE)
         f(i)=l
         case(ETS_ADD)
         f(i)=l+phistar*b
         case(ETS_MULT)
         if(b<0.0_dp)then
         f(i)=huge(1.0_dp)
         else
         f(i)=l*b**phistar
         end if
         end select
         if(season/=ETS_NONE)then
            j=m-mod(i-1,m)
            if(season==ETS_ADD)f(i)=f(i)+s(j)
            if(season==ETS_MULT)f(i)=f(i)*s(j)
         end if
         if(i<size(f))then
            if(abs(phi-1.0_dp)<1.0e-10_dp)then
            phistar=phistar+1.0_dp
            else
            phistar=phistar+phi**real(i,dp)
            end if
         end if
      end do
   end subroutine

   subroutine ets_update(oldl,l,oldb,b,olds,s,m,trend,season,alpha,beta,gamma,phi,y)
      real(dp),intent(in)::oldl,oldb,olds(:),alpha,beta,gamma,phi,y
      real(dp),intent(out)::l,b,s(:)
      integer,intent(in)::m,trend,season
      real(dp)::q,p,r,t,phib
      integer::j
      select case(trend)
      case(ETS_NONE)
      q=oldl
      phib=0.0_dp
      case(ETS_ADD)
      phib=phi*oldb
      q=oldl+phib
      case(ETS_MULT)
      phib=oldb**phi
      q=oldl*phib
      end select
      select case(season)
      case(ETS_NONE)
      p=y
      case(ETS_ADD)
      p=y-olds(m)
      case(ETS_MULT)
      if(abs(olds(m))<1.0e-12_dp)then
      p=huge(1.0_dp)
      else
      p=y/olds(m)
      end if
      end select
      l=q+alpha*(p-q)
      b=oldb
      if(trend/=ETS_NONE)then
         if(trend==ETS_ADD)then
         r=l-oldl
         else
         if(abs(oldl)<1.0e-12_dp)then
         r=1.0_dp
         else
         r=l/oldl
         end if
         end if
         if(alpha>1.0e-12_dp)b=phib+(beta/alpha)*(r-phib)
      end if
      if(season/=ETS_NONE)then
         if(season==ETS_ADD)then
         t=y-q
         else
         if(abs(q)<1.0e-12_dp)then
         t=1.0_dp
         else
         t=y/q
         end if
         end if
         s(1)=olds(m)+gamma*(t-olds(m))
         do j=2,m
         s(j)=olds(j-1)
         end do
      end if
   end subroutine

   subroutine ets_calc(y,init,m,error_type,trend_type,season_type,alpha,beta,gamma,phi, &
                       residuals,fitted,states,likelihood,amse)
      real(dp),intent(in)::y(:),init(:),alpha,beta,gamma,phi
      integer,intent(in)::m,error_type,trend_type,season_type
      real(dp),allocatable,intent(out)::residuals(:),fitted(:),states(:,:),amse(:)
      real(dp),intent(out)::likelihood
      real(dp),allocatable::s(:),olds(:),f(:),den(:)
      real(dp)::l,b,oldl,oldb,tmp,lik2
      integer::n,nstate,i,j,nmse,idx
      n=size(y)
      nmse=min(30,n)
      nstate=1+merge(1,0,trend_type/=ETS_NONE)+merge(m,0,season_type/=ETS_NONE)
      if(size(init)<nstate)error stop 'ets_calc: initial state too short'
      allocate(residuals(n),fitted(n),states(nstate,n+1),amse(nmse),den(nmse),f(nmse))
      allocate(s(max(1,m)),olds(max(1,m)))
      s=1.0_dp
      olds=s
      l=init(1)
      idx=2
      b=0.0_dp
      if(trend_type/=ETS_NONE)then
      b=init(idx)
      idx=idx+1
      end if
      if(season_type/=ETS_NONE)s(1:m)=init(idx:idx+m-1)
      states(:,1)=init(1:nstate)
      likelihood=0.0_dp
      lik2=0.0_dp
      amse=0.0_dp
      den=0.0_dp
      do i=1,n
         oldl=l
         oldb=b
         olds=s
         call ets_state_forecast(oldl,oldb,olds,m,trend_type,season_type,phi,f)
         fitted(i)=f(1)
         if(error_type==ETS_ADD)then
         residuals(i)=y(i)-fitted(i)
         else
         if(abs(fitted(i))<1.0e-12_dp)then
         residuals(i)=1.0e10_dp
         else
         residuals(i)=(y(i)-fitted(i))/fitted(i)
         end if
         end if
         do j=1,nmse
            if(i+j-1<=n)then
            den(j)=den(j)+1.0_dp
            tmp=y(i+j-1)-f(j)
            amse(j)=amse(j)+(tmp*tmp-amse(j))/den(j)
            end if
         end do
         call ets_update(oldl,l,oldb,b,olds,s,m,trend_type,season_type,alpha,beta,gamma,phi,y(i))
         idx=1
         states(idx,i+1)=l
         idx=idx+1
         if(trend_type/=ETS_NONE)then
         states(idx,i+1)=b
         idx=idx+1
         end if
         if(season_type/=ETS_NONE)states(idx:idx+m-1,i+1)=s(1:m)
         likelihood=likelihood+residuals(i)**2
         lik2=lik2+log(max(abs(fitted(i)),1.0e-300_dp))
      end do
      likelihood=real(n,dp)*log(max(likelihood,1.0e-300_dp))
      if(error_type==ETS_MULT)likelihood=likelihood+2.0_dp*lik2
   end subroutine

   function initial_state(y,m,trend,season) result(init)
      real(dp),intent(in)::y(:)
      integer,intent(in)::m,trend,season
      real(dp),allocatable::init(:)
      real(dp)::l,b
      integer::nstate,j,idx,mm
      mm=max(1,m)
      nstate=1+merge(1,0,trend/=ETS_NONE)+merge(mm,0,season/=ETS_NONE)
      allocate(init(nstate))
      init=0.0_dp
      if(season/=ETS_NONE .and. size(y)>=mm)then
      l=mean_value(y(1:mm))
      else
      l=y(1)
      end if
      init(1)=l
      idx=2
      if(trend/=ETS_NONE)then
         if(size(y)>mm)then
            if(trend==ETS_ADD)b=(mean_value(y(mm+1:min(size(y),2*mm)))-mean_value(y(1:mm)))/real(mm,dp)
            if(trend==ETS_MULT)b=max(1.0e-4_dp,(mean_value(y(mm+1:min(size(y),2*mm)))/max(abs(mean_value(y(1:mm))), &
               & 1.0e-8_dp))**(1.0_dp/real(mm,dp)))
         else
         b=merge(y(min(2,size(y)))-y(1),1.0_dp,trend==ETS_ADD)
         end if
         init(idx)=b
         idx=idx+1
      end if
      if(season/=ETS_NONE)then
         do j=1,mm
            if(season==ETS_ADD)then
            init(idx+j-1)=y(j)-l
            else
            init(idx+j-1)=y(j)/max(abs(l),1.0e-8_dp)
            end if
         end do
         if(season==ETS_ADD)init(idx:idx+mm-1)=init(idx:idx+mm-1)-mean_value(init(idx:idx+mm-1))
         if(season==ETS_MULT)init(idx:idx+mm-1)=init(idx:idx+mm-1)/max(mean_value(init(idx:idx+mm-1)),1.0e-8_dp)
      end if
   end function


   logical function ets_admissible(alpha,beta,gamma,phi,m,has_trend,has_season) result(ok)
      real(dp),intent(in)::alpha,beta,gamma,phi
      integer,intent(in)::m
      logical,intent(in)::has_trend,has_season
      real(dp),allocatable::pcoef(:)
      real(dp)::b,rmax
      integer::j,info
      ok=.false.
      if(phi<=tiny(1.0_dp) .or. phi>1.0_dp+1.0e-8_dp)return
      if(.not.has_season)then
         if(alpha<1.0_dp-1.0_dp/phi .or. alpha>1.0_dp+1.0_dp/phi)return
         if(has_trend)then
            if(beta<alpha*(phi-1.0_dp) .or. beta>(1.0_dp+phi)*(2.0_dp-alpha))return
         end if
      else if(m>1)then
         b=merge(beta,0.0_dp,has_trend)
         if(gamma<max(1.0_dp-1.0_dp/phi-alpha,0.0_dp) .or. gamma>1.0_dp+1.0_dp/phi-alpha)return
         if(alpha<1.0_dp-1.0_dp/phi-gamma*(1.0_dp-real(m,dp)+phi+phi*real(m,dp))/(2.0_dp*phi*real(m,dp)))return
         if(b<-(1.0_dp-phi)*(gamma/real(m,dp)+alpha))return
         allocate(pcoef(m+2))
         pcoef=0.0_dp
         pcoef(1)=phi*(1.0_dp-alpha-gamma)
         pcoef(2)=alpha+b-alpha*phi+gamma-1.0_dp
         if(m>2)then
            do j=3,m
               pcoef(j)=alpha+b-alpha*phi
            end do
         end if
         pcoef(m+1)=alpha+b-phi
         pcoef(m+2)=1.0_dp
         rmax=polynomial_max_root_mod(pcoef,info)
         if(info/=0 .or. rmax>1.0_dp+1.0e-8_dp)return
      end if
      ok=.true.
   end function ets_admissible

   subroutine unpack_ets_vector(x,a,b,g,p,init,valid)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::a,b,g,p
      real(dp),allocatable,intent(out)::init(:)
      logical,intent(out)::valid
      integer::k,idx,j,nstate,nfree
      real(dp)::last
      k=1
      a=x(k)
      k=k+1
      b=0.0_dp
      g=0.0_dp
      p=1.0_dp
      if(ctx_trend/=ETS_NONE)then
         b=x(k)
         k=k+1
      end if
      if(ctx_season/=ETS_NONE)then
         g=x(k)
         k=k+1
      end if
      if(ctx_damped)then
         p=x(k)
         k=k+1
      end if
      valid=.true.
      if(b>a .and. ctx_trend/=ETS_NONE)valid=.false.
      if(g>1.0_dp-a .and. ctx_season/=ETS_NONE)valid=.false.
      if(valid)valid=ets_admissible(a,b,g,p,ctx_m,ctx_trend/=ETS_NONE,ctx_season/=ETS_NONE)
      nstate=1+merge(1,0,ctx_trend/=ETS_NONE)+merge(ctx_m,0,ctx_season/=ETS_NONE)
      allocate(init(nstate))
      if(.not.ctx_opt_init)then
         init=ctx_init
         return
      end if
      idx=1
      init(idx)=x(k)
      k=k+1
      idx=idx+1
      if(ctx_trend/=ETS_NONE)then
         init(idx)=x(k)
         k=k+1
         idx=idx+1
         if(ctx_trend==ETS_MULT .and. init(idx-1)<=0.0_dp)valid=.false.
      end if
      if(ctx_season/=ETS_NONE)then
         nfree=ctx_m-1
         if(nfree>0)then
            init(idx:idx+nfree-1)=x(k:k+nfree-1)
            k=k+nfree
         end if
         if(ctx_season==ETS_ADD)then
            last=-sum(init(idx:idx+nfree-1))
         else
            last=real(ctx_m,dp)-sum(init(idx:idx+nfree-1))
            if(last<=0.0_dp .or. any(init(idx:idx+nfree-1)<=0.0_dp))valid=.false.
         end if
         init(idx+ctx_m-1)=last
      end if
   end subroutine unpack_ets_vector

   function ets_fit(y,m,error_type,trend_type,season_type,damped,optimize,optimize_initial,alpha,beta,gamma,phi) result(model)
      real(dp),intent(in)::y(:)
      integer,intent(in),optional::m,error_type,trend_type,season_type
      logical,intent(in),optional::damped,optimize,optimize_initial
      real(dp),intent(in),optional::alpha,beta,gamma,phi
      type(ets_model)::model
      integer::mm,er,tr,se,k,npar,nsmooth,nfree,idx,j
      logical::dam,opt,optinit,valid
      real(dp),allocatable::x(:),lo(:),hi(:),res(:),fit(:),states(:,:),amse(:),init0(:),init(:)
      real(dp)::lik,obj,rng,a,b,g,pv
      mm=1
      if(present(m))mm=max(1,m)
      er=ETS_ADD
      if(present(error_type))er=error_type
      tr=ETS_NONE
      if(present(trend_type))tr=trend_type
      se=ETS_NONE
      if(present(season_type))se=season_type
      dam=.false.
      if(present(damped))dam=damped
      opt=.true.
      if(present(optimize))opt=optimize
      optinit=.true.
      if(present(optimize_initial))optinit=optimize_initial
      if(se/=ETS_NONE .and. size(y)<2*mm)se=ETS_NONE
      ctx_y=y
      ctx_m=mm
      ctx_err=er
      ctx_trend=tr
      ctx_season=se
      ctx_damped=dam
      ctx_init=initial_state(y,mm,tr,se)
      init0=ctx_init
      ctx_opt_init=optinit
      nsmooth=1+merge(1,0,tr/=ETS_NONE)+merge(1,0,se/=ETS_NONE)+merge(1,0,dam)
      ctx_nsmooth=nsmooth
      nfree=merge(1+merge(1,0,tr/=ETS_NONE)+merge(mm-1,0,se/=ETS_NONE),0,optinit)
      npar=nsmooth+nfree
      allocate(x(npar),lo(npar),hi(npar))
      k=1
      x(k)=0.2_dp
      lo(k)=1.0e-4_dp
      hi(k)=0.9999_dp
      if(present(alpha))then
         x(k)=alpha
         lo(k)=alpha
         hi(k)=alpha
      end if
      k=k+1
      if(tr/=ETS_NONE)then
         x(k)=0.05_dp
         lo(k)=0.0_dp
         hi(k)=0.9999_dp
         if(present(beta))then
            x(k)=beta
            lo(k)=beta
            hi(k)=beta
         end if
         k=k+1
      end if
      if(se/=ETS_NONE)then
         x(k)=0.05_dp
         lo(k)=0.0_dp
         hi(k)=0.9999_dp
         if(present(gamma))then
            x(k)=gamma
            lo(k)=gamma
            hi(k)=gamma
         end if
         k=k+1
      end if
      if(dam)then
         x(k)=0.98_dp
         lo(k)=0.80_dp
         hi(k)=0.98_dp
         if(present(phi))then
            x(k)=phi
            lo(k)=phi
            hi(k)=phi
         end if
         k=k+1
      end if
      if(optinit)then
         rng=max(maxval(y)-minval(y),1.0_dp)
         idx=1
         x(k)=init0(idx)
         if(er==ETS_MULT .or. tr==ETS_MULT .or. se==ETS_MULT)then
            lo(k)=1.0e-6_dp
            hi(k)=max(10.0_dp,maxval(abs(y))*4.0_dp)
         else
            lo(k)=minval(y)-4.0_dp*rng
            hi(k)=maxval(y)+4.0_dp*rng
         end if
         k=k+1
         idx=idx+1
         if(tr/=ETS_NONE)then
            x(k)=init0(idx)
            if(tr==ETS_MULT)then
               lo(k)=1.0e-4_dp
               hi(k)=4.0_dp
            else
               lo(k)=-4.0_dp*rng
               hi(k)=4.0_dp*rng
            end if
            k=k+1
            idx=idx+1
         end if
         if(se/=ETS_NONE .and. mm>1)then
            do j=1,mm-1
               x(k)=init0(idx+j-1)
               if(se==ETS_MULT)then
                  lo(k)=1.0e-4_dp
                  hi(k)=real(mm,dp)
               else
                  lo(k)=-4.0_dp*rng
                  hi(k)=4.0_dp*rng
               end if
               k=k+1
            end do
         end if
      end if
      if(opt .and. npar>0)call pattern_search(ets_objective,x,lo,hi,maxit=260,tol=2.0e-5_dp,fval=obj)
      call unpack_ets_vector(x,a,b,g,pv,init,valid)
      if(.not.valid)init=init0
      model%alpha=a
      model%beta=b
      model%gamma=g
      model%phi=pv
      call ets_calc(y,init,mm,er,tr,se,model%alpha,model%beta,model%gamma,model%phi,res,fit,states,lik,amse)
      model%error_type=er
      model%trend_type=tr
      model%season_type=se
      model%m=mm
      model%damped=dam
      model%state=states(:,size(states,2))
      model%states=states
      model%fitted=fit
      model%residuals=res
      npar=nfree
      if(.not.present(alpha))npar=npar+1
      if(tr/=ETS_NONE .and. .not.present(beta))npar=npar+1
      if(se/=ETS_NONE .and. .not.present(gamma))npar=npar+1
      if(dam .and. .not.present(phi))npar=npar+1
      model%sigma2=sum(res*res)/real(max(1,size(y)-npar),dp)
      ! forecast::ets uses loglik = -0.5 * the native etscalc likelihood criterion.
      model%loglik=-0.5_dp*lik
      npar=npar+1
      model%aic=-2.0_dp*model%loglik+2.0_dp*real(npar,dp)
      if(size(y)>npar+1)then
         model%aicc=model%aic+2.0_dp*real(npar*(npar+1),dp)/real(size(y)-npar-1,dp)
      else
         model%aicc=huge(1.0_dp)
      end if
      model%bic=-2.0_dp*model%loglik+log(real(size(y),dp))*real(npar,dp)
   end function ets_fit

   function ets_objective(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v,a,b,g,pv
      real(dp),allocatable::r(:),f(:),st(:,:),am(:),init(:)
      logical::valid
      call unpack_ets_vector(x,a,b,g,pv,init,valid)
      if(.not.valid)then
         v=1.0e50_dp
         return
      end if
      call ets_calc(ctx_y,init,ctx_m,ctx_err,ctx_trend,ctx_season,a,b,g,pv,r,f,st,v,am)
      if(any(abs(f)>huge(1.0_dp)/100.0_dp))v=1.0e50_dp
      if(.not.(v<huge(1.0_dp)))v=1.0e50_dp
   end function ets_objective

   function ets_auto(y,m,allow_multiplicative,damped,allow_multiplicative_trend,restrict) result(best)
      real(dp),intent(in)::y(:)
      integer,intent(in),optional::m
      logical,intent(in),optional::allow_multiplicative,damped,allow_multiplicative_trend,restrict
      type(ets_model)::best,cand
      integer::mm,er,tr,se,emax,smax,tmax,idamp
      logical::mult,trydamp,allowmt,res
      mm=1
      if(present(m))mm=max(1,m)
      mult=.true.
      if(present(allow_multiplicative))mult=allow_multiplicative
      trydamp=.true.
      if(present(damped))trydamp=damped
      allowmt=.false.
      if(present(allow_multiplicative_trend))allowmt=allow_multiplicative_trend
      res=.true.
      if(present(restrict))res=restrict
      emax=merge(ETS_MULT,ETS_ADD,mult.and.all(y>0.0_dp))
      smax=merge(ETS_MULT,ETS_ADD,mult.and.all(y>0.0_dp))
      tmax=merge(ETS_MULT,ETS_ADD,allowmt.and.all(y>0.0_dp))
      best%aicc=huge(1.0_dp)
      do er=ETS_ADD,emax
         do tr=ETS_NONE,tmax
            do se=ETS_NONE,merge(smax,ETS_NONE,mm>1)
               if(se==ETS_MULT .and. any(y<=0.0_dp))cycle
               if(res)then
                  if(er==ETS_ADD .and. (tr==ETS_MULT .or. se==ETS_MULT))cycle
                  if(er==ETS_MULT .and. tr==ETS_MULT .and. se==ETS_ADD)cycle
               end if
               do idamp=0,merge(1,0,trydamp.and.tr/=ETS_NONE)
                  cand=ets_fit(y,mm,er,tr,se,idamp==1,.true.,.true.)
                  if(cand%aicc<best%aicc)best=cand
               end do
            end do
         end do
      end do
   end function ets_auto

   subroutine ets_class12_variance(model,h,mu,var,ok)
      type(ets_model),intent(in)::model
      integer,intent(in)::h
      real(dp),intent(in)::mu(:)
      real(dp),intent(out)::var(:)
      logical,intent(out)::ok
      real(dp),allocatable::sb(:),ss(:),oldsb(:),oldss(:),c(:),theta(:),f0(:),f1(:)
      real(dp)::lb,bb,ls,bs,oldlb,oldbb,oldls,oldbs,csum
      integer::idx,j,k
      ok=model%trend_type/=ETS_MULT .and. model%season_type/=ETS_MULT
      if(.not.ok)return
      if(size(mu)/=h .or. size(var)/=h)error stop 'ets_class12_variance: shape mismatch'
      allocate(sb(max(1,model%m)),ss(max(1,model%m)),oldsb(max(1,model%m)),oldss(max(1,model%m)))
      allocate(f0(1),f1(1),c(max(0,h-1)))
      idx=1
      lb=model%state(idx)
      ls=lb
      idx=idx+1
      bb=0.0_dp
      bs=0.0_dp
      if(model%trend_type/=ETS_NONE)then
         bb=model%state(idx)
         bs=bb
         idx=idx+1
      end if
      sb=1.0_dp
      ss=1.0_dp
      if(model%season_type/=ETS_NONE)then
         sb(1:model%m)=model%state(idx:idx+model%m-1)
         ss=sb
      end if
      if(h>1)then
         oldlb=lb
         oldbb=bb
         oldsb=sb
         oldls=ls
         oldbs=bs
         oldss=ss
         call ets_state_forecast(oldlb,oldbb,oldsb,model%m,model%trend_type,model%season_type,model%phi,f0)
         call ets_update(oldlb,lb,oldbb,bb,oldsb,sb,model%m,model%trend_type,model%season_type, &
            & model%alpha,model%beta,model%gamma,model%phi,f0(1))
         call ets_update(oldls,ls,oldbs,bs,oldss,ss,model%m,model%trend_type,model%season_type, &
            & model%alpha,model%beta,model%gamma,model%phi,f0(1)+1.0_dp)
         do k=1,h-1
            call ets_state_forecast(lb,bb,sb,model%m,model%trend_type,model%season_type,model%phi,f0)
            call ets_state_forecast(ls,bs,ss,model%m,model%trend_type,model%season_type,model%phi,f1)
            c(k)=f1(1)-f0(1)
            if(k<h-1)then
               oldlb=lb
               oldbb=bb
               oldsb=sb
               oldls=ls
               oldbs=bs
               oldss=ss
               call ets_update(oldlb,lb,oldbb,bb,oldsb,sb,model%m,model%trend_type,model%season_type, &
                  & model%alpha,model%beta,model%gamma,model%phi,f0(1))
               call ets_update(oldls,ls,oldbs,bs,oldss,ss,model%m,model%trend_type,model%season_type, &
                  & model%alpha,model%beta,model%gamma,model%phi,f1(1))
            end if
         end do
      end if
      if(model%error_type==ETS_ADD)then
         var(1)=max(model%sigma2,0.0_dp)
         csum=0.0_dp
         do j=2,h
            csum=csum+c(j-1)**2
            var(j)=max(model%sigma2,0.0_dp)*(1.0_dp+csum)
         end do
      else if(model%error_type==ETS_MULT)then
         allocate(theta(h))
         theta(1)=mu(1)**2
         do j=2,h
            theta(j)=mu(j)**2
            do k=1,j-1
               theta(j)=theta(j)+max(model%sigma2,0.0_dp)*c(k)**2*theta(j-k)
            end do
         end do
         var=max(0.0_dp,(1.0_dp+max(model%sigma2,0.0_dp))*theta-mu**2)
      else
         ok=.false.
      end if
   end subroutine ets_class12_variance

   function kron_matrix(a,b) result(c)
      real(dp),intent(in)::a(:,:),b(:,:)
      real(dp),allocatable::c(:,:)
      integer::i,j,br,bc
      br=size(b,1)
      bc=size(b,2)
      allocate(c(size(a,1)*br,size(a,2)*bc))
      do i=1,size(a,1)
         do j=1,size(a,2)
            c((i-1)*br+1:i*br,(j-1)*bc+1:j*bc)=a(i,j)*b
         end do
      end do
   end function kron_matrix

   subroutine ets_class3_variance(model,h,mu,var,ok)
      type(ets_model),intent(in)::model
      integer,intent(in)::h
      real(dp),intent(out)::mu(:),var(:)
      logical,intent(out)::ok
      real(dp),allocatable::h1(:,:),h2(:,:),f1(:,:),g1(:,:),f2(:,:),g2(:,:)
      real(dp),allocatable::mh(:,:),vh(:,:),h21(:,:),f21(:,:),g21(:,:),km(:,:),vecmh(:),outer(:,:)
      real(dp),allocatable::tmp(:,:),tmp2(:,:)
      integer::n1,m,p,idx,i
      real(dp)::sig
      ok=.false.
      if(model%error_type/=ETS_MULT .or. model%season_type/=ETS_MULT)return
      if(model%trend_type==ETS_MULT)return
      if(size(mu)/=h.or.size(var)/=h)return
      m=model%m
      if(m<2)return
      n1=merge(2,1,model%trend_type/=ETS_NONE)
      p=n1+m
      if(size(model%state)<p)return
      allocate(h1(1,n1),h2(1,m),f1(n1,n1),g1(n1,n1),f2(m,m),g2(m,m))
      h1=1.0_dp
      h2=0.0_dp
      h2(1,m)=1.0_dp
      f1=0.0_dp
      g1=0.0_dp
      if(n1==1)then
         f1(1,1)=1.0_dp
         g1(1,1)=model%alpha
      else
         f1(1,:)=[1.0_dp,1.0_dp]
         f1(2,2)=merge(model%phi,1.0_dp,model%damped)
         g1(1,:)=model%alpha
         g1(2,:)=model%beta
      end if
      f2=0.0_dp
      f2(1,m)=1.0_dp
      do i=2,m
         f2(i,i-1)=1.0_dp
      end do
      g2=0.0_dp
      g2(1,m)=model%gamma
      allocate(mh(n1,m))
      idx=1
      mh=matmul(reshape(model%state(idx:idx+n1-1),[n1,1]), &
         reshape(model%state(idx+n1:idx+n1+m-1),[1,m]))
      allocate(vh(n1*m,n1*m))
      vh=0.0_dp
      h21=kron_matrix(h2,h1)
      f21=kron_matrix(f2,f1)
      g21=kron_matrix(g2,g1)
      km=kron_matrix(g2,f1)+kron_matrix(f2,g1)
      sig=max(model%sigma2,0.0_dp)
      do i=1,h
         mu(i)=dot_product(h1(1,:),matmul(mh,h2(1,:)))
         var(i)=(1.0_dp+sig)*dot_product(h21(1,:),matmul(vh,h21(1,:)))+sig*mu(i)**2
         vecmh=reshape(mh,[n1*m])
         outer=matmul(reshape(vecmh,[n1*m,1]),reshape(vecmh,[1,n1*m]))
         tmp=matmul(matmul(f21,vh),transpose(f21))
         tmp2=matmul(matmul(f21,vh),transpose(g21))+matmul(matmul(g21,vh),transpose(f21))
         vh=tmp+sig*(tmp2+matmul(matmul(km,vh+outer),transpose(km))+ &
            sig*matmul(matmul(g21,3.0_dp*vh+2.0_dp*outer),transpose(g21)))
         mh=matmul(matmul(f1,mh),transpose(f2))+sig*matmul(matmul(g1,mh),transpose(g2))
      end do
      var=max(var,0.0_dp)
      ok=.true.
   end subroutine ets_class3_variance

   function ets_forecast(model,h,levels) result(fc)
      type(ets_model),intent(in)::model
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      real(dp),allocatable::lev(:),s(:),fvar(:)
      real(dp)::l,b,z
      integer::idx,j
      logical::exact_var
      if(present(levels))then
      lev=levels
      else
      lev=[80.0_dp,95.0_dp]
      end if
      idx=1
      l=model%state(idx)
      idx=idx+1
      b=0.0_dp
      if(model%trend_type/=ETS_NONE)then
      b=model%state(idx)
      idx=idx+1
      end if
      allocate(s(max(1,model%m)))
      s=1.0_dp
      if(model%season_type/=ETS_NONE)s(1:model%m)=model%state(idx:idx+model%m-1)
      allocate(fc%mean(h),fc%se(h),fc%level(size(lev)),fc%lower(h,size(lev)),fc%upper(h,size(lev)))
      fc%level=lev
      call ets_state_forecast(l,b,s,model%m,model%trend_type,model%season_type,model%phi,fc%mean)
      allocate(fvar(h))
      call ets_class12_variance(model,h,fc%mean,fvar,exact_var)
      if(.not.exact_var)call ets_class3_variance(model,h,fc%mean,fvar,exact_var)
      if(exact_var)then
         fc%se=sqrt(max(fvar,0.0_dp))
      else
         fc%se=sqrt(max(model%sigma2,0.0_dp))*sqrt([(real(j,dp),j=1,h)])
         if(model%error_type==ETS_MULT)then
            fc%se=abs(fc%mean)*sqrt(max(model%sigma2,0.0_dp))*sqrt([(real(j,dp),j=1,h)])
         end if
      end if
      do j=1,size(lev)
      z=normal_quantile(0.5_dp+lev(j)/200.0_dp)
      fc%lower(:,j)=fc%mean-z*fc%se
      fc%upper(:,j)=fc%mean+z*fc%se
      end do
   end function

   function ets_forecast_simulated(model,h,levels,npaths,bootstrap,innov) result(fc)
      type(ets_model),intent(in)::model
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:),innov(:,:)
      integer,intent(in),optional::npaths
      logical,intent(in),optional::bootstrap
      type(forecast_result)::fc
      real(dp),allocatable::lev(:),paths(:,:),eps(:),sim(:),vals(:)
      real(dp)::u1,u2,mbar,prob
      integer::np,j,k,idx,nres
      logical::boot
      np=5000
      if(present(npaths))np=max(2,npaths)
      boot=.false.
      if(present(bootstrap))boot=bootstrap
      if(present(levels))then
         lev=levels
      else
         lev=[80.0_dp,95.0_dp]
      end if
      if(present(innov))then
         if(size(innov,1)/=h.or.size(innov,2)/=np)error stop 'ets_forecast_simulated: innov shape mismatch'
      end if
      fc=ets_forecast(model,h,lev)
      allocate(paths(h,np),eps(h))
      nres=size(model%residuals)
      do k=1,np
         if(present(innov))then
            eps=innov(:,k)
         else if(boot)then
            do j=1,h
               call random_number(u1)
               idx=1+min(nres-1,int(u1*real(nres,dp)))
               eps(j)=model%residuals(idx)
            end do
         else
            j=1
            do while(j<=h)
               call random_number(u1)
               call random_number(u2)
               u1=max(u1,1.0e-12_dp)
               eps(j)=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)*sqrt(max(model%sigma2,0.0_dp))
               if(j+1<=h)then
                  eps(j+1)=sqrt(-2.0_dp*log(u1))*sin(2.0_dp*acos(-1.0_dp)*u2)*sqrt(max(model%sigma2,0.0_dp))
               end if
               j=j+2
            end do
         end if
         sim=ets_simulate(model,eps)
         paths(:,k)=sim
      end do
      do j=1,h
         mbar=sum(paths(j,:))/real(np,dp)
         fc%se(j)=sqrt(sum((paths(j,:)-mbar)**2)/real(np-1,dp))
      end do
      do k=1,size(lev)
         prob=(100.0_dp-lev(k))/200.0_dp
         do j=1,h
            vals=paths(j,:)
            fc%lower(j,k)=quantile_type8(vals,prob)
            fc%upper(j,k)=quantile_type8(vals,1.0_dp-prob)
         end do
      end do
   end function ets_forecast_simulated

   function ets_simulate(model,innov) result(y)
      type(ets_model),intent(in)::model
      real(dp),intent(in)::innov(:)
      real(dp),allocatable::y(:)
      real(dp),allocatable::s(:),olds(:),f(:)
      real(dp)::l,b,oldl,oldb
      integer::i,idx
      allocate(y(size(innov)),s(max(1,model%m)),olds(max(1,model%m)),f(1))
      idx=1
      l=model%state(idx)
      idx=idx+1
      b=0.0_dp
      if(model%trend_type/=ETS_NONE)then
      b=model%state(idx)
      idx=idx+1
      end if
      s=1.0_dp
      if(model%season_type/=ETS_NONE)s(1:model%m)=model%state(idx:idx+model%m-1)
      do i=1,size(innov)
         oldl=l
         oldb=b
         olds=s
         call ets_state_forecast(oldl,oldb,olds,model%m,model%trend_type,model%season_type,model%phi,f)
         if(model%error_type==ETS_ADD)then
         y(i)=f(1)+innov(i)
         else
         y(i)=f(1)*(1.0_dp+innov(i))
         end if
         call ets_update(oldl,l,oldb,b,olds,s,model%m,model%trend_type,model%season_type,model%alpha,model%beta,model%gamma, &
            & model%phi,y(i))
      end do
   end function
   function ses_fit(y,alpha) result(model)
      real(dp),intent(in)::y(:)
      real(dp),intent(in),optional::alpha
      type(ets_model)::model
      if(present(alpha))then
         model=ets_fit(y,1,ETS_ADD,ETS_NONE,ETS_NONE,.false.,.true.,.true.,alpha=alpha)
      else
         model=ets_fit(y,1,ETS_ADD,ETS_NONE,ETS_NONE,.false.,.true.,.true.)
      end if
   end function
   function holt_fit(y,damped) result(model)
      real(dp),intent(in)::y(:)
      logical,intent(in),optional::damped
      type(ets_model)::model
      logical::d
      d=.false.
      if(present(damped))d=damped
      model=ets_fit(y,1,ETS_ADD,ETS_ADD,ETS_NONE,d,.true.)
   end function
   function hw_fit(y,m,multiplicative,damped) result(model)
      real(dp),intent(in)::y(:)
      integer,intent(in)::m
      logical,intent(in),optional::multiplicative,damped
      type(ets_model)::model
      logical::mul,d
      mul=.false.
      if(present(multiplicative))mul=multiplicative
      d=.false.
      if(present(damped))d=damped
      model=ets_fit(y,m,merge(ETS_MULT,ETS_ADD,mul),ETS_ADD,merge(ETS_MULT,ETS_ADD,mul),d,.true.)
   end function
end module forecast_ets
