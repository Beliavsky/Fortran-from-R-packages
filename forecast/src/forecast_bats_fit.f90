module forecast_bats_fit
   use forecast_kinds, only : dp
   use forecast_types, only : bats_model, forecast_result, arima_model
   use forecast_bats, only : bats_make_w, bats_make_g, bats_make_f, tbats_make_w, tbats_make_g, &
      tbats_make_f, state_space_filter, state_space_forecast
   use forecast_linalg, only : least_squares
   use forecast_optimize, only : pattern_search
   use forecast_transforms, only : boxcox, inv_boxcox, boxcox_lambda_guerrero
   use forecast_stats, only : normal_quantile
   use forecast_arima, only : arima_fit
   implicit none
   private
   public :: bats_fit, tbats_fit, tbats_fit_real, bats_forecast, bats_auto, tbats_auto, bats_refit, tbats_refit

   real(dp), allocatable, save :: fit_y(:), fit_ar0(:), fit_ma0(:), fit_periods_real(:)
   integer, allocatable, save :: fit_periods(:), fit_k(:)
   logical, save :: fit_trend, fit_damped, fit_tbats, fit_boxcox
contains
   subroutine seed_states(y,w,f,g,periods,is_tbats,p,q,x0)
      real(dp),intent(in)::y(:),w(:),f(:,:),g(:)
      integer,intent(in)::periods(:),p,q
      logical,intent(in)::is_tbats
      real(dp),allocatable,intent(out)::x0(:)
      real(dp),allocatable::zero(:),yh(:),e(:),states(:,:),dmat(:,:),wt(:,:),xred(:,:),coef(:),resid(:)
      logical,allocatable::keep(:)
      integer::nstate,n,t,j,pos,adj,info,rank,nkeep
      n=size(y)
      nstate=size(w)
      allocate(zero(nstate))
      zero=0.0_dp
      call state_space_filter(y,w,f,g,zero,yh,e,states)
      allocate(dmat(nstate,nstate))
      dmat=f
      do j=1,nstate
         dmat(:,j)=dmat(:,j)-g*w(j)
      end do
      allocate(wt(n,nstate))
      wt(1,:)=w
      do t=2,n
         wt(t,:)=matmul(wt(t-1,:),dmat)
      end do
      allocate(keep(nstate))
      keep=.true.
      if(p+q>0)keep(nstate-p-q+1:nstate)=.false.
      if(.not.is_tbats .and. size(periods)>0)then
         adj=nstate-sum(periods)-p-q-1
         pos=1+adj
         do j=1,size(periods)
            pos=pos+periods(j)
            keep(pos)=.false.
         end do
      end if
      nkeep=count(keep)
      allocate(xred(n,nkeep))
      nkeep=0
      do j=1,nstate
         if(keep(j))then
            nkeep=nkeep+1
            xred(:,nkeep)=wt(:,j)
         end if
      end do
      call least_squares(xred,e,coef,resid,rank,info)
      allocate(x0(nstate))
      x0=0.0_dp
      if(info/=0)return
      nkeep=0
      do j=1,nstate
         if(keep(j))then
            nkeep=nkeep+1
            x0(j)=coef(nkeep)
         end if
      end do
      if(.not.is_tbats .and. size(periods)>0)then
         adj=nstate-sum(periods)-p-q-1
         pos=1+adj
         do j=1,size(periods)
            x0(pos+periods(j))=-sum(x0(pos+1:pos+periods(j)-1))
            pos=pos+periods(j)
         end do
      end if
   end subroutine seed_states

   logical function arma_coefficients_ok(ar,ma) result(ok)
      real(dp),intent(in)::ar(:),ma(:)
      real(dp),allocatable::a(:),b(:)
      integer::k,j,ncheck
      ok=.true.
      ncheck=256
      if(size(ar)>0)then
         allocate(a(0:ncheck))
         a=0.0_dp
         a(0)=1.0_dp
         do k=1,ncheck
            do j=1,min(k,size(ar))
               a(k)=a(k)+ar(j)*a(k-j)
            end do
            if(abs(a(k))>1.0e6_dp)then
               ok=.false.
               return
            end if
         end do
         if(sum(a*a)>1.0e8_dp)then
            ok=.false.
            return
         end if
      end if
      if(size(ma)>0)then
         allocate(b(0:ncheck))
         b=0.0_dp
         b(0)=1.0_dp
         do k=1,ncheck
            do j=1,min(k,size(ma))
               b(k)=b(k)-ma(j)*b(k-j)
            end do
            if(abs(b(k))>1.0e6_dp)then
               ok=.false.
               return
            end if
         end do
      end if
   end function arma_coefficients_ok

   subroutine unpack_params(x,alpha,beta,phi,g1,g2,ar,ma,lambda)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::alpha,beta,phi,lambda
      real(dp),allocatable,intent(out)::g1(:),g2(:),ar(:),ma(:)
      integer::pos,ns,p,q
      pos=1
      alpha=x(pos)
      pos=pos+1
      beta=0.0_dp
      phi=1.0_dp
      if(fit_trend)then
         beta=x(pos)
         pos=pos+1
         if(fit_damped)then
            phi=x(pos)
            pos=pos+1
         end if
      end if
      ns=size(fit_periods)
      allocate(g1(ns))
      g1=0.0_dp
      if(ns>0)then
         g1=x(pos:pos+ns-1)
         pos=pos+ns
      end if
      allocate(g2(ns))
      g2=0.0_dp
      if(fit_tbats .and. ns>0)then
         g2=x(pos:pos+ns-1)
         pos=pos+ns
      end if
      p=size(fit_ar0)
      q=size(fit_ma0)
      allocate(ar(p),ma(q))
      if(p>0)then
         ar=x(pos:pos+p-1)
         pos=pos+p
      end if
      if(q>0)then
         ma=x(pos:pos+q-1)
         pos=pos+q
      end if
      lambda=1.0_dp
      if(fit_boxcox)lambda=x(pos)
   end subroutine unpack_params

   function bats_objective(x) result(value)
      real(dp),intent(in)::x(:)
      real(dp)::value,alpha,beta,phi,lambda,rss
      real(dp),allocatable::g1(:),g2(:),ar(:),ma(:),w(:),g(:),f(:,:),x0(:),yh(:),e(:),states(:,:),yy(:)
      integer::p,q,n
      call unpack_params(x,alpha,beta,phi,g1,g2,ar,ma,lambda)
      if(.not.arma_coefficients_ok(ar,ma))then
         value=huge(1.0_dp)
         return
      end if
      p=size(ar)
      q=size(ma)
      n=size(fit_y)
      if(fit_boxcox)then
         if(any(fit_y<=0.0_dp))then
            value=huge(1.0_dp)
            return
         end if
         yy=boxcox(fit_y,lambda)
      else
         yy=fit_y
      end if
      if(fit_tbats)then
         w=tbats_make_w(fit_trend,phi,fit_k,ar,ma)
         g=tbats_make_g(alpha,fit_trend,beta,fit_k,g1,g2,p,q)
         f=tbats_make_f(alpha,fit_trend,beta,phi,fit_periods_real,fit_k,g1,g2,ar,ma)
      else
         w=bats_make_w(fit_trend,phi,fit_periods,ar,ma)
         g=bats_make_g(alpha,fit_trend,beta,g1,fit_periods,p,q)
         f=bats_make_f(alpha,fit_trend,beta,phi,fit_periods,g1,ar,ma)
      end if
      call seed_states(yy,w,f,g,fit_periods,fit_tbats,p,q,x0)
      call state_space_filter(yy,w,f,g,x0,yh,e,states)
      rss=sum(e*e)
      if(.not.(rss>0.0_dp .and. rss<huge(1.0_dp)))then
         value=huge(1.0_dp)
         return
      end if
      value=real(n,dp)*log(rss)
      if(fit_boxcox)value=value-2.0_dp*(lambda-1.0_dp)*sum(log(fit_y))
   end function bats_objective

   function fit_common(y,periods,periods_real,kvec,trig,use_trend,damped,boxcox_lambda,ar,ma,optimize) result(model)
      real(dp),intent(in)::y(:),periods_real(:)
      integer,intent(in)::periods(:),kvec(:)
      logical,intent(in)::trig,use_trend,damped
      real(dp),intent(in),optional::boxcox_lambda,ar(:),ma(:)
      logical,intent(in),optional::optimize
      type(bats_model)::model
      real(dp),allocatable::x(:),lo(:),hi(:),g1(:),g2(:),arl(:),mal(:),yy(:),x0(:),states(:,:)
      real(dp)::alpha,beta,phi,lambda,criterion,rss
      logical::opt
      integer::ns,pos,p,q,npar,kcount
      fit_y=y
      fit_periods=periods
      fit_periods_real=periods_real
      fit_k=kvec
      fit_tbats=trig
      fit_trend=use_trend
      fit_damped=damped
      if(allocated(fit_ar0))deallocate(fit_ar0)
      if(allocated(fit_ma0))deallocate(fit_ma0)
      if(present(ar))then
         fit_ar0=ar
      else
         allocate(fit_ar0(0))
      end if
      if(present(ma))then
         fit_ma0=ma
      else
         allocate(fit_ma0(0))
      end if
      fit_boxcox=present(boxcox_lambda)
      ns=size(periods)
      p=size(fit_ar0)
      q=size(fit_ma0)
      npar=1+merge(1,0,use_trend)+merge(1,0,use_trend.and.damped)+ns+merge(ns,0,trig)+p+q+merge(1,0,fit_boxcox)
      allocate(x(npar),lo(npar),hi(npar))
      pos=1
      x(pos)=0.09_dp
      lo(pos)=1.0e-5_dp
      hi(pos)=0.99_dp
      pos=pos+1
      if(use_trend)then
         x(pos)=0.05_dp
         lo(pos)=1.0e-6_dp
         hi(pos)=0.99_dp
         pos=pos+1
         if(damped)then
            x(pos)=0.98_dp
            lo(pos)=0.80_dp
            hi(pos)=0.9999_dp
            pos=pos+1
         end if
      end if
      if(ns>0)then
         x(pos:pos+ns-1)=merge(0.0_dp,0.001_dp,trig)
         lo(pos:pos+ns-1)=-0.5_dp
         hi(pos:pos+ns-1)=0.5_dp
         pos=pos+ns
      end if
      if(trig.and.ns>0)then
         x(pos:pos+ns-1)=0.0_dp
         lo(pos:pos+ns-1)=-0.5_dp
         hi(pos:pos+ns-1)=0.5_dp
         pos=pos+ns
      end if
      if(p>0)then
         x(pos:pos+p-1)=fit_ar0
         lo(pos:pos+p-1)=-0.98_dp
         hi(pos:pos+p-1)=0.98_dp
         pos=pos+p
      end if
      if(q>0)then
         x(pos:pos+q-1)=fit_ma0
         lo(pos:pos+q-1)=-0.98_dp
         hi(pos:pos+q-1)=0.98_dp
         pos=pos+q
      end if
      if(fit_boxcox)then
         x(pos)=boxcox_lambda
         lo(pos)=0.0_dp
         hi(pos)=1.0_dp
      end if
      opt=.true.
      if(present(optimize))opt=optimize
      if(opt)call pattern_search(bats_objective,x,lo,hi,maxit=420,tol=1.0e-5_dp, fval=criterion)
      call unpack_params(x,alpha,beta,phi,g1,g2,arl,mal,lambda)
      if(fit_boxcox)then
         yy=boxcox(y,lambda)
      else
         yy=y
      end if
      model%has_trend=use_trend
      model%trigonometric=trig
      model%use_boxcox=fit_boxcox
      model%lambda=lambda
      model%alpha=alpha
      model%beta=beta
      model%phi=phi
      model%periods=periods
      model%periods_real=periods_real
      model%k=kvec
      model%gamma1=g1
      model%gamma2=g2
      model%ar=arl
      model%ma=mal
      model%p=p
      model%q=q
      if(trig)then
         model%w=tbats_make_w(use_trend,phi,kvec,arl,mal)
         model%g=tbats_make_g(alpha,use_trend,beta,kvec,g1,g2,p,q)
         model%F=tbats_make_f(alpha,use_trend,beta,phi,periods_real,kvec,g1,g2,arl,mal)
      else
         model%w=bats_make_w(use_trend,phi,periods,arl,mal)
         model%g=bats_make_g(alpha,use_trend,beta,g1,periods,p,q)
         model%F=bats_make_f(alpha,use_trend,beta,phi,periods,g1,arl,mal)
      end if
      call seed_states(yy,model%w,model%F,model%g,x0=x0,periods=periods,is_tbats=trig,p=p,q=q)
      model%seed_state=x0
      model%y=y
      call state_space_filter(yy,model%w,model%F,model%g,x0,model%fitted,model%residuals,states)
      model%state=states(:,size(y)+1)
      rss=sum(model%residuals**2)
      model%sigma2=rss/real(size(y),dp)
      criterion=real(size(y),dp)*log(max(rss,tiny(1.0_dp)))
      if(fit_boxcox)criterion=criterion-2.0_dp*(lambda-1.0_dp)*sum(log(y))
      model%loglik=-0.5_dp*criterion
      kcount=npar+size(x0)
      model%aic=criterion+2.0_dp*real(kcount,dp)
      model%bic=criterion+log(real(size(y),dp))*real(kcount,dp)
      if(size(y)>kcount+1)then
         model%aicc=model%aic+2.0_dp*real(kcount*(kcount+1),dp)/real(size(y)-kcount-1,dp)
      else
         model%aicc=huge(1.0_dp)
      end if
      if(fit_boxcox)model%fitted=inv_boxcox(model%fitted,lambda)
   end function fit_common

   function bats_fit(y,periods,use_trend,damped,boxcox_lambda,ar,ma,optimize) result(model)
      real(dp),intent(in)::y(:)
      integer,intent(in),optional::periods(:)
      logical,intent(in),optional::use_trend,damped,optimize
      real(dp),intent(in),optional::boxcox_lambda,ar(:),ma(:)
      type(bats_model)::model
      integer,allocatable::per(:),kv(:)
      real(dp),allocatable::pr(:)
      logical::tr,dm
      if(present(periods))then
         per=periods
      else
         allocate(per(0))
      end if
      allocate(kv(size(per)))
      kv=0
      pr=real(per,dp)
      tr=.true.
      if(present(use_trend))tr=use_trend
      dm=.false.
      if(present(damped))dm=damped
      model=fit_common(y,per,pr,kv,.false.,tr,dm,boxcox_lambda,ar,ma,optimize)
   end function bats_fit

   function tbats_fit(y,periods,kvec,use_trend,damped,boxcox_lambda,ar,ma,optimize) result(model)
      real(dp),intent(in)::y(:)
      integer,intent(in)::periods(:),kvec(:)
      logical,intent(in),optional::use_trend,damped,optimize
      real(dp),intent(in),optional::boxcox_lambda,ar(:),ma(:)
      type(bats_model)::model
      logical::tr,dm
      if(size(periods)/=size(kvec))error stop 'tbats_fit: periods/kvec mismatch'
      tr=.true.
      if(present(use_trend))tr=use_trend
      dm=.false.
      if(present(damped))dm=damped
      model=fit_common(y,periods,real(periods,dp),kvec,.true.,tr,dm,boxcox_lambda,ar,ma,optimize)
   end function tbats_fit

   function tbats_fit_real(y,periods,kvec,use_trend,damped,boxcox_lambda,ar,ma,optimize) result(model)
      real(dp),intent(in)::y(:),periods(:)
      integer,intent(in)::kvec(:)
      logical,intent(in),optional::use_trend,damped,optimize
      real(dp),intent(in),optional::boxcox_lambda,ar(:),ma(:)
      type(bats_model)::model
      integer,allocatable::iper(:)
      logical::tr,dm
      if(size(periods)/=size(kvec))error stop 'tbats_fit_real: periods/kvec mismatch'
      iper=nint(periods)
      tr=.true.
      if(present(use_trend))tr=use_trend
      dm=.false.
      if(present(damped))dm=damped
      model=fit_common(y,iper,periods,kvec,.true.,tr,dm,boxcox_lambda,ar,ma,optimize)
   end function tbats_fit_real

   subroutine select_arma_order(e,max_p,max_q,pbest,qbest)
      real(dp),intent(in)::e(:)
      integer,intent(in)::max_p,max_q
      integer,intent(out)::pbest,qbest
      type(arima_model)::a
      real(dp)::best
      integer::p,q
      best=huge(1.0_dp)
      pbest=0
      qbest=0
      do p=0,max_p
         do q=0,max_q
            if(p+q==0)cycle
            a=arima_fit(e,p,0,q,include_mean=.false.,optimize=.true.,method='ml')
            if(a%aicc<best)then
               best=a%aicc
               pbest=p
               qbest=q
            end if
         end do
      end do
   end subroutine select_arma_order

   function bats_auto(y,periods,use_arma,consider_boxcox,consider_trend,consider_damped,max_p,max_q, &
      boxcox_choice,trend_choice,damped_choice) result(best)
      real(dp),intent(in)::y(:)
      integer,intent(in),optional::periods(:),max_p,max_q,boxcox_choice,trend_choice,damped_choice
      logical,intent(in),optional::use_arma,consider_boxcox,consider_trend,consider_damped
      type(bats_model)::best,cand,first
      integer,allocatable::per(:),empty(:)
      real(dp),allocatable::ar0(:),ma0(:)
      integer::ib,it,id,pb,qb,mp,mq,season_case,bc,trc,dc,ib0,ib1,it0,it1,id0,id1
      logical::ua,cb,ct,cd,tr,dm
      real(dp)::lam
      if(present(periods))then
         per=periods
      else
         allocate(per(0))
      end if
      allocate(empty(0))
      ua=.true.
      if(present(use_arma))ua=use_arma
      cb=.true.
      if(present(consider_boxcox))cb=consider_boxcox
      ct=.true.
      if(present(consider_trend))ct=consider_trend
      cd=.true.
      if(present(consider_damped))cd=consider_damped
      mp=2
      if(present(max_p))mp=max_p
      mq=2
      if(present(max_q))mq=max_q
      bc=0
      if(present(boxcox_choice))bc=max(-1,min(1,boxcox_choice))
      trc=0
      if(present(trend_choice))trc=max(-1,min(1,trend_choice))
      dc=0
      if(present(damped_choice))dc=max(-1,min(1,damped_choice))
      if(bc==1 .and. any(y<=0.0_dp))error stop 'bats_auto: Box-Cox requires positive data'
      ib0=0
      ib1=merge(1,0,cb.and.all(y>0.0_dp))
      if(bc<0)ib1=0
      if(bc>0)then
      ib0=1
      ib1=1
      end if
      it0=0
      it1=merge(1,0,ct)
      if(trc<0)it1=0
      if(trc>0)then
      it0=1
      it1=1
      end if
      best%aic=huge(1.0_dp)
      do ib=ib0,ib1
         if(ib==1)lam=boxcox_lambda_guerrero(y,max(1,merge(maxval(per),1,size(per)>0)),0.0_dp,1.0_dp)
         do it=it0,it1
            tr=(it==1)
            id0=0
            id1=merge(1,0,cd.and.tr)
            if(dc<0)id1=0
            if(dc>0 .and. tr)then
            id0=1
            id1=1
            end if
            if(dc>0 .and. .not.tr)cycle
            do id=id0,id1
               dm=(id==1)
               do season_case=0,merge(1,0,size(per)>0)
                  if(season_case==0)then
                     if(ib==1)then
                        first=bats_fit(y,empty,tr,dm,lam,optimize=.true.)
                     else
                        first=bats_fit(y,empty,tr,dm,optimize=.true.)
                     end if
                  else
                     if(ib==1)then
                        first=bats_fit(y,per,tr,dm,lam,optimize=.true.)
                     else
                        first=bats_fit(y,per,tr,dm,optimize=.true.)
                     end if
                  end if
                  cand=first
                  if(ua)then
                     call select_arma_order(first%residuals,mp,mq,pb,qb)
                     if(pb+qb>0)then
                        allocate(ar0(pb),ma0(qb))
                        ar0=0.0_dp
                        ma0=0.0_dp
                        if(season_case==0)then
                           if(ib==1)then
                              cand=bats_fit(y,empty,tr,dm,lam,ar0,ma0,.true.)
                           else
                              cand=bats_fit(y,empty,tr,dm,ar=ar0,ma=ma0,optimize=.true.)
                           end if
                        else
                           if(ib==1)then
                              cand=bats_fit(y,per,tr,dm,lam,ar0,ma0,.true.)
                           else
                              cand=bats_fit(y,per,tr,dm,ar=ar0,ma=ma0,optimize=.true.)
                           end if
                        end if
                        if(cand%aic>first%aic)cand=first
                        deallocate(ar0,ma0)
                     end if
                  end if
                  if(cand%aic<best%aic)best=cand
               end do
            end do
         end do
      end do
   end function bats_auto

   function tbats_auto(y,periods,use_arma,consider_boxcox,consider_trend,consider_damped,max_p,max_q, &
      boxcox_choice,trend_choice,damped_choice) result(best)
      real(dp),intent(in)::y(:)
      integer,intent(in)::periods(:)
      logical,intent(in),optional::use_arma,consider_boxcox,consider_trend,consider_damped
      integer,intent(in),optional::max_p,max_q,boxcox_choice,trend_choice,damped_choice
      type(bats_model)::best,cand,work,trial
      integer,allocatable::kv(:)
      real(dp),allocatable::ar0(:),ma0(:)
      integer::i,k,maxk,ib,it,id,pb,qb,mp,mq,mref,bc,trc,dc,ib0,ib1,it0,it1,id0,id1
      logical::ua,cb,ct,cd,tr,dm,positive,improved
      real(dp)::lam

      ua=.true.
      if(present(use_arma))ua=use_arma
      cb=.true.
      if(present(consider_boxcox))cb=consider_boxcox
      ct=.true.
      if(present(consider_trend))ct=consider_trend
      cd=.true.
      if(present(consider_damped))cd=consider_damped
      mp=2
      if(present(max_p))mp=max(0,max_p)
      mq=2
      if(present(max_q))mq=max(0,max_q)
      positive=all(y>0.0_dp)
      bc=0
      if(present(boxcox_choice))bc=max(-1,min(1,boxcox_choice))
      trc=0
      if(present(trend_choice))trc=max(-1,min(1,trend_choice))
      dc=0
      if(present(damped_choice))dc=max(-1,min(1,damped_choice))
      if(bc==1 .and. .not.positive)error stop 'tbats_auto: Box-Cox requires positive data'
      best=bats_auto(y,periods,ua,cb,ct,cd,mp,mq,bc,trc,dc)
      if(size(periods)==0)return
      mref=max(2,maxval(periods))
      ib0=0
      ib1=merge(1,0,cb.and.positive)
      if(bc<0)ib1=0
      if(bc>0)then
      ib0=1
      ib1=1
      end if
      it0=0
      it1=merge(1,0,ct)
      if(trc<0)it1=0
      if(trc>0)then
      it0=1
      it1=1
      end if

      do ib=ib0,ib1
         lam=1.0_dp
         if(ib==1)lam=boxcox_lambda_guerrero(y,mref,0.0_dp,1.0_dp)
         do it=it0,it1
            tr=(it==1)
            id0=0
            id1=merge(1,0,cd.and.tr)
            if(dc<0)id1=0
            if(dc>0 .and. tr)then
            id0=1
            id1=1
            end if
            if(dc>0 .and. .not.tr)cycle
            do id=id0,id1
               dm=(id==1)
               allocate(kv(size(periods)))
               kv=1
               if(ib==1)then
                  work=tbats_fit(y,periods,kv,tr,dm,lam,optimize=.true.)
               else
                  work=tbats_fit(y,periods,kv,tr,dm,optimize=.true.)
               end if

               ! Hyndman et al.'s TBATS search expands Fourier harmonics one
               ! seasonal period at a time. Stop when the next harmonic no
               ! longer improves AIC; this avoids an exponential grid while
               ! retaining the upstream nested-model search structure.
               do i=1,size(periods)
                  maxk=max(1,(periods(i)-1)/2)
                  k=kv(i)+1
                  improved=.true.
                  do while(k<=maxk .and. improved)
                     kv(i)=k
                     if(ib==1)then
                        trial=tbats_fit(y,periods,kv,tr,dm,lam,optimize=.true.)
                     else
                        trial=tbats_fit(y,periods,kv,tr,dm,optimize=.true.)
                     end if
                     if(trial%aic+1.0e-9_dp<work%aic)then
                        work=trial
                        k=k+1
                     else
                        kv(i)=k-1
                        improved=.false.
                     end if
                  end do
                  kv=work%k
               end do

               cand=work
               if(ua .and. (mp+mq)>0)then
                  call select_arma_order(work%residuals,mp,mq,pb,qb)
                  if(pb+qb>0)then
                     allocate(ar0(pb),ma0(qb))
                     ar0=0.0_dp
                     ma0=0.0_dp
                     if(ib==1)then
                        trial=tbats_fit(y,periods,work%k,tr,dm,lam,ar0,ma0,.true.)
                     else
                        trial=tbats_fit(y,periods,work%k,tr,dm,ar=ar0,ma=ma0,optimize=.true.)
                     end if
                     if(trial%aic<cand%aic)cand=trial
                     deallocate(ar0,ma0)
                  end if
               end if
               if(cand%aic<best%aic)best=cand
               deallocate(kv)
            end do
         end do
      end do
   end function tbats_auto

   function bats_refit(y,model) result(refit)
      real(dp),intent(in)::y(:)
      type(bats_model),intent(in)::model
      type(bats_model)::refit
      real(dp),allocatable::yy(:),x0(:),states(:,:)
      real(dp)::rss
      refit=model
      if(model%use_boxcox)then
         if(any(y<=0.0_dp))error stop 'bats_refit: Box-Cox model requires positive data'
         yy=boxcox(y,model%lambda)
      else
         yy=y
      end if
      if(allocated(model%seed_state))then
         x0=model%seed_state
      else
         call seed_states(yy,model%w,model%F,model%g,model%periods,model%trigonometric,model%p,model%q,x0)
      end if
      call state_space_filter(yy,model%w,model%F,model%g,x0,refit%fitted,refit%residuals,states)
      refit%seed_state=x0
      refit%state=states(:,size(y)+1)
      refit%y=y
      rss=sum(refit%residuals**2)
      refit%sigma2=rss/real(max(1,size(y)),dp)
      if(model%use_boxcox)refit%fitted=inv_boxcox(refit%fitted,model%lambda)
   end function bats_refit

   function tbats_refit(y,model) result(refit)
      real(dp),intent(in)::y(:)
      type(bats_model),intent(in)::model
      type(bats_model)::refit
      if(.not.model%trigonometric)error stop 'tbats_refit: supplied model is not TBATS'
      refit=bats_refit(y,model)
   end function tbats_refit

   function bats_forecast(model,h,levels) result(fc)
      type(bats_model),intent(in)::model
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      real(dp),allocatable::lev(:),vm(:),frun(:,:),tmp(:,:),lower_t(:,:),upper_t(:,:)
      real(dp)::cj,z
      integer::j,k,nstate
      fc%mean=state_space_forecast(model%w,model%F,model%state,h)
      allocate(fc%se(h),vm(h))
      vm=1.0_dp
      nstate=size(model%w)
      if(h>1)then
         allocate(frun(nstate,nstate))
         frun=0.0_dp
         do j=1,nstate
            frun(j,j)=1.0_dp
         end do
         do j=1,h-1
            cj=dot_product(model%w,matmul(frun,model%g))
            vm(j+1)=vm(j)+cj*cj
            frun=matmul(frun,model%F)
         end do
      end if
      fc%se=sqrt(max(model%sigma2,0.0_dp)*vm)
      if(present(levels))then
         lev=levels
      else
         lev=[80.0_dp,95.0_dp]
      end if
      fc%level=lev
      allocate(lower_t(h,size(lev)),upper_t(h,size(lev)))
      do k=1,size(lev)
         z=normal_quantile(0.5_dp+lev(k)/200.0_dp)
         lower_t(:,k)=fc%mean-z*fc%se
         upper_t(:,k)=fc%mean+z*fc%se
      end do
      if(model%use_boxcox)then
         fc%mean=inv_boxcox(fc%mean,model%lambda,biasadj=.true.,fvar=fc%se**2)
         allocate(fc%lower(h,size(lev)),fc%upper(h,size(lev)))
         do k=1,size(lev)
            fc%lower(:,k)=inv_boxcox(lower_t(:,k),model%lambda)
            fc%upper(:,k)=inv_boxcox(upper_t(:,k),model%lambda)
         end do
      else
         fc%lower=lower_t
         fc%upper=upper_t
      end if
   end function bats_forecast
end module forecast_bats_fit
