module rfast_selection
   use rfast_special, only : dp, pi
   use rfast_arrays, only : standardise_cols, standardise_vector, variance_r
   use rfast_regression, only : regression_result, lmfit, glm_logistic, glm_poisson
   use rfast_regression_v02, only : quasipoisson_regression, proportion_regression, gamma_regression, &
                                    multinomial_regression, multinomial_result
   use rfast_regression_v03, only : normlog_regression, weibull_regression, weibull_regression_result
   use rfast_extra_mle, only : weibull_mle
   use rfast_mle, only : mle_result
   use rfast_linalg, only : logdet_spd
   implicit none
   private

   integer, parameter, public :: OMP_BIC=1, OMP_SSE=2, OMP_ADJR2=3
   integer, parameter, public :: OMP_LOGISTIC=1, OMP_POISSON=2, OMP_QUASIPOISSON=3
   integer, parameter, public :: OMP_QUASIBINOMIAL=4, OMP_NORMLOG=5, OMP_GAMMA=6, OMP_WEIBULL=7, OMP_MV=8, OMP_MULTINOMIAL_FAMILY=9
   type, public :: selection_result
      integer, allocatable :: selected(:)
      real(dp), allocatable :: criterion(:)
      integer :: steps = 0
      integer :: status = 0
   end type selection_result

   public :: ompr, bic_corfsreg, omp_glm, omp_multivariate, omp_multinomial

contains

   function ompr(y,x,method,tol,standardize) result(res)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(in),optional::method
      real(dp),intent(in),optional::tol
      logical,intent(in),optional::standardize
      type(selection_result)::res
      integer::n,p,meth,maxk,k,j,sel
      real(dp)::thr,prev,cur,down,r2,con,best
      logical::st
      real(dp),allocatable::xx(:,:),yy(:),z(:,:),scores(:),crit(:),resid(:)
      integer,allocatable::picked(:)
      logical,allocatable::used(:)
      type(regression_result)::fit
      n=size(y);p=size(x,2);if(size(x,1)/=n)then;res%status=1;return;end if
      meth=OMP_BIC
      if(present(method))meth=method
      thr=2.0_dp
      if(present(tol))thr=tol
      st=.true.
      if(present(standardize))st=standardize
      allocate(xx(n,p),yy(n),scores(p),used(p));xx=x;yy=y;if(st)then;xx=standardise_cols(xx);yy=standardise_vector(yy);end if
      maxk=min(p,max(1,n-2));allocate(picked(maxk),crit(maxk+1));used=.false.;down=sum((yy-sum(yy)/real(n,dp))**2)
      con=real(n,dp)*(log(2.0_dp*pi)+1.0_dp)
      select case(meth)
      case(OMP_BIC);crit(1)=real(n,dp)*log(max(tiny(1.0_dp),down/real(n,dp)))+2.0_dp*log(real(n,dp))+con
      case(OMP_SSE);crit(1)=down
      case(OMP_ADJR2);crit(1)=0.0_dp
      case default;res%status=2;return
      end select
      resid=yy
      do k=1,maxk
         scores=0.0_dp
         do j=1,p;if(.not.used(j))scores(j)=abs(dot_product(xx(:,j),resid));end do
         best=-1.0_dp;sel=0
         do j=1,p;if(.not.used(j).and.scores(j)>best)then;best=scores(j);sel=j;end if;end do
         if(sel==0)exit;used(sel)=.true.;picked(k)=sel;allocate(z(n,k));do j=1,k;z(:,j)=xx(:,picked(j));end do
         fit=lmfit(z,yy,.false.);resid=fit%residuals
         select case(meth)
         case(OMP_BIC)
            cur=real(n,dp)*log(max(tiny(1.0_dp),sum(resid**2)/real(n,dp)))+real(k+2,dp)*log(real(n,dp))+con
            crit(k+1)=cur;prev=crit(k);if(prev-cur<=thr.and.k>1)then;deallocate(z);exit;end if
         case(OMP_SSE)
            cur=sum(resid**2);crit(k+1)=cur;prev=crit(k)
            if(k>1.and.(prev-cur)/max(tiny(1.0_dp),prev)<=thr)then;deallocate(z);exit;end if
         case(OMP_ADJR2)
            r2=1.0_dp-sum(resid**2)/max(tiny(1.0_dp),down)
            cur=1.0_dp-(1.0_dp-r2)*real(n-1,dp)/real(max(1,n-k-1),dp);crit(k+1)=cur;prev=crit(k)
            if(k>1.and.cur-prev<=thr)then;deallocate(z);exit;end if
         end select
         deallocate(z)
      end do
      res%steps=k; if(k>maxk)res%steps=maxk
      ! If the terminating step failed the threshold, keep the upstream-style tested model in the trace.
      allocate(res%selected(res%steps),res%criterion(res%steps+1))
      res%selected=picked(1:res%steps)
      res%criterion=crit(1:res%steps+1)
   end function ompr

   function bic_corfsreg(y,x,tol) result(res)
      real(dp),intent(in)::y(:),x(:,:)
      real(dp),intent(in),optional::tol
      type(selection_result)::res
      integer::n,p,maxk,k,j,sel
      real(dp)::thr,con,logn,best,score,my,mx,sy,sx,corr,cur
      real(dp),allocatable::xx(:,:),yy(:),z(:,:),ery(:),erx(:),crit(:)
      integer,allocatable::picked(:)
      logical,allocatable::used(:)
      type(regression_result)::fy,fx
      n=size(y);p=size(x,2);thr=2.0_dp;if(present(tol))thr=tol
      allocate(xx(n,p),yy(n),used(p));xx=x;yy=y;xx=standardise_cols(xx);yy=standardise_vector(yy);used=.false.
      maxk=min(p,max(1,n-20));allocate(picked(maxk),crit(maxk+1));con=real(n,dp)*(log(2.0_dp*pi)+1.0_dp);logn=log(real(n,dp))
      crit(1)=real(n,dp)*log(max(tiny(1.0_dp),sum(yy**2)/real(n,dp)))+2.0_dp*logn+con;ery=yy
      do k=1,maxk
         best=-1.0_dp;sel=0
         do j=1,p
            if(used(j))cycle
            if(k==1)then;erx=xx(:,j)
            else;fx=lmfit(z,xx(:,j),.false.);erx=fx%residuals;end if
            my=sum(ery)/real(n,dp);mx=sum(erx)/real(n,dp);sy=sqrt(sum((ery-my)**2));sx=sqrt(sum((erx-mx)**2))
            if(sy>0.0_dp.and.sx>0.0_dp)then;corr=dot_product(ery-my,erx-mx)/(sy*sx);else;corr=0.0_dp;end if
            score=abs(corr);if(score>best)then;best=score;sel=j;end if
         end do
         if(sel==0)exit;used(sel)=.true.;picked(k)=sel
         if(allocated(z))deallocate(z);allocate(z(n,k));do j=1,k;z(:,j)=xx(:,picked(j));end do
         fy=lmfit(z,yy,.false.);ery=fy%residuals
         cur=real(n,dp)*log(max(tiny(1.0_dp),sum(ery**2)/real(n,dp)))+real(k+2,dp)*logn+con;crit(k+1)=cur
         if(k>1.and.crit(k)-cur<=thr)exit
      end do
      res%steps=min(k,maxk);allocate(res%selected(res%steps),res%criterion(res%steps+1))
      res%selected=picked(1:res%steps);res%criterion=crit(1:res%steps+1)
   end function bic_corfsreg

   function omp_glm(y,x,family,tol,standardize) result(res)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(in)::family
      real(dp),intent(in),optional::tol
      logical,intent(in),optional::standardize
      type(selection_result)::res
      integer::n,p,maxk,k,j,sel,i
      real(dp)::thr,best,score,m,prev,cur,scale
      logical::st
      real(dp),allocatable::xx(:,:),design(:,:),resid(:),crit(:),phi(:)
      integer,allocatable::picked(:)
      logical,allocatable::used(:)
      type(regression_result)::fit
      type(weibull_regression_result)::wfit
      type(mle_result)::wini
      n=size(y);p=size(x,2);if(size(x,1)/=n)then;res%status=1;return;end if
      thr=log(real(n,dp))+3.841458820694124_dp;if(present(tol))thr=tol
      st=.true.;if(present(standardize))st=standardize
      allocate(xx(n,p),used(p));xx=x;if(st)xx=standardise_cols(xx);used=.false.
      maxk=min(p,max(1,n-3));allocate(picked(maxk),crit(maxk+1),phi(maxk+1),resid(n));phi=1.0_dp
      m=sum(y)/real(n,dp)
      select case(family)
      case(OMP_LOGISTIC,OMP_QUASIBINOMIAL)
         if(m<=0.0_dp.or.m>=1.0_dp)then;res%status=2;return;end if
         crit(1)=0.0_dp
         do i=1,n
            if(y(i)>0.0_dp)crit(1)=crit(1)+2.0_dp*y(i)*log(y(i)/m)
            if(y(i)<1.0_dp)crit(1)=crit(1)+2.0_dp*(1.0_dp-y(i))*log((1.0_dp-y(i))/(1.0_dp-m))
         end do
         resid=y-m
      case(OMP_POISSON,OMP_QUASIPOISSON)
         if(m<=0.0_dp)then;res%status=2;return;end if
         crit(1)=0.0_dp
         do i=1,n
            if(y(i)>0.0_dp)then;crit(1)=crit(1)+2.0_dp*(y(i)*log(y(i)/m)-(y(i)-m))
            else;crit(1)=crit(1)+2.0_dp*m;end if
         end do
         resid=y-m
      case(OMP_NORMLOG)
         crit(1)=sum((y-m)**2);resid=y-m
      case(OMP_GAMMA)
         if(any(y<=0.0_dp))then;res%status=2;return;end if
         crit(1)=2.0_dp*sum(y/m-1.0_dp-log(y/m));resid=y-m
      case(OMP_WEIBULL)
         if(any(y<=0.0_dp))then;res%status=2;return;end if
         wini=weibull_mle(y);if(wini%status/=0)then;res%status=wini%status;return;end if
         crit(1)=2.0_dp*wini%loglik;resid=y-m
      case default;res%status=3;return
      end select
      do k=1,maxk
         best=-1.0_dp;sel=0
         do j=1,p
            if(.not.used(j))then
               score=abs(dot_product(xx(:,j),resid))
               if(score>best)then;best=score;sel=j;end if
            end if
         end do
         if(sel==0)exit;used(sel)=.true.;picked(k)=sel;allocate(design(n,k))
         do j=1,k;design(:,j)=xx(:,picked(j));end do
         select case(family)
         case(OMP_LOGISTIC)
            call fit_with_intercept_logistic(design,y,fit)
            if(fit%status==0)then;cur=fit%deviance;resid=fit%residuals;end if
         case(OMP_POISSON)
            call fit_with_intercept_poisson(design,y,fit)
            if(fit%status==0)then;cur=fit%deviance;resid=fit%residuals;end if
         case(OMP_QUASIPOISSON)
            fit=quasipoisson_regression(y,design)
            if(fit%status==0)then;cur=fit%deviance;phi(k+1)=fit%dispersion;resid=fit%residuals;end if
         case(OMP_QUASIBINOMIAL)
            fit=proportion_regression(y,design,.true.)
            if(fit%status==0)then
               cur=binomial_deviance(y,fit%fitted);phi(k+1)=fit%dispersion;resid=fit%residuals
            end if
         case(OMP_NORMLOG)
            fit=normlog_regression(y,design)
            if(fit%status==0)then
               cur=fit%deviance;phi(k+1)=fit%deviance/real(max(1,n-size(fit%beta)),dp);resid=fit%residuals
            end if
         case(OMP_GAMMA)
            fit=gamma_regression(y,design)
            if(fit%status==0)then;cur=fit%deviance;phi(k+1)=fit%dispersion;resid=fit%residuals;end if
         case(OMP_WEIBULL)
            wfit=weibull_regression(y,design)
            if(wfit%status==0)then;cur=2.0_dp*wfit%loglik;resid=y-wfit%fitted;end if
         end select
         if((family/=OMP_WEIBULL.and.fit%status/=0).or.(family==OMP_WEIBULL.and.wfit%status/=0))then
            res%status=4;deallocate(design);exit
         end if
         crit(k+1)=cur;prev=crit(k);deallocate(design)
         select case(family)
         case(OMP_QUASIPOISSON,OMP_QUASIBINOMIAL,OMP_NORMLOG,OMP_GAMMA)
            scale=max(tiny(1.0_dp),phi(k+1));if((prev-cur)/scale<=thr)exit
         case(OMP_WEIBULL)
            if(cur-prev<=thr)exit
         case default
            if(prev-cur<=thr)exit
         end select
      end do
      res%steps=min(k,maxk);allocate(res%selected(res%steps),res%criterion(res%steps+1))
      res%selected=picked(1:res%steps);res%criterion=crit(1:res%steps+1)
   end function omp_glm

   function omp_multivariate(y,x,tol,standardize) result(res)
      real(dp),intent(in)::y(:,:),x(:,:)
      real(dp),intent(in),optional::tol
      logical,intent(in),optional::standardize
      type(selection_result)::res
      real(dp),allocatable::xx(:,:),resid(:,:),design(:,:),crit(:),scores(:),cov(:,:)
      integer,allocatable::picked(:)
      logical,allocatable::used(:)
      real(dp)::thr,best,cur,prev,con,ld
      integer::n,p,d,maxk,k,j,c,sel,info
      type(regression_result)::fit
      n=size(y,1);d=size(y,2);p=size(x,2);if(size(x,1)/=n)then;res%status=1;return;end if
      thr=log(real(n,dp))+3.841458820694124_dp;if(present(tol))thr=tol
      allocate(xx(n,p));xx=x
      if(present(standardize))then
         if(standardize)xx=standardise_cols(xx)
      else
         xx=standardise_cols(xx)
      end if
      allocate(resid(n,d));resid=y-spread(sum(y,dim=1)/real(n,dp),1,n);allocate(cov(d,d))
      cov=matmul(transpose(resid),resid)/real(max(1,n-1),dp);ld=logdet_spd(cov,info)
      if(info/=0)then;res%status=2;return;end if
      con=real(n*d,dp)*(log(2.0_dp*pi)+1.0_dp);maxk=min(p,max(1,n-d-2))
      allocate(picked(maxk),crit(maxk+1),scores(p),used(p));crit(1)=con+real(n,dp)*ld;used=.false.
      do k=1,maxk
         scores=0.0_dp
         do j=1,p
            if(.not.used(j))then
               do c=1,d;scores(j)=scores(j)+dot_product(xx(:,j),resid(:,c))**2;end do
            end if
         end do
         best=-1.0_dp;sel=0
         do j=1,p;if(.not.used(j).and.scores(j)>best)then;best=scores(j);sel=j;end if;end do
         if(sel==0)exit;used(sel)=.true.;picked(k)=sel;allocate(design(n,k+1));design(:,1)=1.0_dp
         do j=1,k;design(:,j+1)=xx(:,picked(j));end do
         do c=1,d
            fit=lmfit(design,y(:,c),.false.);if(fit%status/=0)then;res%status=fit%status;exit;end if
            resid(:,c)=fit%residuals
         end do
         if(res%status/=0)then;deallocate(design);exit;end if
         cov=matmul(transpose(resid),resid)/real(max(1,n-1),dp);ld=logdet_spd(cov,info)
         if(info/=0)then;res%status=2;deallocate(design);exit;end if
         cur=con+real(n,dp)*ld;crit(k+1)=cur;prev=crit(k);deallocate(design)
         if(prev-cur<=thr)exit
      end do
      res%steps=min(k,maxk);allocate(res%selected(res%steps),res%criterion(res%steps+1))
      res%selected=picked(1:res%steps);res%criterion=crit(1:res%steps+1)
   end function omp_multivariate

   function omp_multinomial(y,x,tol,standardize) result(res)
      integer,intent(in)::y(:)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(in),optional::tol
      logical,intent(in),optional::standardize
      type(selection_result)::res
      real(dp),allocatable::xx(:,:),resid(:,:),design(:,:),crit(:),scores(:),eta(:,:),prob(:,:)
      integer,allocatable::picked(:),counts(:)
      logical,allocatable::used(:)
      real(dp)::thr,best,cur,prev,den,shift,null_ll
      integer::n,p,nclass,d,maxk,k,j,c,i,sel
      type(multinomial_result)::fit
      n=size(y);p=size(x,2);if(size(x,1)/=n.or.n==0.or.minval(y)<1)then;res%status=1;return;end if
      nclass=maxval(y);d=nclass-1;if(d<1)then;res%status=2;return;end if
      thr=log(real(n,dp))+3.841458820694124_dp;if(present(tol))thr=tol
      allocate(xx(n,p));xx=x
      if(present(standardize))then
         if(standardize)xx=standardise_cols(xx)
      else
         xx=standardise_cols(xx)
      end if
      allocate(counts(nclass));do c=1,nclass;counts(c)=count(y==c);end do
      if(any(counts==0))then;res%status=3;return;end if
      null_ll=0.0_dp;do c=1,nclass;null_ll=null_ll+real(counts(c),dp)*log(real(counts(c),dp)/real(n,dp));end do
      allocate(resid(n,d));do c=1,d;resid(:,c)=merge(1.0_dp,0.0_dp,y==c+1)-real(counts(c+1),dp)/real(n,dp);end do
      maxk=min(p,max(1,n-d-2));allocate(picked(maxk),crit(maxk+1),scores(p),used(p));crit(1)=-2.0_dp*null_ll;used=.false.
      do k=1,maxk
         scores=0.0_dp
         do j=1,p
            if(.not.used(j))then
               do c=1,d;scores(j)=scores(j)+dot_product(xx(:,j),resid(:,c))**2;end do
            end if
         end do
         best=-1.0_dp;sel=0
         do j=1,p;if(.not.used(j).and.scores(j)>best)then;best=scores(j);sel=j;end if;end do
         if(sel==0)exit;used(sel)=.true.;picked(k)=sel;allocate(design(n,k));do j=1,k;design(:,j)=xx(:,picked(j));end do
         fit=multinomial_regression(y,design);if(fit%status/=0)then;res%status=fit%status;deallocate(design);exit;end if
         cur=-2.0_dp*fit%loglik;crit(k+1)=cur;prev=crit(k)
         allocate(eta(n,d),prob(n,d));eta=matmul(augment_intercept(design),fit%beta)
         do i=1,n
            shift=max(0.0_dp,maxval(eta(i,:)));den=exp(-shift)+sum(exp(eta(i,:)-shift));prob(i,:)=exp(eta(i,:)-shift)/den
         end do
         do c=1,d;resid(:,c)=merge(1.0_dp,0.0_dp,y==c+1)-prob(:,c);end do
         deallocate(eta,prob,design)
         if(prev-cur<=thr)exit
      end do
      res%steps=min(k,maxk);allocate(res%selected(res%steps),res%criterion(res%steps+1))
      res%selected=picked(1:res%steps);res%criterion=crit(1:res%steps+1)
   end function omp_multinomial

   function augment_intercept(x) result(xx)
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable::xx(:,:)
      allocate(xx(size(x,1),size(x,2)+1));xx(:,1)=1.0_dp;xx(:,2:)=x
   end function augment_intercept

   subroutine fit_with_intercept_logistic(x,y,fit)
      real(dp),intent(in)::x(:,:),y(:);type(regression_result),intent(out)::fit
      real(dp),allocatable::xx(:,:)
      allocate(xx(size(x,1),size(x,2)+1));xx(:,1)=1.0_dp;xx(:,2:)=x;fit=glm_logistic(xx,y)
   end subroutine fit_with_intercept_logistic

   subroutine fit_with_intercept_poisson(x,y,fit)
      real(dp),intent(in)::x(:,:),y(:);type(regression_result),intent(out)::fit
      real(dp),allocatable::xx(:,:)
      allocate(xx(size(x,1),size(x,2)+1));xx(:,1)=1.0_dp;xx(:,2:)=x;fit=glm_poisson(xx,y)
   end subroutine fit_with_intercept_poisson

   real(dp) function binomial_deviance(y,p) result(v)
      real(dp),intent(in)::y(:),p(:);integer::i
      v=0.0_dp
      do i=1,size(y)
         if(y(i)>0.0_dp)v=v+2.0_dp*y(i)*log(y(i)/max(tiny(1.0_dp),p(i)))
         if(y(i)<1.0_dp)v=v+2.0_dp*(1.0_dp-y(i))*log((1.0_dp-y(i))/max(tiny(1.0_dp),1.0_dp-p(i)))
      end do
   end function binomial_deviance


end module rfast_selection
