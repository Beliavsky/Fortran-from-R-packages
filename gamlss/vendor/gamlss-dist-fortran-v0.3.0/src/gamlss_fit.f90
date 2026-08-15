! Generic maximum-likelihood regression layer for translated gamlss.dist families.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_fit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gamlss_kinds, only : dp
   use gamlss_optim, only : bfgs_minimize, numerical_hessian
   use gamlss_linalg, only : invert_matrix
   use gamlss_continuous
   use gamlss_discrete
   use gamlss_boxcox
   use gamlss_continuous_v02
   use gamlss_discrete_v02
   use gamlss_continuous_v03
   use gamlss_discrete_v03
   use gamlss_flexible_v03
   implicit none
   private
   integer, parameter, public :: GAMLSS_NO=1, GAMLSS_GA=2, GAMLSS_BE=3
   integer, parameter, public :: GAMLSS_NBI=4, GAMLSS_NBII=5, GAMLSS_ZIP=6
   integer, parameter, public :: GAMLSS_GG=7, GAMLSS_EGB2=8, GAMLSS_GB2=9
   integer, parameter, public :: GAMLSS_JSUO=10, GAMLSS_TF=11, GAMLSS_PIG=12
   integer, parameter, public :: GAMLSS_BEINF=13, GAMLSS_WEI=14, GAMLSS_LNO=15
   integer, parameter, public :: GAMLSS_BCCG=16, GAMLSS_BCT=17, GAMLSS_BCPE=18
   integer, parameter, public :: GAMLSS_GIG=19, GAMLSS_SHASHO=20, GAMLSS_SHASH=21
   integer, parameter, public :: GAMLSS_SIMPLEX=22, GAMLSS_SEP=23, GAMLSS_SEP1=24
   integer, parameter, public :: GAMLSS_SEP2=25, GAMLSS_ST1=26, GAMLSS_ST2=27
   integer, parameter, public :: GAMLSS_ST5=28, GAMLSS_NET=29, GAMLSS_GPO=30
   integer, parameter, public :: GAMLSS_DPO=31, GAMLSS_DEL=32, GAMLSS_SI=33
   integer, parameter, public :: GAMLSS_SICHEL=34, GAMLSS_YULE=35, GAMLSS_WARING=36
   integer, parameter, public :: GAMLSS_ZIPF=37
   integer, parameter, public :: GAMLSS_ST3=38, GAMLSS_ST4=39
   integer, parameter, public :: GAMLSS_SEP3=40, GAMLSS_SEP4=41
   integer, parameter, public :: GAMLSS_ST3C=42, GAMLSS_SN1=43, GAMLSS_SN2=44
   integer, parameter, public :: GAMLSS_SST=45, GAMLSS_GT=46, GAMLSS_EXGAUS=47
   integer, parameter, public :: GAMLSS_PARETO=48, GAMLSS_PARETO1=49, GAMLSS_PARETO2=50, GAMLSS_PARETO2O=51
   integer, parameter, public :: GAMLSS_PIG2=52, GAMLSS_ZIPIG=53, GAMLSS_ZAPIG=54
   integer, parameter, public :: GAMLSS_ZISICHEL=55, GAMLSS_ZASICHEL=56
   integer, parameter, public :: GAMLSS_ZIBNB=57, GAMLSS_ZABNB=58, GAMLSS_ZAZIPF=59
   integer, parameter, public :: GAMLSS_GAF=60, GAMLSS_NBF=61, GAMLSS_ZINBF=62

   type, public :: gamlss_fit_result_t
      real(dp), allocatable :: beta_mu(:), beta_sigma(:), beta_nu(:), beta_tau(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: fitted_mu(:), fitted_sigma(:), fitted_nu(:), fitted_tau(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: iterations_status = 0
      integer :: status = 0
      logical :: converged = .false.
   end type

   type :: fit_context_t
      real(dp),allocatable :: y(:),x1(:,:),x2(:,:),x3(:,:),x4(:,:),w(:)
      integer :: family=0,p1=0,p2=0,p3=0,p4=0,np=0,n=0
   end type fit_context_t
   type(fit_context_t),save :: fit_ctx

   public :: fit_gamlss, family_npar, map_parameters, inverse_link, family_logpdf, default_parameters
contains

   subroutine fit_gamlss(y,x_mu,family,result,x_sigma,x_nu,x_tau,weights,start,max_iter,tol)
      real(dp),intent(in)::y(:),x_mu(:,:)
      integer,intent(in)::family
      type(gamlss_fit_result_t),intent(out)::result
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:),start(:)
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tol
      integer::n,p1,p2,p3,p4,np,npar,status_inv,i
      real(dp),allocatable::theta(:),hess(:,:),cov(:,:),w(:)
      real(dp)::fval,tolerance

      n=size(y)
      p1=size(x_mu,2)
      npar=family_npar(family)
      p2=0
      p3=0
      p4=0
      if(npar>=2)then
         if(.not.present(x_sigma))then
         result%status=11
         return
         end if
         if(size(x_sigma,1)/=n)then
         result%status=12
         return
         end if
         p2=size(x_sigma,2)
      end if
      if(npar>=3)then
         if(.not.present(x_nu))then
         result%status=13
         return
         end if
         if(size(x_nu,1)/=n)then
         result%status=14
         return
         end if
         p3=size(x_nu,2)
      end if
      if(npar>=4)then
         if(.not.present(x_tau))then
         result%status=15
         return
         end if
         if(size(x_tau,1)/=n)then
         result%status=16
         return
         end if
         p4=size(x_tau,2)
      end if
      np=p1+p2+p3+p4
      if(np<=0.or.size(x_mu,1)/=n)then
      result%status=10
      return
      end if
      allocate(w(n))
      w=1.0_dp
      if(present(weights))then
      if(size(weights)/=n)then
      result%status=17
      return
      end if
      w=weights
      end if
      allocate(theta(np))
      theta=0.0_dp
      if(present(start))then
         if(size(start)/=np)then
         result%status=18
         return
         end if
         theta=start
      else
         call initialize_theta(theta,y,x_mu,family,p1,p2,p3,p4)
      end if
      tolerance=1e-7_dp
      if(present(tol))tolerance=tol
      call set_fit_context(y,x_mu,family,w,p1,p2,p3,p4,x_sigma,x_nu,x_tau)
      call bfgs_minimize(fit_objective,theta,fval,result%iterations_status,max_iter=max_iter,tol=tolerance)
      result%converged=(result%iterations_status==0)
      result%status=result%iterations_status
      result%loglik=-fval
      result%aic=2.0_dp*real(np,dp)-2.0_dp*result%loglik
      call unpack_result(theta,p1,p2,p3,p4,result)
      call fill_fitted(theta,p1,p2,p3,p4,result)
      allocate(hess(np,np))
      call numerical_hessian(fit_objective,theta,hess)
      call invert_matrix(hess,cov,status_inv)
      if(status_inv==0)then
      result%covariance=cov
      else
      allocate(result%covariance(0,0))
      end if
      call clear_fit_context()
   contains
      subroutine fill_fitted(par,pp1,pp2,pp3,pp4,res)
         real(dp),intent(in)::par(:)
         integer,intent(in)::pp1,pp2,pp3,pp4
         type(gamlss_fit_result_t),intent(inout)::res
         integer::ii,o2,o3,o4
         real(dp)::e1,e2,e3,e4,a,b,c,d
         o2=pp1
         o3=pp1+pp2
         o4=pp1+pp2+pp3
         allocate(res%fitted_mu(n))
         if(pp2>0)allocate(res%fitted_sigma(n))
         if(pp3>0)allocate(res%fitted_nu(n))
         if(pp4>0)allocate(res%fitted_tau(n))
         do ii=1,n
            e1=dot_product(x_mu(ii,:),par(1:pp1))
            e2=0
            e3=0
            e4=0
            if(pp2>0)e2=dot_product(x_sigma(ii,:),par(o2+1:o3))
            if(pp3>0)e3=dot_product(x_nu(ii,:),par(o3+1:o4))
            if(pp4>0)e4=dot_product(x_tau(ii,:),par(o4+1:np))
            call map_parameters(family,e1,e2,e3,e4,a,b,c,d)
            res%fitted_mu(ii)=a
            if(pp2>0)res%fitted_sigma(ii)=b
            if(pp3>0)res%fitted_nu(ii)=c
            if(pp4>0)res%fitted_tau(ii)=d
         end do
      end subroutine fill_fitted

   end subroutine

   subroutine set_fit_context(y,x1,family,w,p1,p2,p3,p4,x2,x3,x4)
      real(dp),intent(in)::y(:),x1(:,:),w(:)
      integer,intent(in)::family,p1,p2,p3,p4
      real(dp),intent(in),optional::x2(:,:),x3(:,:),x4(:,:)
      call clear_fit_context()
      fit_ctx%y=y;fit_ctx%x1=x1;fit_ctx%w=w;fit_ctx%family=family
      fit_ctx%n=size(y);fit_ctx%p1=p1;fit_ctx%p2=p2;fit_ctx%p3=p3;fit_ctx%p4=p4
      fit_ctx%np=p1+p2+p3+p4
      if(present(x2))then;fit_ctx%x2=x2;else;allocate(fit_ctx%x2(fit_ctx%n,0));end if
      if(present(x3))then;fit_ctx%x3=x3;else;allocate(fit_ctx%x3(fit_ctx%n,0));end if
      if(present(x4))then;fit_ctx%x4=x4;else;allocate(fit_ctx%x4(fit_ctx%n,0));end if
   end subroutine set_fit_context

   subroutine clear_fit_context()
      if(allocated(fit_ctx%y))deallocate(fit_ctx%y)
      if(allocated(fit_ctx%x1))deallocate(fit_ctx%x1)
      if(allocated(fit_ctx%x2))deallocate(fit_ctx%x2)
      if(allocated(fit_ctx%x3))deallocate(fit_ctx%x3)
      if(allocated(fit_ctx%x4))deallocate(fit_ctx%x4)
      if(allocated(fit_ctx%w))deallocate(fit_ctx%w)
      fit_ctx%n=0;fit_ctx%family=0
   end subroutine clear_fit_context

   real(dp) function fit_objective(par) result(nll)
      real(dp),intent(in)::par(:)
      real(dp)::eta1,eta2,eta3,eta4,a,b,c,d,lp
      integer::i,o2,o3,o4
      nll=0.0_dp;o2=fit_ctx%p1;o3=o2+fit_ctx%p2;o4=o3+fit_ctx%p3
      do i=1,fit_ctx%n
         eta1=dot_product(fit_ctx%x1(i,:),par(1:fit_ctx%p1))
         eta2=0.0_dp;eta3=0.0_dp;eta4=0.0_dp
         if(fit_ctx%p2>0)eta2=dot_product(fit_ctx%x2(i,:),par(o2+1:o3))
         if(fit_ctx%p3>0)eta3=dot_product(fit_ctx%x3(i,:),par(o3+1:o4))
         if(fit_ctx%p4>0)eta4=dot_product(fit_ctx%x4(i,:),par(o4+1:fit_ctx%np))
         call map_parameters(fit_ctx%family,eta1,eta2,eta3,eta4,a,b,c,d)
         lp=family_logpdf(fit_ctx%family,fit_ctx%y(i),a,b,c,d)
         if(.not.ieee_is_finite(lp))then
            nll=nll+1.0e12_dp+1.0e6_dp*sum(par*par)
         else
            nll=nll-fit_ctx%w(i)*lp
         end if
      end do
   end function fit_objective

   integer pure function family_npar(family) result(n)
      integer,intent(in)::family
      select case(family)
      case(GAMLSS_YULE,GAMLSS_ZIPF,GAMLSS_PARETO,GAMLSS_PARETO1)
      n=1
      case(GAMLSS_NO,GAMLSS_GA,GAMLSS_BE,GAMLSS_NBI,GAMLSS_NBII,GAMLSS_ZIP,GAMLSS_PIG,GAMLSS_WEI, &
           GAMLSS_SIMPLEX,GAMLSS_GPO,GAMLSS_DPO,GAMLSS_WARING,GAMLSS_PARETO2,GAMLSS_PARETO2O, &
           GAMLSS_PIG2,GAMLSS_ZAZIPF)
      n=2
      case(GAMLSS_GG,GAMLSS_TF,GAMLSS_LNO,GAMLSS_BCCG,GAMLSS_GIG,GAMLSS_DEL,GAMLSS_SI,GAMLSS_SICHEL, &
           GAMLSS_GAF,GAMLSS_NBF, &
           GAMLSS_SN1,GAMLSS_SN2,GAMLSS_EXGAUS,GAMLSS_ZIPIG,GAMLSS_ZAPIG)
      n=3
      case(GAMLSS_EGB2,GAMLSS_GB2,GAMLSS_BEINF,GAMLSS_JSUO,GAMLSS_BCT,GAMLSS_BCPE, &
           GAMLSS_SHASHO,GAMLSS_SHASH,GAMLSS_SEP,GAMLSS_SEP1,GAMLSS_SEP2, &
           GAMLSS_ST1,GAMLSS_ST2,GAMLSS_ST3,GAMLSS_ST4,GAMLSS_ST5, &
           GAMLSS_SEP3,GAMLSS_SEP4,GAMLSS_NET,GAMLSS_ST3C,GAMLSS_SST,GAMLSS_GT, &
           GAMLSS_ZISICHEL,GAMLSS_ZASICHEL,GAMLSS_ZIBNB,GAMLSS_ZABNB,GAMLSS_ZINBF)
      n=4
      case default
      n=0
      end select
   end function

   subroutine initialize_theta(theta,y,xmu,family,p1,p2,p3,p4)
      real(dp),intent(out)::theta(:)
      real(dp),intent(in)::y(:),xmu(:,:)
      integer,intent(in)::family,p1,p2,p3,p4
      real(dp)::m,s,pmu,psig,pnu,ptau
      integer::o2,o3,o4
      theta=0
      m=sum(y)/real(size(y),dp)
      s=sqrt(max(sum((y-m)**2)/real(max(1,size(y)-1),dp),1e-8_dp))
      call default_parameters(family,m,s,pmu,psig,pnu,ptau)
      o2=p1
      o3=p1+p2
      o4=p1+p2+p3
      if(p1>0.and.is_intercept(xmu(:,1)))theta(1)=inverse_link(family,1,pmu)
      if(p2>0)theta(o2+1)=inverse_link(family,2,psig)
      if(p3>0)theta(o3+1)=inverse_link(family,3,pnu)
      if(p4>0)theta(o4+1)=inverse_link(family,4,ptau)
   end subroutine

   logical pure function is_intercept(x) result(ok)
      real(dp),intent(in)::x(:)
      ok=maxval(abs(x-1.0_dp))<1e-12_dp
   end function

   subroutine default_parameters(family,m,s,a,b,c,d)
      integer,intent(in)::family
      real(dp),intent(in)::m,s
      real(dp),intent(out)::a,b,c,d
      a=m
      b=max(s,0.2_dp)
      c=1.0_dp
      d=1.0_dp
      select case(family)
      case(GAMLSS_GA,GAMLSS_WEI,GAMLSS_NBI,GAMLSS_NBII,GAMLSS_PIG)
      a=max(abs(m),0.5_dp)
      b=max(s/max(abs(m),.5_dp),.2_dp)
      case(GAMLSS_BE)
      a=min(.9_dp,max(.1_dp,m))
      b=.2_dp
      case(GAMLSS_ZIP)
      a=max(abs(m),.5_dp)
      b=.1_dp
      case(GAMLSS_GG)
      a=max(abs(m),.5_dp)
      b=.5_dp
      c=1.0_dp
      case(GAMLSS_TF)
      a=m
      b=max(s,.2_dp)
      c=8.0_dp
      case(GAMLSS_JSUO)
      a=m
      b=max(s,.2_dp)
      c=0.0_dp
      d=1.0_dp
      case(GAMLSS_EGB2)
      a=m
      b=max(s,.2_dp)
      c=1.0_dp
      d=1.0_dp
      case(GAMLSS_GB2)
      a=max(abs(m),.5_dp)
      b=1.0_dp
      c=1.0_dp
      d=1.0_dp
      case(GAMLSS_BEINF)
      a=min(.9_dp,max(.1_dp,m))
      b=.2_dp
      c=.05_dp
      d=.05_dp
      case(GAMLSS_LNO)
      a=log(max(abs(m),.5_dp))
      b=max(s/max(abs(m),.5_dp),.2_dp)
      c=0.0_dp
      case(GAMLSS_BCCG)
      a=max(abs(m),.5_dp)
      b=max(s/max(abs(m),.5_dp),.1_dp)
      c=1.0_dp
      case(GAMLSS_BCT,GAMLSS_BCPE)
         a=max(abs(m),.5_dp)
         b=max(s/max(abs(m),.5_dp),.1_dp)
         c=1.0_dp
         d=2.0_dp
      case(GAMLSS_GIG)
         a=max(abs(m),0.5_dp)
         b=0.7_dp
         c=0.0_dp
      case(GAMLSS_SHASHO)
         a=m
         b=max(s,0.2_dp)
         c=0.0_dp
         d=1.0_dp
      case(GAMLSS_SHASH)
         a=m
         b=max(s,0.2_dp)
         c=1.0_dp
         d=1.0_dp
      case(GAMLSS_SIMPLEX)
         a=min(0.9_dp,max(0.1_dp,m))
         b=0.5_dp
      case(GAMLSS_SEP,GAMLSS_SEP1,GAMLSS_SEP2,GAMLSS_ST1,GAMLSS_ST2)
         a=m
         b=max(s,0.2_dp)
         c=0.1_dp
         d=2.0_dp
      case(GAMLSS_ST5)
         a=m
         b=max(s,0.2_dp)
         c=0.0_dp
         d=1.0_dp
      case(GAMLSS_ST3,GAMLSS_ST4)
         a=m
         b=max(s,0.2_dp)
         c=1.0_dp
         d=10.0_dp
      case(GAMLSS_SEP3,GAMLSS_SEP4)
         a=m
         b=max(s,0.2_dp)
         c=2.0_dp
         d=2.0_dp
      case(GAMLSS_NET)
         a=m
         b=max(s,0.2_dp)
         c=1.5_dp
         d=0.7_dp
      case(GAMLSS_GPO,GAMLSS_DPO)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
      case(GAMLSS_DEL)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
         c=0.5_dp
      case(GAMLSS_SI,GAMLSS_SICHEL)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
         c=-0.5_dp
      case(GAMLSS_YULE,GAMLSS_ZIPF)
         a=max(abs(m),0.5_dp)
      case(GAMLSS_WARING)
         a=max(abs(m),0.5_dp)
         b=1.0_dp
      case(GAMLSS_ST3C)
         a=m
         b=max(s,0.2_dp)
         c=1.0_dp
         d=10.0_dp
      case(GAMLSS_SST)
         a=m
         b=max(s,0.2_dp)
         c=1.0_dp
         d=7.0_dp
      case(GAMLSS_GT)
         a=m
         b=max(s,0.2_dp)
         c=3.0_dp
         d=1.5_dp
      case(GAMLSS_SN1)
         a=m
         b=max(s,0.2_dp)
         c=0.0_dp
      case(GAMLSS_SN2)
         a=m
         b=max(s,0.2_dp)
         c=1.0_dp
      case(GAMLSS_EXGAUS)
         a=m
         b=max(0.7_dp*s,0.2_dp)
         c=max(0.3_dp*s,0.2_dp)
      case(GAMLSS_PARETO,GAMLSS_PARETO1)
         a=2.0_dp
      case(GAMLSS_PARETO2,GAMLSS_PARETO2O)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
      case(GAMLSS_PIG2)
         a=max(abs(m),0.5_dp)
         b=max(s,0.2_dp)
      case(GAMLSS_ZIPIG,GAMLSS_ZAPIG)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
         c=0.1_dp
      case(GAMLSS_ZISICHEL,GAMLSS_ZASICHEL)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
         c=-0.5_dp
         d=0.1_dp
      case(GAMLSS_ZIBNB,GAMLSS_ZABNB)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
         c=1.0_dp
         d=0.1_dp
      case(GAMLSS_ZAZIPF)
         a=1.0_dp
         b=0.1_dp
      case(GAMLSS_GAF)
         a=max(abs(m),0.5_dp)
         b=max(s/max(abs(m),0.5_dp),0.2_dp)
         c=2.0_dp
      case(GAMLSS_NBF)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
         c=2.0_dp
      case(GAMLSS_ZINBF)
         a=max(abs(m),0.5_dp)
         b=0.5_dp
         c=2.0_dp
         d=0.1_dp
      end select
   end subroutine

   elemental real(dp) function logistic_local(x) result(v)
      real(dp),intent(in)::x
      if(x>=0)then
      v=1/(1+exp(-min(x,700.0_dp)))
      else
      v=exp(max(x,-700.0_dp))/(1+exp(max(x,-700.0_dp)))
      end if
   end function

   subroutine map_parameters(family,e1,e2,e3,e4,a,b,c,d)
      integer,intent(in)::family
      real(dp),intent(in)::e1,e2,e3,e4
      real(dp),intent(out)::a,b,c,d
      a=e1
      b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      c=e3
      d=e4
      select case(family)
      case(GAMLSS_GA,GAMLSS_NBI,GAMLSS_NBII,GAMLSS_PIG,GAMLSS_WEI)
      a=exp(max(-50.0_dp,min(50.0_dp,e1)))
      b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      case(GAMLSS_BE)
      a=logistic_local(e1)
      b=logistic_local(e2)
      case(GAMLSS_ZIP)
      a=exp(max(-50.0_dp,min(50.0_dp,e1)))
      b=logistic_local(e2)
      case(GAMLSS_GG)
      a=exp(max(-50.0_dp,min(50.0_dp,e1)))
      b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      c=e3
      case(GAMLSS_EGB2)
      a=e1
      b=e2
      if(abs(b)<1e-7_dp)b=sign(1e-7_dp,b+1e-20_dp)
      c=exp(max(-50.0_dp,min(50.0_dp,e3)))
      d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_GB2)
      a=exp(max(-50.0_dp,min(50.0_dp,e1)))
      b=e2
      if(abs(b)<1e-7_dp)b=sign(1e-7_dp,b+1e-20_dp)
      c=exp(max(-50.0_dp,min(50.0_dp,e3)))
      d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_JSUO)
      a=e1
      b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      c=e3
      d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_TF)
      a=e1
      b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      c=exp(max(-50.0_dp,min(50.0_dp,e3)))
      case(GAMLSS_BEINF)
      a=logistic_local(e1)
      b=logistic_local(e2)
      c=exp(max(-50.0_dp,min(50.0_dp,e3)))
      d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_LNO)
      a=e1
      b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      c=e3
      case(GAMLSS_BCCG)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
      case(GAMLSS_BCT,GAMLSS_BCPE)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
         d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_GIG)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
      case(GAMLSS_SHASHO)
         a=e1
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
         d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_SHASH)
         a=e1
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=exp(max(-50.0_dp,min(50.0_dp,e3)))
         d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_SIMPLEX)
         a=logistic_local(e1)
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      case(GAMLSS_SEP,GAMLSS_SEP1,GAMLSS_SEP2,GAMLSS_ST1,GAMLSS_ST2,GAMLSS_ST5)
         a=e1
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
         d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_ST3,GAMLSS_ST4,GAMLSS_SEP3,GAMLSS_SEP4)
         a=e1
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=exp(max(-50.0_dp,min(50.0_dp,e3)))
         d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_NET)
         a=e1
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=exp(max(-50.0_dp,min(50.0_dp,e3)))
         d=c+exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_GPO,GAMLSS_DPO,GAMLSS_WARING)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      case(GAMLSS_DEL)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=logistic_local(e3)
      case(GAMLSS_SI,GAMLSS_SICHEL)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
      case(GAMLSS_YULE,GAMLSS_ZIPF,GAMLSS_PARETO,GAMLSS_PARETO1)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
      case(GAMLSS_ST3C,GAMLSS_SST,GAMLSS_GT)
         a=e1
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=exp(max(-50.0_dp,min(50.0_dp,e3)))
         d=exp(max(-50.0_dp,min(50.0_dp,e4)))
      case(GAMLSS_SN1)
         a=e1
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
      case(GAMLSS_SN2,GAMLSS_EXGAUS)
         a=e1
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=exp(max(-50.0_dp,min(50.0_dp,e3)))
      case(GAMLSS_PARETO2,GAMLSS_PARETO2O,GAMLSS_PIG2)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
      case(GAMLSS_ZIPIG,GAMLSS_ZAPIG)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=logistic_local(e3)
      case(GAMLSS_ZISICHEL,GAMLSS_ZASICHEL)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
         d=logistic_local(e4)
      case(GAMLSS_ZIBNB,GAMLSS_ZABNB)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=exp(max(-50.0_dp,min(50.0_dp,e3)))
         d=logistic_local(e4)
      case(GAMLSS_ZAZIPF)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=logistic_local(e2)
      case(GAMLSS_GAF)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=e3
      case(GAMLSS_NBF)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=exp(max(-50.0_dp,min(50.0_dp,e3)))
      case(GAMLSS_ZINBF)
         a=exp(max(-50.0_dp,min(50.0_dp,e1)))
         b=exp(max(-50.0_dp,min(50.0_dp,e2)))
         c=exp(max(-50.0_dp,min(50.0_dp,e3)))
         d=logistic_local(e4)
      end select
   end subroutine

   real(dp) function inverse_link(family,j,p) result(e)
      integer,intent(in)::family,j
      real(dp),intent(in)::p
      e=p
      select case(family)
      case(GAMLSS_GA,GAMLSS_NBI,GAMLSS_NBII,GAMLSS_PIG,GAMLSS_WEI)
         e=log(max(p,1e-10_dp))
      case(GAMLSS_BE)
         e=log(max(p,1e-10_dp)/max(1-p,1e-10_dp))
      case(GAMLSS_ZIP)
         if(j==1)e=log(max(p,1e-10_dp))
         if(j==2)e=log(max(p,1e-10_dp)/max(1-p,1e-10_dp))
      case(GAMLSS_GG)
         if(j<=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_EGB2)
         if(j>=3)e=log(max(p,1e-10_dp))
      case(GAMLSS_GB2)
         if(j==1.or.j>=3)e=log(max(p,1e-10_dp))
      case(GAMLSS_JSUO)
         if(j==2.or.j==4)e=log(max(p,1e-10_dp))
      case(GAMLSS_TF)
         if(j>=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_BEINF)
         if(j<=2)e=log(max(p,1e-10_dp)/max(1-p,1e-10_dp))
         if(j>=3)e=log(max(p,1e-10_dp))
      case(GAMLSS_LNO)
         if(j==2)e=log(max(p,1e-10_dp))
      case(GAMLSS_BCCG)
         if(j<=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_BCT,GAMLSS_BCPE)
         if(j<=2.or.j==4)e=log(max(p,1e-10_dp))
      case(GAMLSS_GIG)
         if(j<=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_SHASHO)
         if(j==2.or.j==4)e=log(max(p,1e-10_dp))
      case(GAMLSS_SHASH)
         if(j>=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_SIMPLEX)
         if(j==1)e=log(max(p,1e-10_dp)/max(1.0_dp-p,1e-10_dp))
         if(j==2)e=log(max(p,1e-10_dp))
      case(GAMLSS_SEP,GAMLSS_SEP1,GAMLSS_SEP2,GAMLSS_ST1,GAMLSS_ST2,GAMLSS_ST5)
         if(j==2.or.j==4)e=log(max(p,1e-10_dp))
      case(GAMLSS_ST3,GAMLSS_ST4,GAMLSS_SEP3,GAMLSS_SEP4)
         if(j>=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_NET)
         if(j==2.or.j==3.or.j==4)e=log(max(p,1e-10_dp))
      case(GAMLSS_GPO,GAMLSS_DPO,GAMLSS_WARING)
         e=log(max(p,1e-10_dp))
      case(GAMLSS_DEL)
         if(j<=2)e=log(max(p,1e-10_dp))
         if(j==3)e=log(max(p,1e-10_dp)/max(1.0_dp-p,1e-10_dp))
      case(GAMLSS_SI,GAMLSS_SICHEL)
         if(j<=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_YULE,GAMLSS_ZIPF,GAMLSS_PARETO,GAMLSS_PARETO1)
         e=log(max(p,1e-10_dp))
      case(GAMLSS_ST3C,GAMLSS_SST,GAMLSS_GT)
         if(j>=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_SN1)
         if(j==2)e=log(max(p,1e-10_dp))
      case(GAMLSS_SN2,GAMLSS_EXGAUS)
         if(j>=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_PARETO2,GAMLSS_PARETO2O,GAMLSS_PIG2)
         e=log(max(p,1e-10_dp))
      case(GAMLSS_ZIPIG,GAMLSS_ZAPIG)
         if(j<=2)e=log(max(p,1e-10_dp))
         if(j==3)e=log(max(p,1e-10_dp)/max(1.0_dp-p,1e-10_dp))
      case(GAMLSS_ZISICHEL,GAMLSS_ZASICHEL)
         if(j<=2)e=log(max(p,1e-10_dp))
         if(j==4)e=log(max(p,1e-10_dp)/max(1.0_dp-p,1e-10_dp))
      case(GAMLSS_ZIBNB,GAMLSS_ZABNB)
         if(j<=3)e=log(max(p,1e-10_dp))
         if(j==4)e=log(max(p,1e-10_dp)/max(1.0_dp-p,1e-10_dp))
      case(GAMLSS_ZAZIPF)
         if(j==1)e=log(max(p,1e-10_dp))
         if(j==2)e=log(max(p,1e-10_dp)/max(1.0_dp-p,1e-10_dp))
      case(GAMLSS_GAF)
         if(j<=2)e=log(max(p,1e-10_dp))
      case(GAMLSS_NBF)
         e=log(max(p,1e-10_dp))
      case(GAMLSS_ZINBF)
         if(j<=3)e=log(max(p,1e-10_dp))
         if(j==4)e=log(max(p,1e-10_dp)/max(1.0_dp-p,1e-10_dp))
      end select
   end function

   real(dp) function family_logpdf(family,y,a,b,c,d) result(lp)
      integer,intent(in)::family
      real(dp),intent(in)::y,a,b,c,d
      select case(family)
      case(GAMLSS_NO)
      lp=dNO(y,a,b,.true.)
      case(GAMLSS_GA)
      lp=dGA(y,a,b,.true.)
      case(GAMLSS_BE)
      lp=dBE(y,a,b,.true.)
      case(GAMLSS_NBI)
      lp=dNBI(y,a,b,.true.)
      case(GAMLSS_NBII)
      lp=dNBII(y,a,b,.true.)
      case(GAMLSS_ZIP)
      lp=dZIP(y,a,b,.true.)
      case(GAMLSS_GG)
      lp=dGG(y,a,b,c,.true.)
      case(GAMLSS_EGB2)
      lp=dEGB2(y,a,b,c,d,.true.)
      case(GAMLSS_GB2)
      lp=dGB2(y,a,b,c,d,.true.)
      case(GAMLSS_JSUO)
      lp=dJSUo(y,a,b,c,d,.true.)
      case(GAMLSS_TF)
      lp=dTF(y,a,b,c,.true.)
      case(GAMLSS_PIG)
      lp=dPIG(y,a,b,.true.)
      case(GAMLSS_BEINF)
      lp=dBEINF(y,a,b,c,d,.true.)
      case(GAMLSS_WEI)
      lp=dWEI(y,a,b,.true.)
      case(GAMLSS_LNO)
      lp=dLNO(y,a,b,c,.true.)
      case(GAMLSS_BCCG)
      lp=dBCCG(y,a,b,c,.true.)
      case(GAMLSS_BCT)
      lp=dBCT(y,a,b,c,d,.true.)
      case(GAMLSS_BCPE)
      lp=dBCPE(y,a,b,c,d,.true.)
      case(GAMLSS_GIG)
      lp=dGIG(y,a,b,c,.true.)
      case(GAMLSS_SHASHO)
      lp=dSHASHo(y,a,b,c,d,.true.)
      case(GAMLSS_SHASH)
      lp=dSHASH(y,a,b,c,d,.true.)
      case(GAMLSS_SIMPLEX)
      lp=dSIMPLEX(y,a,b,.true.)
      case(GAMLSS_SEP)
      lp=dSEP(y,a,b,c,d,.true.)
      case(GAMLSS_SEP1)
      lp=dSEP1(y,a,b,c,d,.true.)
      case(GAMLSS_SEP2)
      lp=dSEP2(y,a,b,c,d,.true.)
      case(GAMLSS_ST1)
      lp=dST1(y,a,b,c,d,.true.)
      case(GAMLSS_ST2)
      lp=dST2(y,a,b,c,d,.true.)
      case(GAMLSS_ST3)
      lp=dST3(y,a,b,c,d,.true.)
      case(GAMLSS_ST4)
      lp=dST4(y,a,b,c,d,.true.)
      case(GAMLSS_SEP3)
      lp=dSEP3(y,a,b,c,d,.true.)
      case(GAMLSS_SEP4)
      lp=dSEP4(y,a,b,c,d,.true.)
      case(GAMLSS_ST5)
      lp=dST5(y,a,b,c,d,.true.)
      case(GAMLSS_NET)
      lp=dNET(y,a,b,c,d,.true.)
      case(GAMLSS_GPO)
      lp=dGPO(y,a,b,.true.)
      case(GAMLSS_DPO)
      lp=dDPO(y,a,b,.true.)
      case(GAMLSS_DEL)
      lp=dDEL(y,a,b,c,.true.)
      case(GAMLSS_SI)
      lp=dSI(y,a,b,c,.true.)
      case(GAMLSS_SICHEL)
      lp=dSICHEL(y,a,b,c,.true.)
      case(GAMLSS_YULE)
      lp=dYULE(y,a,.true.)
      case(GAMLSS_WARING)
      lp=dWARING(y,a,b,.true.)
      case(GAMLSS_ZIPF)
      lp=dZIPF(y,a,.true.)
      case(GAMLSS_ST3C)
      lp=dST3C(y,a,b,c,d,.true.)
      case(GAMLSS_SN1)
      lp=dSN1(y,a,b,c,.true.)
      case(GAMLSS_SN2)
      lp=dSN2(y,a,b,c,.true.)
      case(GAMLSS_SST)
      lp=dSST(y,a,b,c,d,.true.)
      case(GAMLSS_GT)
      lp=dGT(y,a,b,c,d,.true.)
      case(GAMLSS_EXGAUS)
      lp=dexGAUS(y,a,b,c,.true.)
      case(GAMLSS_PARETO)
      lp=dPARETO(y,a,.true.)
      case(GAMLSS_PARETO1)
      lp=dPARETO1(y,a,.true.)
      case(GAMLSS_PARETO2)
      lp=dPARETO2(y,a,b,.true.)
      case(GAMLSS_PARETO2O)
      lp=dPARETO2o(y,a,b,.true.)
      case(GAMLSS_PIG2)
      lp=dPIG2(y,a,b,.true.)
      case(GAMLSS_ZIPIG)
      lp=dZIPIG(y,a,b,c,.true.)
      case(GAMLSS_ZAPIG)
      lp=dZAPIG(y,a,b,c,.true.)
      case(GAMLSS_ZISICHEL)
      lp=dZISICHEL(y,a,b,c,d,.true.)
      case(GAMLSS_ZASICHEL)
      lp=dZASICHEL(y,a,b,c,d,.true.)
      case(GAMLSS_ZIBNB)
      lp=dZIBNB(y,a,b,c,d,.true.)
      case(GAMLSS_ZABNB)
      lp=dZABNB(y,a,b,c,d,.true.)
      case(GAMLSS_ZAZIPF)
      lp=dZAZIPF(y,a,b,.true.)
      case(GAMLSS_GAF)
      lp=dGAF(y,a,b,c,.true.)
      case(GAMLSS_NBF)
      lp=dNBF(y,a,b,c,.true.)
      case(GAMLSS_ZINBF)
      lp=dZINBF(y,a,b,c,d,.true.)
      case default
      lp=-huge(1.0_dp)
      end select
   end function

   subroutine unpack_result(theta,p1,p2,p3,p4,result)
      real(dp),intent(in)::theta(:)
      integer,intent(in)::p1,p2,p3,p4
      type(gamlss_fit_result_t),intent(inout)::result
      integer::o2,o3,o4
      o2=p1
      o3=p1+p2
      o4=p1+p2+p3
      result%beta_mu=theta(1:p1)
      if(p2>0)result%beta_sigma=theta(o2+1:o3)
      if(p3>0)result%beta_nu=theta(o3+1:o4)
      if(p4>0)result%beta_tau=theta(o4+1:size(theta))
   end subroutine


end module gamlss_fit
