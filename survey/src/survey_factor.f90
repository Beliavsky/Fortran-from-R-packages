! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_factor
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, factor_result_t, FACT_N_NONE, FACT_N_SAMPLE, &
    FACT_N_DEGF, FACT_N_EFFECTIVE, FACT_N_MIN_EFFECTIVE
  use survey_design, only : design_degf
  use survey_estimators, only : svy_covariance
  use survey_linalg, only : inverse_matrix, symmetric_eigen
  use survey_special, only : chisq_survival
  use minqa_module, only : minqa_result_t, minqa_control_t, bobyqa
  implicit none
  private
  public :: svy_factanal, factor_ml_cov, varimax_rotate

  real(dp), allocatable :: active_corr(:,:)
  integer :: active_nf = 0
contains

  subroutine svy_factanal(x,design,factors,result,n_mode,rotate,maxfun)
    real(dp),intent(in)::x(:,:)
    type(survey_design_t),intent(in)::design
    integer,intent(in)::factors
    type(factor_result_t),intent(out)::result
    integer,intent(in),optional::n_mode,maxfun
    logical,intent(in),optional::rotate
    real(dp),allocatable::cov(:,:),varcov(:,:),ses2(:),neff(:)
    real(dp)::nobs
    integer::p,i,k,mode
    if(size(x,1)/=design%n)error stop 'svy_factanal: row mismatch'
    p=size(x,2);if(factors<1.or.factors>=(p-1))error stop 'svy_factanal: invalid number of factors'
    allocate(cov(p,p),varcov(p*p,p*p),ses2(p),neff(p))
    call svy_covariance(x,design,cov,varcov)
    do i=1,p
      k=(i-1)*p+i
      ses2(i)=sqrt(max(0.0_dp,varcov(k,k)))
      if(ses2(i)>tiny(1.0_dp))then
        neff(i)=2.0_dp*(cov(i,i)/ses2(i))**2
      else
        neff(i)=huge(1.0_dp)
      end if
    end do
    mode=FACT_N_NONE;if(present(n_mode))mode=n_mode
    select case(mode)
    case(FACT_N_NONE);nobs=0.0_dp
    case(FACT_N_SAMPLE);nobs=real(design%n,dp)
    case(FACT_N_DEGF);nobs=real(design_degf(design)+1,dp)
    case(FACT_N_EFFECTIVE);nobs=1.0_dp/(sum(1.0_dp/max(neff,1.0_dp))/real(p,dp))+1.0_dp
    case(FACT_N_MIN_EFFECTIVE);nobs=minval(neff)+1.0_dp
    case default;error stop 'svy_factanal: invalid n_mode'
    end select
    call factor_ml_cov(cov,factors,result,nobs,rotate,maxfun)
  end subroutine svy_factanal

  subroutine factor_ml_cov(cov,factors,result,nobs,rotate,maxfun)
    real(dp),intent(in)::cov(:,:)
    integer,intent(in)::factors
    type(factor_result_t),intent(out)::result
    real(dp),intent(in),optional::nobs
    logical,intent(in),optional::rotate
    integer,intent(in),optional::maxfun
    type(minqa_result_t)::opt
    type(minqa_control_t)::ctrl
    real(dp),allocatable::corr(:,:),sd(:),psi(:),lo(:),hi(:),cinv(:,:),loads(:,:)
    real(dp)::nn,bc
    integer::p,i,j,info
    logical::dorot
    p=size(cov,1);if(size(cov,2)/=p.or.factors<1.or.factors>=p-1)error stop 'factor_ml_cov: invalid shape/factors'
    allocate(corr(p,p),sd(p),psi(p),lo(p),hi(p),cinv(p,p))
    do i=1,p
      if(cov(i,i)<=0.0_dp)error stop 'factor_ml_cov: nonpositive variance'
      sd(i)=sqrt(cov(i,i))
    end do
    do j=1,p;do i=1,p;corr(i,j)=cov(i,j)/(sd(i)*sd(j));end do;end do
    do i=1,p;corr(i,i)=1.0_dp;end do
    call inverse_matrix(corr,cinv,info)
    if(info==0)then
      do i=1,p;psi(i)=min(0.99_dp,max(0.05_dp,1.0_dp/max(cinv(i,i),1.0_dp)));end do
    else
      psi=0.5_dp
    end if
    lo=0.005_dp;hi=1.0_dp
    if(allocated(active_corr))error stop 'factor_ml_cov: nested/concurrent calls are not supported'
    allocate(active_corr,source=corr);active_nf=factors
    ctrl%maxfun=10000;if(present(maxfun))ctrl%maxfun=maxfun
    call bobyqa(factor_objective,psi,opt,lo,hi,ctrl);psi=opt%x
    call factor_loadings(psi,corr,factors,loads)
    dorot=.true.;if(present(rotate))dorot=rotate
    if(dorot.and.factors>1)call varimax_rotate(loads)
    allocate(result%loadings(p,factors),result%uniqueness(p),result%communalities(p))
    result%loadings=loads;result%uniqueness=psi;result%communalities=sum(loads*loads,dim=2)
    result%criterion=factor_objective(psi);result%factors=factors;result%iterations=opt%evaluations;result%converged=(opt%status==0)
    nn=0.0_dp;if(present(nobs))nn=nobs;result%effective_n=nn
    result%dof=((p-factors)*(p-factors)-p-factors)/2
    if(nn>0.0_dp.and.result%dof>0)then
      bc=nn-1.0_dp-(2.0_dp*real(p,dp)+5.0_dp)/6.0_dp-2.0_dp*real(factors,dp)/3.0_dp
      result%statistic=max(0.0_dp,bc*result%criterion)
      result%p_value=chisq_survival(result%statistic,real(result%dof,dp))
    else
      result%statistic=0.0_dp;result%p_value=1.0_dp
    end if
    deallocate(active_corr);active_nf=0
  end subroutine factor_ml_cov

  real(dp) function factor_objective(psi) result(value)
    real(dp),intent(in)::psi(:)
    real(dp),allocatable::sstar(:,:),eval(:),evec(:,:)
    integer::p,i,j,info
    if(.not.allocated(active_corr))error stop 'factor_objective: no active covariance'
    p=size(psi);allocate(sstar(p,p),eval(p),evec(p,p))
    do j=1,p;do i=1,p;sstar(i,j)=active_corr(i,j)/sqrt(psi(i)*psi(j));end do;end do
    call symmetric_eigen(sstar,eval,evec,info)
    if(info/=0.or.any(eval(active_nf+1:p)<=tiny(1.0_dp)))then
      value=huge(1.0_dp)/100.0_dp;return
    end if
    value=-sum(log(eval(active_nf+1:p))-eval(active_nf+1:p))+real(active_nf-p,dp)
  end function factor_objective

  subroutine factor_loadings(psi,corr,nf,loadings)
    real(dp),intent(in)::psi(:),corr(:,:)
    integer,intent(in)::nf
    real(dp),allocatable,intent(out)::loadings(:,:)
    real(dp),allocatable::sstar(:,:),eval(:),evec(:,:)
    integer::p,i,j,k,info
    p=size(psi);allocate(sstar(p,p),eval(p),evec(p,p),loadings(p,nf))
    do j=1,p;do i=1,p;sstar(i,j)=corr(i,j)/sqrt(psi(i)*psi(j));end do;end do
    call symmetric_eigen(sstar,eval,evec,info);if(info/=0)error stop 'factor_loadings: eigensolver failed'
    do k=1,nf
      loadings(:,k)=sqrt(psi)*evec(:,k)*sqrt(max(eval(k)-1.0_dp,0.0_dp))
    end do
  end subroutine factor_loadings

  subroutine varimax_rotate(loadings,maxit,tol)
    real(dp),intent(inout)::loadings(:,:)
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    integer::p,q,mi,it,a,b,i
    real(dp)::eps,u,v,aa,bb,cc,dd,theta,c,s,maxang,la,lb
    p=size(loadings,1);q=size(loadings,2);if(q<2)return
    mi=100;if(present(maxit))mi=maxit;eps=1e-8_dp;if(present(tol))eps=tol
    do it=1,mi
      maxang=0.0_dp
      do a=1,q-1
        do b=a+1,q
          aa=0.0_dp;bb=0.0_dp;cc=0.0_dp;dd=0.0_dp
          do i=1,p
            u=loadings(i,a)**2-loadings(i,b)**2;v=2.0_dp*loadings(i,a)*loadings(i,b)
            aa=aa+u;bb=bb+v;cc=cc+u*u-v*v;dd=dd+2.0_dp*u*v
          end do
          cc=cc-(aa*aa-bb*bb)/real(p,dp);dd=dd-2.0_dp*aa*bb/real(p,dp)
          theta=0.25_dp*atan2(dd,cc);maxang=max(maxang,abs(theta));c=cos(theta);s=sin(theta)
          do i=1,p
            la=loadings(i,a);lb=loadings(i,b);loadings(i,a)=c*la+s*lb;loadings(i,b)=-s*la+c*lb
          end do
        end do
      end do
      if(maxang<eps)exit
    end do
  end subroutine varimax_rotate

end module survey_factor
