! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_inference
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, svystat_t, glm_result_t, FAMILY_GAUSSIAN, LINK_IDENTITY
  use survey_estimators, only : svy_mean
  use survey_glm, only : svy_glm
  use survey_design, only : design_degf
  use survey_quantiles, only : student_t_quantile, normal_quantile
  implicit none
  private
  integer,parameter,public::RANK_WILCOXON=1,RANK_VANDERWAERDEN=2,RANK_MEDIAN=3
  public :: svy_ttest_one, svy_ttest_two, svy_ci_prop, svy_ranktest, linear_contrast
contains

  subroutine svy_ttest_one(y,design,mean0,estimate,se,tstat,df,lower,upper,level)
    real(dp),intent(in)::y(:);type(survey_design_t),intent(in)::design;real(dp),intent(in),optional::mean0,level
    real(dp),intent(out)::estimate,se,tstat,lower,upper;integer,intent(out)::df
    real(dp),allocatable::x(:,:);type(svystat_t)::m;real(dp)::m0,lev,crit
    allocate(x(size(y),1));x(:,1)=y;m=svy_mean(x,design);estimate=m%estimate(1);se=sqrt(max(0.0_dp,m%variance(1,1)));m0=0;if(present(mean0))m0=mean0
    df=max(1,design_degf(design)-1);tstat=(estimate-m0)/max(se,tiny(1.0_dp));lev=0.95_dp;if(present(level))lev=level;crit=student_t_quantile((1+lev)/2,df);lower=estimate-crit*se;upper=estimate+crit*se
  end subroutine svy_ttest_one

  subroutine svy_ttest_two(y,group,design,difference,se,tstat,df,lower,upper,level)
    real(dp),intent(in)::y(:);integer,intent(in)::group(:);type(survey_design_t),intent(in)::design
    real(dp),intent(out)::difference,se,tstat,lower,upper;integer,intent(out)::df;real(dp),intent(in),optional::level
    real(dp),allocatable::x(:,:);type(glm_result_t)::fit;real(dp)::lev,crit;integer::g0,g1
    if(size(group)/=size(y))error stop 'svy_ttest_two: size mismatch';g0=minval(group);g1=maxval(group);if(g0==g1)error stop 'svy_ttest_two: group must be binary'
    allocate(x(size(y),2));x(:,1)=1;x(:,2)=merge(1.0_dp,0.0_dp,group==g1);call svy_glm(x,y,design,fit,FAMILY_GAUSSIAN,LINK_IDENTITY)
    difference=fit%coef(2);se=sqrt(max(0.0_dp,fit%vcov(2,2)));tstat=difference/max(se,tiny(1.0_dp));df=max(1,fit%df_residual);lev=.95_dp;if(present(level))lev=level;crit=student_t_quantile((1+lev)/2,df);lower=difference-crit*se;upper=difference+crit*se
  end subroutine svy_ttest_two

  subroutine svy_ci_prop(y,design,estimate,se,lower,upper,method,level)
    real(dp),intent(in)::y(:);type(survey_design_t),intent(in)::design;real(dp),intent(out)::estimate,se,lower,upper
    character(len=*),intent(in),optional::method;real(dp),intent(in),optional::level
    real(dp),allocatable::x(:,:);type(svystat_t)::m;real(dp)::lev,alpha,crit,eta,seeta,theta,neff,den,sq;integer::df;character(len=16)::meth
    allocate(x(size(y),1));x(:,1)=y;m=svy_mean(x,design);estimate=m%estimate(1);se=sqrt(max(0.0_dp,m%variance(1,1)));lev=.95_dp;if(present(level))lev=level;alpha=1-lev;df=max(1,design_degf(design));crit=student_t_quantile(1-alpha/2,df);meth='mean';if(present(method))meth=adjustl(method)
    select case(trim(meth))
    case('mean');lower=max(0.0_dp,estimate-crit*se);upper=min(1.0_dp,estimate+crit*se)
    case('logit','xlogit')
      theta=max(1e-12_dp,min(1-1e-12_dp,estimate));eta=log(theta/(1-theta));seeta=se/(theta*(1-theta));lower=1/(1+exp(-(eta-crit*seeta)));upper=1/(1+exp(-(eta+crit*seeta)))
    case('asin')
      theta=max(1e-12_dp,min(1-1e-12_dp,estimate));eta=asin(sqrt(theta));seeta=se/(2*sqrt(theta*(1-theta)));lower=sin(max(0.0_dp,eta-crit*seeta))**2;upper=sin(min(acos(-1.0_dp)/2,eta+crit*seeta))**2
    case('wilson')
      if(se>0)then;neff=max(1.0_dp,estimate*(1-estimate)/(se*se));else;neff=huge(1.0_dp);end if;den=1+crit*crit/neff;sq=sqrt(max(0.0_dp,4*neff*estimate*(1-estimate)+crit*crit));lower=(estimate+crit*crit/(2*neff)-crit*sq/(2*neff))/den;upper=(estimate+crit*crit/(2*neff)+crit*sq/(2*neff))/den
    case default;error stop 'svy_ci_prop: supported methods are mean, logit/xlogit, asin, wilson'
    end select
  end subroutine svy_ci_prop

  subroutine svy_ranktest(y,group,design,statistic,estimate,se,df,test)
    real(dp),intent(in)::y(:);integer,intent(in)::group(:);type(survey_design_t),intent(in)::design
    real(dp),intent(out)::statistic,estimate,se;integer,intent(out)::df;integer,intent(in),optional::test
    integer::rt,n,i,j,k,g1;integer,allocatable::ord(:);real(dp),allocatable::rankhat(:),score(:),x(:,:),ys(:),ws(:);real(dp)::Ntot,mid
    type(glm_result_t)::fit
    n=size(y);if(size(group)/=n.or.design%n/=n)error stop 'svy_ranktest: size mismatch';rt=RANK_WILCOXON;if(present(test))rt=test
    allocate(ord(n),ys(n),ws(n),rankhat(n),score(n),x(n,2));ord=[(i,i=1,n)];call sort_order(y,ord);ys=y(ord);ws=design%weight(ord);Ntot=sum(ws);rankhat=0
    i=1;do while(i<=n);j=i;do while(j<n.and.abs(ys(j+1)-ys(i))<=epsilon(1.0_dp)*max(1.0_dp,abs(ys(i))));j=j+1;end do;mid=0;do k=i,j;mid=mid+(sum(ws(1:k))-ws(k)/2);end do;mid=mid/real(j-i+1,dp);do k=i,j;rankhat(ord(k))=mid;end do;i=j+1;end do
    select case(rt);case(RANK_WILCOXON);score=rankhat/Ntot;case(RANK_VANDERWAERDEN);do i=1,n;score(i)=normal_quantile(max(1e-10_dp,min(1-1e-10_dp,rankhat(i)/Ntot)));end do;case(RANK_MEDIAN);score=merge(1.0_dp,0.0_dp,rankhat>Ntot/2);end select
    g1=maxval(group);x(:,1)=1;x(:,2)=merge(1.0_dp,0.0_dp,group==g1);call svy_glm(x,score,design,fit,FAMILY_GAUSSIAN,LINK_IDENTITY);estimate=fit%coef(2);se=sqrt(max(0.0_dp,fit%vcov(2,2)));statistic=estimate/max(se,tiny(1.0_dp));df=max(1,design_degf(design)-1)
  end subroutine svy_ranktest

  subroutine linear_contrast(coef,vcov,L,estimate,covariance)
    real(dp),intent(in)::coef(:),vcov(:,:),L(:,:);real(dp),intent(out)::estimate(:),covariance(:,:)
    if(size(L,2)/=size(coef).or.any(shape(vcov)/=[size(coef),size(coef)]))error stop 'linear_contrast: shape mismatch';estimate=matmul(L,coef);covariance=matmul(L,matmul(vcov,transpose(L)))
  end subroutine linear_contrast

  subroutine sort_order(y,ord);real(dp),intent(in)::y(:);integer,intent(inout)::ord(:);integer::i,j,k;do i=2,size(ord);k=ord(i);j=i-1;do while(j>=1.and.y(ord(j))>y(k));ord(j+1)=ord(j);j=j-1;end do;ord(j+1)=k;end do;end subroutine sort_order
end module survey_inference
