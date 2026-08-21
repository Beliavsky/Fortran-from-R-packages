! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_score
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, svystat_t, model_test_t, &
    FAMILY_GAUSSIAN,FAMILY_BINOMIAL,FAMILY_POISSON,LINK_IDENTITY,LINK_LOGIT,LINK_LOG
  use survey_estimators, only : svy_total
  use survey_design, only : design_degf
  use survey_linalg, only : sym_pinv, symmetric_eigen, outer_product
  use survey_special, only : f_survival, chisq_survival, weighted_f_survival, weighted_chisq_survival
  implicit none
  private
  public :: glm_pseudoscore_test, glm_working_score_test, wald_term_test, misspecified_lrt_test
contains

  subroutine glm_pseudoscore_test(x,y,design,beta0,test_index,result,family,link)
    real(dp),intent(in) :: x(:,:),y(:),beta0(:)
    type(survey_design_t),intent(in) :: design
    integer,intent(in) :: test_index(:)
    type(model_test_t),intent(out) :: result
    integer,intent(in),optional :: family,link
    real(dp),allocatable :: u(:,:),info(:,:),vc(:,:),coef(:),vinv(:,:)
    type(svystat_t) :: st
    integer :: fam,lnk,q,rank,istat,ddf,nactive
    call residual_scores(x,y,design,beta0,test_index,u,info,fam,lnk,family,link)
    q=size(u,2);st=svy_total(u,design);vc=st%variance;coef=st%estimate
    allocate(vinv(q,q));call sym_pinv(vc,vinv,rank,info=istat)
    result%statistic=dot_product(coef,matmul(vinv,coef));result%df=rank;nactive=count(abs(beta0)>sqrt(epsilon(1.0_dp)));ddf=max(1,design_degf(design)+1-nactive);result%ddf=ddf
    if(ddf>0) then;result%p_value=f_survival(result%statistic/real(max(1,rank),dp),real(max(1,rank),dp),real(ddf,dp));else;result%p_value=chisq_survival(result%statistic,real(max(1,rank),dp));end if
    allocate(result%lambda(0))
  end subroutine glm_pseudoscore_test

  subroutine glm_working_score_test(x,y,design,beta0,test_index,result,family,link,dispersion,saddlepoint)
    real(dp),intent(in) :: x(:,:),y(:),beta0(:)
    type(survey_design_t),intent(in) :: design
    integer,intent(in) :: test_index(:)
    type(model_test_t),intent(out) :: result
    integer,intent(in),optional :: family,link
    real(dp),intent(in),optional :: dispersion
    logical,intent(in),optional :: saddlepoint
    real(dp),allocatable :: u(:,:),info_schur(:,:),robust(:,:),coef(:),invinfo(:,:),lambda(:)
    real(dp) :: phi,mwt
    type(svystat_t) :: st
    integer :: fam,lnk,q,rank,istat,ddf,nactive
    logical :: sad
    call residual_scores(x,y,design,beta0,test_index,u,info_schur,fam,lnk,family,link)
    q=size(u,2);st=svy_total(u,design);coef=st%estimate;robust=st%variance
    allocate(invinfo(q,q));call sym_pinv(info_schur,invinfo,rank,info=istat)
    phi=1.0_dp;if(present(dispersion))phi=dispersion
    if(phi<=0.0_dp) error stop 'glm_working_score_test: dispersion must be positive'
    mwt=sum(design%weight)/real(design%n,dp)
    result%statistic=dot_product(coef,matmul(invinfo,coef))/(phi*mwt)
    call generalized_eigenvalues(info_schur,robust/(phi*mwt),lambda)
    result%df=size(lambda);nactive=count(abs(beta0)>sqrt(epsilon(1.0_dp)));ddf=max(1,design_degf(design)+1-nactive);result%ddf=ddf;sad=.true.;if(present(saddlepoint))sad=saddlepoint
    result%p_value=weighted_f_survival(result%statistic,lambda,real(ddf,dp),saddlepoint=sad)
    allocate(result%lambda(size(lambda)));result%lambda=lambda
  end subroutine glm_working_score_test

  subroutine wald_term_test(beta,vcov,index,ddf,result,chisq_test)
    real(dp),intent(in) :: beta(:),vcov(:,:)
    integer,intent(in) :: index(:),ddf
    type(model_test_t),intent(out) :: result
    logical,intent(in),optional :: chisq_test
    real(dp),allocatable :: b(:),v(:,:),vinv(:,:)
    integer :: i,j,q,rank,istat
    logical :: ch
    q=size(index);if(q<1)error stop 'wald_term_test: empty index';allocate(b(q),v(q,q),vinv(q,q))
    do i=1,q;b(i)=beta(index(i));do j=1,q;v(i,j)=vcov(index(i),index(j));end do;end do
    call sym_pinv(v,vinv,rank,info=istat);result%statistic=dot_product(b,matmul(vinv,b));result%df=rank;result%ddf=max(1,ddf);ch=.false.;if(present(chisq_test))ch=chisq_test
    if(ch) then;result%p_value=chisq_survival(result%statistic,real(max(1,rank),dp));else;result%p_value=f_survival(result%statistic/real(max(1,rank),dp),real(max(1,rank),dp),real(result%ddf,dp));end if
    allocate(result%lambda(0))
  end subroutine wald_term_test

  subroutine misspecified_lrt_test(chisq,model_vcov,robust_vcov,ddf,result,chisq_test,saddlepoint)
    real(dp),intent(in) :: chisq,model_vcov(:,:),robust_vcov(:,:)
    integer,intent(in) :: ddf
    type(model_test_t),intent(out) :: result
    logical,intent(in),optional :: chisq_test,saddlepoint
    real(dp),allocatable :: lambda(:)
    logical :: ch,sad
    if(any(shape(model_vcov)/=shape(robust_vcov)).or.size(model_vcov,1)/=size(model_vcov,2)) error stop 'misspecified_lrt_test: covariance shape'
    call generalized_eigenvalues(model_vcov,robust_vcov,lambda);ch=.false.;if(present(chisq_test))ch=chisq_test;sad=.true.;if(present(saddlepoint))sad=saddlepoint
    result%statistic=chisq;result%df=size(lambda);result%ddf=max(1,ddf);allocate(result%lambda(size(lambda)));result%lambda=lambda
    if(ch) then;result%p_value=weighted_chisq_survival(chisq,lambda,saddlepoint=sad);else;result%p_value=weighted_f_survival(chisq,lambda,real(result%ddf,dp),saddlepoint=sad);end if
  end subroutine misspecified_lrt_test

  subroutine residual_scores(x,y,design,beta0,test_index,uout_in,info_schur,fam,lnk,family,link)
    real(dp),intent(in) :: x(:,:),y(:),beta0(:)
    type(survey_design_t),intent(in) :: design
    integer,intent(in) :: test_index(:)
    real(dp),allocatable,intent(out) :: uout_in(:,:),info_schur(:,:)
    integer,intent(out) :: fam,lnk
    integer,intent(in),optional :: family,link
    logical,allocatable :: out(:)
    integer,allocatable :: ii(:),oo(:)
    real(dp),allocatable :: eta(:),mu(:),dmu(:),varmu(:),uin(:,:),uout(:,:),inf(:,:),ai(:,:),cross(:,:)
    integer :: n,p,ni,no,i,j,k,rank,istat
    n=size(x,1);p=size(x,2);if(size(y)/=n.or.size(beta0)/=p.or.design%n/=n)error stop 'score test: shape mismatch'
    allocate(out(p));out=.false.;do i=1,size(test_index);if(test_index(i)<1.or.test_index(i)>p)error stop 'score test: index';out(test_index(i))=.true.;end do
    no=count(out);ni=p-no;if(no<1)error stop 'score test: no tested coefficients';allocate(oo(no),ii(ni));oo=pack([(i,i=1,p)],out);ii=pack([(i,i=1,p)],.not.out)
    fam=FAMILY_GAUSSIAN;if(present(family))fam=family;lnk=default_link(fam);if(present(link))lnk=link
    allocate(eta(n),mu(n),dmu(n),varmu(n),uout(n,no),inf(p,p));eta=matmul(x,beta0);call mean_components(eta,fam,lnk,mu,dmu,varmu)
    inf=0.0_dp
    do i=1,n
      do j=1,p
        do k=1,p;inf(j,k)=inf(j,k)+design%weight(i)*x(i,j)*x(i,k)*dmu(i)*dmu(i)/max(varmu(i),tiny(1.0_dp));end do
      end do
      do j=1,no;uout(i,j)=x(i,oo(j))*dmu(i)*(y(i)-mu(i))/max(varmu(i),tiny(1.0_dp));end do
    end do
    if(ni==0) then
      allocate(uout_in(n,no),info_schur(no,no));uout_in=uout;info_schur=inf(oo,oo);return
    end if
    allocate(uin(n,ni),ai(ni,ni),cross(no,ni),uout_in(n,no),info_schur(no,no))
    do i=1,n;do j=1,ni;uin(i,j)=x(i,ii(j))*dmu(i)*(y(i)-mu(i))/max(varmu(i),tiny(1.0_dp));end do;end do
    call sym_pinv(inf(ii,ii),ai,rank,info=istat);cross=matmul(inf(oo,ii),ai)
    uout_in=uout-matmul(uin,transpose(cross))
    info_schur=inf(oo,oo)-matmul(cross,inf(ii,oo))
  end subroutine residual_scores

  subroutine generalized_eigenvalues(a,b,lambda)
    real(dp),intent(in) :: a(:,:),b(:,:)
    real(dp),allocatable,intent(out) :: lambda(:)
    real(dp),allocatable :: eval(:),evec(:,:),invsqrt(:,:),c(:,:),ceval(:),cevec(:,:)
    real(dp) :: tol
    integer :: n,i,info,r
    n=size(a,1);if(size(a,2)/=n.or.any(shape(b)/=[n,n]))error stop 'generalized_eigenvalues: shape'
    allocate(eval(n),evec(n,n),invsqrt(n,n));call symmetric_eigen(a,eval,evec,info);tol=256.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(eval)));invsqrt=0.0_dp;r=0
    do i=1,n;if(eval(i)>tol)then;invsqrt=invsqrt+outer_product(evec(:,i),evec(:,i))/sqrt(eval(i));r=r+1;end if;end do
    if(r<1)error stop 'generalized_eigenvalues: singular model covariance';allocate(c(n,n),ceval(n),cevec(n,n));c=matmul(invsqrt,matmul(b,invsqrt));call symmetric_eigen(c,ceval,cevec,info)
    allocate(lambda(count(ceval>tol)));lambda=pack(ceval,ceval>tol)
    if(size(lambda)<1)error stop 'generalized_eigenvalues: no positive eigenvalues'
  end subroutine generalized_eigenvalues

  subroutine mean_components(eta,fam,lnk,mu,dmu,varmu)
    real(dp),intent(in)::eta(:);integer,intent(in)::fam,lnk;real(dp),intent(out)::mu(:),dmu(:),varmu(:);integer::i
    select case(lnk)
    case(LINK_IDENTITY);mu=eta;dmu=1.0_dp
    case(LINK_LOGIT);do i=1,size(eta);mu(i)=1.0_dp/(1.0_dp+exp(-max(-35.0_dp,min(35.0_dp,eta(i)))));dmu(i)=mu(i)*(1.0_dp-mu(i));end do
    case(LINK_LOG);mu=exp(max(-35.0_dp,min(35.0_dp,eta)));dmu=mu
    case default;error stop 'score test: unsupported link'
    end select
    select case(fam);case(FAMILY_GAUSSIAN);varmu=1.0_dp;case(FAMILY_BINOMIAL);varmu=max(mu*(1.0_dp-mu),tiny(1.0_dp));case(FAMILY_POISSON);varmu=max(mu,tiny(1.0_dp));case default;error stop 'score test: unsupported family';end select
  end subroutine mean_components

  integer pure function default_link(fam) result(lnk)
    integer,intent(in)::fam;select case(fam);case(FAMILY_GAUSSIAN);lnk=LINK_IDENTITY;case(FAMILY_BINOMIAL);lnk=LINK_LOGIT;case(FAMILY_POISSON);lnk=LINK_LOG;case default;lnk=LINK_IDENTITY;end select
  end function default_link
end module survey_score
