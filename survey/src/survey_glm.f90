! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_glm
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, rep_design_t, glm_result_t, &
    FAMILY_GAUSSIAN,FAMILY_BINOMIAL,FAMILY_POISSON,LINK_IDENTITY,LINK_LOGIT,LINK_LOG
  use survey_linalg, only : sym_pinv
  use survey_taylor, only : svyrecvar
  use survey_replicates, only : svr_var
  use survey_design, only : design_degf
  implicit none
  private
  public :: svy_glm, rep_glm, glm_predict, glm_pseudo_r2
contains

  subroutine svy_glm(x,y,design,result,family,link,start,maxit,tol,rescale_weights)
    real(dp), intent(in) :: x(:,:),y(:)
    type(survey_design_t), intent(in) :: design
    type(glm_result_t), intent(out) :: result
    integer, intent(in), optional :: family,link,maxit
    real(dp), intent(in), optional :: start(:),tol
    logical, intent(in), optional :: rescale_weights
    real(dp), allocatable :: beta(:),beta_new(:),eta(:),mu(:),dmu(:),varmu(:),z(:),ww(:),prior(:), &
      xtwx(:,:),xtwz(:),ainv(:,:),score(:,:),infl(:,:)
    real(dp) :: eps,dev,devold,disp,sw
    integer :: n,p,fam,lnk,mi,it,i,j,k,rank,info
    logical :: resc
    n=size(x,1);p=size(x,2)
    if(size(y)/=n .or. design%n/=n) error stop 'svy_glm: shape mismatch'
    fam=FAMILY_GAUSSIAN; if(present(family)) fam=family
    lnk=LINK_IDENTITY
    if(present(link)) then; lnk=link; else; select case(fam);case(FAMILY_GAUSSIAN);lnk=LINK_IDENTITY;case(FAMILY_BINOMIAL);lnk=LINK_LOGIT;case(FAMILY_POISSON);lnk=LINK_LOG;end select;end if
    mi=50;if(present(maxit))mi=maxit;eps=1e-8_dp;if(present(tol))eps=tol;resc=.true.;if(present(rescale_weights))resc=rescale_weights
    allocate(beta(p),beta_new(p),eta(n),mu(n),dmu(n),varmu(n),z(n),ww(n),prior(n),xtwx(p,p),xtwz(p),ainv(p,p),score(n,p),infl(n,p))
    prior=design%weight
    if(resc) prior=prior/(sum(prior)/real(n,dp))
    if(present(start)) then;if(size(start)/=p)error stop 'svy_glm: start size mismatch';beta=start;else;beta=0;call initial_beta(x,y,prior,lnk,beta);end if
    devold=huge(1.0_dp)
    do it=1,mi
      eta=matmul(x,beta);call mean_components(eta,fam,lnk,mu,dmu,varmu)
      ww=prior*dmu*dmu/max(varmu,tiny(1.0_dp));z=eta+(y-mu)/merge(dmu,tiny(1.0_dp),abs(dmu)>tiny(1.0_dp))
      xtwx=0;xtwz=0
      do i=1,n
        do j=1,p
          xtwz(j)=xtwz(j)+ww(i)*x(i,j)*z(i)
          do k=1,p;xtwx(j,k)=xtwx(j,k)+ww(i)*x(i,j)*x(i,k);end do
        end do
      end do
      call sym_pinv(xtwx,ainv,rank,info=info);beta_new=matmul(ainv,xtwz)
      eta=matmul(x,beta_new);call mean_components(eta,fam,lnk,mu,dmu,varmu);dev=glm_deviance(y,mu,prior,fam)
      if(maxval(abs(beta_new-beta)/(1+abs(beta)))<eps .and. abs(dev-devold)/(1+abs(dev))<sqrt(eps)) then; beta=beta_new; exit;end if
      beta=beta_new;devold=dev
    end do
    eta=matmul(x,beta);call mean_components(eta,fam,lnk,mu,dmu,varmu)
    ww=prior*dmu*dmu/max(varmu,tiny(1.0_dp));xtwx=0
    do i=1,n;do j=1,p;do k=1,p;xtwx(j,k)=xtwx(j,k)+ww(i)*x(i,j)*x(i,k);end do;end do;end do
    call sym_pinv(xtwx,ainv,rank,info=info)
    do i=1,n;do j=1,p;score(i,j)=x(i,j)*prior(i)*(y(i)-mu(i))*dmu(i)/max(varmu(i),tiny(1.0_dp));end do;end do
    infl=matmul(score,ainv)
    allocate(result%coef(p),result%vcov(p,p),result%naive_vcov(p,p),result%fitted(n),result%residual(n))
    result%coef=beta;result%vcov=svyrecvar(infl,design);result%fitted=mu;result%residual=y-mu;result%rank=rank
    result%iterations=min(it,mi);result%converged=(it<=mi);result%deviance=glm_deviance(y,mu,prior,fam)
    select case(fam)
    case(FAMILY_GAUSSIAN)
      sw=sum(prior);disp=sum(prior*(y-mu)**2)/max(1.0_dp,sw-real(rank,dp));result%naive_vcov=ainv*disp
    case default
      result%naive_vcov=ainv
    end select
    result%df_residual=max(0,design_degf(design)+1-rank)
  end subroutine svy_glm

  subroutine rep_glm(x,y,design,result,family,link,maxit,tol)
    real(dp),intent(in)::x(:,:),y(:);type(rep_design_t),intent(in)::design;type(glm_result_t),intent(out)::result
    integer,intent(in),optional::family,link,maxit;real(dp),intent(in),optional::tol
    real(dp),allocatable::coefrep(:,:),coef(:),vv(:,:);integer::r,p,fam,lnk,mi;real(dp)::eps
    type(glm_result_t)::tmp;type(survey_design_t)::d
    p=size(x,2);fam=FAMILY_GAUSSIAN;if(present(family))fam=family;lnk=default_link(fam);if(present(link))lnk=link;mi=50;if(present(maxit))mi=maxit;eps=1e-8_dp;if(present(tol))eps=tol
    call pseudo_design(design%weight,d);call svy_glm(x,y,d,tmp,fam,lnk,maxit=mi,tol=eps)
    allocate(coef(p),coefrep(design%r,p),vv(p,p));coef=tmp%coef
    do r=1,design%r;call pseudo_design(design%repweights(:,r),d);call svy_glm(x,y,d,tmp,fam,lnk,maxit=mi,tol=eps,rescale_weights=.false.);coefrep(r,:)=tmp%coef;end do
    vv=svr_var(coefrep,design%scale,design%rscales,design%mse,coef)
    allocate(result%coef(p),result%vcov(p,p),result%naive_vcov(p,p),result%fitted(size(y)),result%residual(size(y)))
    call pseudo_design(design%weight,d);call svy_glm(x,y,d,tmp,fam,lnk,maxit=mi,tol=eps)
    result=tmp;result%vcov=vv
  end subroutine rep_glm

  function glm_predict(xnew,fit,link_scale) result(pred)
    real(dp),intent(in)::xnew(:,:);type(glm_result_t),intent(in)::fit;logical,intent(in),optional::link_scale
    real(dp)::pred(size(xnew,1));logical::ls
    ls=.true.;if(present(link_scale))ls=link_scale;pred=matmul(xnew,fit%coef)
    ! Generic fit does not retain family/link. Caller can request link scale; response-scale
    ! prediction should use the corresponding inverse link explicitly.
  end function glm_predict

  real(dp) function glm_pseudo_r2(y,mu,weight,family) result(r2)
    real(dp),intent(in)::y(:),mu(:),weight(:);integer,intent(in)::family
    real(dp)::m0,d0,d1
    m0=dot_product(weight,y)/sum(weight);d1=glm_deviance(y,mu,weight,family);d0=glm_deviance(y,spread_scalar(m0,size(y)),weight,family);r2=1-d1/d0
  end function glm_pseudo_r2

  subroutine mean_components(eta,fam,lnk,mu,dmu,varmu)
    real(dp),intent(in)::eta(:);integer,intent(in)::fam,lnk;real(dp),intent(out)::mu(:),dmu(:),varmu(:);integer::i
    select case(lnk)
    case(LINK_IDENTITY);mu=eta;dmu=1
    case(LINK_LOGIT);do i=1,size(eta);mu(i)=1/(1+exp(-max(-35.0_dp,min(35.0_dp,eta(i)))));dmu(i)=mu(i)*(1-mu(i));end do
    case(LINK_LOG);mu=exp(max(-35.0_dp,min(35.0_dp,eta)));dmu=mu
    case default;error stop 'svy_glm: unsupported link'
    end select
    select case(fam)
    case(FAMILY_GAUSSIAN);varmu=1
    case(FAMILY_BINOMIAL);varmu=max(mu*(1-mu),tiny(1.0_dp))
    case(FAMILY_POISSON);varmu=max(mu,tiny(1.0_dp))
    case default;error stop 'svy_glm: unsupported family'
    end select
  end subroutine mean_components

  real(dp) function glm_deviance(y,mu,w,fam) result(d)
    real(dp),intent(in)::y(:),mu(:),w(:);integer,intent(in)::fam;integer::i;real(dp)::term
    d=0
    select case(fam)
    case(FAMILY_GAUSSIAN);d=sum(w*(y-mu)**2)
    case(FAMILY_BINOMIAL)
      do i=1,size(y);term=0;if(y(i)>0)term=term+y(i)*log(y(i)/max(mu(i),tiny(1.0_dp)));if(y(i)<1)term=term+(1-y(i))*log((1-y(i))/max(1-mu(i),tiny(1.0_dp)));d=d+2*w(i)*term;end do
    case(FAMILY_POISSON)
      do i=1,size(y);if(y(i)>0)then;term=y(i)*log(y(i)/max(mu(i),tiny(1.0_dp)))-(y(i)-mu(i));else;term=mu(i);end if;d=d+2*w(i)*term;end do
    end select
  end function glm_deviance

  subroutine initial_beta(x,y,w,lnk,beta)
    real(dp),intent(in)::x(:,:),y(:),w(:);integer,intent(in)::lnk;real(dp),intent(out)::beta(:)
    real(dp),allocatable::z(:),a(:,:),rhs(:),ainv(:,:);integer::i,j,k,p,rank,info
    p=size(x,2);allocate(z(size(y)),a(p,p),rhs(p),ainv(p,p));select case(lnk);case(LINK_IDENTITY);z=y;case(LINK_LOGIT);z=log(max(1e-4_dp,min(1-1e-4_dp,(y+0.5_dp)/2)));case(LINK_LOG);z=log(max(y,0.1_dp));end select
    a=0;rhs=0;do i=1,size(y);do j=1,p;rhs(j)=rhs(j)+w(i)*x(i,j)*z(i);do k=1,p;a(j,k)=a(j,k)+w(i)*x(i,j)*x(i,k);end do;end do;end do;call sym_pinv(a,ainv,rank,info=info);beta=matmul(ainv,rhs)
  end subroutine initial_beta

  integer pure function default_link(fam) result(lnk)
    integer,intent(in)::fam;select case(fam);case(FAMILY_GAUSSIAN);lnk=LINK_IDENTITY;case(FAMILY_BINOMIAL);lnk=LINK_LOGIT;case(FAMILY_POISSON);lnk=LINK_LOG;case default;lnk=LINK_IDENTITY;end select
  end function default_link
  function spread_scalar(x,n) result(v);real(dp),intent(in)::x;integer,intent(in)::n;real(dp)::v(n);v=x;end function spread_scalar
  subroutine pseudo_design(w,d)
    real(dp),intent(in)::w(:);type(survey_design_t),intent(out)::d;integer::n,i
    n=size(w);d%n=n;d%stages=1;allocate(d%weight(n),d%cluster(n,1),d%strata(n,1),d%samp_size(n,1),d%pop_size(n,1));d%weight=w;d%strata=1;d%samp_size=real(n,dp);d%pop_size=huge(1.0_dp);do i=1,n;d%cluster(i,1)=i;end do
  end subroutine pseudo_design
end module survey_glm
