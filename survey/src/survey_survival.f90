! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_survival
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, survival_curve_t, cox_result_t, aft_result_t, svystat_t
  use survey_estimators, only : svy_total
  use survey_taylor, only : svyrecvar
  use survival, only : coxph_fit, coxph_result, survreg_fit, survreg_loglik, surv_aft_result => aft_result
  implicit none
  private
  public :: svy_km, svy_coxph, svy_logrank, svy_survreg
contains

  subroutine svy_km(time,status,design,curve,se)
    real(dp), intent(in) :: time(:)
    integer, intent(in) :: status(:)
    type(survey_design_t), intent(in) :: design
    type(survival_curve_t), intent(out) :: curve
    logical, intent(in), optional :: se
    logical :: dose
    real(dp), allocatable :: et(:), risk(:), deaths(:), z(:,:),cov(:,:),dvn(:,:),dvy(:,:),dcv(:,:),dv(:,:)
    type(svystat_t) :: totals
    real(dp) :: sw,surv,h,riskj,deathj
    integer :: m,j,k,i
    if(size(time)/=design%n .or. size(status)/=design%n) error stop 'svy_km: size mismatch'
    dose=.false.;if(present(se))dose=se
    call unique_event_times(time,status,design%weight,et);m=size(et)
    if(.not.dose) then
      allocate(curve%time(m+1),curve%survival(m+1),curve%hazard(m+1),curve%variance(m+1))
      curve%time(1)=0;curve%survival(1)=1;curve%hazard(1)=0;curve%variance=0;sw=sum(design%weight);surv=1;h=0
      do j=1,m
        riskj=sum(design%weight,mask=time>=et(j));deathj=sum(design%weight,mask=(status/=0 .and. same_time_vec(time,et(j))))
        if(riskj>0)then;surv=max(0.0_dp,surv*(1-deathj/riskj));h=h+deathj/riskj;end if
        curve%time(j+1)=et(j);curve%survival(j+1)=surv;curve%hazard(j+1)=h
      end do
      return
    end if
    allocate(z(design%n,2*m));z=0
    do j=1,m
      do i=1,design%n
        if(status(i)/=0 .and. abs(time(i)-et(j))<=time_tol(et(j))) z(i,j)=1
        if(time(i)>=et(j)) z(i,m+j)=1
      end do
    end do
    totals=svy_total(z,design);allocate(risk(m),deaths(m));deaths=totals%estimate(1:m);risk=totals%estimate(m+1:2*m);cov=totals%variance
    allocate(dvn(m,m),dvy(m,m),dcv(m,m),dv(m,m));dvn=0;dvy=0;dcv=0
    do j=1,m;do k=1,m
      dvn(j,k)=cov(j,k)/(risk(j)*risk(k))
      dvy(j,k)=cov(m+j,m+k)*(deaths(j)/risk(j)**2)*(deaths(k)/risk(k)**2)
      dcv(j,k)=-cov(j,m+k)*(1/risk(j))*(deaths(k)/risk(k)**2)
    end do;end do
    dv=dvn+dvy+dcv+transpose(dcv)
    allocate(curve%time(m),curve%survival(m),curve%hazard(m),curve%variance(m));h=0
    do j=1,m
      h=h+deaths(j)/risk(j);curve%time(j)=et(j);curve%hazard(j)=h;curve%survival(j)=exp(-h)
      curve%variance(j)=sum(dv(1:j,1:j))
    end do
  end subroutine svy_km

  subroutine svy_coxph(time,status,x,design,result,method,maxiter,eps)
    real(dp),intent(in)::time(:),x(:,:);integer,intent(in)::status(:);type(survey_design_t),intent(in)::design;type(cox_result_t),intent(out)::result
    character(len=*),intent(in),optional::method;integer,intent(in),optional::maxiter;real(dp),intent(in),optional::eps
    type(coxph_result)::fit;real(dp),allocatable::score(:,:),infl(:,:);character(len=12)::meth
    meth='efron';if(present(method))meth=method
    call coxph_fit(time,status,x,fit,method=meth,weights=design%weight,maxiter=maxiter,eps=eps)
    call cox_score_residuals(time,status,x,design%weight,fit%coef,meth,score)
    allocate(infl(size(x,1),size(x,2)));infl=matmul(score,fit%var)
    allocate(result%coef(size(fit%coef)),result%naive_vcov(size(fit%coef),size(fit%coef)),result%vcov(size(fit%coef),size(fit%coef)))
    result%coef=fit%coef;result%naive_vcov=fit%var;result%vcov=svyrecvar(infl,design);result%loglik=fit%loglik;result%iterations=fit%iterations;result%rank=fit%rank;result%converged=fit%converged
  end subroutine svy_coxph

  subroutine svy_survreg(time,status,x,design,dist,result,maxiter,eps)
    real(dp),intent(in)::time(:),x(:,:)
    integer,intent(in)::status(:)
    type(survey_design_t),intent(in)::design
    character(len=*),intent(in)::dist
    type(aft_result_t),intent(out)::result
    integer,intent(in),optional::maxiter
    real(dp),intent(in),optional::eps
    type(surv_aft_result)::fit
    real(dp),allocatable::w(:),theta(:),score(:,:),infl(:,:),trial(:)
    real(dp)::h,fp,fm
    integer::n,p,q,i,j
    logical::fixed_scale
    n=size(time);p=size(x,2)
    if(size(status)/=n.or.size(x,1)/=n.or.design%n/=n)error stop 'svy_survreg: size mismatch'
    if(any(time<=0.0_dp).and.is_log_time_dist(dist)) error stop 'svy_survreg: log-time distributions require positive times'
    allocate(w(n));w=design%weight/(sum(design%weight)/real(n,dp))
    call survreg_fit(time,status,x,dist,fit,weights=w,maxiter=maxiter,eps=eps)
    fixed_scale=(trim(dist)=='exponential'.or.trim(dist)=='rayleigh');if(fixed_scale)then;q=p;else;q=p+1;end if
    allocate(theta(q));theta(1:p)=fit%coef;if(.not.fixed_scale)theta(q)=log(fit%scale)
    allocate(score(n,q),trial(q));score=0.0_dp
    do j=1,q
      h=1.0e-5_dp*max(1.0_dp,abs(theta(j)))
      do i=1,n
        trial=theta;trial(j)=theta(j)+h;fp=survreg_loglik(time(i:i),status(i:i),x(i:i,:),trial,dist,w(i:i))
        trial=theta;trial(j)=theta(j)-h;fm=survreg_loglik(time(i:i),status(i:i),x(i:i,:),trial,dist,w(i:i))
        score(i,j)=(fp-fm)/(2.0_dp*h)
      end do
    end do
    allocate(infl(n,q));infl=matmul(score,fit%var)
    allocate(result%coef(p),result%naive_vcov(q,q),result%vcov(q,q))
    result%coef=fit%coef;result%naive_vcov=fit%var;result%vcov=svyrecvar(infl,design);result%scale=fit%scale
    result%loglik=fit%loglik;result%iterations=fit%iterations;result%converged=fit%converged
  end subroutine svy_survreg

  logical pure function is_log_time_dist(dist) result(ans)
    character(len=*),intent(in)::dist
    ans=trim(dist)=='weibull'.or.trim(dist)=='exponential'.or.trim(dist)=='rayleigh'.or.trim(dist)=='lognormal'.or.trim(dist)=='loggaussian'.or.trim(dist)=='loglogistic'
  end function is_log_time_dist

  subroutine svy_logrank(time,status,group,design,score,se,z,rho,gamma)
    real(dp),intent(in)::time(:);integer,intent(in)::status(:),group(:);type(survey_design_t),intent(in)::design
    real(dp),intent(out)::score,se,z;real(dp),intent(in),optional::rho,gamma
    real(dp),allocatable::et(:),u(:,:),surv_before(:),risk(:),death(:),emean(:),wf(:);real(dp)::rh,ga,surv
    integer::i,j,m;type(svystat_t)::tot
    rh=0;if(present(rho))rh=rho;ga=0;if(present(gamma))ga=gamma
    if(size(time)/=design%n.or.size(status)/=design%n.or.size(group)/=design%n)error stop 'svy_logrank: size mismatch'
    call unique_event_times(time,status,design%weight,et);m=size(et);allocate(surv_before(m),risk(m),death(m),emean(m),wf(m));surv=1
    do j=1,m
      surv_before(j)=surv;risk(j)=sum(design%weight,mask=time>=et(j));death(j)=sum(design%weight,mask=(status/=0.and.same_time_vec(time,et(j))))
      if(risk(j)>0)then;emean(j)=sum(design%weight*real(group,dp),mask=time>=et(j))/risk(j);surv=surv*(1-death(j)/risk(j));else;emean(j)=0;end if
      wf(j)=surv_before(j)**rh*(1-surv_before(j))**ga
    end do
    allocate(u(design%n,1));u=0
    do i=1,design%n
      if(status(i)/=0)then;j=find_time(et,time(i));if(j>0)u(i,1)=u(i,1)+(real(group(i),dp)-emean(j))*wf(j);end if
      do j=1,m
        if(time(i)>=et(j).and.risk(j)>0)u(i,1)=u(i,1)-(real(group(i),dp)-emean(j))*wf(j)*death(j)/risk(j)
      end do
    end do
    tot=svy_total(u,design);score=tot%estimate(1);se=sqrt(max(0.0_dp,tot%variance(1,1)));if(se>0)then;z=score/se;else;z=0;end if
  end subroutine svy_logrank

  subroutine cox_score_residuals(time,status,x,w,beta,method,score)
    real(dp),intent(in)::time(:),x(:,:),w(:),beta(:);integer,intent(in)::status(:);character(len=*),intent(in)::method
    real(dp),allocatable,intent(out)::score(:,:)
    real(dp),allocatable::et(:),rr(:);real(dp)::denom,eventden,eventw,frac,part,coef
    integer::i,j,l,m,nd,p
    p=size(x,2);allocate(score(size(time),p));score=0;rr=w*exp(min(matmul(x,beta),700.0_dp));call unique_event_times(time,status,w,et)
    do j=1,size(et)
      denom=sum(rr,mask=time>=et(j));eventden=sum(rr,mask=(status/=0.and.same_time_vec(time,et(j))));eventw=sum(w,mask=(status/=0.and.same_time_vec(time,et(j))));nd=count(status/=0.and.same_time_vec(time,et(j)))
      do i=1,size(time);if(status(i)/=0.and.abs(time(i)-et(j))<=time_tol(et(j)))score(i,:)=score(i,:)+w(i)*x(i,:);end do
      if(index(method,'efron')==1.and.nd>1)then;m=nd;else;m=1;end if
      do l=0,m-1
        if(m>1)then;frac=real(l,dp)/real(m,dp);else;frac=0;end if;part=denom-frac*eventden;if(part<=0)cycle
        do i=1,size(time)
          if(time(i)>=et(j))then
            coef=(eventw/real(m,dp))*rr(i)/part
            if(status(i)/=0.and.abs(time(i)-et(j))<=time_tol(et(j)))coef=coef*(1-frac)
            score(i,:)=score(i,:)-coef*x(i,:)
          end if
        end do
      end do
    end do
  end subroutine cox_score_residuals

  subroutine unique_event_times(time,status,w,et)
    real(dp),intent(in)::time(:),w(:);integer,intent(in)::status(:);real(dp),allocatable,intent(out)::et(:);real(dp),allocatable::tmp(:);integer::i,m
    allocate(tmp(size(time)));m=0;do i=1,size(time);if(status(i)/=0.and.w(i)>0)then;if(m==0.or..not.any(abs(tmp(1:m)-time(i))<=time_tol(time(i))))then;m=m+1;tmp(m)=time(i);end if;end if;end do
    call sort_real(tmp,m);allocate(et(m));if(m>0)et=tmp(1:m)
  end subroutine unique_event_times
  subroutine sort_real(x,n);real(dp),intent(inout)::x(:);integer,intent(in)::n;integer::i,j;real(dp)::v;do i=2,n;v=x(i);j=i-1;do while(j>=1.and.x(j)>v);x(j+1)=x(j);j=j-1;end do;x(j+1)=v;end do;end subroutine sort_real
  pure real(dp) function time_tol(t) result(e);real(dp),intent(in)::t;e=64*epsilon(1.0_dp)*max(1.0_dp,abs(t));end function time_tol
  elemental logical function same_time_vec(t,v) result(ok);real(dp),intent(in)::t,v;ok=abs(t-v)<=time_tol(v);end function same_time_vec
  integer function find_time(et,t) result(k);real(dp),intent(in)::et(:),t;integer::i;k=0;do i=1,size(et);if(abs(et(i)-t)<=time_tol(t))then;k=i;return;end if;end do;end function find_time
end module survey_survival
