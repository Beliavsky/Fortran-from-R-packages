! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_nls
  use survey_kinds, only : dp
  use survey_types, only : survey_design_t, nls_result_t
  use survey_linalg, only : sym_pinv
  use survey_taylor, only : svyrecvar
  implicit none
  private
  public :: svy_nls, observation_model
  abstract interface
    function observation_model(theta,i) result(mu)
      import dp
      real(dp),intent(in)::theta(:)
      integer,intent(in)::i
      real(dp)::mu
    end function observation_model
  end interface
contains
  subroutine svy_nls(model,y,design,start,result,prior_weight,maxit,tol,variance_power,precision_iterations)
    procedure(observation_model)::model
    real(dp),intent(in)::y(:),start(:)
    type(survey_design_t),intent(in)::design
    type(nls_result_t),intent(out)::result
    real(dp),intent(in),optional::prior_weight(:),tol,variance_power
    integer,intent(in),optional::maxit,precision_iterations
    real(dp),allocatable::theta(:),newtheta(:),mu(:),res(:),jac(:,:),w(:),pw(:),a(:,:),ainv(:,:),rhs(:),score(:,:),infl(:,:)
    real(dp)::eps,oldrss,newrss,alpha,disp,vpow
    integer::n,p,mi,it,i,j,k,rank,info,piter,pit
    n=size(y);p=size(start);if(design%n/=n)error stop 'svy_nls: size mismatch'
    mi=100;if(present(maxit))mi=maxit;eps=1e-8_dp;if(present(tol))eps=tol;vpow=0.0_dp;if(present(variance_power))vpow=variance_power
    piter=0;if(present(precision_iterations))piter=precision_iterations;it=0
    allocate(theta(p),newtheta(p),mu(n),res(n),jac(n,p),w(n),pw(n),a(p,p),ainv(p,p),rhs(p),score(n,p),infl(n,p));theta=start
    pw=1.0_dp;if(present(prior_weight))then;if(size(prior_weight)/=n)error stop 'svy_nls: prior_weight size';pw=prior_weight;end if
    w=pw*design%weight/(sum(design%weight)/real(n,dp));disp=1.0_dp
    do pit=0,piter
      call eval_model(model,theta,mu);oldrss=sum(w*(y-mu)**2)
      do it=1,mi
        call numeric_jacobian(model,theta,mu,jac)
        a=0.0_dp;rhs=0.0_dp;res=y-mu
        do i=1,n;do j=1,p;rhs(j)=rhs(j)+w(i)*jac(i,j)*res(i);do k=1,p;a(j,k)=a(j,k)+w(i)*jac(i,j)*jac(i,k);end do;end do;end do
        call sym_pinv(a,ainv,rank,info=info);newtheta=theta+matmul(ainv,rhs);alpha=1.0_dp
        do
          call eval_model(model,theta+alpha*(newtheta-theta),mu);newrss=sum(w*(y-mu)**2)
          if(newrss<=oldrss.or.alpha<1e-7_dp)exit;alpha=alpha/2.0_dp
        end do
        newtheta=theta+alpha*(newtheta-theta)
        if(maxval(abs(newtheta-theta)/(1+abs(theta)))<eps)then;theta=newtheta;exit;end if
        theta=newtheta;oldrss=newrss
      end do
      call eval_model(model,theta,mu);res=y-mu
      if(pit<piter.and.abs(vpow)>tiny(1.0_dp))then
        disp=sum((res*res)**max(0.0_dp,vpow/2.0_dp))/max(sum(abs(mu)**vpow),tiny(1.0_dp));w=pw*design%weight/(sum(design%weight)/real(n,dp))/max(disp*abs(mu)**vpow,tiny(1.0_dp))
      end if
    end do
    call eval_model(model,theta,mu);res=y-mu;call numeric_jacobian(model,theta,mu,jac);a=0.0_dp
    do i=1,n;do j=1,p;do k=1,p;a(j,k)=a(j,k)+w(i)*jac(i,j)*jac(i,k);end do;end do;end do
    call sym_pinv(a,ainv,rank,info=info)
    do i=1,n;do j=1,p;score(i,j)=w(i)*res(i)*jac(i,j);end do;end do;infl=matmul(score,ainv)
    allocate(result%coef(p),result%vcov(p,p),result%naive_vcov(p,p),result%fitted(n),result%residual(n));result%coef=theta;result%fitted=mu;result%residual=res
    result%rss=sum(w*res*res);result%vcov=svyrecvar(infl,design);result%naive_vcov=ainv*result%rss/max(1.0_dp,real(count(w>0)-rank,dp));result%iterations=min(it,mi);result%converged=(it<=mi)
  end subroutine svy_nls

  subroutine eval_model(model,theta,mu)
    procedure(observation_model)::model;real(dp),intent(in)::theta(:);real(dp),intent(out)::mu(:);integer::i
    do i=1,size(mu);mu(i)=model(theta,i);end do
  end subroutine eval_model
  subroutine numeric_jacobian(model,theta,mu,jac)
    procedure(observation_model)::model;real(dp),intent(in)::theta(:),mu(:);real(dp),intent(out)::jac(:,:)
    real(dp),allocatable::trial(:);real(dp)::h,fp,fm;integer::i,j
    if(size(mu)/=size(jac,1))error stop 'numeric_jacobian: shape';allocate(trial(size(theta)))
    do j=1,size(theta);h=1e-5_dp*max(1.0_dp,abs(theta(j)));do i=1,size(mu);trial=theta;trial(j)=theta(j)+h;fp=model(trial,i);trial(j)=theta(j)-h;fm=model(trial,i);jac(i,j)=(fp-fm)/(2*h);end do;end do
  end subroutine numeric_jacobian
end module survey_nls
