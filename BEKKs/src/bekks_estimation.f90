! SPDX-License-Identifier: MIT
module bekks_estimation
  use bekks_kinds, only: dp
  use bekks_types
  use bekks_model
  use bekks_linalg, only: general_inverse
  use bekks_rng, only: rng_state
  implicit none
  private
  public :: score_bekk, hessian_bekk, bhhh_fit, bekk_fit
  public :: qml_covariance, rmse_parameters

contains

  subroutine score_bekk(theta,data,model_type,asymmetric,signs,score,status,relative_step)
    real(dp), intent(in) :: theta(:),data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    real(dp), intent(in), optional :: signs(:)
    real(dp), allocatable, intent(out) :: score(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: relative_step
    real(dp), allocatable :: tp(:),tm(:),cp(:),cm(:),c0(:)
    real(dp) :: h,rel
    integer :: p,j,sp,sm,s0
    p=size(theta);allocate(score(size(data,1),p),tp(p),tm(p))
    rel=epsilon(1.0_dp)**(1.0_dp/3.0_dp)
    if(present(relative_step))rel=relative_step
    c0=log_likelihood_contributions(theta,data,model_type,asymmetric,signs,s0)
    if(s0/=bekk_ok)then;status=s0;score=0.0_dp;return;end if
    do j=1,p
      h=rel*max(1.0_dp,abs(theta(j)))
      tp=theta;tm=theta;tp(j)=tp(j)+h;tm(j)=tm(j)-h
      cp=log_likelihood_contributions(tp,data,model_type,asymmetric,signs,sp)
      cm=log_likelihood_contributions(tm,data,model_type,asymmetric,signs,sm)
      if(sp==bekk_ok .and. sm==bekk_ok)then
        score(:,j)=(cp-cm)/(2.0_dp*h)
      else if(sp==bekk_ok)then
        score(:,j)=(cp-c0)/h
      else if(sm==bekk_ok)then
        score(:,j)=(c0-cm)/h
      else
        score(:,j)=0.0_dp
      end if
    end do
    status=bekk_ok
  end subroutine score_bekk

  subroutine hessian_bekk(theta,data,model_type,asymmetric,signs,hessian,status,relative_step)
    real(dp), intent(in) :: theta(:),data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    real(dp), intent(in), optional :: signs(:)
    real(dp), allocatable, intent(out) :: hessian(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: relative_step
    real(dp), allocatable :: tpp(:),tpm(:),tmp(:),tmm(:)
    real(dp) :: hi,hj,fpp,fpm,fmp,fmm,f0,fi1,fi2,rel
    integer :: p,i,j,st
    p=size(theta);allocate(hessian(p,p),tpp(p),tpm(p),tmp(p),tmm(p));hessian=0.0_dp
    rel=epsilon(1.0_dp)**0.25_dp
    if(present(relative_step))rel=relative_step
    f0=log_likelihood(theta,data,model_type,asymmetric,signs,st)
    if(st/=bekk_ok)then;status=st;return;end if
    do i=1,p
      hi=rel*max(1.0_dp,abs(theta(i)))
      tpp=theta;tmm=theta;tpp(i)=tpp(i)+hi;tmm(i)=tmm(i)-hi
      fi1=log_likelihood(tpp,data,model_type,asymmetric,signs,st)
      fi2=log_likelihood(tmm,data,model_type,asymmetric,signs,st)
      if(fi1>-1.0e24_dp .and. fi2>-1.0e24_dp)hessian(i,i)=(fi1-2.0_dp*f0+fi2)/(hi*hi)
      do j=i+1,p
        hj=rel*max(1.0_dp,abs(theta(j)))
        tpp=theta;tpm=theta;tmp=theta;tmm=theta
        tpp(i)=tpp(i)+hi;tpp(j)=tpp(j)+hj
        tpm(i)=tpm(i)+hi;tpm(j)=tpm(j)-hj
        tmp(i)=tmp(i)-hi;tmp(j)=tmp(j)+hj
        tmm(i)=tmm(i)-hi;tmm(j)=tmm(j)-hj
        fpp=log_likelihood(tpp,data,model_type,asymmetric,signs)
        fpm=log_likelihood(tpm,data,model_type,asymmetric,signs)
        fmp=log_likelihood(tmp,data,model_type,asymmetric,signs)
        fmm=log_likelihood(tmm,data,model_type,asymmetric,signs)
        if(min(fpp,fpm,fmp,fmm)>-1.0e24_dp)then
          hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj);hessian(j,i)=hessian(i,j)
        end if
      end do
    end do
    status=bekk_ok
  end subroutine hessian_bekk

  subroutine qml_covariance(theta,data,model_type,asymmetric,signs,covariance,robust,status)
    real(dp), intent(in) :: theta(:),data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    real(dp), intent(in), optional :: signs(:)
    real(dp), allocatable, intent(out) :: covariance(:,:),robust(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: score(:,:),hess(:,:),opg(:,:),hinv(:,:),opginv(:,:)
    integer :: p,info,st
    p=size(theta);allocate(covariance(p,p),robust(p,p),opg(p,p),hinv(p,p),opginv(p,p))
    call score_bekk(theta,data,model_type,asymmetric,signs,score,st);if(st/=bekk_ok)then;status=st;return;end if
    opg=matmul(transpose(score),score)
    call general_inverse(opg,opginv,info);if(info/=0)then;status=bekk_linalg_failure;return;end if
    covariance=opginv
    call hessian_bekk(theta,data,model_type,asymmetric,signs,hess,st)
    if(st/=bekk_ok)then;robust=covariance;status=bekk_ok;return;end if
    call general_inverse(-hess,hinv,info)
    if(info/=0)then
      robust=covariance
    else
      robust=matmul(hinv,matmul(opg,hinv))
    end if
    status=bekk_ok
  end subroutine qml_covariance

  subroutine bhhh_fit(data,spec,result,max_iter,criterion,use_qml)
    real(dp), intent(in) :: data(:,:)
    type(bekk_spec_type), intent(in) :: spec
    type(bekk_fit_result), intent(out) :: result
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: criterion
    logical, intent(in), optional :: use_qml
    real(dp), parameter :: steps(21)=[9.9_dp,9.0_dp,8.0_dp,7.0_dp,6.0_dp,5.0_dp,4.0_dp,3.0_dp,2.0_dp,1.0_dp, &
      0.5_dp,0.25_dp,0.1_dp,0.01_dp,0.005_dp,0.001_dp,0.0005_dp,0.0001_dp,0.00005_dp,0.00001_dp,0.0_dp]
    type(bekk_parameters) :: par0
    real(dp), allocatable :: theta(:),candidate(:),best_theta(:),score(:,:),opg(:,:),opginv(:,:),gradient(:)
    real(dp), allocatable :: covariance(:,:),robust(:,:),signs(:)
    real(dp) :: ll,best_ll,new_ll,crit,scale,expected
    integer :: it,maxit,j,info,st,p,n
    logical :: qml
    n=size(data,2)
    maxit=50
    if(present(max_iter))maxit=max_iter
    crit=1.0e-9_dp
    if(present(criterion))crit=criterion
    qml=.false.
    if(present(use_qml))qml=use_qml
    allocate(signs(n));signs=-1.0_dp;if(allocated(spec%signs))signs=spec%signs
    if(allocated(spec%initial_theta))then
      theta=spec%initial_theta
    else
      call initial_parameters(data,spec%model_type,spec%asymmetric,par0,st)
      if(st/=bekk_ok)then;result%status=st;return;end if
      theta=pack_parameters(par0)
    end if
    p=size(theta);allocate(candidate(p),best_theta(p),opg(p,p),opginv(p,p),gradient(p),result%likelihood_path(maxit+1))
    result%likelihood_path=0.0_dp;ll=log_likelihood(theta,data,spec%model_type,spec%asymmetric,signs,st)
    if(st/=bekk_ok)then;result%status=st;return;end if
    result%likelihood_path(1)=ll;result%converged=.false.
    do it=1,maxit
      call score_bekk(theta,data,spec%model_type,spec%asymmetric,signs,score,st)
      if(st/=bekk_ok)exit
      opg=matmul(transpose(score),score)
      call general_inverse(opg,opginv,info)
      if(info/=0)then;st=bekk_linalg_failure;exit;end if
      gradient=sum(score,dim=1);best_ll=ll;best_theta=theta
      do j=1,size(steps)
        scale=0.1_dp*steps(j)
        candidate=theta+scale*matmul(opginv,gradient)
        new_ll=log_likelihood(candidate,data,spec%model_type,spec%asymmetric,signs)
        if(new_ll>best_ll)then;best_ll=new_ll;best_theta=candidate;end if
      end do
      result%likelihood_path(it+1)=best_ll
      if(best_ll<ll)exit
      theta=best_theta
      if((best_ll-ll)**2/max(1.0_dp,abs(ll))<crit)then
        ll=best_ll;result%converged=.true.;exit
      end if
      ll=best_ll
    end do
    result%iterations=min(it,maxit);result%theta=theta;result%log_likelihood=ll;result%spec=spec;result%data=data;result%signs=signs
    call unpack_parameters(theta,n,spec%model_type,spec%asymmetric,result%parameters,st)
    expected=expected_indicator_value(data,signs);result%expected_indicator=expected
    result%stationary=valid_parameters(result%parameters,expected)
    call filter_bekk(theta,data,spec%model_type,spec%asymmetric,signs,result%h,result%residuals,st)
    if(st/=bekk_ok)then;result%status=st;return;end if
    call score_bekk(theta,data,spec%model_type,spec%asymmetric,signs,result%score,st)
    call qml_covariance(theta,data,spec%model_type,spec%asymmetric,signs,covariance,robust,st)
    if(st==bekk_ok)then
      result%covariance=covariance;result%robust_covariance=robust
      allocate(result%standard_error(p),result%t_value(p))
      if(qml)then
        result%standard_error=sqrt(max(diagonal(robust),0.0_dp))
      else
        result%standard_error=sqrt(max(diagonal(covariance),0.0_dp))
      end if
      where(result%standard_error>0.0_dp)
        result%t_value=theta/result%standard_error
      elsewhere
        result%t_value=0.0_dp
      end where
    end if
    result%aic=-2.0_dp*ll+2.0_dp*real(p,dp)
    result%bic=-2.0_dp*ll+log(real(size(data,1),dp))*real(p,dp)
    result%status=merge(bekk_ok,bekk_no_convergence,result%converged)
  contains
    function diagonal(a) result(d)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: d(min(size(a,1),size(a,2)))
      integer :: k
      do k=1,size(d);d(k)=a(k,k);end do
    end function diagonal
  end subroutine bhhh_fit

  subroutine bekk_fit(spec,data,result,max_iter,criterion,use_qml)
    type(bekk_spec_type), intent(in) :: spec
    real(dp), intent(in) :: data(:,:)
    type(bekk_fit_result), intent(out) :: result
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: criterion
    logical, intent(in), optional :: use_qml
    call bhhh_fit(data,spec,result,max_iter,criterion,use_qml)
  end subroutine bekk_fit

  real(dp) function rmse_parameters(theta_estimated,theta_true) result(v)
    real(dp), intent(in) :: theta_estimated(:),theta_true(:)
    real(dp), allocatable :: scale(:)
    allocate(scale(size(theta_true)));scale=max(abs(theta_true),sqrt(epsilon(1.0_dp)))
    v=sum(abs((theta_estimated-theta_true)/scale))/real(size(theta_true),dp)
  end function rmse_parameters

end module bekks_estimation
