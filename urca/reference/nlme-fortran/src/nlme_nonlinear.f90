! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_nonlinear
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp, pi_dp
  use nlme_status, only : NLME_SUCCESS, NLME_INVALID_ARGUMENT, NLME_DIMENSION_ERROR, &
       NLME_MAX_ITER, NLME_CALLBACK_ERROR
  use nlme_types, only : nonlinear_model, nonlinear_result, correlation_spec, variance_spec, &
       pd_spec, nlme_control, lme_result
  use nlme_covariance, only : build_residual_covariance, default_group_vector, group_level_indices
  use nlme_linalg, only : solve_spd, inverse_spd, find_group_indices, vector_norm2
  use nlme_lme, only : fit_lme
  implicit none
  private
  public :: fit_gnls, fit_nlme, numerical_model_jacobian
contains
  subroutine numerical_model_jacobian(model,theta,x,yhat,jacobian,status,step)
    procedure(nonlinear_model) :: model
    real(dp), intent(in) :: theta(:),x(:,:)
    real(dp), allocatable, intent(out) :: yhat(:),jacobian(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: step
    real(dp), allocatable :: tp(:),tm(:),yp(:),ym(:)
    real(dp) :: h
    integer :: j,n,p,st
    n=size(x,1)
    p=size(theta)
    h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)
    if(present(step))h=step
    allocate(yhat(n),jacobian(n,p),tp(p),tm(p),yp(n),ym(n))
    call model(theta,x,yhat,st)
    if(st/=NLME_SUCCESS .or. any(.not.ieee_is_finite(yhat)))then
    status=NLME_CALLBACK_ERROR
    return
    end if
    do j=1,p
      tp=theta
      tm=theta
      tp(j)=theta(j)+h*max(1.0_dp,abs(theta(j)))
      tm(j)=theta(j)-h*max(1.0_dp,abs(theta(j)))
      call model(tp,x,yp,st)
      if(st/=NLME_SUCCESS)then
      status=NLME_CALLBACK_ERROR
      return
      end if
      call model(tm,x,ym,st)
      if(st/=NLME_SUCCESS)then
      status=NLME_CALLBACK_ERROR
      return
      end if
      if(any(.not.ieee_is_finite(yp)) .or. any(.not.ieee_is_finite(ym)))then
      status=NLME_CALLBACK_ERROR
      return
      end if
      jacobian(:,j)=(yp-ym)/(tp(j)-tm(j))
    end do
    status=NLME_SUCCESS
  end subroutine numerical_model_jacobian

  subroutine fit_gnls(model,y,xdata,theta0,result,correlation,variance,time,group,var_covariate,var_group,&
       coordinates,control)
    procedure(nonlinear_model) :: model
    real(dp), intent(in) :: y(:),xdata(:,:),theta0(:)
    type(nonlinear_result), intent(out) :: result
    type(correlation_spec), intent(in), optional :: correlation
    type(variance_spec), intent(in), optional :: variance
    real(dp), intent(in), optional :: time(:),var_covariate(:),coordinates(:,:)
    integer, intent(in), optional :: group(:),var_group(:)
    type(nlme_control), intent(in), optional :: control
    type(correlation_spec) :: corr
    type(variance_spec) :: var
    type(nlme_control) :: ctl
    real(dp), allocatable :: theta(:),candidate(:),f(:),fc(:),jmat(:,:),v(:,:),vinvj(:,:),vinvr(:),&
         normal(:,:),rhs(:),delta(:),covp(:,:),r(:),rc(:),tvec(:),vcov(:),coords(:,:)
    integer, allocatable :: grp(:),vgrp(:)
    real(dp) :: lambda,obj,newobj,ldv,s2,logdetn
    integer :: n,p,it,status,st

    n=size(y)
    p=size(theta0)
    if(size(xdata,1)/=n .or. n<=p .or. p<1)then
    result%status=NLME_DIMENSION_ERROR
    return
    end if
    ctl=nlme_control()
    if(present(control))ctl=control
    corr=correlation_spec()
    if(present(correlation))corr=correlation
    var=variance_spec()
    if(present(variance))var=variance
    allocate(tvec(n),vcov(n),theta(p),candidate(p),fc(n))
    theta=theta0
    if(present(time))then
    if(size(time)/=n)then
    result%status=NLME_DIMENSION_ERROR
    return
    end if
    tvec=time
    else
    tvec=[(real(it,dp),it=1,n)]
    end if
    if(present(group))then
    if(size(group)/=n)then
    result%status=NLME_DIMENSION_ERROR
    return
    end if
    allocate(grp(n))
    grp=group
    else
    call default_group_vector(n,grp)
    end if
    if(present(var_covariate))then
      if(size(var_covariate)/=n)then
        result%status=NLME_DIMENSION_ERROR
        return
      end if
      vcov=var_covariate
    else
      vcov=abs(y)
    end if
    if(present(var_group))then
      if(size(var_group)/=n)then
        result%status=NLME_DIMENSION_ERROR
        return
      end if
      allocate(vgrp(n))
      vgrp=var_group
    else
      allocate(vgrp(n))
      vgrp=1
    end if
    if(present(coordinates))then
      if(size(coordinates,1)/=n)then
      result%status=NLME_DIMENSION_ERROR
      return
      end if
      allocate(coords(n,size(coordinates,2)))
      coords=coordinates
      call build_residual_covariance(corr,var,tvec,grp,vcov,vgrp,v,status,coords)
    else
      call build_residual_covariance(corr,var,tvec,grp,vcov,vgrp,v,status)
    end if
    if(status/=NLME_SUCCESS)then
    result%status=status
    return
    end if
    lambda=1.0e-3_dp
    result%converged=.false.
    do it=1,ctl%max_iter
      call numerical_model_jacobian(model,theta,xdata,f,jmat,status,ctl%fd_step)
      if(status/=NLME_SUCCESS)then
      result%status=status
      return
      end if
      r=y-f
      call solve_spd(v,jmat,vinvj,status,ldv)
      if(status/=NLME_SUCCESS)then
      result%status=status
      return
      end if
      call solve_spd(v,r,vinvr,status)
      if(status/=NLME_SUCCESS)then
      result%status=status
      return
      end if
      normal=matmul(transpose(jmat),vinvj)
      rhs=matmul(transpose(jmat),vinvr)
      obj=dot_product(r,vinvr)
      do st=1,p
      normal(st,st)=normal(st,st)+lambda*max(1.0_dp,normal(st,st))
      end do
      call solve_spd(normal,rhs,delta,status)
      if(status/=NLME_SUCCESS)then
      lambda=lambda*10.0_dp
      cycle
      end if
      candidate=theta+delta
      call model(candidate,xdata,fc,st)
      if(st/=NLME_SUCCESS .or. any(.not.ieee_is_finite(fc)))then
      lambda=lambda*10.0_dp
      cycle
      end if
      rc=y-fc
      call solve_spd(v,rc,vinvr,status)
      if(status/=NLME_SUCCESS)then
      result%status=status
      return
      end if
      newobj=dot_product(rc,vinvr)
      if(newobj<obj)then
        theta=candidate
        lambda=max(1.0e-12_dp,lambda/3.0_dp)
        if(vector_norm2(delta)<=ctl%tolerance*(1.0_dp+vector_norm2(theta)))then
        result%converged=.true.
        exit
        end if
      else
        lambda=min(1.0e12_dp,lambda*10.0_dp)
      end if
    end do
    call numerical_model_jacobian(model,theta,xdata,f,jmat,status,ctl%fd_step)
    if(status/=NLME_SUCCESS)then
    result%status=status
    return
    end if
    r=y-f
    call solve_spd(v,jmat,vinvj,status,ldv)
    if(status/=NLME_SUCCESS)then
    result%status=status
    return
    end if
    call solve_spd(v,r,vinvr,status)
    if(status/=NLME_SUCCESS)then
    result%status=status
    return
    end if
    normal=matmul(transpose(jmat),vinvj)
    call inverse_spd(normal,covp,status,logdetn)
    if(status/=NLME_SUCCESS)then
    result%status=status
    return
    end if
    s2=dot_product(r,vinvr)/real(n-p,dp)
    allocate(result%parameters(p),result%parameter_cov(p,p),result%fitted(n),result%residuals(n),&
         result%random_effects(0,0),result%group_levels(0))
    result%parameters=theta
    result%parameter_cov=s2*covp
    result%fitted=f
    result%residuals=r
    result%sigma=sqrt(s2)
    result%log_likelihood=-0.5_dp*(real(n,dp)*(log(2.0_dp*pi_dp)+log(s2))+ldv+dot_product(r,vinvr)/s2)
    result%iterations=min(it,ctl%max_iter)
    result%status=merge(NLME_SUCCESS,NLME_MAX_ITER,result%converged)
  end subroutine fit_gnls

  subroutine fit_nlme(model,y,xdata,group,theta0,random_index,result,random,correlation,variance,time,&
       var_covariate,var_group,coordinates,control,fixed_sigma)
    procedure(nonlinear_model) :: model
    real(dp), intent(in) :: y(:),xdata(:,:),theta0(:)
    integer, intent(in) :: group(:),random_index(:)
    type(nonlinear_result), intent(out) :: result
    type(pd_spec), intent(in), optional :: random
    type(correlation_spec), intent(in), optional :: correlation
    type(variance_spec), intent(in), optional :: variance
    real(dp), intent(in), optional :: time(:),var_covariate(:),coordinates(:,:)
    integer, intent(in), optional :: var_group(:)
    type(nlme_control), intent(in), optional :: control
    real(dp), intent(in), optional :: fixed_sigma
    type(nlme_control) :: ctl,lctl
    type(lme_result) :: linear_fit
    real(dp), allocatable :: theta(:),theta_old(:),bre(:,:),f(:),jall(:,:),xlin(:,:),zlin(:,:),ystar(:),&
         local_theta(:),flocal(:),jlocal(:,:),final_fit(:)
    integer, allocatable :: levels(:),gindex(:),idx(:)
    integer :: n,p,q,ng,g,j,status,it,st
    real(dp) :: change

    n=size(y)
    p=size(theta0)
    q=size(random_index)
    if(size(xdata,1)/=n .or. size(group)/=n .or. p<1 .or. q<1 .or. any(random_index<1) .or. &
       any(random_index>p))then
       result%status=NLME_DIMENSION_ERROR
       return
       end if
    ctl=nlme_control()
    if(present(control))ctl=control
    call group_level_indices(group,levels,gindex,status)
    if(status/=NLME_SUCCESS)then
    result%status=status
    return
    end if
    ng=size(levels)
    allocate(theta(p),theta_old(p),bre(ng,q),f(n),jall(n,p),xlin(n,p),zlin(n,q),ystar(n))
    theta=theta0
    bre=0.0_dp
    result%converged=.false.
    do it=1,ctl%max_outer
      theta_old=theta
      do g=1,ng
        call find_group_indices(group,levels(g),idx)
        allocate(local_theta(p))
        local_theta=theta
        do j=1,q
        local_theta(random_index(j))=local_theta(random_index(j))+bre(g,j)
        end do
        call numerical_model_jacobian(model,local_theta,xdata(idx,:),flocal,jlocal,status,ctl%fd_step)
        deallocate(local_theta)
        if(status/=NLME_SUCCESS)then
        result%status=status
        return
        end if
        f(idx)=flocal
        jall(idx,:)=jlocal
      end do
      xlin=jall
      do j=1,q
      zlin(:,j)=jall(:,random_index(j))
      end do
      ystar=y-f+sum(jall*spread(theta,1,n),dim=2)
      do g=1,ng
        call find_group_indices(group,levels(g),idx)
        ystar(idx)=ystar(idx)+matmul(zlin(idx,:),bre(g,:))
      end do
      lctl=ctl
      lctl%max_iter=max(80,ctl%max_iter/2)
      if(present(fixed_sigma))then
        call fit_lme(y=ystar,x=xlin,z=zlin,group=group,result=linear_fit,random=random,correlation=correlation,&
             variance=variance,time=time,var_covariate=var_covariate,var_group=var_group,coordinates=coordinates,&
             control=lctl,fixed_sigma=fixed_sigma)
      else
        call fit_lme(y=ystar,x=xlin,z=zlin,group=group,result=linear_fit,random=random,correlation=correlation,&
             variance=variance,time=time,var_covariate=var_covariate,var_group=var_group,coordinates=coordinates,control=lctl)
      end if
      if(linear_fit%status/=NLME_SUCCESS .and. linear_fit%status/=NLME_MAX_ITER)then
        result%status=linear_fit%status
        return
      end if
      theta=linear_fit%beta
      bre=linear_fit%random_effects
      change=maxval(abs(theta-theta_old))/max(1.0_dp,maxval(abs(theta_old)))
      if(change<=ctl%tolerance)then
      result%converged=.true.
      exit
      end if
    end do
    allocate(final_fit(n))
    final_fit=0.0_dp
    if(allocated(flocal))deallocate(flocal)
    do g=1,ng
      call find_group_indices(group,levels(g),idx)
      allocate(local_theta(p))
      local_theta=theta
      do j=1,q
      local_theta(random_index(j))=local_theta(random_index(j))+bre(g,j)
      end do
      allocate(flocal(size(idx)))
      call model(local_theta,xdata(idx,:),flocal,st)
      deallocate(local_theta)
      if(st/=NLME_SUCCESS)then
      result%status=NLME_CALLBACK_ERROR
      return
      end if
      final_fit(idx)=flocal
      deallocate(flocal)
    end do
    allocate(result%parameters(p),result%parameter_cov(p,p),result%random_effects(ng,q),result%group_levels(ng),&
         result%fitted(n),result%residuals(n))
    result%parameters=theta
    result%parameter_cov=linear_fit%beta_cov
    result%random_effects=bre
    result%group_levels=levels
    result%fitted=final_fit
    result%residuals=y-final_fit
    result%sigma=linear_fit%sigma
    result%log_likelihood=linear_fit%log_likelihood
    result%iterations=min(it,ctl%max_outer)
    result%status=merge(NLME_SUCCESS,NLME_MAX_ITER,result%converged)
  end subroutine fit_nlme
end module nlme_nonlinear
