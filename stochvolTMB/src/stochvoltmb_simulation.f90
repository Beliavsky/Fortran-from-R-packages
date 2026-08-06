! SPDX-License-Identifier: GPL-3.0-only
module stochvoltmb_simulation
  use stochvoltmb_kinds, only : dp, pi, tiny_dp
  use stochvoltmb_status, only : sv_ok, sv_invalid_argument, sv_singular
  use stochvoltmb_math, only : finite_real, normal_cdf, normal_quantile, &
                              student_t_cdf, skew_normal_cdf, quantile_type7
  use stochvoltmb_linalg, only : cholesky_lower
  use stochvoltmb_rng, only : sv_rng_state
  use stochvoltmb_types
  implicit none
  private

  public :: sim_sv, simulate_parameters, predict_sv, summarize_prediction
  public :: standardized_residuals, one_step_residuals

contains

  real(dp) function skew_standard_draw(rng, alpha) result(z)
    type(sv_rng_state), intent(inout) :: rng
    real(dp), intent(in) :: alpha
    real(dp) :: delta, omega, xi, u, v
    delta=alpha/sqrt(1.0_dp+alpha*alpha)
    omega=1.0_dp/sqrt(max(tiny_dp,1.0_dp-2.0_dp*delta*delta/pi))
    xi=-omega*delta*sqrt(2.0_dp/pi)
    u=rng%normal(); v=rng%normal()
    z=xi+omega*(delta*abs(u)+sqrt(max(0.0_dp,1.0_dp-delta*delta))*v)
  end function skew_standard_draw

  subroutine sim_sv(params,nobs,model,rng,sim)
    type(sv_parameters), intent(in) :: params
    integer, intent(in) :: nobs, model
    type(sv_rng_state), intent(inout) :: rng
    type(sv_simulation), intent(out) :: sim
    real(dp), allocatable :: hall(:)
    real(dp) :: eta, eps, z
    integer :: nlatent, i
    sim%model=model; sim%params=params
    if (nobs<2 .or. params%sigma_y<0.0_dp .or. params%sigma_h<0.0_dp .or. &
        abs(params%phi)>=1.0_dp) then
      sim%status=sv_invalid_argument; sim%message='invalid simulation parameters'; return
    end if
    if (model==sv_student_t .and. params%df<=2.0_dp) then
      sim%status=sv_invalid_argument; sim%message='df must be greater than 2'; return
    end if
    if ((model==sv_leverage .or. model==sv_skew_gaussian_leverage) .and. abs(params%rho)>=1.0_dp) then
      sim%status=sv_invalid_argument; sim%message='rho must be between -1 and 1'; return
    end if
    nlatent=nobs
    if (model==sv_leverage .or. model==sv_skew_gaussian_leverage) nlatent=nobs+1
    allocate(hall(nlatent),sim%y(nobs),sim%h(nobs))
    hall(1)=params%sigma_h/sqrt(1.0_dp-params%phi*params%phi)*rng%normal()
    do i=2,nlatent
      hall(i)=params%phi*hall(i-1)+params%sigma_h*rng%normal()
    end do
    select case(model)
    case(sv_gaussian)
      do i=1,nobs
        sim%y(i)=params%sigma_y*exp(0.5_dp*hall(i))*rng%normal()
      end do
    case(sv_student_t)
      do i=1,nobs
        sim%y(i)=params%sigma_y*exp(0.5_dp*hall(i))*sqrt((params%df-2.0_dp)/params%df)* &
          rng%student_t(params%df)
      end do
    case(sv_skew_gaussian)
      do i=1,nobs
        sim%y(i)=params%sigma_y*exp(0.5_dp*hall(i))*skew_standard_draw(rng,params%alpha)
      end do
    case(sv_leverage)
      do i=1,nobs
        eta=(hall(i+1)-params%phi*hall(i))/params%sigma_h
        eps=rng%normal()
        sim%y(i)=params%sigma_y*exp(0.5_dp*hall(i))* &
          (params%rho*eta+sqrt(1.0_dp-params%rho*params%rho)*eps)
      end do
    case(sv_skew_gaussian_leverage)
      do i=1,nobs
        eta=(hall(i+1)-params%phi*hall(i))/params%sigma_h
        z=skew_standard_draw(rng,params%alpha)
        sim%y(i)=params%sigma_y*exp(0.5_dp*hall(i))* &
          (params%rho*eta+sqrt(1.0_dp-params%rho*params%rho)*z)
      end do
    case default
      sim%status=sv_invalid_argument; sim%message='unknown model'; return
    end select
    sim%h=hall(1:nobs)
    if (any(.not. finite_real(sim%y))) then
      sim%status=sv_invalid_argument; sim%message='nonfinite simulated observations'; return
    end if
    sim%status=sv_ok; sim%message='ok'
  end subroutine sim_sv

  subroutine draw_theta_parameters(fit,nsim,rng,theta_draws,info)
    type(sv_fit_result), intent(in) :: fit
    integer, intent(in) :: nsim
    type(sv_rng_state), intent(inout) :: rng
    real(dp), allocatable, intent(out) :: theta_draws(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:), z(:)
    integer :: p, i, j, stat
    p=size(fit%theta); info=sv_ok
    if (nsim<1 .or. .not. allocated(fit%theta_cov)) then
      info=sv_invalid_argument; return
    end if
    allocate(theta_draws(p,nsim),l(p,p),z(p))
    call cholesky_lower(fit%theta_cov,l,stat,jitter=1.0e-10_dp)
    if (stat/=sv_ok) then
      info=sv_singular; return
    end if
    do j=1,nsim
      do i=1,p
        z(i)=rng%normal()
      end do
      theta_draws(:,j)=fit%theta+matmul(l,z)
    end do
  end subroutine draw_theta_parameters

  subroutine simulate_parameters(fit,nsim,rng,draws,info)
    type(sv_fit_result), intent(in) :: fit
    integer, intent(in) :: nsim
    type(sv_rng_state), intent(inout) :: rng
    real(dp), allocatable, intent(out) :: draws(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: theta_draws(:,:)
    type(sv_parameters) :: par
    integer :: p,j,stat
    p=size(fit%theta)
    call draw_theta_parameters(fit,nsim,rng,theta_draws,info)
    if (info/=sv_ok) return
    allocate(draws(p,nsim))
    do j=1,nsim
      call local_theta_to_params(theta_draws(:,j),fit%model,par,stat)
      if (stat/=sv_ok) then
        info=stat; return
      end if
      draws(1,j)=par%sigma_y
      draws(2,j)=par%sigma_h
      draws(3,j)=par%phi
      select case(fit%model)
      case(sv_student_t)
        draws(4,j)=par%df
      case(sv_skew_gaussian)
        draws(4,j)=par%alpha
      case(sv_leverage)
        draws(4,j)=par%rho
      case(sv_skew_gaussian_leverage)
        draws(4,j)=par%alpha
        draws(5,j)=par%rho
      end select
    end do
    info=sv_ok
  contains
    subroutine local_theta_to_params(theta,model,params,istat)
      use stochvoltmb_model, only : theta_to_parameters
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: model
      type(sv_parameters), intent(out) :: params
      integer, intent(out) :: istat
      call theta_to_parameters(theta,model,params,istat)
    end subroutine local_theta_to_params
  end subroutine simulate_parameters

  subroutine predict_sv(fit,steps,nsim,rng,pred,include_parameters)
    type(sv_fit_result), intent(in) :: fit
    integer, intent(in) :: steps,nsim
    type(sv_rng_state), intent(inout) :: rng
    type(sv_prediction), intent(out) :: pred
    logical, intent(in), optional :: include_parameters
    logical :: include_par
    real(dp), allocatable :: theta_draws(:,:), hwork(:)
    type(sv_parameters) :: par
    real(dp) :: hprev, hnext, eta, z, eps
    integer :: j,i,p,stat,extra
    include_par=.true.; if (present(include_parameters)) include_par=include_parameters
    if (steps<1 .or. nsim<1 .or. .not. allocated(fit%h)) then
      pred%status=sv_invalid_argument; pred%message='invalid prediction request'; return
    end if
    pred%steps=steps; pred%nsim=nsim
    allocate(pred%y(steps,nsim),pred%h(steps,nsim),pred%h_exp(steps,nsim))
    p=size(fit%theta)
    if (include_par .and. allocated(fit%theta_cov)) then
      call draw_theta_parameters(fit,nsim,rng,theta_draws,stat)
      if (stat/=sv_ok) include_par=.false.
    end if
    extra=0
    if (fit%model==sv_leverage .or. fit%model==sv_skew_gaussian_leverage) extra=1
    allocate(hwork(steps+extra))
    do j=1,nsim
      if (include_par) then
        call local_theta_to_params(theta_draws(:,j),fit%model,par,stat)
        if (stat/=sv_ok) par=fit%params
      else
        par=fit%params
      end if
      hprev=fit%h(size(fit%h))+fit%h_se(size(fit%h))*rng%normal()
      do i=1,steps+extra
        hnext=par%phi*hprev+par%sigma_h*rng%normal()
        hwork(i)=hnext
        hprev=hnext
      end do
      do i=1,steps
        select case(fit%model)
        case(sv_gaussian)
          pred%y(i,j)=par%sigma_y*exp(0.5_dp*hwork(i))*rng%normal()
        case(sv_student_t)
          pred%y(i,j)=par%sigma_y*exp(0.5_dp*hwork(i))*sqrt((par%df-2.0_dp)/par%df)*rng%student_t(par%df)
        case(sv_skew_gaussian)
          pred%y(i,j)=par%sigma_y*exp(0.5_dp*hwork(i))*skew_standard_draw(rng,par%alpha)
        case(sv_leverage)
          eta=(hwork(i+1)-par%phi*hwork(i))/par%sigma_h
          eps=rng%normal()
          pred%y(i,j)=par%sigma_y*exp(0.5_dp*hwork(i))* &
            (par%rho*eta+sqrt(1.0_dp-par%rho*par%rho)*eps)
        case(sv_skew_gaussian_leverage)
          eta=(hwork(i+1)-par%phi*hwork(i))/par%sigma_h
          z=skew_standard_draw(rng,par%alpha)
          pred%y(i,j)=par%sigma_y*exp(0.5_dp*hwork(i))* &
            (par%rho*eta+sqrt(1.0_dp-par%rho*par%rho)*z)
        end select
        pred%h(i,j)=hwork(i)
        pred%h_exp(i,j)=par%sigma_y*exp(0.5_dp*hwork(i))
      end do
    end do
    pred%status=sv_ok; pred%message='ok'
  contains
    subroutine local_theta_to_params(theta,model,params,istat)
      use stochvoltmb_model, only : theta_to_parameters
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: model
      type(sv_parameters), intent(out) :: params
      integer, intent(out) :: istat
      call theta_to_parameters(theta,model,params,istat)
    end subroutine local_theta_to_params
  end subroutine predict_sv

  subroutine summarize_prediction(pred,probabilities,summary,include_mean)
    type(sv_prediction), intent(in) :: pred
    real(dp), intent(in) :: probabilities(:)
    type(sv_prediction_summary), intent(out) :: summary
    logical, intent(in), optional :: include_mean
    logical :: means
    integer :: i,j
    means=.true.; if (present(include_mean)) means=include_mean
    if (pred%status/=sv_ok .or. size(probabilities)<1 .or. any(probabilities<0.0_dp) .or. &
        any(probabilities>1.0_dp)) then
      summary%status=sv_invalid_argument; summary%message='invalid prediction summary request'; return
    end if
    allocate(summary%probabilities(size(probabilities)))
    summary%probabilities=probabilities
    allocate(summary%y_quantiles(pred%steps,size(probabilities)))
    allocate(summary%h_quantiles(pred%steps,size(probabilities)))
    allocate(summary%h_exp_quantiles(pred%steps,size(probabilities)))
    do i=1,pred%steps
      do j=1,size(probabilities)
        summary%y_quantiles(i,j)=quantile_type7(pred%y(i,:),probabilities(j))
        summary%h_quantiles(i,j)=quantile_type7(pred%h(i,:),probabilities(j))
        summary%h_exp_quantiles(i,j)=quantile_type7(pred%h_exp(i,:),probabilities(j))
      end do
    end do
    if (means) then
      allocate(summary%y_mean(pred%steps),summary%h_mean(pred%steps),summary%h_exp_mean(pred%steps))
      do i=1,pred%steps
        summary%y_mean(i)=sum(pred%y(i,:))/real(pred%nsim,dp)
        summary%h_mean(i)=sum(pred%h(i,:))/real(pred%nsim,dp)
        summary%h_exp_mean(i)=sum(pred%h_exp(i,:))/real(pred%nsim,dp)
      end do
    end if
    summary%status=sv_ok; summary%message='ok'
  end subroutine summarize_prediction

  subroutine standardized_residuals(y,fit,residuals,info)
    real(dp), intent(in) :: y(:)
    type(sv_fit_result), intent(in) :: fit
    real(dp), intent(out) :: residuals(:)
    integer, intent(out) :: info
    real(dp) :: scale,eta,mu,sd,delta,omega,xi,p,z
    integer :: i,n
    n=size(y); info=sv_ok
    if (size(residuals)/=n .or. .not. allocated(fit%h) .or. size(fit%h)/=n) then
      info=sv_invalid_argument; return
    end if
    do i=1,n
      scale=fit%params%sigma_y*exp(0.5_dp*fit%h(i))
      select case(fit%model)
      case(sv_gaussian)
        p=normal_cdf(y(i)/scale)
      case(sv_student_t)
        z=y(i)/(scale*sqrt((fit%params%df-2.0_dp)/fit%params%df))
        p=student_t_cdf(z,fit%params%df)
      case(sv_skew_gaussian)
        delta=fit%params%alpha/sqrt(1.0_dp+fit%params%alpha**2)
        omega=scale/sqrt(1.0_dp-2.0_dp*delta*delta/pi)
        xi=-omega*delta*sqrt(2.0_dp/pi)
        p=skew_normal_cdf(y(i),xi,omega,fit%params%alpha)
      case(sv_leverage)
        if (i<n) then
          eta=(fit%h(i+1)-fit%params%phi*fit%h(i))/fit%params%sigma_h
        else
          eta=0.0_dp
        end if
        mu=scale*fit%params%rho*eta
        sd=scale*sqrt(1.0_dp-fit%params%rho**2)
        p=normal_cdf((y(i)-mu)/sd)
      case(sv_skew_gaussian_leverage)
        if (i<n) then
          eta=(fit%h(i+1)-fit%params%phi*fit%h(i))/fit%params%sigma_h
        else
          eta=0.0_dp
        end if
        mu=scale*fit%params%rho*eta
        sd=scale*sqrt(1.0_dp-fit%params%rho**2)
        delta=fit%params%alpha/sqrt(1.0_dp+fit%params%alpha**2)
        omega=sd/sqrt(1.0_dp-2.0_dp*delta*delta/pi)
        xi=mu-omega*delta*sqrt(2.0_dp/pi)
        p=skew_normal_cdf(y(i),xi,omega,fit%params%alpha)
      end select
      p=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,p))
      residuals(i)=normal_quantile(p)
    end do
  end subroutine standardized_residuals

  subroutine one_step_residuals(y,fit,residuals,info)
    real(dp), intent(in) :: y(:)
    type(sv_fit_result), intent(in) :: fit
    real(dp), intent(out) :: residuals(:)
    integer, intent(out) :: info
    call standardized_residuals(y,fit,residuals,info)
  end subroutine one_step_residuals

end module stochvoltmb_simulation
