! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_profiles
  use acdm_kinds, only : dp, tiny_pos, huge_penalty, ACDM_SUCCESS, ACDM_BAD_INPUT
  use acdm_math, only : empirical_quantile
  use acdm_distributions, only : distribution_pdf, distribution_cdf, &
                                 distribution_quantile
  use acdm_models, only : acd_order, model_parameter_count, acd_loglik
  implicit none
  private

  type, public :: hazard_result
    real(dp), allocatable :: empirical_x(:), empirical_hazard(:)
    real(dp), allocatable :: implied_x(:), implied_hazard(:)
    integer :: status = ACDM_BAD_INPUT
  end type hazard_result

  type, public :: likelihood_profile_result
    real(dp), allocatable :: parameter1(:), parameter2(:)
    real(dp), allocatable :: loglik(:, :)
    integer :: status = ACDM_BAD_INPUT
  end type likelihood_profile_result

  public :: hazard_diagnostics, likelihood_profile

contains

  subroutine hazard_diagnostics(residuals, dist, dist_parameters, force_mean, &
                                result, breaks, implied_points)
    real(dp), intent(in) :: residuals(:), dist_parameters(:)
    integer, intent(in) :: dist
    logical, intent(in) :: force_mean
    type(hazard_result), intent(out) :: result
    integer, intent(in), optional :: breaks, implied_points
    integer :: b, m, i
    real(dp), allocatable :: edges(:)
    real(dp) :: probability, survivor, dx

    result%status = ACDM_BAD_INPUT
    if (size(residuals) < 10 .or. any(residuals <= 0.0_dp)) return
    b = 20
    if (present(breaks)) b = max(4, breaks)
    m = 201
    if (present(implied_points)) m = max(20, implied_points)
    allocate(edges(b), result%empirical_x(b-1), result%empirical_hazard(b-1))
    do i = 1, b
      probability = real(i-1, dp) / real(b, dp)
      edges(i) = empirical_quantile(residuals, probability)
    end do
    do i = 1, b-1
      probability = (real(i, dp)-0.5_dp) / real(b, dp)
      result%empirical_x(i) = empirical_quantile(residuals, probability)
      survivor = 1.0_dp - real(i, dp)/real(b, dp) + 0.5_dp/real(b,dp)
      dx = edges(i+1)-edges(i)
      if (dx <= tiny_pos) then
        result%empirical_hazard(i) = huge(1.0_dp)
      else
        result%empirical_hazard(i) = (1.0_dp/real(b,dp)) / survivor / dx
      end if
    end do
    allocate(result%implied_x(m), result%implied_hazard(m))
    do i=1,m
      probability = 0.0025_dp + real(i-1,dp)*(0.9925_dp-0.0025_dp)/real(m-1,dp)
      result%implied_x(i) = distribution_quantile(probability,dist,dist_parameters,force_mean)
      survivor = 1.0_dp-distribution_cdf(result%implied_x(i),dist,dist_parameters,force_mean)
      result%implied_hazard(i) = distribution_pdf(result%implied_x(i),dist,dist_parameters,force_mean) / &
                                 max(tiny_pos,survivor)
    end do
    result%status = ACDM_SUCCESS
  end subroutine hazard_diagnostics

  subroutine likelihood_profile(x, model, order, parameters, dist, &
                                dist_parameters, force_mean, parameter1_index, &
                                parameter1_values, result, parameter2_index, &
                                parameter2_values, breakpoints, exogenous, &
                                new_day)
    real(dp), intent(in) :: x(:), parameters(:), dist_parameters(:)
    integer, intent(in) :: model, dist, parameter1_index
    type(acd_order), intent(in) :: order
    logical, intent(in) :: force_mean
    real(dp), intent(in) :: parameter1_values(:)
    type(likelihood_profile_result), intent(out) :: result
    integer, intent(in), optional :: parameter2_index
    real(dp), intent(in), optional :: parameter2_values(:), breakpoints(:)
    real(dp), intent(in), optional :: exogenous(:, :)
    integer, intent(in), optional :: new_day(:)

    real(dp), allocatable :: theta(:), mu(:), residual(:)
    integer :: i,j,n2,idx2,st,base_n,n_exo,nb

    result%status=ACDM_BAD_INPUT
    n_exo=0;if(present(exogenous))n_exo=size(exogenous,2)
    nb=0;if(present(breakpoints))nb=size(breakpoints)
    base_n=model_parameter_count(model,order,nb)+n_exo
    if(size(parameters)/=base_n .or. parameter1_index<1 .or. parameter1_index>base_n) return
    if(size(x)<1 .or. size(parameter1_values)<1)return
    n2=1;idx2=0
    if(present(parameter2_index))then
      idx2=parameter2_index
      if(idx2<1.or.idx2>base_n.or..not.present(parameter2_values))return
      n2=size(parameter2_values);if(n2<1)return
    end if
    allocate(result%parameter1(size(parameter1_values)))
    result%parameter1=parameter1_values
    if(n2>1.or.idx2>0)then
      allocate(result%parameter2(n2));result%parameter2=parameter2_values
    else
      allocate(result%parameter2(0))
    end if
    allocate(result%loglik(size(parameter1_values),n2),theta(base_n),mu(size(x)),residual(size(x)))
    do i=1,size(parameter1_values)
      do j=1,n2
        theta=parameters;theta(parameter1_index)=parameter1_values(i)
        if(idx2>0)theta(idx2)=parameter2_values(j)
        call profile_ll(x,model,order,theta,dist,dist_parameters,force_mean,mu,residual,st, &
                        result%loglik(i,j),breakpoints,exogenous,new_day)
        if(st/=ACDM_SUCCESS)result%loglik(i,j)=-huge_penalty
      end do
    end do
    result%status=ACDM_SUCCESS
  end subroutine likelihood_profile

  subroutine profile_ll(x,model,order,theta,dist,dpar,force_mean,mu,residual,status,ll, &
                        breakpoints,exogenous,new_day)
    real(dp),intent(in)::x(:),theta(:),dpar(:)
    integer,intent(in)::model,dist
    type(acd_order),intent(in)::order
    logical,intent(in)::force_mean
    real(dp),intent(out)::mu(:),residual(:),ll
    integer,intent(out)::status
    real(dp),intent(in),optional::breakpoints(:),exogenous(:,:)
    integer,intent(in),optional::new_day(:)
    if(present(breakpoints))then
      if(present(exogenous))then
        if(present(new_day))then
          ll=acd_loglik(x,model,order,theta,dist,dpar,force_mean,mu,residual,status, &
                        breakpoints=breakpoints,exogenous=exogenous,new_day=new_day)
        else
          ll=acd_loglik(x,model,order,theta,dist,dpar,force_mean,mu,residual,status, &
                        breakpoints=breakpoints,exogenous=exogenous)
        end if
      else if(present(new_day))then
        ll=acd_loglik(x,model,order,theta,dist,dpar,force_mean,mu,residual,status, &
                      breakpoints=breakpoints,new_day=new_day)
      else
        ll=acd_loglik(x,model,order,theta,dist,dpar,force_mean,mu,residual,status,breakpoints=breakpoints)
      end if
    else if(present(exogenous))then
      if(present(new_day))then
        ll=acd_loglik(x,model,order,theta,dist,dpar,force_mean,mu,residual,status, &
                      exogenous=exogenous,new_day=new_day)
      else
        ll=acd_loglik(x,model,order,theta,dist,dpar,force_mean,mu,residual,status,exogenous=exogenous)
      end if
    else if(present(new_day))then
      ll=acd_loglik(x,model,order,theta,dist,dpar,force_mean,mu,residual,status,new_day=new_day)
    else
      ll=acd_loglik(x,model,order,theta,dist,dpar,force_mean,mu,residual,status)
    end if
  end subroutine profile_ll

end module acdm_profiles
