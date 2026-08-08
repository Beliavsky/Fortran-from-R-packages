module nfcp_mle
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nfcp_types, only : dp, nfcp_model_t, nfcp_parameterization_t, nfcp_filter_result_t, &
                         nfcp_mle_control_t, nfcp_mle_result_t, nfcp_ok, nfcp_invalid_input, &
                         nfcp_not_converged
  use nfcp_math, only : nfcp_rng_t, finite_difference_hessian, inverse_spd
  use nfcp_parameters, only : unpack_parameters, pack_parameters
  use nfcp_kalman, only : nfcp_kalman_filter
  use lbfgsb3_mod, only : lbfgsb_control_t, lbfgsb_result_t, lbfgsb_minimize_fd
  implicit none
  private
  public :: nfcp_fit_mle

  type :: mle_data_t
    type(nfcp_parameterization_t) :: spec
    real(dp), allocatable :: log_futures(:,:)
    real(dp), allocatable :: futures_ttm(:,:)
    real(dp), allocatable :: me_ttm(:)
    real(dp), allocatable :: fixed_initial_state(:)
    real(dp) :: dt = 0.0_dp
    real(dp) :: seasonal_trend = 0.0_dp
    logical :: has_me_ttm = .false.
    logical :: has_fixed_initial = .false.
    integer :: evaluations = 0
  end type mle_data_t

contains

  subroutine nfcp_fit_mle(spec, log_futures, dt, futures_ttm, lower, upper, result, &
                          initial_state, me_ttm, seasonal_trend, starting_values, control)
    type(nfcp_parameterization_t), intent(in) :: spec
    real(dp), intent(in) :: log_futures(:,:), dt, futures_ttm(:,:), lower(:), upper(:)
    type(nfcp_mle_result_t), intent(out) :: result
    real(dp), intent(in), optional :: initial_state(:), me_ttm(:), seasonal_trend, starting_values(:)
    type(nfcp_mle_control_t), intent(in), optional :: control

    type(nfcp_mle_control_t) :: ctrl
    type(mle_data_t) :: data
    type(nfcp_rng_t) :: rng
    type(lbfgsb_control_t) :: local_ctrl
    type(lbfgsb_result_t) :: local_result
    type(nfcp_model_t) :: fitted_model
    type(nfcp_filter_result_t) :: filter
    real(dp), allocatable :: pop(:,:), value(:), trial(:), mutant(:), best(:), theta(:), x0(:)
    real(dp), allocatable :: hessian(:,:), inverse_hessian(:,:)
    integer :: n, npop, i, j, g, a, b, c, forced, status, no_improve, best_index
    real(dp) :: trial_value, best_value, previous_best

    result%status=nfcp_invalid_input;result%message='invalid input';result%converged=.false.
    ctrl=nfcp_mle_control_t();if(present(control))ctrl=control
    n=spec%count()
    if(n<1 .or. size(lower)/=n .or. size(upper)/=n .or. any(lower>=upper) .or. &
       any(shape(log_futures)/=shape(futures_ttm)) .or. dt<=0.0_dp) return
    if(present(starting_values)) then
      if(size(starting_values)/=n) then;result%message='starting_values has wrong length';return;end if
    end if
    if(.not.spec%estimate_initial_state .and. present(initial_state)) then
      if(size(initial_state)/=spec%n_factors) then;result%message='initial_state has wrong length';return;end if
    end if

    data%spec=spec;data%log_futures=log_futures;data%futures_ttm=futures_ttm;data%dt=dt
    if(present(seasonal_trend))data%seasonal_trend=seasonal_trend
    if(present(me_ttm))then;data%me_ttm=me_ttm;data%has_me_ttm=.true.;end if
    if(present(initial_state))then;data%fixed_initial_state=initial_state;data%has_fixed_initial=.true.;end if

    npop=ctrl%population_size
    if(npop<=0)npop=max(20,10*n)
    npop=max(npop,4)
    allocate(pop(n,npop),value(npop),trial(n),mutant(n),best(n))
    call rng%seed(ctrl%seed)
    do j=1,npop
      do i=1,n
        pop(i,j)=lower(i)+(upper(i)-lower(i))*rng%uniform()
      end do
      value(j)=mle_objective(pop(:,j),data)
    end do
    if(present(starting_values))then
      pop(:,1)=max(lower,min(upper,starting_values));value(1)=mle_objective(pop(:,1),data)
    end if
    best_index=minloc(value,dim=1);best=pop(:,best_index);best_value=value(best_index)
    no_improve=0

    result%generations=0
    do g=1,ctrl%generations
      previous_best=best_value
      do j=1,npop
        call distinct_indices(rng,npop,j,a,b,c)
        mutant=pop(:,a)+ctrl%differential_weight*(pop(:,b)-pop(:,c))
        mutant=max(lower,min(upper,mutant))
        forced=1+int(rng%uniform()*real(n,dp));forced=min(n,max(1,forced))
        do i=1,n
          if(rng%uniform()<ctrl%crossover_probability .or. i==forced)then
            trial(i)=mutant(i)
          else
            trial(i)=pop(i,j)
          end if
        end do
        trial_value=mle_objective(trial,data)
        if(trial_value<=value(j))then
          pop(:,j)=trial;value(j)=trial_value
          if(trial_value<best_value)then;best_value=trial_value;best=trial;end if
        end if
      end do
      if(ctrl%trace>0)write(*,'(a,i0,a,es16.8)')'generation ',g,': negative log likelihood = ',best_value
      if(abs(previous_best-best_value)<=ctrl%solution_tolerance*max(1.0_dp,abs(best_value)))then
        no_improve=no_improve+1
      else
        no_improve=0
      end if
      result%generations=g
      if(no_improve>=ctrl%wait_generations)exit
    end do

    theta=best
    if(ctrl%local_refinement)then
      local_ctrl=lbfgsb_control_t()
      local_ctrl%max_evaluations=ctrl%local_max_evaluations
      local_ctrl%pgtol=1.0e-6_dp
      local_ctrl%reltol=1.0e-9_dp
      local_ctrl%trace=max(0,ctrl%trace-1)
      call lbfgsb_minimize_fd(theta,mle_objective_callback,local_result,lower,upper,local_ctrl,data)
      if(ieee_is_finite(local_result%value) .and. local_result%value<best_value)then
        best_value=local_result%value;best=theta
      end if
    end if

    call unpack_parameters(best,spec,fitted_model,x0,status)
    if(status/=nfcp_ok)then;result%message='failed to unpack fitted parameters';return;end if
    call order_model_factors(fitted_model,x0,spec)
    call pack_parameters(fitted_model,spec,x0,best,status)
    if(spec%estimate_initial_state)then
      call filter_for_data(fitted_model,data,filter,x0)
    else if(data%has_fixed_initial)then
      call filter_for_data(fitted_model,data,filter,data%fixed_initial_state)
    else
      call filter_for_data(fitted_model,data,filter)
    end if

    result%parameters=best;result%log_likelihood=filter%log_likelihood;result%filter=filter
    result%evaluations=data%evaluations
    allocate(hessian(n,n),inverse_hessian(n,n),result%standard_errors(n))
    call finite_difference_hessian(objective_at_solution,best,hessian,status)
    result%hessian=hessian
    if(status==nfcp_ok)then
      call inverse_spd(hessian,inverse_hessian,status)
    end if
    if(status==nfcp_ok)then
      do i=1,n;result%standard_errors(i)=sqrt(max(0.0_dp,inverse_hessian(i,i)));end do
    else
      result%standard_errors=huge(1.0_dp)
    end if
    result%converged=filter%status==nfcp_ok
    result%status=merge(nfcp_ok,nfcp_not_converged,result%converged)
    result%message=merge('ok           ','not converged',result%converged)

  contains
    real(dp) function objective_at_solution(x) result(f)
      real(dp),intent(in)::x(:)
      f=mle_objective(x,data)
    end function objective_at_solution
  end subroutine nfcp_fit_mle

  real(dp) function mle_objective_callback(x,user_data) result(f)
    real(dp),intent(in)::x(:)
    class(*),intent(inout),optional::user_data
    if(.not.present(user_data))then;f=huge(1.0_dp);return;end if
    select type(user_data)
    type is(mle_data_t)
      f=mle_objective(x,user_data)
    class default
      f=huge(1.0_dp)
    end select
  end function mle_objective_callback

  real(dp) function mle_objective(theta,data) result(value)
    real(dp),intent(in)::theta(:)
    type(mle_data_t),intent(inout)::data
    type(nfcp_model_t)::model
    type(nfcp_filter_result_t)::filter
    real(dp),allocatable::x0(:)
    integer::status
    data%evaluations=data%evaluations+1
    call unpack_parameters(theta,data%spec,model,x0,status)
    if (status /= nfcp_ok) then
      value = 1.0e100_dp
      return
    end if
    if (.not. model%valid()) then
      value = 1.0e100_dp
      return
    end if
    if(data%spec%estimate_initial_state)then
      call filter_for_data(model,data,filter,x0)
    else if(data%has_fixed_initial)then
      call filter_for_data(model,data,filter,data%fixed_initial_state)
    else
      call filter_for_data(model,data,filter)
    end if
    if(filter%status/=nfcp_ok .or. .not.ieee_is_finite(filter%log_likelihood))then
      value=1.0e100_dp
    else
      value=-filter%log_likelihood
    end if
  end function mle_objective

  subroutine filter_for_data(model,data,filter,initial_state)
    type(nfcp_model_t),intent(in)::model
    type(mle_data_t),intent(in)::data
    type(nfcp_filter_result_t),intent(out)::filter
    real(dp),intent(in),optional::initial_state(:)
    if(present(initial_state))then
      if(data%has_me_ttm)then
        call nfcp_kalman_filter(model,data%log_futures,data%dt,data%futures_ttm,filter, &
          initial_state=initial_state,me_ttm=data%me_ttm,seasonal_trend=data%seasonal_trend, &
          n_parameters=data%spec%count())
      else
        call nfcp_kalman_filter(model,data%log_futures,data%dt,data%futures_ttm,filter, &
          initial_state=initial_state,seasonal_trend=data%seasonal_trend,n_parameters=data%spec%count())
      end if
    else
      if(data%has_me_ttm)then
        call nfcp_kalman_filter(model,data%log_futures,data%dt,data%futures_ttm,filter, &
          me_ttm=data%me_ttm,seasonal_trend=data%seasonal_trend,n_parameters=data%spec%count())
      else
        call nfcp_kalman_filter(model,data%log_futures,data%dt,data%futures_ttm,filter, &
          seasonal_trend=data%seasonal_trend,n_parameters=data%spec%count())
      end if
    end if
  end subroutine filter_for_data

  subroutine distinct_indices(rng,n,excluded,a,b,c)
    type(nfcp_rng_t),intent(inout)::rng
    integer,intent(in)::n,excluded
    integer,intent(out)::a,b,c
    do;a=1+int(rng%uniform()*n);a=min(n,a);if(a/=excluded)exit;end do
    do;b=1+int(rng%uniform()*n);b=min(n,b);if(b/=excluded.and.b/=a)exit;end do
    do;c=1+int(rng%uniform()*n);c=min(n,c);if(c/=excluded.and.c/=a.and.c/=b)exit;end do
  end subroutine distinct_indices

  subroutine order_model_factors(model,x0,spec)
    type(nfcp_model_t),intent(inout)::model
    real(dp),intent(inout)::x0(:)
    type(nfcp_parameterization_t),intent(in)::spec
    integer::first,i,j,tmp,n
    integer,allocatable::order(:),full_order(:)
    real(dp),allocatable::lambda(:),kappa(:),sigma(:),rho(:,:),xcopy(:)
    first=merge(2,1,model%gbm);n=model%n_factors-first+1
    if(n<=1)return
    allocate(order(n));do i=1,n;order(i)=first+i-1;end do
    do i=1,n-1;do j=i+1,n
      if(model%kappa(order(j))<model%kappa(order(i)))then;tmp=order(i);order(i)=order(j);order(j)=tmp;end if
    end do;end do
    allocate(full_order(model%n_factors))
    if(model%gbm)then;full_order(1)=1;full_order(2:)=order;else;full_order=order;end if
    lambda=model%lambda;kappa=model%kappa;sigma=model%sigma;rho=model%rho;xcopy=x0
    model%lambda=lambda(full_order);model%kappa=kappa(full_order);model%sigma=sigma(full_order)
    do i=1,model%n_factors;do j=1,model%n_factors
      model%rho(i,j)=rho(full_order(i),full_order(j))
    end do;end do
    if(spec%estimate_initial_state)x0=xcopy(full_order)
  end subroutine order_model_factors

end module nfcp_mle
