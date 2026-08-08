module nfcp_simulation
  use nfcp_types, only : dp, nfcp_model_t, nfcp_simulation_result_t, &
                         nfcp_futures_simulation_result_t, nfcp_ok, nfcp_invalid_input, nfcp_singular
  use nfcp_math, only : nfcp_rng_t, nfcp_covariance, nfcp_a_t, nfcp_seasonality, &
                        correlated_normals
  implicit none
  private
  public :: spot_price_simulate, futures_price_simulate

contains

  subroutine spot_price_simulate(model, initial_state, horizon, dt, n_simulations, result, &
                                 antithetic, seasonal_trend, seed)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: initial_state(:), horizon, dt
    integer, intent(in) :: n_simulations
    type(nfcp_simulation_result_t), intent(out) :: result
    logical, intent(in), optional :: antithetic
    real(dp), intent(in), optional :: seasonal_trend
    integer, intent(in), optional :: seed

    logical :: anti
    integer :: nsteps, nsim_work, nhalf, t, s, i, f, status, seed_value
    real(dp) :: trend, sign_value, mean_level
    real(dp), allocatable :: correlation(:,:), draws(:,:)
    type(nfcp_rng_t) :: rng

    result%status = nfcp_invalid_input
    result%message = 'invalid input'
    if (size(initial_state) /= model%n_factors .or. horizon < 0.0_dp .or. dt <= 0.0_dp .or. &
        n_simulations < 1) return
    anti = .true.
    if (present(antithetic)) anti = antithetic
    trend = 0.0_dp
    if (present(seasonal_trend)) trend = seasonal_trend
    seed_value = 12345
    if (present(seed)) seed_value = seed
    call rng%seed(seed_value)

    nsteps = nint(horizon/dt)
    if (abs(real(nsteps,dp)*dt-horizon) > 1.0e-10_dp*max(1.0_dp,horizon)) then
      result%message = 'horizon must be an integer multiple of dt'
      return
    end if
    if (anti) then
      nsim_work = n_simulations + mod(n_simulations,2)
      nhalf = nsim_work/2
    else
      nsim_work = n_simulations
      nhalf = n_simulations
    end if

    allocate(result%times(nsteps+1), result%states(nsteps+1,n_simulations,model%n_factors))
    allocate(result%spot_prices(nsteps+1,n_simulations))
    allocate(correlation(model%n_factors,model%n_factors), draws(model%n_factors,nsteps*nhalf))
    correlation = model%rho
    call correlated_normals(rng,correlation,draws,status)
    if (status /= nfcp_ok) then
      result%status = nfcp_singular
      result%message = 'factor correlation matrix is not positive definite'
      return
    end if
    draws = sqrt(dt)*draws

    do t = 1, nsteps+1
      result%times(t)=real(t-1,dp)*dt
    end do
    do s = 1, n_simulations
      result%states(1,s,:) = initial_state
    end do

    do t = 2, nsteps+1
      do s = 1, n_simulations
        if (anti) then
          i = (s+1)/2
          sign_value = merge(1.0_dp,-1.0_dp,mod(s,2)==1)
        else
          i = s
          sign_value = 1.0_dp
        end if
        do f=1,model%n_factors
          if (model%gbm .and. f == 1) then
            result%states(t,s,f) = result%states(t-1,s,f) + model%mu_rn*dt + &
                                   sign_value*model%sigma(f)*draws(f,(t-2)*nhalf+i)
          else
            if (model%kappa(f) > 0.0_dp) then
              mean_level = -model%lambda(f)/model%kappa(f)
              result%states(t,s,f) = result%states(t-1,s,f) + &
                model%kappa(f)*(mean_level-result%states(t-1,s,f))*dt + &
                sign_value*model%sigma(f)*draws(f,(t-2)*nhalf+i)
            else
              result%states(t,s,f) = result%states(t-1,s,f) + &
                sign_value*model%sigma(f)*draws(f,(t-2)*nhalf+i)
            end if
          end if
        end do
      end do
    end do

    do t = 1, nsteps+1
      do s = 1, n_simulations
        result%spot_prices(t,s) = exp(model%equilibrium + sum(result%states(t,s,:)) + &
          nfcp_seasonality(model,result%times(t)+trend))
      end do
    end do
    result%status=nfcp_ok
    result%message='ok'
  end subroutine spot_price_simulate

  subroutine futures_price_simulate(model, initial_state, dt, futures_ttm, result, &
                                    me_ttm, seasonal_trend, seed)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: initial_state(:), dt
    real(dp), intent(in) :: futures_ttm(:,:)
    type(nfcp_futures_simulation_result_t), intent(out) :: result
    real(dp), intent(in), optional :: me_ttm(:), seasonal_trend
    integer, intent(in), optional :: seed

    integer :: nt, nc, nf, t, j, status, seed_value
    real(dp) :: trend, d, me_sd
    real(dp), allocatable :: q(:,:), omega(:,:), state(:), gdiag(:)
    type(nfcp_rng_t) :: rng

    result%status=nfcp_invalid_input
    result%message='invalid input'
    nt=size(futures_ttm,1); nc=size(futures_ttm,2); nf=model%n_factors
    if (size(initial_state)/=nf .or. nt<1 .or. nc<1 .or. dt<=0.0_dp) return
    if (model%n_me>1 .and. model%n_me<nc) then
      if (.not.present(me_ttm) .or. size(me_ttm)/=model%n_me) then
        result%message='me_ttm required for grouped measurement errors'; return
      end if
    end if
    trend=0.0_dp
    if(present(seasonal_trend)) trend=seasonal_trend
    seed_value=12345
    if(present(seed)) seed_value=seed
    call rng%seed(seed_value)

    allocate(result%states(nt,nf),result%futures_prices(nt,nc),result%spot_prices(nt))
    allocate(q(nf,nf),omega(nf,max(1,nt-1)),state(nf),gdiag(nf))
    call nfcp_covariance(model,dt,q)
    if(nt>1) then
      call correlated_normals(rng,q,omega(:,1:nt-1),status)
      if(status/=nfcp_ok) then
        result%status=nfcp_singular; result%message='state covariance is not positive definite'; return
      end if
    end if
    state=initial_state
    gdiag=exp(-model%kappa*dt)
    do t=1,nt
      result%states(t,:)=state
      result%spot_prices(t)=exp(model%equilibrium+sum(state)+ &
        nfcp_seasonality(model,real(t-1,dp)*dt+trend))
      do j=1,nc
        d=model%equilibrium+nfcp_seasonality(model,futures_ttm(t,j)+real(t-1,dp)*dt+trend)+ &
          nfcp_a_t(model,futures_ttm(t,j))
        me_sd=measurement_sd(model,j,futures_ttm(t,j),nc,me_ttm)
        result%futures_prices(t,j)=exp(d+dot_product(exp(-model%kappa*futures_ttm(t,j)),state)+ &
          me_sd*rng%normal())
      end do
      if(t<nt) then
        state=gdiag*state+omega(:,t)
        if(model%gbm) state(1)=state(1)+model%mu*dt
      end if
    end do
    result%status=nfcp_ok; result%message='ok'
  end subroutine futures_price_simulate

  real(dp) function measurement_sd(model,contract,maturity,n_contracts,me_ttm) result(sd)
    type(nfcp_model_t),intent(in)::model
    integer,intent(in)::contract,n_contracts
    real(dp),intent(in)::maturity
    real(dp),intent(in),optional::me_ttm(:)
    integer::k
    if(model%n_me<=0) then
      sd=0.0_dp
    else if(model%n_me==1) then
      sd=model%measurement_error(1)
    else if(model%n_me==n_contracts) then
      sd=model%measurement_error(contract)
    else
      k=model%n_me
      if(present(me_ttm)) then
        do k=1,model%n_me
          if(maturity<me_ttm(k)) exit
        end do
        k=min(k,model%n_me)
      end if
      sd=model%measurement_error(k)
    end if
  end function measurement_sd

end module nfcp_simulation
