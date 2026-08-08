module nfcp_options
  use nfcp_types, only : dp, nfcp_model_t, nfcp_option_result_t, nfcp_simulation_result_t, &
                         nfcp_ok, nfcp_invalid_input, nfcp_singular
  use nfcp_math, only : nfcp_a_t, nfcp_covariance, nfcp_seasonality, normal_cdf, nfcp_pi, &
                        inverse_spd, sample_sd
  use nfcp_forecast, only : futures_price_forecast
  use nfcp_simulation, only : spot_price_simulate
  implicit none
  private
  public :: european_option_value, american_option_value

contains

  subroutine european_option_value(model, initial_state, futures_maturity, option_maturity, strike, &
                                   risk_free_rate, is_call, result, seasonal_trend)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: initial_state(:), futures_maturity, option_maturity, strike, risk_free_rate
    logical, intent(in) :: is_call
    type(nfcp_option_result_t), intent(out) :: result
    real(dp), intent(in), optional :: seasonal_trend
    real(dp), allocatable :: fmat(:,:), covariance(:,:)
    real(dp) :: f0, total_sd, d1, d2, phi, h, vp, vm
    integer :: i,j,status

    result%status=nfcp_invalid_input;result%message='invalid input'
    if(size(initial_state)/=model%n_factors .or. futures_maturity<option_maturity .or. &
       option_maturity<0.0_dp .or. strike<=0.0_dp) return
    if(present(seasonal_trend)) then
      call futures_price_forecast(model,initial_state,0.0_dp,[futures_maturity],fmat, &
                                  seasonal_trend=seasonal_trend,status=status)
    else
      call futures_price_forecast(model,initial_state,0.0_dp,[futures_maturity],fmat,status=status)
    end if
    if(status/=nfcp_ok) return
    f0=fmat(1,1)
    allocate(covariance(model%n_factors,model%n_factors))
    call nfcp_covariance(model,option_maturity,covariance)
    total_sd=0.0_dp
    do i=1,model%n_factors
      do j=1,model%n_factors
        total_sd=total_sd+covariance(i,j)*exp(-(model%kappa(i)+model%kappa(j))* &
          (futures_maturity-option_maturity))
      end do
    end do
    total_sd=sqrt(max(0.0_dp,total_sd))
    result%value=black_futures_value(f0,strike,total_sd,risk_free_rate,option_maturity,is_call)
    if (option_maturity > 0.0_dp) then
      result%annualized_volatility = total_sd / sqrt(option_maturity)
    else
      result%annualized_volatility = 0.0_dp
    end if
    if(total_sd>sqrt(epsilon(1.0_dp))) then
      d1=(log(f0/strike)+0.5_dp*total_sd**2)/total_sd
      d2=d1-total_sd
      phi=exp(-0.5_dp*d1*d1)/sqrt(2.0_dp*nfcp_pi)
      if(is_call) then
        result%delta=exp(-risk_free_rate*option_maturity)*normal_cdf(d1)
      else
        result%delta=exp(-risk_free_rate*option_maturity)*(normal_cdf(d1)-1.0_dp)
      end if
      result%gamma=exp(-risk_free_rate*option_maturity)*phi/(f0*total_sd)
      result%vega=exp(-risk_free_rate*option_maturity)*f0*phi
    else
      result%delta=0.0_dp;result%gamma=0.0_dp;result%vega=0.0_dp
    end if
    result%rho=-option_maturity*result%value
    if(option_maturity>1.0e-6_dp) then
      h=min(1.0e-5_dp,max(1.0e-8_dp,0.1_dp*option_maturity))
      vp=european_value_at_maturity(model,initial_state,futures_maturity,option_maturity+h, &
                                    strike,risk_free_rate,is_call,seasonal_trend)
      vm=european_value_at_maturity(model,initial_state,futures_maturity,max(0.0_dp,option_maturity-h), &
                                    strike,risk_free_rate,is_call,seasonal_trend)
      result%theta=-(vp-vm)/(option_maturity+h-max(0.0_dp,option_maturity-h))
    end if
    result%status=nfcp_ok;result%message='ok'
  end subroutine european_option_value

  real(dp) function european_value_at_maturity(model,x0,futures_maturity,option_maturity,strike,r,is_call,trend) result(v)
    type(nfcp_model_t),intent(in)::model
    real(dp),intent(in)::x0(:),futures_maturity,option_maturity,strike,r
    logical,intent(in)::is_call
    real(dp),intent(in),optional::trend
    real(dp),allocatable::fmat(:,:),cov(:,:)
    real(dp)::sd
    integer::i,j,status
    if(present(trend)) then
      call futures_price_forecast(model,x0,0.0_dp,[futures_maturity],fmat,seasonal_trend=trend,status=status)
    else
      call futures_price_forecast(model,x0,0.0_dp,[futures_maturity],fmat,status=status)
    end if
    if(status/=nfcp_ok) then;v=huge(1.0_dp);return;end if
    allocate(cov(model%n_factors,model%n_factors));call nfcp_covariance(model,option_maturity,cov)
    sd=0.0_dp
    do i=1,model%n_factors;do j=1,model%n_factors
      sd=sd+cov(i,j)*exp(-(model%kappa(i)+model%kappa(j))*(futures_maturity-option_maturity))
    end do;end do
    v=black_futures_value(fmat(1,1),strike,sqrt(max(0.0_dp,sd)),r,option_maturity,is_call)
  end function european_value_at_maturity

  pure real(dp) function black_futures_value(f,k,sd,r,t,is_call) result(v)
    real(dp),intent(in)::f,k,sd,r,t
    logical,intent(in)::is_call
    real(dp)::d1,d2,disc
    disc=exp(-r*t)
    if(sd<=sqrt(epsilon(1.0_dp))) then
      if(is_call) then;v=disc*max(f-k,0.0_dp);else;v=disc*max(k-f,0.0_dp);end if
      return
    end if
    d1=(log(f/k)+0.5_dp*sd*sd)/sd;d2=d1-sd
    if(is_call) then
      v=disc*(f*normal_cdf(d1)-k*normal_cdf(d2))
    else
      v=disc*(k*normal_cdf(-d2)-f*normal_cdf(-d1))
    end if
  end function black_futures_value

  subroutine american_option_value(model,initial_state,futures_maturity,option_maturity,strike, &
                                   risk_free_rate,is_call,n_simulations,dt,result,degree,basis, &
                                   seasonal_trend,seed,antithetic)
    type(nfcp_model_t),intent(in)::model
    real(dp),intent(in)::initial_state(:),futures_maturity,option_maturity,strike,risk_free_rate,dt
    logical,intent(in)::is_call
    integer,intent(in)::n_simulations
    type(nfcp_option_result_t),intent(out)::result
    integer,intent(in),optional::degree,seed
    character(len=*),intent(in),optional::basis
    real(dp),intent(in),optional::seasonal_trend
    logical,intent(in),optional::antithetic

    type(nfcp_simulation_result_t)::sim
    integer::deg,nsteps,t,s,i,nitm,status,seed_value
    logical::anti
    character(len=16)::basis_name
    real(dp)::trend,remaining,intrinsic0,continuation0
    real(dp),allocatable::futures(:,:),intrinsic(:,:),cashflow(:),y(:),x(:),cont(:)
    integer,allocatable::exercise_time(:),idx(:),counts(:)

    result%status=nfcp_invalid_input;result%message='invalid input'
    if(futures_maturity<option_maturity .or. option_maturity<=0.0_dp .or. dt<=0.0_dp .or. &
       n_simulations<2 .or. strike<=0.0_dp) return
    nsteps=nint(option_maturity/dt)
    if(abs(real(nsteps,dp)*dt-option_maturity)>1.0e-10_dp) then
      result%message='option_maturity must be an integer multiple of dt';return
    end if
    deg=2;if(present(degree))deg=max(1,degree)
    basis_name='power';if(present(basis))basis_name=adjustl(basis)
    trend=0.0_dp;if(present(seasonal_trend))trend=seasonal_trend
    seed_value=12345;if(present(seed))seed_value=seed
    anti=.true.;if(present(antithetic))anti=antithetic
    call spot_price_simulate(model,initial_state,option_maturity,dt,n_simulations,sim, &
                             anti,trend,seed_value)
    if(sim%status/=nfcp_ok) then;result%message=sim%message;return;end if
    allocate(futures(nsteps+1,n_simulations),intrinsic(nsteps+1,n_simulations))
    do t=1,nsteps+1
      remaining=futures_maturity-real(t-1,dp)*dt
      do s=1,n_simulations
        futures(t,s)=exp(model%equilibrium+nfcp_seasonality(model,remaining+trend)+ &
          nfcp_a_t(model,remaining)+dot_product(sim%states(t,s,:),exp(-model%kappa*remaining)))
        if(is_call) then;intrinsic(t,s)=max(futures(t,s)-strike,0.0_dp)
        else;intrinsic(t,s)=max(strike-futures(t,s),0.0_dp);end if
      end do
    end do
    allocate(cashflow(n_simulations),exercise_time(n_simulations))
    cashflow=intrinsic(nsteps+1,:);exercise_time=nsteps+1
    do t=nsteps,2,-1
      nitm=count(intrinsic(t,:)>0.0_dp)
      if(nitm<deg+1) cycle
      allocate(idx(nitm),x(nitm),y(nitm),cont(nitm))
      i=0
      do s=1,n_simulations
        if(intrinsic(t,s)>0.0_dp) then
          i=i+1;idx(i)=s;x(i)=futures(t,s)
          y(i)=cashflow(s)*exp(-risk_free_rate*dt*real(exercise_time(s)-t,dp))
        end if
      end do
      call continuation_regression(x,y,deg,basis_name,cont,status)
      if(status==nfcp_ok) then
        do i=1,nitm
          s=idx(i)
          if(intrinsic(t,s)>cont(i)) then
            cashflow(s)=intrinsic(t,s);exercise_time(s)=t
          end if
        end do
      end if
      deallocate(idx,x,y,cont)
    end do
    continuation0=sum(cashflow*exp(-risk_free_rate*dt*real(exercise_time-1,dp)))/real(n_simulations,dp)
    intrinsic0=intrinsic(1,1)
    result%value=max(intrinsic0,continuation0)
    result%standard_error=sample_sd(cashflow*exp(-risk_free_rate*dt*real(exercise_time-1,dp)))/ &
                          sqrt(real(n_simulations,dp))
    allocate(counts(nsteps+1),result%exercise_probability(nsteps+1));counts=0
    do s=1,n_simulations;counts(exercise_time(s))=counts(exercise_time(s))+1;end do
    result%exercise_probability=real(counts,dp)/real(n_simulations,dp)
    result%status=nfcp_ok;result%message='ok'
  end subroutine american_option_value

  subroutine continuation_regression(x,y,degree,basis_name,prediction,status)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in)::degree
    character(len=*),intent(in)::basis_name
    real(dp),intent(out)::prediction(:)
    integer,intent(out)::status
    real(dp),allocatable::design(:,:),xtx(:,:),xty(:),beta(:),inv(:,:)
    real(dp)::mean_x,scale_x
    integer::i
    allocate(design(size(x),degree+1),xtx(degree+1,degree+1),xty(degree+1), &
             beta(degree+1),inv(degree+1,degree+1))
    mean_x=sum(x)/real(size(x),dp);scale_x=max(sqrt(sum((x-mean_x)**2)/real(max(1,size(x)-1),dp)),1.0e-10_dp)
    do i=1,size(x);call basis_values((x(i)-mean_x)/scale_x,degree,basis_name,design(i,:));end do
    xtx=matmul(transpose(design),design)
    do i=1,degree+1;xtx(i,i)=xtx(i,i)+1.0e-10_dp;end do
    xty=matmul(transpose(design),y)
    call inverse_spd(xtx,inv,status)
    if(status/=nfcp_ok) then;prediction=0.0_dp;return;end if
    beta=matmul(inv,xty);prediction=matmul(design,beta)
  end subroutine continuation_regression

  subroutine basis_values(x,degree,name,v)
    real(dp),intent(in)::x
    integer,intent(in)::degree
    character(len=*),intent(in)::name
    real(dp),intent(out)::v(:)
    integer::k
    v=0.0_dp;v(1)=1.0_dp
    if(degree==0)return
    select case(trim(lowercase(name)))
    case('laguerre')
      v(2)=1.0_dp-x
      do k=2,degree;v(k+1)=((2.0_dp*k-1.0_dp-x)*v(k)-(k-1.0_dp)*v(k-1))/real(k,dp);end do
    case('hermite')
      v(2)=x
      do k=2,degree;v(k+1)=x*v(k)-real(k-1,dp)*v(k-1);end do
    case('legendre')
      v(2)=x
      do k=2,degree;v(k+1)=((2.0_dp*k-1.0_dp)*x*v(k)-(k-1.0_dp)*v(k-1))/real(k,dp);end do
    case('chebyshev')
      v(2)=x
      do k=2,degree;v(k+1)=2.0_dp*x*v(k)-v(k-1);end do
    case default
      do k=1,degree;v(k+1)=x**k;end do
    end select
  end subroutine basis_values

  pure function lowercase(s) result(out)
    character(len=*),intent(in)::s
    character(len=len(s))::out
    integer::i,c
    do i=1,len(s);c=iachar(s(i:i));if(c>=65.and.c<=90)then;out(i:i)=achar(c+32);else;out(i:i)=s(i:i);end if;end do
  end function lowercase

end module nfcp_options
