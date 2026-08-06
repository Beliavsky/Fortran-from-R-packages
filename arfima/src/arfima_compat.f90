module arfima_compat
  use arfima_kinds, only : dp
  use arfima_types
  use arfima_status, only : arfima_invalid_input
  use arfima_polynomial, only : ar_to_pacf, pacf_to_ar, integrate_series, psi_weights, exact_integration_weights
  use arfima_autocov, only : tacvf_arma, tacvf_fdwn, tacvf_fgn, tacvf_pla, tacvf_arfima
  use arfima_durbin, only : dl_loglikelihood
  use arfima_transfer, only : apply_transfer_function
  use arfima_information_mod, only : arfima_information, identifiable_invertible
  use arfima_simulation, only : arfima_simulate
  use arfima_forecast_mod, only : arfima_forecast
  use arfima_fit, only : fit_arfima, fit_arfima_modes, arfima0_fit
  implicit none
  private
  public :: ARToPacf, PacfToAR, IdentInvertQ, iARFIMA
  public :: lARFIMA, tacvfARFIMA, tacvfARMA, tacvfFDWN, tacvfFGN, tacvfHD
  public :: arfima_sim, sim_from_fitted, predict_ARFIMA, predict_from_fitted, arfima_estimate, arfima0
  public :: integ, psiwts, wtsforexact, tacvf_from_fitted, lARFIMAwTF
contains

  function ARToPacf(phi) result(pacf)
    real(dp),intent(in)::phi(:)
    real(dp),allocatable::pacf(:)
    pacf=ar_to_pacf(phi)
  end function ARToPacf

  function PacfToAR(pacf) result(phi)
    real(dp),intent(in)::pacf(:)
    real(dp),allocatable::phi(:)
    phi=pacf_to_ar(pacf)
  end function PacfToAR

  logical function IdentInvertQ(spec,params,ident,error)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    logical,intent(in),optional::ident
    type(arfima_error),intent(out),optional::error
    IdentInvertQ=identifiable_invertible(spec,params,ident,error)
  end function IdentInvertQ

  subroutine iARFIMA(spec,params,information,error)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    real(dp),allocatable,intent(out)::information(:,:)
    type(arfima_error),intent(out)::error
    call arfima_information(spec,params,information,error)
  end subroutine iARFIMA

  real(dp) function lARFIMA(z,spec,params,error) result(loglik)
    real(dp),intent(in)::z(:)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    type(arfima_error),intent(out)::error
    real(dp),allocatable::r(:)
    call tacvf_arfima(spec,params,size(z)-1,1.0_dp,r,error)
    if(error%code/=0) then; loglik=-1.0e100_dp; return; end if
    loglik=dl_loglikelihood(r,z-params%mean,error)
  end function lARFIMA

  real(dp) function lARFIMAwTF(z,spec,params,error) result(loglik)
    real(dp),intent(in)::z(:)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    type(arfima_error),intent(out)::error
    type(transfer_spec)::tr
    type(arfima_spec)::stationary_spec
    real(dp),allocatable::adjusted(:),effect(:),r(:)
    if(.not.spec%use_transfer) then
      loglik=lARFIMA(z,spec,params,error)
      return
    end if
    tr=spec%transfer
    tr%delta=params%delta
    tr%omega=params%omega
    call apply_transfer_function(z,tr,params%mean,adjusted,effect,error)
    if(error%code/=0) then; loglik=-1.0e100_dp; return; end if
    stationary_spec=spec
    stationary_spec%use_transfer=.false.
    call tacvf_arfima(stationary_spec,params,size(z)-1,1.0_dp,r,error)
    if(error%code/=0) then; loglik=-1.0e100_dp; return; end if
    loglik=dl_loglikelihood(r,adjusted,error)
  end function lARFIMAwTF

  subroutine tacvfARFIMA(spec,params,maxlag,gamma_out,error,sigma2)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    integer,intent(in)::maxlag
    real(dp),allocatable,intent(out)::gamma_out(:)
    type(arfima_error),intent(out)::error
    real(dp),intent(in),optional::sigma2
    real(dp)::s2
    s2=1.0_dp; if(present(sigma2)) s2=sigma2
    call tacvf_arfima(spec,params,maxlag,s2,gamma_out,error)
  end subroutine tacvfARFIMA

  subroutine tacvfARMA(phi,theta,maxlag,gamma_out,error,sigma2)
    real(dp),intent(in)::phi(:),theta(:)
    integer,intent(in)::maxlag
    real(dp),allocatable,intent(out)::gamma_out(:)
    type(arfima_error),intent(out)::error
    real(dp),intent(in),optional::sigma2
    call tacvf_arma(phi,theta,maxlag,sigma2,gamma_out,error)
  end subroutine tacvfARMA

  subroutine tacvfFDWN(dfrac,maxlag,gamma_out,error)
    real(dp),intent(in)::dfrac
    integer,intent(in)::maxlag
    real(dp),allocatable,intent(out)::gamma_out(:)
    type(arfima_error),intent(out)::error
    call tacvf_fdwn(dfrac,maxlag,gamma_out,error)
  end subroutine tacvfFDWN

  subroutine tacvfFGN(hurst,maxlag,gamma_out,error)
    real(dp),intent(in)::hurst
    integer,intent(in)::maxlag
    real(dp),allocatable,intent(out)::gamma_out(:)
    type(arfima_error),intent(out)::error
    call tacvf_fgn(hurst,maxlag,gamma_out,error)
  end subroutine tacvfFGN

  subroutine tacvfHD(alpha,maxlag,gamma_out,error)
    real(dp),intent(in)::alpha
    integer,intent(in)::maxlag
    real(dp),allocatable,intent(out)::gamma_out(:)
    type(arfima_error),intent(out)::error
    call tacvf_pla(alpha,maxlag,gamma_out,error)
  end subroutine tacvfHD

  subroutine arfima_sim(n,spec,params,sigma2,z,error,innov,zinit)
    integer,intent(in)::n
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    real(dp),intent(in)::sigma2
    real(dp),allocatable,intent(out)::z(:)
    type(arfima_error),intent(out)::error
    real(dp),intent(in),optional::innov(:),zinit(:)
    call arfima_simulate(spec,params,n,sigma2,z,error,innov,zinit)
  end subroutine arfima_sim

  subroutine sim_from_fitted(n,fit,z,error,xreg)
    integer,intent(in)::n
    type(arfima_fit_result),intent(in)::fit
    real(dp),allocatable,intent(out)::z(:)
    type(arfima_error),intent(out)::error
    real(dp),intent(in),optional::xreg(:,:)
    type(arfima_spec)::spec
    spec=fit%spec
    if(spec%use_regression) then
      if(.not.present(xreg)) then
        call set_error(error,arfima_invalid_input,'xreg is required for a fitted regression model')
        allocate(z(0)); return
      end if
      spec%xreg=xreg
    end if
    call arfima_simulate(spec,fit%parameters,n,fit%sigma2,z,error)
  end subroutine sim_from_fitted

  subroutine predict_ARFIMA(spec,params,z,h,forecast,error,xreg_future)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    real(dp),intent(in)::z(:)
    integer,intent(in)::h
    type(arfima_forecast_result),intent(out)::forecast
    type(arfima_error),intent(out)::error
    real(dp),intent(in),optional::xreg_future(:,:)
    call arfima_forecast(spec,params,z,h,forecast,error,xreg_future)
  end subroutine predict_ARFIMA

  subroutine predict_from_fitted(fit,z,h,forecast,error,xreg_future)
    type(arfima_fit_result),intent(in)::fit
    real(dp),intent(in)::z(:)
    integer,intent(in)::h
    type(arfima_forecast_result),intent(out)::forecast
    type(arfima_error),intent(out)::error
    real(dp),intent(in),optional::xreg_future(:,:)
    call arfima_forecast(fit%spec,fit%parameters,z,h,forecast,error,xreg_future)
    if(error%code==0) then
      forecast%covariance=fit%sigma2*forecast%covariance
      forecast%standard_error=sqrt(fit%sigma2)*forecast%standard_error
    end if
  end subroutine predict_from_fitted

  subroutine tacvf_from_fitted(fit,maxlag,gamma_out,error)
    type(arfima_fit_result),intent(in)::fit
    integer,intent(in)::maxlag
    real(dp),allocatable,intent(out)::gamma_out(:)
    type(arfima_error),intent(out)::error
    call tacvf_arfima(fit%spec,fit%parameters,maxlag,fit%sigma2,gamma_out,error)
  end subroutine tacvf_from_fitted

  subroutine arfima_estimate(spec,z,result,start)
    type(arfima_spec),intent(in)::spec
    real(dp),intent(in)::z(:)
    type(arfima_fit_result),intent(out)::result
    type(arfima_parameters),intent(in),optional::start
    call fit_arfima(spec,z,result,start)
  end subroutine arfima_estimate

  subroutine arfima0(z,p,dint,q,lmodel,result)
    real(dp),intent(in)::z(:)
    integer,intent(in)::p,dint,q,lmodel
    type(arfima_fit_result),intent(out)::result
    call arfima0_fit(z,p,dint,q,lmodel,result)
  end subroutine arfima0

  subroutine integ(y,zinit,dint,dseas,period,z,error)
    real(dp),intent(in)::y(:),zinit(:)
    integer,intent(in)::dint,dseas,period
    real(dp),allocatable,intent(out)::z(:)
    type(arfima_error),intent(out)::error
    call integrate_series(y,zinit,dint,dseas,period,z,error)
  end subroutine integ

  function psiwts(phi,theta,phiseas,thetaseas,dfrac,dfs,dint,dseas,period,n) result(w)
    real(dp),intent(in)::phi(:),theta(:),phiseas(:),thetaseas(:),dfrac,dfs
    integer,intent(in)::dint,dseas,period,n
    real(dp),allocatable::w(:)
    w=psi_weights(phi,theta,phiseas,thetaseas,dfrac,dfs,dint,dseas,period,n)
  end function psiwts

  function wtsforexact(dint,dseas,period,n) result(w)
    integer,intent(in)::dint,dseas,period,n
    real(dp),allocatable::w(:)
    w=exact_integration_weights(dint,dseas,period,n)
  end function wtsforexact
end module arfima_compat
