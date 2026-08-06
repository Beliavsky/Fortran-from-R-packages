module arfima_simulation
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok, arfima_invalid_input
  use arfima_types, only : arfima_spec, arfima_parameters, arfima_error, set_error
  use arfima_autocov, only : tacvf_arfima
  use arfima_durbin, only : dl_simulate
  use arfima_polynomial, only : integrate_series
  use arfima_random, only : random_normal_vector
  implicit none
  private
  public :: arfima_simulate
contains

  subroutine arfima_simulate(spec,params,n,sigma2,z,error,innov,zinit,center_sample)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    integer,intent(in)::n
    real(dp),intent(in)::sigma2
    real(dp),allocatable,intent(out)::z(:)
    type(arfima_error),intent(out)::error
    real(dp),intent(in),optional::innov(:),zinit(:)
    logical,intent(in),optional::center_sample
    real(dp),allocatable::r(:),e(:),stationary(:),integrated(:),zi(:)
    type(arfima_error)::err
    logical::center
    integer::lag

    call set_error(error,arfima_ok,'')
    if(n<1 .or. sigma2<0.0_dp) then
      allocate(z(0)); call set_error(error,arfima_invalid_input,'n must be positive and sigma2 nonnegative'); return
    end if
    call tacvf_arfima(spec,params,n-1,sigma2,r,err)
    if(err%code/=arfima_ok) then; allocate(z(0)); error=err; return; end if
    allocate(e(n))
    if(present(innov)) then
      if(size(innov)/=n) then; allocate(z(0)); call set_error(error,arfima_invalid_input,'innov must have length n'); return; end if
      e=innov
    else
      call random_normal_vector(e)
    end if
    call dl_simulate(r,e,stationary,err)
    if(err%code/=arfima_ok) then; allocate(z(0)); error=err; return; end if
    center=.true.; if(present(center_sample)) center=center_sample
    if(center) stationary=stationary-sum(stationary)/real(n,dp)+params%mean
    if(.not.center) stationary=stationary+params%mean
    lag=spec%dint+spec%dseas*spec%period
    if(lag>0) then
      allocate(zi(lag))
      if(present(zinit)) then
        if(size(zinit)/=lag) then; allocate(z(0)); call set_error(error,arfima_invalid_input,'zinit has wrong length'); return; end if
        zi=zinit
      else
        zi=0.0_dp
      end if
      call integrate_series(stationary,zi,spec%dint,spec%dseas,spec%period,integrated,err)
      if(err%code/=arfima_ok) then; allocate(z(0)); error=err; return; end if
      allocate(z(n)); z=integrated(lag+1:lag+n)
    else
      allocate(z(n)); z=stationary
    end if
    if(spec%use_regression) then
      if(.not.allocated(spec%xreg) .or. .not.allocated(params%beta) .or. size(spec%xreg,1)/=n .or. &
         size(spec%xreg,2)/=size(params%beta)) then
        call set_error(error,arfima_invalid_input,'regression matrix or coefficients have incompatible dimensions'); return
      end if
      z=z+matmul(spec%xreg,params%beta)
    end if
    if(spec%use_transfer) then
      call set_error(error,arfima_invalid_input,'simulation from dynamic transfer-function fits is not implemented')
    end if
  end subroutine arfima_simulate
end module arfima_simulation
