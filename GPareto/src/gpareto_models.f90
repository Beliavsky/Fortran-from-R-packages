! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_models
  use gpareto_kinds, only : dp
  use dk_model, only : km_model, km_prediction, km_fit, km_predict, km_update
  use dk_model, only : trend_constant, trend_linear, trend_linear_interactions, trend_quadratic
  use dk_model, only : km_control
  implicit none
  private
  integer, parameter, public :: trend_const=1, trend_lin=2, trend_inter=3, trend_quad=4

  type, public :: gp_model
    type(km_model) :: km
    integer :: trend_kind = trend_const
  end type gp_model
  type, public :: gp_model_set
    type(gp_model), allocatable :: model(:)
  contains
    procedure :: nobj => modelset_nobj
    procedure :: dim => modelset_dim
  end type gp_model_set
  public :: fit_gp_model, predict_gp, predict_gps, update_gp, make_trend
contains
  integer function modelset_nobj(self) result(n)
    class(gp_model_set), intent(in) :: self
    if(allocated(self%model)) then
    n=size(self%model)
    else
    n=0
    end if
  end function modelset_nobj

  integer function modelset_dim(self) result(d)
    class(gp_model_set), intent(in) :: self
    if(allocated(self%model)) then
    d=self%model(1)%km%d
    else
    d=0
    end if
  end function modelset_dim

  subroutine make_trend(kind,x,f)
    integer,intent(in)::kind
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::f(:,:)
    select case(kind)
    case(trend_const)
    call trend_constant(x,f)
    case(trend_lin)
    call trend_linear(x,f)
    case(trend_inter)
    call trend_linear_interactions(x,f)
    case(trend_quad)
    call trend_quadratic(x,f)
    case default
    error stop 'make_trend: unsupported trend kind'
    end select
  end subroutine make_trend

  subroutine fit_gp_model(model,x,y,covtype,trend_kind,iso,noise_var,nugget,nugget_estim,estim_method,control)
    type(gp_model),intent(out)::model
    real(dp),intent(in)::x(:,:),y(:)
    character(len=*),intent(in),optional::covtype,estim_method
    integer,intent(in),optional::trend_kind
    logical,intent(in),optional::iso,nugget_estim
    real(dp),intent(in),optional::noise_var(:),nugget
    type(km_control),intent(in),optional::control
    real(dp),allocatable::f(:,:)
    character(len=16)::ct,em
    logical::isiso,ne
    integer::tk
    ct='matern5_2'
    if(present(covtype))ct=covtype
    em='MLE'
    if(present(estim_method))em=estim_method
    tk=trend_const
    if(present(trend_kind))tk=trend_kind
    isiso=.false.
    if(present(iso))isiso=iso
    ne=.false.
    if(present(nugget_estim))ne=nugget_estim
    model%trend_kind=tk
    call make_trend(tk,x,f)
    if(present(noise_var)) then
      call km_fit(model%km,x,y,f,trim(ct),noise_var=noise_var,estim_method=trim(em),iso=isiso,control=control)
    else if(present(nugget)) then
      call km_fit(model%km,x,y,f,trim(ct),nugget=nugget,nugget_estim=ne,estim_method=trim(em),iso=isiso,control=control)
    else
      call km_fit(model%km,x,y,f,trim(ct),estim_method=trim(em),iso=isiso,control=control)
    end if
  end subroutine fit_gp_model

  subroutine predict_gp(model,x,mean,sd,cov,kind)
    type(gp_model),intent(in)::model
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::mean(:),sd(:)
    real(dp),allocatable,intent(out),optional::cov(:,:)
    character(len=*),intent(in),optional::kind
    real(dp),allocatable::f(:,:)
    type(km_prediction)::p
    character(len=2)::kt
    kt='UK'
    if(present(kind))kt=kind
    call make_trend(model%trend_kind,x,f)
    if(present(cov)) then
      call km_predict(model%km,x,f,kt,p,se_compute=.true.,cov_compute=.true.)
      cov=p%cov
    else
      call km_predict(model%km,x,f,kt,p,se_compute=.true.,cov_compute=.false.)
    end if
    mean=p%mean
    sd=p%sd
  end subroutine predict_gp

  subroutine predict_gps(models,x,mean,sd,kind)
    type(gp_model_set),intent(in)::models
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::mean(:,:),sd(:,:)
    character(len=*),intent(in),optional::kind
    real(dp),allocatable::m(:),s(:)
    integer::j
    allocate(mean(size(x,1),models%nobj()),sd(size(x,1),models%nobj()))
    do j=1,models%nobj()
      call predict_gp(models%model(j),x,m,s,kind=kind)
      mean(:,j)=m
      sd(:,j)=s
    end do
  end subroutine predict_gps

  subroutine update_gp(model,newx,newy,new_noise,cov_reestimate)
    type(gp_model),intent(inout)::model
    real(dp),intent(in)::newx(:,:),newy(:)
    real(dp),intent(in),optional::new_noise(:)
    logical,intent(in),optional::cov_reestimate
    real(dp),allocatable::f(:,:)
    logical::cr
    cr=model%km%param_estim
    if(present(cov_reestimate))cr=cov_reestimate
    call make_trend(model%trend_kind,newx,f)
    if(present(new_noise)) then
      call km_update(model%km,newx,newy,f,newnoise_var=new_noise,cov_reestimate=cr)
    else
      call km_update(model%km,newx,newy,f,cov_reestimate=cr)
    end if
  end subroutine update_gp
end module gpareto_models
