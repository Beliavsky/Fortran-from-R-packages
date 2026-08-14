module mlrmbo_gp
  use mlrmbo_kinds, only : dp
  use mlrmbo_types, only : mbo_control, mbo_path
  use dicekriging, only : km_model, km_control, km_prediction, km_fit, km_predict, km_update, trend_constant
  implicit none
  private
  type, public :: mbo_surrogates
    type(km_model), allocatable :: model(:)
    integer :: m=0
  end type mbo_surrogates
  public :: fit_surrogates, predict_surrogate, predict_surrogates, update_surrogate
contains
  subroutine fit_surrogates(path,control,sur)
    type(mbo_path), intent(in) :: path
    type(mbo_control), intent(in) :: control
    type(mbo_surrogates), intent(out) :: sur
    type(km_control) :: kc
    real(dp), allocatable :: f(:,:)
    integer :: j
    if(path%n<2) error stop 'fit_surrogates: at least two observations required'
    if(size(path%y,2)/=control%n_objectives) error stop 'fit_surrogates: objective mismatch'
    call trend_constant(path%x,f)
    kc%multistart=max(1,control%km_multistart)
    kc%pop_size=max(kc%multistart,control%km_pop_size)
    kc%max_iter=control%km_max_iter; kc%tol=control%km_tol; kc%use_gradient=.true.
    sur%m=control%n_objectives; allocate(sur%model(sur%m))
    do j=1,sur%m
      call km_fit(sur%model(j),path%x,path%y(:,j),f,trim(control%covariance), &
        nugget=1.0e-10_dp,control=kc)
    end do
  end subroutine fit_surrogates

  subroutine predict_surrogate(model,x,mean,sd)
    type(km_model), intent(in) :: model
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mean(:),sd(:)
    real(dp), allocatable :: f(:,:)
    type(km_prediction) :: p
    call trend_constant(x,f)
    call km_predict(model,x,f,'UK',p,se_compute=.true.)
    mean=p%mean; sd=p%sd
  end subroutine predict_surrogate

  subroutine predict_surrogates(sur,x,mean,sd)
    type(mbo_surrogates), intent(in) :: sur
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mean(:,:),sd(:,:)
    real(dp), allocatable :: a(:),b(:)
    integer :: j
    allocate(mean(size(x,1),sur%m),sd(size(x,1),sur%m))
    do j=1,sur%m
      call predict_surrogate(sur%model(j),x,a,b)
      mean(:,j)=a; sd(:,j)=b
    end do
  end subroutine predict_surrogates

  subroutine update_surrogate(model,x,y,reestimate)
    type(km_model), intent(inout) :: model
    real(dp), intent(in) :: x(:,:),y(:)
    logical, intent(in) :: reestimate
    real(dp), allocatable :: f(:,:)
    call trend_constant(x,f)
    call km_update(model,x,y,f,cov_reestimate=reestimate,trend_reestimate=reestimate)
  end subroutine update_surrogate
end module mlrmbo_gp
