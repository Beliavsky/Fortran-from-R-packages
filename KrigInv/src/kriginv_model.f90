module kriginv_model
  use kriginv_kinds, only : dp
  use kriginv_linalg, only : cholesky_lower, solve_lower, solve_chol, invert_spd, symmetrize
  use dk_model, only : km_model, km_control, km_prediction, dk_km_fit => km_fit, &
                       dk_km_predict => km_predict, dk_km_update => km_update, &
                       dk_trend_constant => trend_constant, dk_trend_linear => trend_linear, &
                       dk_trend_interactions => trend_linear_interactions, &
                       dk_trend_quadratic => trend_quadratic
  use dk_covariance, only : scaling_axis, covariance_name, dk_covariance_matrix => covariance_matrix, &
                            dk_covariance_cross => covariance_cross
  implicit none
  private

  type, public :: krig_model
    integer :: n=0, d=0, p=0
    integer :: trend_order=0
    character(len=16) :: covariance='gauss'
    real(dp) :: variance=1.0_dp
    real(dp) :: nugget=0.0_dp
    real(dp), allocatable :: x(:,:), y(:), noise(:), lengthscale(:)
    real(dp), allocatable :: l(:,:), f(:,:), m(:,:), beta(:), z(:), ginv(:,:)
    logical :: ready=.false.
    logical :: uses_dicekriging=.false.
    logical :: param_estim=.false.
    logical :: cov_reestimate_default=.false.
    type(km_model), allocatable :: dice_model
  end type krig_model

  type, public :: krig_prediction
    real(dp), allocatable :: mean(:), sd(:), covariance(:,:)
    real(dp), allocatable :: c(:,:), linv_c(:,:), fnew(:,:)
    logical :: ok=.true.
  end type krig_prediction

  public :: km_control, scaling_axis
  public :: init_krig_model, init_dice_krig_model, fit_krig_model, update_krig_model
  public :: predict_nobias_km, covariance_matrix, covariance_cross, posterior_covariance

contains
  subroutine init_krig_model(model,x,y,lengthscale,variance,nugget,noise,covariance,trend_order,ok)
    type(krig_model), intent(out) :: model
    real(dp), intent(in) :: x(:,:),y(:),lengthscale(:)
    real(dp), intent(in), optional :: variance,nugget,noise(:)
    character(len=*), intent(in), optional :: covariance
    integer, intent(in), optional :: trend_order
    logical, intent(out), optional :: ok
    logical :: good
    model%n=size(x,1); model%d=size(x,2)
    if(size(y)/=model%n .or. size(lengthscale)/=model%d .or. any(lengthscale<=0.0_dp)) then
      if(present(ok)) ok=.false.
      return
    end if
    model%x=x; model%y=y; model%lengthscale=lengthscale
    if(present(variance)) model%variance=variance
    if(present(nugget)) model%nugget=max(0.0_dp,nugget)
    if(present(covariance)) model%covariance=adjustl(covariance)
    if(present(trend_order)) model%trend_order=trend_order
    allocate(model%noise(model%n)); model%noise=0.0_dp
    if(present(noise)) then
      if(size(noise)/=model%n) then
        if(present(ok)) ok=.false.
        return
      end if
      model%noise=max(0.0_dp,noise)
    end if
    model%uses_dicekriging=.false.
    model%param_estim=.false.
    model%cov_reestimate_default=.false.
    call prepare_model(model,good)
    if(present(ok)) ok=good
  end subroutine init_krig_model

  subroutine init_dice_krig_model(model,x,y,lengthscale,variance,nugget,noise,covariance,trend_order,iso,ok)
    type(krig_model), intent(out) :: model
    real(dp), intent(in) :: x(:,:),y(:),lengthscale(:)
    real(dp), intent(in), optional :: variance,nugget,noise(:)
    character(len=*), intent(in), optional :: covariance
    integer, intent(in), optional :: trend_order
    logical, intent(in), optional :: iso
    logical, intent(out), optional :: ok
    character(len=16) :: covname
    integer :: tr
    real(dp) :: var
    logical :: good

    covname='gauss'; if(present(covariance)) covname=adjustl(covariance)
    tr=0; if(present(trend_order)) tr=trend_order
    var=1.0_dp; if(present(variance)) var=variance
    if(present(nugget) .and. present(noise)) then
      if(any(noise>0.0_dp) .and. nugget>0.0_dp) then
        if(present(ok)) ok=.false.
        return
      end if
    end if
    call fit_krig_model(model,x,y,covariance=covname,trend_order=tr,coef_cov=lengthscale, &
                        coef_var=var,nugget=nugget,noise=noise,iso=iso,ok=good)
    if(present(ok)) ok=good
  end subroutine init_dice_krig_model

  subroutine fit_krig_model(model,x,y,covariance,trend_order,coef_cov,coef_var,coef_trend,nugget, &
                            nugget_estim,noise,estim_method,iso,lower,upper,parinit,control, &
                            scaling_axes,scad_lambda,ok)
    type(krig_model), intent(out) :: model
    real(dp), intent(in) :: x(:,:),y(:)
    character(len=*), intent(in), optional :: covariance,estim_method
    integer, intent(in), optional :: trend_order
    real(dp), intent(in), optional :: coef_cov(:),coef_var,coef_trend(:),nugget,noise(:)
    logical, intent(in), optional :: nugget_estim,iso
    real(dp), intent(in), optional :: lower(:),upper(:),parinit(:),scad_lambda
    type(km_control), intent(in), optional :: control
    type(scaling_axis), intent(in), optional :: scaling_axes(:)
    logical, intent(out), optional :: ok
    character(len=16) :: covname
    real(dp), allocatable :: f(:,:)
    integer :: tr
    logical :: good,ne

    good=.false.
    if(size(y)/=size(x,1) .or. size(x,1)<1 .or. size(x,2)<1) then
      if(present(ok)) ok=.false.
      return
    end if
    tr=0; if(present(trend_order)) tr=trend_order
    if(tr<0 .or. tr>3) then
      if(present(ok)) ok=.false.
      return
    end if
    if(present(noise)) then
      if(size(noise)/=size(y) .or. any(noise<0.0_dp)) then
        if(present(ok)) ok=.false.
        return
      end if
    end if
    ne=.false.; if(present(nugget_estim)) ne=nugget_estim
    if(present(noise)) then
      if(any(noise>0.0_dp)) then
        if((present(nugget) .and. nugget>0.0_dp) .or. ne) then
          if(present(ok)) ok=.false.
          return
        end if
      end if
    end if

    covname='matern5_2'; if(present(covariance)) covname=adjustl(covariance)
    call make_trend(x,tr,f)
    allocate(model%dice_model)
    call dk_km_fit(model%dice_model,x,y,f,trim(covname),coef_cov=coef_cov,coef_var=coef_var, &
                   coef_trend=coef_trend,nugget=nugget,nugget_estim=nugget_estim,noise_var=noise, &
                   estim_method=estim_method,iso=iso,lower=lower,upper=upper,parinit=parinit, &
                   control=control,scaling_axes=scaling_axes,scad_lambda=scad_lambda)
    model%trend_order=tr
    model%uses_dicekriging=.true.
    call sync_from_dice(model,good)
    if(present(ok)) ok=good
  end subroutine fit_krig_model

  subroutine sync_from_dice(model,ok)
    type(krig_model), intent(inout) :: model
    logical, intent(out) :: ok
    real(dp), allocatable :: g(:,:)
    logical :: good

    if(.not.allocated(model%dice_model)) then
      ok=.false.; return
    end if
    model%n=model%dice_model%n; model%d=model%dice_model%d; model%p=model%dice_model%p
    model%x=model%dice_model%x; model%y=model%dice_model%y; model%f=model%dice_model%f
    model%l=model%dice_model%l; model%m=model%dice_model%m
    model%beta=model%dice_model%trend_coef; model%z=model%dice_model%z
    model%variance=model%dice_model%covariance%sd2
    model%nugget=0.0_dp
    if(model%dice_model%covariance%nugget_flag) model%nugget=model%dice_model%covariance%nugget
    model%covariance=covariance_name(model%dice_model%covariance%kind)
    if(allocated(model%dice_model%covariance%range)) then
      model%lengthscale=model%dice_model%covariance%range
    else
      allocate(model%lengthscale(model%d)); model%lengthscale=1.0_dp
    end if
    if(allocated(model%noise)) deallocate(model%noise)
    allocate(model%noise(model%n)); model%noise=0.0_dp
    if(model%dice_model%noise_flag) model%noise=model%dice_model%noise_var
    g=matmul(transpose(model%m),model%m)
    call invert_spd(g,model%ginv,good)
    if(.not.good) then
      model%ready=.false.; ok=.false.; return
    end if
    model%param_estim=model%dice_model%param_estim
    model%cov_reestimate_default=model%dice_model%param_estim
    model%uses_dicekriging=.true.
    model%ready=.true.; ok=.true.
  end subroutine sync_from_dice

  subroutine prepare_model(model,ok)
    type(krig_model), intent(inout) :: model
    logical, intent(out) :: ok
    real(dp), allocatable :: k(:,:),rhs(:,:),kinvf(:,:),kinvy(:,:),g(:,:),ginv(:,:),tmp(:,:)
    integer :: i
    logical :: good
    k=covariance_matrix(model,model%x,.true.)
    do i=1,model%n
      k(i,i)=k(i,i)+model%noise(i)
    end do
    call cholesky_lower(k,model%l,good)
    if(.not.good) then
      do i=1,model%n
        k(i,i)=k(i,i)+max(1.0e-12_dp,1.0e-10_dp*model%variance)
      end do
      call cholesky_lower(k,model%l,good)
    end if
    if(.not.good) then
      model%ready=.false.; ok=.false.; return
    end if
    model%f=trend_matrix(model,model%x); model%p=size(model%f,2)
    call solve_lower(model%l,model%f,model%m)
    call solve_chol(model%l,model%f,kinvf)
    allocate(g(model%p,model%p)); g=matmul(transpose(model%f),kinvf)
    call invert_spd(g,ginv,good)
    if(.not.good) then
      model%ready=.false.; ok=.false.; return
    end if
    model%ginv=ginv
    allocate(rhs(model%n,1)); rhs(:,1)=model%y
    call solve_chol(model%l,rhs,kinvy)
    allocate(tmp(model%p,1)); tmp=matmul(transpose(model%f),kinvy)
    model%beta=matmul(ginv,tmp(:,1))
    rhs(:,1)=model%y-matmul(model%f,model%beta)
    call solve_lower(model%l,rhs,tmp)
    model%z=tmp(:,1)
    model%ready=.true.; ok=.true.
  end subroutine prepare_model

  subroutine make_trend(x,order,f)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: order
    real(dp), allocatable, intent(out) :: f(:,:)
    select case(order)
    case(0)
      call dk_trend_constant(x,f)
    case(1)
      call dk_trend_linear(x,f)
    case(2)
      call dk_trend_interactions(x,f)
    case(3)
      call dk_trend_quadratic(x,f)
    case default
      error stop 'kriginv: trend_order must be 0, 1, 2, or 3'
    end select
  end subroutine make_trend

  function trend_matrix(model,x) result(f)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: f(:,:)
    call make_trend(x,model%trend_order,f)
  end function trend_matrix

  real(dp) function kernel_corr(model,a,b) result(r)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: a(:),b(:)
    real(dp) :: q,rr,sq5,sq3
    integer :: j
    q=0.0_dp
    do j=1,size(a)
      q=q+((a(j)-b(j))/model%lengthscale(j))**2
    end do
    rr=sqrt(max(0.0_dp,q))
    select case(trim(model%covariance))
    case('gauss','gaussian')
      r=exp(-q)
    case('exp','exponential')
      r=exp(-rr)
    case('matern3_2','matern32')
      sq3=sqrt(3.0_dp); r=(1.0_dp+sq3*rr)*exp(-sq3*rr)
    case('matern5_2','matern52')
      sq5=sqrt(5.0_dp); r=(1.0_dp+sq5*rr+5.0_dp*q/3.0_dp)*exp(-sq5*rr)
    case default
      r=exp(-q)
    end select
  end function kernel_corr

  function covariance_cross(model,x1,x2) result(c)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: x1(:,:),x2(:,:)
    real(dp), allocatable :: c(:,:)
    integer :: i,j
    if(model%uses_dicekriging .and. allocated(model%dice_model)) then
      call dk_covariance_cross(model%dice_model%covariance,x1,x2,c,include_nugget=.false.)
      return
    end if
    allocate(c(size(x1,1),size(x2,1)))
    do j=1,size(x2,1)
      do i=1,size(x1,1)
        c(i,j)=model%variance*kernel_corr(model,x1(i,:),x2(j,:))
      end do
    end do
  end function covariance_cross

  function covariance_matrix(model,x,include_nugget) result(c)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: include_nugget
    real(dp), allocatable :: c(:,:)
    logical :: nug
    integer :: i,j
    nug=.true.; if(present(include_nugget)) nug=include_nugget
    if(model%uses_dicekriging .and. allocated(model%dice_model)) then
      call dk_covariance_matrix(model%dice_model%covariance,x,c,include_nugget=nug)
      return
    end if
    allocate(c(size(x,1),size(x,1)))
    do j=1,size(x,1)
      do i=j,size(x,1)
        c(i,j)=model%variance*kernel_corr(model,x(i,:),x(j,:))
        c(j,i)=c(i,j)
      end do
    end do
    if(nug) then
      do i=1,size(x,1)
        c(i,i)=c(i,i)+model%nugget
      end do
    end if
  end function covariance_matrix

  function predict_nobias_km(model,newdata,type,cov_compute) result(pred)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: newdata(:,:)
    character(len=*), intent(in), optional :: type
    logical, intent(in), optional :: cov_compute
    type(krig_prediction) :: pred
    type(km_prediction) :: dkp
    real(dp), allocatable :: linvc(:,:),fnew(:,:),delta(:,:),base(:,:),cnew(:,:)
    real(dp) :: total_var
    character(len=2) :: kind
    integer :: i,m
    logical :: want_cov,is_uk
    if(.not.model%ready .or. size(newdata,2)/=model%d) then
      pred%ok=.false.; return
    end if
    m=size(newdata,1); want_cov=.false.; if(present(cov_compute)) want_cov=cov_compute
    kind='UK'; if(present(type)) kind=adjustl(type(1:min(2,len_trim(type))))
    is_uk=(trim(kind)/='SK')
    fnew=trend_matrix(model,newdata); pred%fnew=fnew

    if(model%uses_dicekriging .and. allocated(model%dice_model)) then
      call dk_km_predict(model%dice_model,newdata,fnew,kind,dkp,se_compute=.true.,cov_compute=want_cov)
      pred%mean=dkp%mean; pred%sd=dkp%sd
      if(want_cov) pred%covariance=dkp%cov
      call dk_covariance_cross(model%dice_model%covariance,model%x,newdata,pred%c, &
                               include_nugget=model%dice_model%covariance%nugget_flag)
      call solve_lower(model%l,pred%c,linvc); pred%linv_c=linvc
      return
    end if

    pred%c=transpose(covariance_cross(model,model%x,newdata))
    pred%c=transpose(pred%c)
    call solve_lower(model%l,pred%c,linvc); pred%linv_c=linvc
    allocate(pred%mean(m)); pred%mean=matmul(fnew,model%beta)+matmul(transpose(linvc),model%z)
    total_var=model%variance+model%nugget
    allocate(pred%sd(m)); pred%sd=0.0_dp
    if(is_uk) then
      delta=fnew-matmul(transpose(linvc),model%m)
      do i=1,m
        pred%sd(i)=sqrt(max(0.0_dp,total_var-dot_product(linvc(:,i),linvc(:,i))+ &
                       dot_product(delta(i,:),matmul(model%ginv,delta(i,:)))))
      end do
    else
      do i=1,m
        pred%sd(i)=sqrt(max(0.0_dp,total_var-dot_product(linvc(:,i),linvc(:,i))))
      end do
    end if
    if(want_cov) then
      cnew=covariance_matrix(model,newdata,.true.)
      base=cnew-matmul(transpose(linvc),linvc)
      if(is_uk) then
        delta=fnew-matmul(transpose(linvc),model%m)
        base=base+matmul(matmul(delta,model%ginv),transpose(delta))
      end if
      call symmetrize(base)
      do i=1,m
        base(i,i)=max(0.0_dp,base(i,i))
      end do
      pred%covariance=base
    end if
  end function predict_nobias_km

  function posterior_covariance(model,a,b) result(cab)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: a(:,:),b(:,:)
    real(dp), allocatable :: cab(:,:)
    type(km_prediction) :: dkp
    real(dp), allocatable :: cxa(:,:),cxb(:,:),la(:,:),lb(:,:),fa(:,:),fb(:,:),da(:,:),db(:,:)
    real(dp), allocatable :: ab(:,:),fab(:,:)
    integer :: na,nb
    if(model%uses_dicekriging .and. allocated(model%dice_model)) then
      na=size(a,1); nb=size(b,1)
      allocate(ab(na+nb,model%d)); ab(1:na,:)=a; ab(na+1:,:)=b
      fab=trend_matrix(model,ab)
      call dk_km_predict(model%dice_model,ab,fab,'UK',dkp,se_compute=.false.,cov_compute=.true.)
      cab=dkp%cov(1:na,na+1:na+nb)
      return
    end if
    cxa=covariance_cross(model,model%x,a)
    cxb=covariance_cross(model,model%x,b)
    call solve_lower(model%l,cxa,la)
    call solve_lower(model%l,cxb,lb)
    cab=covariance_cross(model,a,b)-matmul(transpose(la),lb)
    fa=trend_matrix(model,a); fb=trend_matrix(model,b)
    da=fa-matmul(transpose(la),model%m)
    db=fb-matmul(transpose(lb),model%m)
    cab=cab+matmul(matmul(da,model%ginv),transpose(db))
  end function posterior_covariance

  subroutine update_krig_model(model,newx,newy,newnoise,ok,cov_reestimate,trend_reestimate,nugget_reestimate)
    type(krig_model), intent(inout) :: model
    real(dp), intent(in) :: newx(:,:),newy(:)
    real(dp), intent(in), optional :: newnoise(:)
    logical, intent(out), optional :: ok
    logical, intent(in), optional :: cov_reestimate,trend_reestimate,nugget_reestimate
    type(krig_model) :: tmp
    real(dp), allocatable :: xa(:,:),ya(:),na(:),newf(:,:)
    integer :: n0,n1
    logical :: good,cr,tr,nr,pass_noise
    n0=model%n; n1=size(newx,1)
    if(size(newx,2)/=model%d .or. size(newy)/=n1) then
      if(present(ok)) ok=.false.
      return
    end if
    if(present(newnoise)) then
      if(size(newnoise)/=n1 .or. any(newnoise<0.0_dp)) then
        if(present(ok)) ok=.false.
        return
      end if
    end if

    if(model%uses_dicekriging .and. allocated(model%dice_model)) then
      cr=model%cov_reestimate_default; if(present(cov_reestimate)) cr=cov_reestimate
      tr=.true.; if(present(trend_reestimate)) tr=trend_reestimate
      nr=.false.; if(present(nugget_reestimate)) nr=nugget_reestimate
      newf=trend_matrix(model,newx)
      pass_noise=model%dice_model%noise_flag
      if(present(newnoise)) pass_noise=pass_noise .or. any(newnoise>0.0_dp)
      if(pass_noise) then
        if(present(newnoise)) then
          call dk_km_update(model%dice_model,newx,newy,newf,cr,tr,nr,newnoise)
        else
          allocate(na(n1)); na=0.0_dp
          call dk_km_update(model%dice_model,newx,newy,newf,cr,tr,nr,na)
        end if
      else
        call dk_km_update(model%dice_model,newx,newy,newf,cr,tr,nr)
      end if
      call sync_from_dice(model,good)
      if(present(ok)) ok=good
      return
    end if

    allocate(xa(n0+n1,model%d),ya(n0+n1),na(n0+n1))
    xa(1:n0,:)=model%x; xa(n0+1:,:)=newx
    ya(1:n0)=model%y; ya(n0+1:)=newy
    na(1:n0)=model%noise; na(n0+1:)=0.0_dp
    if(present(newnoise)) na(n0+1:)=max(0.0_dp,newnoise)
    call init_krig_model(tmp,xa,ya,model%lengthscale,model%variance,model%nugget,na, &
                         model%covariance,model%trend_order,good)
    if(good) model=tmp
    if(present(ok)) ok=good
  end subroutine update_krig_model
end module kriginv_model
