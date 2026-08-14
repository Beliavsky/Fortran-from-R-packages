! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! Vendored in KrigInv-fortran under the GPL-3 option; see
! licenses/DiceKriging/LICENSE-GPL-3 and UPSTREAM.md.
module dk_model
  use dk_kinds, only : dp, pi_dp
  use dk_linalg, only : chol_lower, solve_lower, solve_chol, invert_spd, least_squares_normal
  use dk_linalg, only : logdet_from_chol, diag_aba, normal_fill
  use dk_covariance, only : covariance_model, scaling_axis, covariance_kind, covariance_param_count
  use dk_covariance, only : covariance_matrix, covariance_cross, covariance_derivative
  use dk_covariance, only : covariance_bounds, covariance_vector_dx, get_cov_params, set_cov_params
  use dk_covariance, only : cov_powexp, cov_gauss
  use dk_optimizer, only : bounded_bfgs
  implicit none
  private

  integer, parameter :: case_default=1, case_nugget=2, case_noisy=3

  type, public :: km_control
    integer :: multistart = 1
    integer :: pop_size = 20
    integer :: max_iter = 300
    real(dp) :: tol = 1.0e-7_dp
    logical :: use_gradient = .true.
    real(dp) :: upper_alpha = 1.0_dp-1.0e-8_dp
  end type km_control

  type, public :: km_prediction
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: sd(:)
    real(dp), allocatable :: cov(:,:)
    real(dp), allocatable :: trend(:)
    real(dp), allocatable :: lower95(:), upper95(:)
  end type km_prediction

  type, public :: km_model
    real(dp), allocatable :: x(:,:), y(:), f(:,:)
    real(dp), allocatable :: noise_var(:)
    type(covariance_model) :: covariance
    real(dp), allocatable :: trend_coef(:)
    real(dp), allocatable :: l(:,:), m(:,:), z(:)
    real(dp), allocatable :: lower(:), upper(:), parinit(:)
    integer :: n=0, d=0, p=0
    logical :: noise_flag=.false.
    logical :: trend_known=.false.
    logical :: param_estim=.true.
    logical :: estimate_cov=.true.
    logical :: estimate_var=.true.
    logical :: estimate_trend=.true.
    logical :: nugget_reestim=.false.
    character(len=8) :: method='MLE'
    character(len=8) :: optim_method='BFGS'
    integer :: fit_case=case_default
    real(dp) :: loglik=-huge(1.0_dp)
    logical :: pmle=.false.
    real(dp) :: scad_lambda=0.0_dp
    type(km_control) :: control
  end type km_model

  public :: km_fit, km_fit_cov, km_estimate, km_recompute, km_predict, km_simulate, km_update, km_update_response
  public :: loglik_fun, loglik_grad, leave_one_out_fun, leave_one_out_grad
  public :: leave_one_out, km_loglik, km_cov_vector_dx
  public :: trend_constant, trend_linear, trend_linear_interactions, trend_quadratic
  public :: trend_gradient_constant, trend_gradient_linear, trend_gradient_linear_interactions, trend_gradient_quadratic
  public :: scad, scad_derivative

contains

  subroutine km_fit(model,x,y,f,covtype,coef_cov,coef_var,coef_trend,nugget,nugget_estim,noise_var, &
                    estim_method,iso,lower,upper,parinit,control,scaling_axes,scad_lambda)
    type(km_model), intent(out) :: model
    real(dp), intent(in) :: x(:,:), y(:), f(:,:)
    character(len=*), intent(in) :: covtype
    real(dp), intent(in), optional :: coef_cov(:), coef_var, coef_trend(:), nugget, noise_var(:)
    logical, intent(in), optional :: nugget_estim, iso
    character(len=*), intent(in), optional :: estim_method
    real(dp), intent(in), optional :: lower(:),upper(:),parinit(:),scad_lambda
    type(km_control), intent(in), optional :: control
    type(scaling_axis), intent(in), optional :: scaling_axes(:)
    type(covariance_model) :: cov
    integer :: d
    logical :: is_iso

    d=size(x,2); cov%kind=covariance_kind(covtype)
    if(cov%kind==0) error stop 'km_fit: unsupported covariance type'
    is_iso=.false.; if(present(iso))is_iso=iso
    cov%iso=is_iso
    if(present(scaling_axes)) then
      cov%scaling=.true.; cov%iso=.false.; cov%axis=scaling_axes
      allocate(cov%range(d)); cov%range=1.0_dp
    else
      if(is_iso) then
        allocate(cov%range(d)); cov%range=1.0_dp
      else
        allocate(cov%range(d)); cov%range=1.0_dp
      end if
      if(cov%kind==cov_powexp) then
        allocate(cov%shape(d)); cov%shape=1.5_dp
      end if
    end if
    cov%sd2=1.0_dp
    if(present(coef_var)) cov%sd2=coef_var
    cov%nugget_flag=.false.; cov%nugget_estim=.false.; cov%nugget=0.0_dp
    if(present(nugget)) then
      cov%nugget_flag=.true.; cov%nugget=nugget
    end if
    if(present(nugget_estim)) then
      cov%nugget_estim=nugget_estim
      if(nugget_estim) cov%nugget_flag=.true.
    end if
    if(present(coef_cov)) call set_cov_params(cov,d,coef_cov)

    call km_fit_cov(model,x,y,f,cov,coef_trend=coef_trend, &
      estimate_cov=.not.present(coef_cov), estimate_var=.not.present(coef_var), &
      estimate_trend=.not.present(coef_trend), noise_var=noise_var, estim_method=estim_method, &
      lower=lower,upper=upper,parinit=parinit,control=control,scad_lambda=scad_lambda)

  end subroutine km_fit

  subroutine km_fit_cov(model,x,y,f,cov,coef_trend,estimate_cov,estimate_var,estimate_trend, &
                        noise_var,estim_method,lower,upper,parinit,control,scad_lambda)
    type(km_model), intent(out) :: model
    real(dp), intent(in) :: x(:,:),y(:),f(:,:)
    type(covariance_model), intent(in) :: cov
    real(dp), intent(in), optional :: coef_trend(:),noise_var(:),lower(:),upper(:),parinit(:),scad_lambda
    logical, intent(in), optional :: estimate_cov,estimate_var,estimate_trend
    character(len=*), intent(in), optional :: estim_method
    type(km_control), intent(in), optional :: control
    logical :: ec,ev,et

    if(size(x,1)/=size(y) .or. size(f,1)/=size(y)) error stop 'km_fit_cov: row mismatch'
    model%x=x; model%y=y; model%f=f; model%n=size(x,1); model%d=size(x,2); model%p=size(f,2)
    model%covariance=cov
    ec=.true.; ev=.true.; et=.true.
    if(present(estimate_cov))ec=estimate_cov
    if(present(estimate_var))ev=estimate_var
    if(present(estimate_trend))et=estimate_trend
    model%estimate_cov=ec; model%estimate_var=ev; model%estimate_trend=et
    model%trend_known=.not.et
    if(present(coef_trend)) then
      model%trend_coef=coef_trend; model%trend_known=.true.; model%estimate_trend=.false.
    end if
    if(present(noise_var)) then
      if(size(noise_var)/=model%n) error stop 'km_fit_cov: noise_var length mismatch'
      model%noise_var=noise_var; model%noise_flag=.true.
    end if
    if(model%noise_flag .and. model%covariance%nugget_flag) error stop 'km_fit_cov: nugget and noise cannot both be used'
    if(present(estim_method)) model%method=adjustl(estim_method)
    if(present(control)) model%control=control
    if(present(scad_lambda)) then
      if(model%covariance%kind/=cov_gauss) error stop 'km_fit_cov: SCAD PMLE requires Gaussian covariance'
      model%pmle=.true.; model%scad_lambda=scad_lambda; model%method='PMLE'
    end if
    if(model%covariance%nugget_estim) then
      model%fit_case=case_nugget
    else if(model%noise_flag .or. model%covariance%nugget_flag) then
      model%fit_case=case_noisy
    else
      model%fit_case=case_default
    end if
    if(trim(model%method)=='LOO' .and. model%fit_case/=case_default) error stop 'km_fit_cov: LOO requires no nugget/noise'

    if(present(lower)) model%lower=lower
    if(present(upper)) model%upper=upper
    if(present(parinit)) model%parinit=parinit
    model%param_estim=ec .or. ev .or. et .or. model%covariance%nugget_estim

    if(.not.ec .and. .not.ev) then
      call km_recompute(model,reestimate_trend=et)
      model%loglik=km_loglik(model)
      return
    end if
    call km_estimate(model)
  end subroutine km_fit_cov

  subroutine km_estimate(model)
    type(km_model), intent(inout), target :: model
    real(dp), allocatable :: covlo(:),covup(:),lo(:),up(:),starts(:,:),vals(:),x0(:),xb(:),best(:)
    real(dp), allocatable :: p0(:)
    real(dp) :: fb,bestf,var0
    integer :: ncp,npop,nstart,i,j,info,besti

    ncp=covariance_param_count(model%covariance,model%d)
    call covariance_bounds(model%covariance,model%x,covlo,covup)
    if(allocated(model%lower)) covlo=model%lower(1:ncp)
    if(allocated(model%upper)) covup=model%upper(1:ncp)
    select case(model%fit_case)
    case(case_default)
      allocate(lo(ncp),up(ncp)); lo=covlo; up=covup
    case(case_nugget)
      allocate(lo(ncp+1),up(ncp+1)); lo(1:ncp)=covlo; up(1:ncp)=covup
      lo(ncp+1)=1.0e-8_dp; up(ncp+1)=model%control%upper_alpha
    case(case_noisy)
      allocate(lo(ncp+1),up(ncp+1)); lo(1:ncp)=covlo; up(1:ncp)=covup
      var0=residual_variance(model%f,model%y)
      lo(ncp+1)=max(1.0e-12_dp,1.0e-3_dp*var0)
      up(ncp+1)=max(lo(ncp+1)*100.0_dp,100.0_dp*max(var0,1.0e-10_dp))
    end select
    if(allocated(model%lower)) then
      if(size(model%lower)==size(lo))lo=model%lower
    end if
    if(allocated(model%upper)) then
      if(size(model%upper)==size(up))up=model%upper
    end if
    model%lower=lo; model%upper=up

    npop=max(model%control%pop_size,model%control%multistart)
    allocate(starts(size(lo),npop),vals(npop))
    call random_number(starts)
    do j=1,npop
      if (model%covariance%scaling) then
        starts(1:ncp,j)=1.0_dp/(1.0_dp/up(1:ncp)+starts(1:ncp,j)*(1.0_dp/lo(1:ncp)-1.0_dp/up(1:ncp)))
        if (size(lo)>ncp) starts(ncp+1:,j)=lo(ncp+1:)+starts(ncp+1:,j)*(up(ncp+1:)-lo(ncp+1:))
      else
        starts(:,j)=lo+starts(:,j)*(up-lo)
      end if
    end do
    call get_cov_params(model%covariance,model%d,p0)
    if(size(p0)==ncp) starts(1:ncp,1)=min(max(p0,covlo),covup)
    if(model%fit_case==case_nugget) starts(ncp+1,1)=0.9_dp
    if(model%fit_case==case_noisy) &
      starts(ncp+1,1)=max(lo(ncp+1),min(up(ncp+1),residual_variance(model%f,model%y)))
    if(allocated(model%parinit)) then
      if(size(model%parinit)==ncp) starts(1:ncp,1)=min(max(model%parinit,covlo),covup)
      if(size(model%parinit)==size(lo)) starts(:,1)=min(max(model%parinit,lo),up)
    end if
    do j=1,npop
      vals(j)=objective(starts(:,j))
    end do
    nstart=min(model%control%multistart,npop)
    bestf=huge(1.0_dp); besti=1
    do i=1,nstart
      j=minloc(vals,dim=1)
      x0=starts(:,j); vals(j)=huge(1.0_dp)
      if(model%control%use_gradient) then
        call bounded_bfgs(objective,x0,lo,up,xb,fb,info,gradient=objective_grad, &
          max_iter=model%control%max_iter,tol=model%control%tol)
      else
        call bounded_bfgs(objective,x0,lo,up,xb,fb,info,max_iter=model%control%max_iter,tol=model%control%tol)
      end if
      if(fb<bestf) then
        bestf=fb; best=xb; besti=i
      end if
    end do
    model%parinit=best

    if(trim(model%method)=='LOO') then
      call finalize_loo(model,best)
      model%loglik=-bestf
    else
      call finalize_mle(model,best)
      model%loglik=-bestf
    end if
    model%param_estim=.true.

  contains
    function objective(q) result(v)
      real(dp), intent(in) :: q(:)
      real(dp) :: v
      if(trim(model%method)=='LOO') then
        v=leave_one_out_fun(model,q)
      else
        v=-loglik_fun(model,q)
      end if
      if(.not.(v<huge(1.0_dp))) v=huge(1.0_dp)/10.0_dp
    end function objective
    subroutine objective_grad(q,g)
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: g(:)
      if(trim(model%method)=='LOO') then
        call leave_one_out_grad(model,q,g)
      else
        call loglik_grad(model,q,g); g=-g
      end if
    end subroutine objective_grad
  end subroutine km_estimate

  function loglik_fun(model,param) result(ll)
    type(km_model), intent(in) :: model
    real(dp), intent(in) :: param(:)
    real(dp) :: ll
    type(covariance_model) :: cov
    real(dp), allocatable :: c(:,:),l(:,:),xx(:,:),mm(:,:),beta(:,:),res(:,:)
    real(dp) :: sigma2,v,alpha
    integer :: ncp,info
    ll=-huge(1.0_dp)
    ncp=covariance_param_count(model%covariance,model%d); cov=model%covariance
    call set_cov_params(cov,model%d,param(1:ncp))
    select case(model%fit_case)
    case(case_default)
      cov%sd2=1.0_dp; cov%nugget_flag=.false.
      call covariance_matrix(cov,model%x,c,include_nugget=.false.)
      call chol_lower(c,l,info)
      if(info/=0) then; ll=-huge(1.0_dp); return; end if
      call profile_whitened(model,l,xx,mm,beta,res)
      sigma2=dot_product(res(:,1),res(:,1))/real(model%n,dp)
      if(sigma2<=tiny(1.0_dp)) then; ll=-huge(1.0_dp); return; end if
      ll=-0.5_dp*(real(model%n,dp)*log(2.0_dp*pi_dp*sigma2)+logdet_from_chol(l)+real(model%n,dp))
    case(case_nugget)
      cov%sd2=1.0_dp; cov%nugget=0.0_dp; cov%nugget_flag=.false.; alpha=param(ncp+1)
      call covariance_matrix(cov,model%x,c,include_nugget=.false.)
      c=alpha*c
      c=c+(1.0_dp-alpha)*identity_matrix(model%n)
      call chol_lower(c,l,info)
      if(info/=0) then; ll=-huge(1.0_dp); return; end if
      call profile_whitened(model,l,xx,mm,beta,res)
      v=dot_product(res(:,1),res(:,1))/real(model%n,dp)
      if(v<=tiny(1.0_dp)) then; ll=-huge(1.0_dp); return; end if
      ll=-0.5_dp*(real(model%n,dp)*log(2.0_dp*pi_dp*v)+logdet_from_chol(l)+real(model%n,dp))
    case(case_noisy)
      cov%sd2=param(ncp+1)
      if(model%noise_flag) then
        call covariance_matrix(cov,model%x,c,noise_var=model%noise_var,include_nugget=.false.)
      else
        call covariance_matrix(cov,model%x,c,include_nugget=.true.)
      end if
      call chol_lower(c,l,info)
      if(info/=0) then; ll=-huge(1.0_dp); return; end if
      call profile_whitened(model,l,xx,mm,beta,res)
      ll=-0.5_dp*(real(model%n,dp)*log(2.0_dp*pi_dp)+logdet_from_chol(l)+dot_product(res(:,1),res(:,1)))
    end select
    if(model%pmle) ll=ll-real(model%n,dp)*sum(scad(1.0_dp/cov%range**2,model%scad_lambda))
  end function loglik_fun

  subroutine loglik_grad(model,param,g)
    type(km_model), intent(in) :: model
    real(dp), intent(in) :: param(:)
    real(dp), intent(out) :: g(:)
    type(covariance_model) :: cov
    real(dp), allocatable :: c(:,:),r0(:,:),l(:,:),cinv(:,:),xx(:,:),mm(:,:),beta(:,:),res(:,:),xinv(:,:)
    real(dp), allocatable :: dc(:,:)
    real(dp) :: sigma2,v,alpha,term1,term2
    integer :: ncp,k,info
    ncp=covariance_param_count(model%covariance,model%d); g=0.0_dp; cov=model%covariance
    call set_cov_params(cov,model%d,param(1:ncp))
    select case(model%fit_case)
    case(case_default)
      cov%sd2=1.0_dp; cov%nugget_flag=.false.; call covariance_matrix(cov,model%x,c,include_nugget=.false.)
      call chol_lower(c,l,info); if(info/=0)return
      call profile_whitened(model,l,xx,mm,beta,res); sigma2=dot_product(res(:,1),res(:,1))/real(model%n,dp)
      call invert_spd(l,cinv); call solve_chol(l,reshape(model%y,[model%n,1])-matmul(model%f,beta),xinv)
      do k=1,ncp
        call covariance_derivative(cov,model%x,c,k,dc)
        term1=dot_product(xinv(:,1),matmul(dc,xinv(:,1)))/sigma2
        term2=sum(cinv*transpose(dc))
        g(k)=0.5_dp*(term1-term2)
      end do
    case(case_nugget)
      cov%sd2=1.0_dp; cov%nugget=0.0_dp; cov%nugget_flag=.false.; alpha=param(ncp+1)
      call covariance_matrix(cov,model%x,r0,include_nugget=.false.); c=alpha*r0+(1.0_dp-alpha)*identity_matrix(model%n)
      call chol_lower(c,l,info); if(info/=0)return
      call profile_whitened(model,l,xx,mm,beta,res); v=dot_product(res(:,1),res(:,1))/real(model%n,dp)
      call invert_spd(l,cinv); call solve_chol(l,reshape(model%y,[model%n,1])-matmul(model%f,beta),xinv)
      do k=1,ncp
        call covariance_derivative(cov,model%x,r0,k,dc); dc=alpha*dc
        g(k)=0.5_dp*(dot_product(xinv(:,1),matmul(dc,xinv(:,1)))/v-sum(cinv*transpose(dc)))
      end do
      dc=r0-identity_matrix(model%n)
      g(ncp+1)=0.5_dp*(dot_product(xinv(:,1),matmul(dc,xinv(:,1)))/v-sum(cinv*transpose(dc)))
    case(case_noisy)
      cov%sd2=param(ncp+1)
      if(model%noise_flag) then
        call covariance_matrix(cov,model%x,c,noise_var=model%noise_var,include_nugget=.false.)
      else
        call covariance_matrix(cov,model%x,c,include_nugget=.true.)
      end if
      call chol_lower(c,l,info); if(info/=0)return
      call profile_whitened(model,l,xx,mm,beta,res); call invert_spd(l,cinv)
      call solve_chol(l,reshape(model%y,[model%n,1])-matmul(model%f,beta),xinv)
      do k=1,ncp+1
        call covariance_derivative(cov,model%x,c-diagonal_noise(model,cov),k,dc)
        g(k)=0.5_dp*(dot_product(xinv(:,1),matmul(dc,xinv(:,1)))-sum(cinv*transpose(dc)))
      end do
    end select
    if(model%pmle) then
      do k=1,min(ncp,size(cov%range))
        g(k)=g(k)-real(model%n,dp)*scad_derivative(1.0_dp/cov%range(k)**2,model%scad_lambda)*(-2.0_dp/cov%range(k)**3)
      end do
    end if
  end subroutine loglik_grad

  function leave_one_out_fun(model,param) result(v)
    type(km_model), intent(in) :: model
    real(dp), intent(in) :: param(:)
    real(dp) :: v
    type(covariance_model) :: cov
    real(dp), allocatable :: r(:,:),l(:,:),rinv(:,:),q(:,:),qy(:,:),a(:,:),ata(:,:),li(:,:),tmp(:,:),err(:),s2(:)
    integer :: info,ncp,i
    ncp=covariance_param_count(model%covariance,model%d); cov=model%covariance
    call set_cov_params(cov,model%d,param(1:ncp)); cov%sd2=1.0_dp; cov%nugget_flag=.false.
    call covariance_matrix(cov,model%x,r,include_nugget=.false.); call chol_lower(r,l,info)
    if(info/=0) then; v=huge(1.0_dp); return; end if
    call invert_spd(l,rinv)
    if(model%trend_known) then
      q=rinv; qy=matmul(q,reshape(model%y,[model%n,1])-matmul(model%f,reshape(model%trend_coef,[model%p,1])))
    else
      a=matmul(rinv,model%f); ata=matmul(transpose(model%f),a)
      call chol_lower(ata,li,info); if(info/=0) then; v=huge(1.0_dp); return; end if
      call solve_chol(li,transpose(a),tmp)
      q=rinv-matmul(a,tmp); qy=matmul(q,reshape(model%y,[model%n,1]))
    end if
    allocate(err(model%n),s2(model%n))
    do i=1,model%n
      if(q(i,i)<=tiny(1.0_dp)) then; v=huge(1.0_dp); return; end if
      s2(i)=1.0_dp/q(i,i); err(i)=s2(i)*qy(i,1)
    end do
    v=dot_product(err,err)/real(model%n,dp)
  end function leave_one_out_fun

  subroutine leave_one_out_grad(model,param,g)
    type(km_model), intent(in) :: model
    real(dp), intent(in) :: param(:)
    real(dp), intent(out) :: g(:)
    type(covariance_model) :: cov
    real(dp), allocatable :: r(:,:),l(:,:),rinv(:,:),q(:,:),qy(:,:),a(:,:),ata(:,:),li(:,:),tmp(:,:),dc(:,:)
    real(dp), allocatable :: err(:),s2(:),dqdiag(:),ds2(:),derr(:)
    integer :: info,ncp,i,k
    ncp=covariance_param_count(model%covariance,model%d); g=0.0_dp; cov=model%covariance
    call set_cov_params(cov,model%d,param(1:ncp)); cov%sd2=1.0_dp; cov%nugget_flag=.false.
    call covariance_matrix(cov,model%x,r,include_nugget=.false.); call chol_lower(r,l,info); if(info/=0)return
    call invert_spd(l,rinv)
    if(model%trend_known) then
      q=rinv; qy=matmul(q,reshape(model%y,[model%n,1])-matmul(model%f,reshape(model%trend_coef,[model%p,1])))
    else
      a=matmul(rinv,model%f); ata=matmul(transpose(model%f),a); call chol_lower(ata,li,info); if(info/=0)return
      call solve_chol(li,transpose(a),tmp); q=rinv-matmul(a,tmp); qy=matmul(q,reshape(model%y,[model%n,1]))
    end if
    allocate(err(model%n),s2(model%n),dqdiag(model%n),ds2(model%n),derr(model%n))
    do i=1,model%n
      s2(i)=1.0_dp/q(i,i); err(i)=s2(i)*qy(i,1)
    end do
    do k=1,ncp
      call covariance_derivative(cov,model%x,r,k,dc)
      dqdiag=-diag_aba(q,dc); ds2=-(s2*s2)*dqdiag
      derr=ds2*qy(:,1)-s2*matmul(q,matmul(dc,qy(:,1)))
      g(k)=2.0_dp*dot_product(err,derr)/real(model%n,dp)
    end do
  end subroutine leave_one_out_grad

  subroutine finalize_mle(model,param)
    type(km_model), intent(inout) :: model
    real(dp), intent(in) :: param(:)
    type(covariance_model) :: cov
    real(dp), allocatable :: c(:,:),l(:,:),xx(:,:),mm(:,:),beta(:,:),res(:,:),r0(:,:)
    real(dp) :: sigma2,v,alpha
    integer :: ncp,info
    ncp=covariance_param_count(model%covariance,model%d); cov=model%covariance
    call set_cov_params(cov,model%d,param(1:ncp))
    select case(model%fit_case)
    case(case_default)
      cov%sd2=1.0_dp; cov%nugget_flag=.false.; call covariance_matrix(cov,model%x,c,include_nugget=.false.)
      call chol_lower(c,l,info); call profile_whitened(model,l,xx,mm,beta,res)
      sigma2=dot_product(res(:,1),res(:,1))/real(model%n,dp); cov%sd2=sigma2
    case(case_nugget)
      cov%sd2=1.0_dp; cov%nugget=0.0_dp; cov%nugget_flag=.false.; alpha=param(ncp+1)
      call covariance_matrix(cov,model%x,r0,include_nugget=.false.); c=alpha*r0+(1.0_dp-alpha)*identity_matrix(model%n)
      call chol_lower(c,l,info); call profile_whitened(model,l,xx,mm,beta,res)
      v=dot_product(res(:,1),res(:,1))/real(model%n,dp); cov%sd2=alpha*v; cov%nugget=(1.0_dp-alpha)*v
      cov%nugget_flag=.true.; cov%nugget_estim=.true.
    case(case_noisy)
      cov%sd2=param(ncp+1)
    end select
    model%covariance=cov
    call km_recompute(model,reestimate_trend=.not.model%trend_known)
    if (allocated(model%lower)) model%lower=model%lower(1:ncp)
    if (allocated(model%upper)) model%upper=model%upper(1:ncp)
    if (allocated(model%parinit)) model%parinit=model%parinit(1:ncp)
  end subroutine finalize_mle

  subroutine finalize_loo(model,param)
    type(km_model), intent(inout) :: model
    real(dp), intent(in) :: param(:)
    type(covariance_model) :: cov
    real(dp), allocatable :: r(:,:),l(:,:),rinv(:,:),q(:,:),qy(:,:),a(:,:),ata(:,:),li(:,:),tmp(:,:),err(:),s2(:)
    real(dp) :: sigma2
    integer :: ncp,info,i
    ncp=covariance_param_count(model%covariance,model%d); cov=model%covariance
    call set_cov_params(cov,model%d,param(1:ncp)); cov%sd2=1.0_dp; cov%nugget_flag=.false.
    call covariance_matrix(cov,model%x,r,include_nugget=.false.); call chol_lower(r,l,info); call invert_spd(l,rinv)
    if(model%trend_known) then
      q=rinv; qy=matmul(q,reshape(model%y,[model%n,1])-matmul(model%f,reshape(model%trend_coef,[model%p,1])))
    else
      a=matmul(rinv,model%f); ata=matmul(transpose(model%f),a); call chol_lower(ata,li,info)
      call solve_chol(li,transpose(a),tmp); q=rinv-matmul(a,tmp); qy=matmul(q,reshape(model%y,[model%n,1]))
    end if
    allocate(err(model%n),s2(model%n))
    do i=1,model%n
      s2(i)=1.0_dp/q(i,i); err(i)=s2(i)*qy(i,1)
    end do
    sigma2=sum(err*err/s2)/real(model%n,dp); cov%sd2=sigma2; model%covariance=cov
    call km_recompute(model,reestimate_trend=.not.model%trend_known)
    if (allocated(model%lower)) model%lower=model%lower(1:ncp)
    if (allocated(model%upper)) model%upper=model%upper(1:ncp)
    if (allocated(model%parinit)) model%parinit=model%parinit(1:ncp)
  end subroutine finalize_loo

  subroutine km_recompute(model,reestimate_trend)
    type(km_model), intent(inout) :: model
    logical, intent(in), optional :: reestimate_trend
    real(dp), allocatable :: c(:,:),ym(:,:),b(:,:)
    integer :: info
    logical :: rt
    rt=model%estimate_trend; if(present(reestimate_trend))rt=reestimate_trend
    if(model%noise_flag) then
      call covariance_matrix(model%covariance,model%x,c,noise_var=model%noise_var,include_nugget=.false.)
    else
      call covariance_matrix(model%covariance,model%x,c,include_nugget=.true.)
    end if
    call chol_lower(c,model%l,info); if(info/=0) error stop 'km_recompute: covariance is not positive definite'
    ym=reshape(model%y,[model%n,1]); call solve_lower(model%l,ym,b); model%z=b(:,1)
    call solve_lower(model%l,model%f,model%m)
    if(rt .or. .not.allocated(model%trend_coef)) then
      call least_squares_normal(model%m,reshape(model%z,[model%n,1]),b,info)
      if(info/=0) error stop 'km_recompute: singular trend matrix'
      model%trend_coef=b(:,1)
    end if
    model%z=model%z-matmul(model%m,model%trend_coef)
  end subroutine km_recompute

  subroutine km_predict(model,newx,newf,kind,pred,se_compute,cov_compute,bias_correct)
    type(km_model), intent(in) :: model
    real(dp), intent(in) :: newx(:,:),newf(:,:)
    character(len=*), intent(in) :: kind
    type(km_prediction), intent(out) :: pred
    logical, intent(in), optional :: se_compute,cov_compute,bias_correct
    real(dp), allocatable :: c(:,:),a(:,:),delta(:,:),tm(:,:),knew(:,:),cc(:,:)
    real(dp) :: q95
    logical :: se,co,bc
    integer :: j,nnew
    se=.true.; if(present(se_compute))se=se_compute
    co=.false.; if(present(cov_compute))co=cov_compute
    bc=.false.; if(present(bias_correct))bc=bias_correct
    nnew=size(newx,1)
    if(size(newx,2)/=model%d .or. size(newf,2)/=model%p) error stop 'km_predict: dimension mismatch'
    call covariance_cross(model%covariance,model%x,newx,c,include_nugget=model%covariance%nugget_flag)
    call solve_lower(model%l,c,a)
    allocate(pred%trend(nnew),pred%mean(nnew)); pred%trend=matmul(newf,model%trend_coef)
    pred%mean=pred%trend+matmul(transpose(a),model%z)
    if(se .or. co) then
      call covariance_matrix(model%covariance,newx,knew,include_nugget=.true.)
      cc=knew-matmul(transpose(a),a)
      if(trim(kind)=='UK') then
        delta=newf-matmul(transpose(a),model%m)
        call uk_correction(model%m,delta,tm)
        cc=cc+matmul(transpose(tm),tm)
        if(bc) cc=cc*real(model%n,dp)/real(model%n-model%p,dp)
      end if
      if(co) pred%cov=cc
      if(se) then
        allocate(pred%sd(nnew),pred%lower95(nnew),pred%upper95(nnew))
        do j=1,nnew
          pred%sd(j)=sqrt(max(cc(j,j),0.0_dp))
        end do
        if(trim(kind)=='UK') then
          q95=t_quantile_975(model%n-model%p)
        else
          q95=1.959963984540054_dp
        end if
        pred%lower95=pred%mean-q95*pred%sd
        pred%upper95=pred%mean+q95*pred%sd
      end if
    end if
  end subroutine km_predict

  subroutine uk_correction(m,delta,t)
    real(dp), intent(in) :: m(:,:),delta(:,:)
    real(dp), allocatable, intent(out) :: t(:,:)
    real(dp), allocatable :: mtm(:,:),l(:,:),rhs(:,:)
    integer :: info
    mtm=matmul(transpose(m),m); call chol_lower(mtm,l,info)
    if(info/=0) error stop 'uk_correction: singular trend information'
    call solve_lower(l,transpose(delta),rhs)
    t=rhs
  end subroutine uk_correction

  subroutine km_simulate(model,nsim,newx,newf,draws,conditional,nugget_sim)
    type(km_model), intent(in) :: model
    integer, intent(in) :: nsim
    real(dp), intent(in) :: newx(:,:),newf(:,:)
    real(dp), allocatable, intent(out) :: draws(:,:)
    logical, intent(in), optional :: conditional
    real(dp), intent(in), optional :: nugget_sim
    logical :: cond
    real(dp) :: ns
    real(dp), allocatable :: k(:,:),l(:,:),zrand(:,:),mean(:),c(:,:),a(:,:),kc(:,:)
    integer :: info,i
    cond=.false.; if(present(conditional))cond=conditional
    ns=0.0_dp; if(present(nugget_sim))ns=nugget_sim
    call covariance_matrix(model%covariance,newx,k,include_nugget=.true.)
    do i=1,size(k,1); k(i,i)=k(i,i)+ns; end do
    if(cond) then
      call covariance_cross(model%covariance,model%x,newx,c,include_nugget=.false.)
      call solve_lower(model%l,c,a)
      mean=matmul(newf,model%trend_coef)+matmul(transpose(a),model%z)
      kc=k-matmul(transpose(a),a)
    else
      mean=matmul(newf,model%trend_coef); kc=k
    end if
    call chol_lower(kc,l,info,jitter=1.0e-12_dp*max(1.0_dp,maxval(abs(kc))))
    if(info/=0) error stop 'km_simulate: conditional covariance is not positive definite'
    allocate(zrand(size(newx,1),nsim)); call normal_fill(zrand)
    allocate(draws(nsim,size(newx,1)))
    draws=transpose(spread(mean,2,nsim)+matmul(l,zrand))
  end subroutine km_simulate

  subroutine km_update(model,newx,newy,newf,cov_reestimate,trend_reestimate,nugget_reestimate,newnoise_var)
    type(km_model), intent(inout) :: model
    real(dp), intent(in) :: newx(:,:),newy(:),newf(:,:)
    logical, intent(in), optional :: cov_reestimate,trend_reestimate,nugget_reestimate
    real(dp), intent(in), optional :: newnoise_var(:)
    logical :: cr,tr,nr
    real(dp), allocatable :: xall(:,:),yall(:),fall(:,:),nall(:)
    type(covariance_model) :: cov0
    real(dp), allocatable :: beta0(:),oldlo(:),oldup(:)
    cr=.true.; tr=.true.; nr=.false.
    if(present(cov_reestimate))cr=cov_reestimate
    if(present(trend_reestimate))tr=trend_reestimate
    if(present(nugget_reestimate))nr=nugget_reestimate
    call append_rows(model%x,newx,xall); yall=[model%y,newy]; call append_rows(model%f,newf,fall)
    if(model%noise_flag .or. present(newnoise_var)) then
      allocate(nall(size(yall))); nall=0.0_dp
      if(model%noise_flag)nall(1:model%n)=model%noise_var
      if(present(newnoise_var))nall(model%n+1:)=newnoise_var
    end if
    cov0=model%covariance
    if (allocated(model%lower)) oldlo=model%lower
    if (allocated(model%upper)) oldup=model%upper
    if(nr) cov0%nugget_estim=.true.
    if((cr .or. nr) .and. model%param_estim) then
      if(tr) then
        if(allocated(nall)) then
          if(allocated(oldlo) .and. allocated(oldup)) then
            call km_fit_cov(model,xall,yall,fall,cov0,estimate_cov=cr,estimate_var=cr,estimate_trend=.true., &
              noise_var=nall,estim_method=model%method,control=model%control,lower=oldlo,upper=oldup)
          else
            call km_fit_cov(model,xall,yall,fall,cov0,estimate_cov=cr,estimate_var=cr,estimate_trend=.true., &
              noise_var=nall,estim_method=model%method,control=model%control)
          end if
        else
          if(allocated(oldlo) .and. allocated(oldup)) then
            call km_fit_cov(model,xall,yall,fall,cov0,estimate_cov=cr,estimate_var=cr,estimate_trend=.true., &
              estim_method=model%method,control=model%control,lower=oldlo,upper=oldup)
          else
            call km_fit_cov(model,xall,yall,fall,cov0,estimate_cov=cr,estimate_var=cr,estimate_trend=.true., &
              estim_method=model%method,control=model%control)
          end if
        end if
      else
        beta0=model%trend_coef
        if(allocated(nall)) then
          if(allocated(oldlo) .and. allocated(oldup)) then
            call km_fit_cov(model,xall,yall,fall,cov0,coef_trend=beta0,estimate_cov=cr,estimate_var=cr, &
              estimate_trend=.false.,noise_var=nall,estim_method=model%method,control=model%control,lower=oldlo,upper=oldup)
          else
            call km_fit_cov(model,xall,yall,fall,cov0,coef_trend=beta0,estimate_cov=cr,estimate_var=cr, &
              estimate_trend=.false.,noise_var=nall,estim_method=model%method,control=model%control)
          end if
        else
          if(allocated(oldlo) .and. allocated(oldup)) then
            call km_fit_cov(model,xall,yall,fall,cov0,coef_trend=beta0,estimate_cov=cr,estimate_var=cr, &
              estimate_trend=.false.,estim_method=model%method,control=model%control,lower=oldlo,upper=oldup)
          else
            call km_fit_cov(model,xall,yall,fall,cov0,coef_trend=beta0,estimate_cov=cr,estimate_var=cr, &
              estimate_trend=.false.,estim_method=model%method,control=model%control)
          end if
        end if
      end if
    else
      model%x=xall; model%y=yall; model%f=fall; model%n=size(yall)
      if(allocated(nall)) then; model%noise_var=nall; model%noise_flag=.true.; end if
      model%estimate_trend=tr
      call km_recompute(model,reestimate_trend=tr)
      model%loglik=km_loglik(model)
    end if
  end subroutine km_update

  subroutine km_update_response(model,newy)
    type(km_model), intent(inout) :: model
    real(dp), intent(in) :: newy(:)
    integer :: k, first
    real(dp), allocatable :: ym(:,:), b(:,:)
    if(size(newy)>model%n) error stop 'km_update_response: too many responses'
    k=size(newy); first=model%n-k+1
    model%y(first:model%n)=newy
    ym=reshape(model%y,[model%n,1])
    call solve_lower(model%l,ym,b)
    model%z=b(:,1)-matmul(model%m,model%trend_coef)
    model%loglik=km_loglik(model)
  end subroutine km_update_response

  subroutine leave_one_out(model,kind,mean,sd,trend_reestimate)
    type(km_model), intent(in) :: model
    character(len=*), intent(in) :: kind
    real(dp), allocatable, intent(out) :: mean(:),sd(:)
    logical, intent(in), optional :: trend_reestimate
    real(dp), allocatable :: c(:,:),cinv(:,:),q(:,:),qy(:,:),a(:,:),ata(:,:),l2(:,:),tmp(:,:)
    logical :: tr
    integer :: i,info
    tr=.false.; if(present(trend_reestimate))tr=trend_reestimate
    if(model%noise_flag) error stop 'leave_one_out: noisy observations are not supported'
    c=matmul(model%l,transpose(model%l)); call invert_spd(model%l,cinv)
    if(trim(kind)=='UK' .and. tr) then
      a=matmul(cinv,model%f); ata=matmul(transpose(model%f),a); call chol_lower(ata,l2,info)
      call solve_chol(l2,transpose(a),tmp); q=cinv-matmul(a,tmp); qy=matmul(q,reshape(model%y,[model%n,1]))
    else if(trim(kind)=='SK' .and. .not.tr) then
      q=cinv; qy=matmul(q,reshape(model%y,[model%n,1])-matmul(model%f,reshape(model%trend_coef,[model%p,1])))
    else
      call loo_slow(model,kind,tr,mean,sd); return
    end if
    allocate(mean(model%n),sd(model%n))
    do i=1,model%n
      sd(i)=sqrt(max(1.0_dp/q(i,i),0.0_dp)); mean(i)=model%y(i)-(1.0_dp/q(i,i))*qy(i,1)
    end do
  end subroutine leave_one_out

  subroutine loo_slow(model,kind,tr,mean,sd)
    type(km_model),intent(in)::model
    character(len=*),intent(in)::kind
    logical,intent(in)::tr
    real(dp),allocatable,intent(out)::mean(:),sd(:)
    type(km_model)::m2 = km_model()
    type(km_prediction)::pr
    real(dp),allocatable::xx(:,:),yy(:),ff(:,:),xt(:,:),ft(:,:)
    integer::i,j,k
    allocate(mean(model%n),sd(model%n))
    do i=1,model%n
      allocate(xx(model%n-1,model%d),yy(model%n-1),ff(model%n-1,model%p),xt(1,model%d),ft(1,model%p))
      k=0
      do j=1,model%n
        if(j==i)cycle
        k=k+1; xx(k,:)=model%x(j,:); yy(k)=model%y(j); ff(k,:)=model%f(j,:)
      end do
      xt(1,:)=model%x(i,:); ft(1,:)=model%f(i,:); m2=model; m2%x=xx;m2%y=yy;m2%f=ff;m2%n=model%n-1
      call km_recompute(m2,reestimate_trend=tr); call km_predict(m2,xt,ft,kind,pr,se_compute=.true.)
      mean(i)=pr%mean(1); sd(i)=pr%sd(1)
      deallocate(xx,yy,ff,xt,ft)
    end do
  end subroutine loo_slow

  function km_loglik(model) result(ll)
    type(km_model), intent(in) :: model
    real(dp) :: ll
    real(dp), allocatable :: c(:,:),res(:,:),xsol(:,:)
    integer :: n
    n=model%n
    if(model%noise_flag) then
      call covariance_matrix(model%covariance,model%x,c,noise_var=model%noise_var,include_nugget=.false.)
    else
      call covariance_matrix(model%covariance,model%x,c,include_nugget=.true.)
    end if
    res=reshape(model%y,[model%n,1])-matmul(model%f,reshape(model%trend_coef,[model%p,1]))
    call solve_lower(model%l,res,xsol)
    ll=-0.5_dp*(real(n,dp)*log(2.0_dp*pi_dp)+logdet_from_chol(model%l)+dot_product(xsol(:,1),xsol(:,1)))
  end function km_loglik

  subroutine km_cov_vector_dx(model,x,c,grad)
    type(km_model), intent(in) :: model
    real(dp), intent(in) :: x(:),c(:)
    real(dp), allocatable, intent(out) :: grad(:,:)
    call covariance_vector_dx(model%covariance,x,model%x,c,grad)
  end subroutine km_cov_vector_dx

  subroutine profile_whitened(model,l,xw,mw,beta,res)
    type(km_model), intent(in) :: model
    real(dp), intent(in) :: l(:,:)
    real(dp), allocatable, intent(out) :: xw(:,:),mw(:,:),beta(:,:),res(:,:)
    integer :: info
    call solve_lower(l,reshape(model%y,[model%n,1]),xw); call solve_lower(l,model%f,mw)
    if(model%trend_known .and. allocated(model%trend_coef)) then
      beta=reshape(model%trend_coef,[model%p,1])
    else
      call least_squares_normal(mw,xw,beta,info)
      if(info/=0) then
        allocate(beta(model%p,1)); beta=0.0_dp
      end if
    end if
    res=xw-matmul(mw,beta)
  end subroutine profile_whitened


  function identity_matrix(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0.0_dp; do i=1,n;a(i,i)=1.0_dp;end do
  end function identity_matrix

  function diagonal_noise(model,cov) result(a)
    type(km_model),intent(in)::model
    type(covariance_model),intent(in)::cov
    real(dp)::a(model%n,model%n)
    integer::i
    a=0.0_dp
    if(model%noise_flag) then
      do i=1,model%n;a(i,i)=model%noise_var(i);end do
    else if(cov%nugget_flag) then
      do i=1,model%n;a(i,i)=cov%nugget;end do
    end if
  end function diagonal_noise

  real(dp) function residual_variance(f,y) result(v)
    real(dp),intent(in)::f(:,:),y(:)
    real(dp),allocatable::b(:,:),r(:,:)
    integer::info
    call least_squares_normal(f,reshape(y,[size(y),1]),b,info)
    if(info==0) then
      r=reshape(y,[size(y),1])-matmul(f,b); v=sum(r*r)/real(max(1,size(y)-1),dp)
    else
      v=sum((y-sum(y)/real(size(y),dp))**2)/real(max(1,size(y)-1),dp)
    end if
    v=max(v,1.0e-10_dp)
  end function residual_variance

  subroutine append_rows(a,b,c)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp),allocatable,intent(out)::c(:,:)
    allocate(c(size(a,1)+size(b,1),size(a,2)))
    c(1:size(a,1),:)=a; c(size(a,1)+1:,:)=b
  end subroutine append_rows

  subroutine trend_gradient_constant(x,g)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: g(:,:)
    allocate(g(1,size(x))); g=0.0_dp
  end subroutine trend_gradient_constant

  subroutine trend_gradient_linear(x,g)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: g(:,:)
    integer :: d,i
    d=size(x); allocate(g(d+1,d)); g=0.0_dp
    do i=1,d; g(i+1,i)=1.0_dp; end do
  end subroutine trend_gradient_linear

  subroutine trend_gradient_linear_interactions(x,g)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: g(:,:)
    integer :: d,p,i,j,k
    d=size(x); p=1+d+d*(d-1)/2; allocate(g(p,d)); g=0.0_dp
    do i=1,d; g(i+1,i)=1.0_dp; end do
    k=d+2
    do i=1,d-1
      do j=i+1,d
        g(k,i)=x(j); g(k,j)=x(i); k=k+1
      end do
    end do
  end subroutine trend_gradient_linear_interactions

  subroutine trend_gradient_quadratic(x,g)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: g(:,:)
    integer :: d,p,i,j,k
    d=size(x); p=1+d+d+d*(d-1)/2; allocate(g(p,d)); g=0.0_dp
    do i=1,d; g(i+1,i)=1.0_dp; end do
    k=d+2
    do i=1,d; g(k,i)=2.0_dp*x(i); k=k+1; end do
    do i=1,d-1
      do j=i+1,d
        g(k,i)=x(j); g(k,j)=x(i); k=k+1
      end do
    end do
  end subroutine trend_gradient_quadratic

  subroutine trend_constant(x,f)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::f(:,:)
    allocate(f(size(x,1),1));f=1.0_dp
  end subroutine trend_constant

  subroutine trend_linear(x,f)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::f(:,:)
    allocate(f(size(x,1),size(x,2)+1));f(:,1)=1.0_dp;f(:,2:)=x
  end subroutine trend_linear


  subroutine trend_linear_interactions(x,f)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::f(:,:)
    integer::d,p,i,j,k
    d=size(x,2);p=1+d+d*(d-1)/2;allocate(f(size(x,1),p));f(:,1)=1.0_dp;f(:,2:d+1)=x;k=d+2
    do i=1,d-1
      do j=i+1,d
        f(:,k)=x(:,i)*x(:,j);k=k+1
      end do
    end do
  end subroutine trend_linear_interactions

  subroutine trend_quadratic(x,f)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::f(:,:)
    integer::d,p,i,j,k
    d=size(x,2);p=1+d+d+d*(d-1)/2;allocate(f(size(x,1),p));f(:,1)=1.0_dp;f(:,2:d+1)=x;k=d+2
    do j=1,d;f(:,k)=x(:,j)*x(:,j);k=k+1;end do
    do i=1,d-1;do j=i+1,d;f(:,k)=x(:,i)*x(:,j);k=k+1;end do;end do
  end subroutine trend_quadratic

  pure real(dp) function t_quantile_975(df) result(q)
    integer, intent(in) :: df
    real(dp), parameter :: tab(30) = [ &
      12.7062047364_dp,4.30265272975_dp,3.18244630528_dp,2.77644510520_dp,2.57058183564_dp, &
      2.44691184879_dp,2.36462425101_dp,2.30600413503_dp,2.26215716285_dp,2.22813885196_dp, &
      2.20098516008_dp,2.17881282966_dp,2.16036865646_dp,2.14478668792_dp,2.13144954556_dp, &
      2.11990529922_dp,2.10981557783_dp,2.10092204024_dp,2.09302405441_dp,2.08596344727_dp, &
      2.07961384473_dp,2.07387306790_dp,2.06865761042_dp,2.06389856163_dp,2.05953855275_dp, &
      2.05552943864_dp,2.05183051648_dp,2.04840714180_dp,2.04522964213_dp,2.04227245630_dp ]
    real(dp), parameter :: z=1.959963984540054_dp
    real(dp) :: nu
    if(df<=0) error stop 't_quantile_975: nonpositive degrees of freedom'
    if(df<=30) then
      q=tab(df)
    else
      nu=real(df,dp)
      q=z+(z**3+z)/(4.0_dp*nu) &
        +(5.0_dp*z**5+16.0_dp*z**3+3.0_dp*z)/(96.0_dp*nu**2) &
        +(3.0_dp*z**7+19.0_dp*z**5+17.0_dp*z**3-15.0_dp*z)/(384.0_dp*nu**3)
    end if
  end function t_quantile_975

  elemental real(dp) function scad(theta,lambda) result(v)
    real(dp),intent(in)::theta,lambda
    real(dp),parameter::a=3.7_dp
    real(dp)::t
    t=abs(theta)
    if(t<=lambda) then;v=lambda*t
    else if(t<=a*lambda) then;v=(-t*t+2.0_dp*a*lambda*t-lambda*lambda)/(2.0_dp*(a-1.0_dp))
    else;v=0.5_dp*(a+1.0_dp)*lambda*lambda
    end if
  end function scad

  elemental real(dp) function scad_derivative(theta,lambda) result(v)
    real(dp),intent(in)::theta,lambda
    real(dp),parameter::a=3.7_dp
    real(dp)::t,s
    t=abs(theta);s=merge(1.0_dp,-1.0_dp,theta>=0.0_dp)
    if(t<=lambda) then;v=lambda*s
    else if(t<=a*lambda) then;v=(a*lambda-t)/(a-1.0_dp)*s
    else;v=0.0_dp
    end if
  end function scad_derivative

end module dk_model
