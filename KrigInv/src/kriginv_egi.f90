module kriginv_egi
  use kriginv_kinds, only : dp
  use kriginv_model, only : krig_model, krig_prediction, predict_nobias_km, update_krig_model
  use kriginv_integration, only : integration_control, integration_result, integration_design
  use kriginv_optimize, only : optimizer_control
  use kriginv_maximize, only : optimization_result, max_infill_criterion, max_sur_parallel, max_timse_parallel, &
                               max_vorob_parallel, max_futurevol_parallel
  use kriginv_sobol, only : sobol_points
  use anmc_conservative, only : conservative_estimate
  use anmc_types, only : conservative_result
  implicit none
  private
  abstract interface
    function kriginv_blackbox(x) result(y)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: y
    end function kriginv_blackbox
  end interface
  type, public :: egi_result
    real(dp), allocatable :: par(:,:), value(:)
    integer :: npoints=0, nsteps=0
    type(krig_model) :: lastmodel
    real(dp) :: lastcriterion=0.0_dp
    logical :: ok=.true.
  end type egi_result
  public :: egi, egi_parallel
contains
  function egi(thresholds,model,method,fun,iter,lower,upper,new_noise,method_param,opt_control,int_control, &
               cov_reestimate,trend_reestimate,nugget_reestimate) result(res)
    real(dp), intent(in) :: thresholds(:),lower(:),upper(:)
    type(krig_model), intent(in) :: model
    character(len=*), intent(in) :: method
    procedure(kriginv_blackbox) :: fun
    integer, intent(in) :: iter
    real(dp), intent(in), optional :: new_noise,method_param
    type(optimizer_control), intent(in), optional :: opt_control
    type(integration_control), intent(in), optional :: int_control
    logical, intent(in), optional :: cov_reestimate,trend_reestimate,nugget_reestimate
    type(egi_result) :: res
    res=egi_parallel(thresholds,model,method,fun,iter,1,lower,upper,new_noise,method_param,opt_control,int_control, &
                     cov_reestimate,trend_reestimate,nugget_reestimate)
  end function egi

  function egi_parallel(thresholds,model,method,fun,iter,batchsize,lower,upper,new_noise,method_param, &
                        opt_control,int_control,cov_reestimate,trend_reestimate,nugget_reestimate) result(res)
    real(dp), intent(in) :: thresholds(:),lower(:),upper(:)
    type(krig_model), intent(in) :: model
    character(len=*), intent(in) :: method
    procedure(kriginv_blackbox) :: fun
    integer, intent(in) :: iter,batchsize
    real(dp), intent(in), optional :: new_noise,method_param
    type(optimizer_control), intent(in), optional :: opt_control
    type(integration_control), intent(in), optional :: int_control
    logical, intent(in), optional :: cov_reestimate,trend_reestimate,nugget_reestimate
    type(egi_result) :: res
    type(krig_model) :: work
    type(optimization_result) :: opt
    type(integration_result) :: integ
    type(integration_control) :: ictl
    type(optimizer_control) :: octl
    real(dp), allocatable :: newy(:),noise(:),allx(:,:),ally(:)
    real(dp) :: nn,mp
    integer :: step,j,b,d,nsel,offset
    logical :: ok,cr,tr,nr
    work=model; d=model%d; b=batchsize; nn=0.0_dp; if(present(new_noise)) nn=new_noise
    cr=work%cov_reestimate_default; if(present(cov_reestimate)) cr=cov_reestimate
    tr=.true.; if(present(trend_reestimate)) tr=trend_reestimate
    nr=.false.; if(present(nugget_reestimate)) nr=nugget_reestimate
    if(trim(method)=='tmse' .or. trim(method)=='timse' .or. trim(method)=='imse') then
      mp=0.0_dp
    else
      mp=1.0_dp
    end if
    if(present(method_param)) mp=method_param
    if(present(opt_control)) octl=opt_control
    if(present(int_control)) ictl=int_control
    if(trim(method)=='tmse' .or. trim(method)=='ranjan' .or. trim(method)=='bichon' .or. trim(method)=='tsee') b=1
    if(iter<1 .or. b<1) then; res%ok=.false.; return; end if
    allocate(allx(iter*b,d),ally(iter*b)); offset=0
    do step=1,iter
      select case(trim(method))
      case('tmse','ranjan','bichon','tsee')
        opt=max_infill_criterion(lower,upper,method,thresholds,work,mp,octl)
      case('sur')
        integ=integration_design(lower,upper,work,thresholds,ictl)
        opt=max_sur_parallel(lower,upper,b,integ,thresholds,work,nn,.false.,octl)
      case('jn')
        ictl%distrib='jn'; integ=integration_design(lower,upper,work,thresholds,ictl)
        opt=max_sur_parallel(lower,upper,b,integ,thresholds,work,nn,.true.,octl)
      case('timse')
        ictl%distrib='timse'; integ=integration_design(lower,upper,work,thresholds,ictl)
        opt=max_timse_parallel(lower,upper,b,integ,thresholds,work,nn,mp,.false.,octl)
      case('imse')
        ictl%distrib='imse'; integ=integration_design(lower,upper,work,thresholds,ictl)
        opt=max_timse_parallel(lower,upper,b,integ,thresholds,work,nn,0.0_dp,.true.,octl)
      case('vorob')
        ictl%distrib='vorob'; integ=integration_design(lower,upper,work,thresholds,ictl)
        opt=max_vorob_parallel(lower,upper,b,integ,thresholds(1),work,nn,mp,'>',octl)
      case('vorobCons','vorobVol')
        ictl%distrib='vorob'; integ=integration_design(lower,upper,work,thresholds,ictl)
        call conservative_level(work,lower,upper,thresholds(1),integ%alpha,ok)
        if(ok) integ%has_alpha=.true.
        if(trim(method)=='vorobCons') then
          opt=max_vorob_parallel(lower,upper,b,integ,thresholds(1),work,nn,mp,'>',octl)
        else
          opt=max_futurevol_parallel(lower,upper,b,integ,thresholds(1),work,nn,'>',octl)
        end if
      case default
        opt=max_infill_criterion(lower,upper,'ranjan',thresholds,work,1.0_dp,octl)
      end select
      if(.not.opt%ok) then; res%ok=.false.; return; end if
      nsel=size(opt%par,1); allocate(newy(nsel),noise(nsel)); noise=nn
      do j=1,nsel; newy(j)=fun(opt%par(j,:)); end do
      allx(offset+1:offset+nsel,:)=opt%par; ally(offset+1:offset+nsel)=newy; offset=offset+nsel
      call update_krig_model(work,opt%par,newy,noise,ok,cov_reestimate=cr, &
                             trend_reestimate=tr,nugget_reestimate=nr)
      if(.not.ok) then; res%ok=.false.; return; end if
      res%lastcriterion=opt%value
      deallocate(newy,noise)
    end do
    res%par=allx(1:offset,:); res%value=ally(1:offset); res%npoints=b; res%nsteps=offset; res%lastmodel=work
  end function egi_parallel

  subroutine conservative_level(model,lower,upper,threshold,level,ok)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: lower(:),upper(:),threshold
    real(dp), intent(out) :: level
    logical, intent(out) :: ok
    type(krig_prediction) :: pr
    type(conservative_result) :: ce
    real(dp), allocatable :: u(:,:),design(:,:)
    integer :: j,n
    n=max(20,500*model%d); u=sobol_points(n,model%d); design=u
    do j=1,model%d; design(:,j)=lower(j)+u(:,j)*(upper(j)-lower(j)); end do
    pr=predict_nobias_km(model,design,'UK',.true.)
    do j=1,n; pr%covariance(j,j)=pr%covariance(j,j)+1.0e-7_dp; end do
    ce=conservative_estimate(alpha=0.95_dp,mean=pr%mean,covariance=pr%covariance,design=design, &
                             threshold=threshold,excursion_type='>',algo='GANMC',verb=0)
    ok=ce%ok; level=ce%level
  end subroutine conservative_level
end module kriginv_egi
