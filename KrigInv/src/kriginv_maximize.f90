module kriginv_maximize
  use kriginv_kinds, only : dp
  use kriginv_sobol, only : sobol_points
  use kriginv_model, only : krig_model, krig_prediction, predict_nobias_km
  use kriginv_math, only : normal_pdf, normal_cdf
  use kriginv_criteria, only : ranjan_optim, bichon_optim, tmse_optim, tsee_optim, excursion_probability, &
    vorob_threshold, sur_optim_parallel, jn_optim_parallel, timse_optim_parallel, vorob_optim_parallel, &
    vorobvol_optim_parallel, compute_real_volume_constant
  use kriginv_integration, only : integration_result
  use kriginv_optimize, only : optimizer_control, bounded_de, discrete_optimum
  use anmc_conservative, only : conservative_estimate
  use anmc_types, only : conservative_result
  implicit none
  private
  abstract interface
    function local_objective(x) result(v)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: v
    end function local_objective
  end interface
  type, public :: optimization_result
    real(dp), allocatable :: par(:,:), allvalues(:)
    real(dp) :: value=0.0_dp
    real(dp) :: alpha=0.5_dp
    logical :: ok=.true.
  end type optimization_result
  public :: max_infill_criterion, max_sur_parallel, max_timse_parallel, max_vorob_parallel, max_futurevol_parallel
contains
  function max_infill_criterion(lower,upper,method,thresholds,model,method_param,control) result(res)
    real(dp), intent(in) :: lower(:),upper(:),thresholds(:)
    character(len=*), intent(in) :: method
    type(krig_model), intent(in) :: model
    real(dp), intent(in), optional :: method_param
    type(optimizer_control), intent(in), optional :: control
    type(optimization_result) :: res
    type(optimizer_control) :: ctl
    real(dp), allocatable :: best(:),pts(:,:),vals(:)
    real(dp) :: bv
    integer :: d
    d=model%d; if(size(lower)/=d .or. size(upper)/=d) then; res%ok=.false.; return; end if
    if(present(control)) ctl=control
    if(trim(ctl%method)=='discrete') then
      if(allocated(ctl%optim_points)) then; pts=ctl%optim_points
      else; pts=scale_points(sobol_points(100*d,d),lower,upper); end if
      call discrete_optimum(obj,pts,.true.,best,bv,vals); res%allvalues=vals
    else
      call bounded_de(obj,lower,upper,.true.,ctl,best,bv)
    end if
    allocate(res%par(1,d)); res%par(1,:)=best; res%value=bv
  contains
    real(dp) function obj(v) result(f)
      real(dp), intent(in) :: v(:)
      real(dp) :: xx(1,size(v)),a
      real(dp), allocatable :: z(:)
      xx(1,:)=v
      if(trim(method)=='tmse') then; a=0.0_dp; else; a=1.0_dp; end if
      if(present(method_param)) a=method_param
      select case(trim(method))
      case('tmse'); z=tmse_optim(xx,model,thresholds,a)
      case('bichon'); z=bichon_optim(xx,model,thresholds(1),a)
      case('tsee'); z=tsee_optim(xx,model,thresholds(1))
      case default; z=ranjan_optim(xx,model,thresholds(1),a)
      end select
      f=z(1)
    end function obj
  end function max_infill_criterion

  function max_sur_parallel(lower,upper,batchsize,integration,thresholds,model,newnoise,real_volume_variance,control) result(res)
    real(dp), intent(in) :: lower(:),upper(:),thresholds(:)
    integer, intent(in) :: batchsize
    type(integration_result), intent(in) :: integration
    type(krig_model), intent(in) :: model
    real(dp), intent(in), optional :: newnoise
    logical, intent(in), optional :: real_volume_variance
    type(optimizer_control), intent(in), optional :: control
    type(optimization_result) :: res
    type(optimizer_control) :: ctl
    type(krig_prediction) :: pr
    real(dp), allocatable :: pn(:),best(:),lo(:),hi(:),vals(:)
    real(dp) :: bv,current,nn,const
    integer :: d,j
    logical :: rv
    d=model%d; rv=.false.; if(present(real_volume_variance)) rv=real_volume_variance
    nn=0.0_dp; if(present(newnoise)) nn=newnoise
    if(present(control)) ctl=control
    pr=predict_nobias_km(model,integration%points,'UK',.false.); pn=excursion_probability(pr%mean,pr%sd,thresholds)
    if(integration%has_weights) then; current=sum(integration%weights*pn*(1.0_dp-pn))
    else; current=sum(pn*(1.0_dp-pn))/real(size(pn),dp); end if
    if(rv) current=0.0_dp
    allocate(lo(d*batchsize),hi(d*batchsize))
    do j=1,batchsize; lo((j-1)*d+1:j*d)=lower; hi((j-1)*d+1:j*d)=upper; end do
    if(trim(ctl%method)=='discrete' .and. allocated(ctl%optim_points) .and. size(ctl%optim_points,2)==d) then
      call greedy_discrete(obj,ctl%optim_points,batchsize,.false.,res%par,bv,vals); res%allvalues=vals
    else if(trim(ctl%method)=='discrete' .and. allocated(ctl%optim_points) .and. size(ctl%optim_points,2)==d*batchsize) then
      call discrete_optimum(obj,ctl%optim_points,.false.,best,bv,vals); res%allvalues=vals; res%par=flat_to_points(best,d,batchsize)
    else
      call bounded_de(obj,lo,hi,.false.,ctl,best,bv); res%par=flat_to_points(best,d,batchsize)
    end if
    const=0.0_dp
    if(rv .and. size(thresholds)==1) then
      if(integration%has_weights .and. size(integration%weights)==size(integration%points,1)) then
        const=compute_real_volume_constant(model,integration%points,thresholds(1),integration%weights)
      else if(.not.integration%has_weights) then
        const=compute_real_volume_constant(model,integration%points,thresholds(1))
      end if
    end if
    res%value=bv+const
  contains
    real(dp) function obj(v) result(f)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: xx(:,:),nv(:)
      xx=flat_to_points(v,d,size(v)/d); allocate(nv(size(v)/d)); nv=nn
      if(rv) then
        if(integration%has_weights) then
          f=jn_optim_parallel(xx,integration%points,integration%weights,pr%mean,pr%sd,model,thresholds,nv,current)
        else
          f=jn_optim_parallel(xx,integration%points,oldmean=pr%mean,oldsd=pr%sd,model=model,thresholds=thresholds, &
                              newnoise=nv,current_sur=current)
        end if
      else
        if(integration%has_weights) then
          f=sur_optim_parallel(xx,integration%points,integration%weights,pr%mean,pr%sd,model,thresholds,nv,current)
        else
          f=sur_optim_parallel(xx,integration%points,oldmean=pr%mean,oldsd=pr%sd,model=model,thresholds=thresholds, &
                               newnoise=nv,current_sur=current)
        end if
      end if
    end function obj
  end function max_sur_parallel

  function max_timse_parallel(lower,upper,batchsize,integration,thresholds,model,newnoise,epsilon,imse,control) result(res)
    real(dp), intent(in) :: lower(:),upper(:),thresholds(:)
    integer, intent(in) :: batchsize
    type(integration_result), intent(in) :: integration
    type(krig_model), intent(in) :: model
    real(dp), intent(in), optional :: newnoise,epsilon
    logical, intent(in), optional :: imse
    type(optimizer_control), intent(in), optional :: control
    type(optimization_result) :: res
    type(optimizer_control) :: ctl
    type(krig_prediction) :: pr
    real(dp), allocatable :: weight(:),best(:),lo(:),hi(:),vals(:)
    real(dp) :: bv,current,nn,eps,v
    integer :: d,j,k
    logical :: do_imse
    d=model%d; nn=0.0_dp; if(present(newnoise)) nn=newnoise; eps=0.0_dp; if(present(epsilon)) eps=epsilon
    do_imse=.false.; if(present(imse)) do_imse=imse; if(present(control)) ctl=control
    pr=predict_nobias_km(model,integration%points,'UK',.false.); allocate(weight(size(pr%mean))); weight=1.0_dp
    if(.not.do_imse) then
      weight=0.0_dp
      do k=1,size(thresholds)
        do j=1,size(weight)
          v=pr%sd(j)**2+eps**2
          if(v>0.0_dp) weight(j)=weight(j)+normal_pdf((pr%mean(j)-thresholds(k))/sqrt(v))/sqrt(v)
        end do
      end do
    end if
    if(integration%has_weights) then; current=sum(integration%weights*weight*pr%sd**2)
    else; current=sum(weight*pr%sd**2)/real(size(weight),dp); end if
    allocate(lo(d*batchsize),hi(d*batchsize))
    do j=1,batchsize; lo((j-1)*d+1:j*d)=lower; hi((j-1)*d+1:j*d)=upper; end do
    if(trim(ctl%method)=='discrete' .and. allocated(ctl%optim_points) .and. size(ctl%optim_points,2)==d) then
      call greedy_discrete(obj,ctl%optim_points,batchsize,.false.,res%par,bv,vals); res%allvalues=vals
    else if(trim(ctl%method)=='discrete' .and. allocated(ctl%optim_points) .and. size(ctl%optim_points,2)==d*batchsize) then
      call discrete_optimum(obj,ctl%optim_points,.false.,best,bv,vals); res%allvalues=vals; res%par=flat_to_points(best,d,batchsize)
    else
      call bounded_de(obj,lo,hi,.false.,ctl,best,bv); res%par=flat_to_points(best,d,batchsize)
    end if
    res%value=bv
  contains
    real(dp) function obj(vv) result(f)
      real(dp), intent(in) :: vv(:)
      real(dp), allocatable :: xx(:,:),nv(:)
      xx=flat_to_points(vv,d,size(vv)/d); allocate(nv(size(vv)/d)); nv=nn
      if(integration%has_weights) then
        f=timse_optim_parallel(xx,integration%points,integration%weights,pr%mean,pr%sd,model,weight,nv,current)
      else
        f=timse_optim_parallel(xx,integration%points,oldmean=pr%mean,oldsd=pr%sd,model=model,weight=weight, &
                               newnoise=nv,current_timse=current)
      end if
    end function obj
  end function max_timse_parallel

  function max_vorob_parallel(lower,upper,batchsize,integration,threshold,model,newnoise,penalisation,type_ex,control) result(res)
    real(dp), intent(in) :: lower(:),upper(:),threshold
    integer, intent(in) :: batchsize
    type(integration_result), intent(in) :: integration
    type(krig_model), intent(in) :: model
    real(dp), intent(in), optional :: newnoise,penalisation
    character(len=*), intent(in), optional :: type_ex
    type(optimizer_control), intent(in), optional :: control
    type(optimization_result) :: res
    type(optimizer_control) :: ctl
    type(krig_prediction) :: pr
    real(dp), allocatable :: pn(:),loss(:),best(:),lo(:),hi(:),vals(:)
    real(dp) :: bv,current,nn,pen,alpha
    integer :: d,j
    logical :: above
    d=model%d; nn=0.0_dp; if(present(newnoise)) nn=newnoise; pen=1.0_dp; if(present(penalisation)) pen=penalisation
    above=.true.; if(present(type_ex)) above=(type_ex(1:1)=='>'); if(present(control)) ctl=control
    pr=predict_nobias_km(model,integration%points,'UK',.false.); allocate(pn(size(pr%mean)))
    do j=1,size(pn)
      if(pr%sd(j)>0.0_dp) then
        if(above) pn(j)=normal_cdf((pr%mean(j)-threshold)/pr%sd(j))
        if(.not.above) pn(j)=normal_cdf((threshold-pr%mean(j))/pr%sd(j))
      else
        if(above) pn(j)=merge(1.0_dp,0.0_dp,pr%mean(j)>threshold)
        if(.not.above) pn(j)=merge(1.0_dp,0.0_dp,pr%mean(j)<threshold)
      end if
    end do
    if(integration%has_alpha) then; alpha=integration%alpha; else; alpha=vorob_threshold(pn); end if
    allocate(loss(size(pn)))
    do j=1,size(pn)
      if(pn(j)>alpha) then; loss(j)=pen*(1.0_dp-pn(j)); else; loss(j)=pn(j); end if
    end do
    if(integration%has_weights) then; current=sum(integration%weights*loss)
    else; current=sum(loss)/real(size(loss),dp); end if
    allocate(lo(d*batchsize),hi(d*batchsize))
    do j=1,batchsize; lo((j-1)*d+1:j*d)=lower; hi((j-1)*d+1:j*d)=upper; end do
    if(trim(ctl%method)=='discrete' .and. allocated(ctl%optim_points) .and. size(ctl%optim_points,2)==d) then
      call greedy_discrete(obj,ctl%optim_points,batchsize,.false.,res%par,bv,vals); res%allvalues=vals
    else if(trim(ctl%method)=='discrete' .and. allocated(ctl%optim_points) .and. size(ctl%optim_points,2)==d*batchsize) then
      call discrete_optimum(obj,ctl%optim_points,.false.,best,bv,vals); res%allvalues=vals; res%par=flat_to_points(best,d,batchsize)
    else
      call bounded_de(obj,lo,hi,.false.,ctl,best,bv); res%par=flat_to_points(best,d,batchsize)
    end if
    res%value=bv; res%alpha=alpha
  contains
    real(dp) function obj(v) result(f)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: xx(:,:),nv(:)
      xx=flat_to_points(v,d,size(v)/d); allocate(nv(size(v)/d)); nv=nn
      if(integration%has_weights) then
        f=vorob_optim_parallel(xx,integration%points,integration%weights,pr%mean,pr%sd,model,threshold,alpha,pen, &
                               merge('>','<',above),nv,current)
      else
        f=vorob_optim_parallel(xx,integration%points,oldmean=pr%mean,oldsd=pr%sd,model=model,threshold=threshold, &
                               alpha=alpha,penalisation=pen,type_ex=merge('>','<',above),newnoise=nv,current_vorob=current)
      end if
    end function obj
  end function max_vorob_parallel

  function max_futurevol_parallel(lower,upper,batchsize,integration,threshold,model,newnoise,type_ex,control) result(res)
    real(dp), intent(in) :: lower(:),upper(:),threshold
    integer, intent(in) :: batchsize
    type(integration_result), intent(in) :: integration
    type(krig_model), intent(in) :: model
    real(dp), intent(in), optional :: newnoise
    character(len=*), intent(in), optional :: type_ex
    type(optimizer_control), intent(in), optional :: control
    type(optimization_result) :: res
    type(optimizer_control) :: ctl
    type(krig_prediction) :: pr
    real(dp), allocatable :: pn(:),loss(:),best(:),lo(:),hi(:),vals(:)
    real(dp) :: bv,current,nn,alpha
    integer :: d,j
    logical :: above
    d=model%d; nn=0.0_dp; if(present(newnoise)) nn=newnoise
    above=.true.; if(present(type_ex)) above=(type_ex(1:1)=='>'); if(present(control)) ctl=control
    pr=predict_nobias_km(model,integration%points,'UK',.false.); allocate(pn(size(pr%mean)))
    do j=1,size(pn)
      if(pr%sd(j)>0.0_dp) then
        if(above) pn(j)=normal_cdf((pr%mean(j)-threshold)/pr%sd(j))
        if(.not.above) pn(j)=normal_cdf((threshold-pr%mean(j))/pr%sd(j))
      else
        pn(j)=0.0_dp
      end if
    end do
    if(integration%has_alpha) then
      alpha=integration%alpha
    else
      call conservative_alpha(model,lower,upper,threshold,above,alpha)
    end if
    allocate(loss(size(pn)))
    do j=1,size(pn); loss(j)=merge(1.0_dp-pn(j),pn(j),pn(j)>alpha); end do
    if(integration%has_weights) then; current=sum(integration%weights*loss); else; current=0.0_dp; end if
    allocate(lo(d*batchsize),hi(d*batchsize))
    do j=1,batchsize; lo((j-1)*d+1:j*d)=lower; hi((j-1)*d+1:j*d)=upper; end do
    if(trim(ctl%method)=='discrete' .and. allocated(ctl%optim_points) .and. size(ctl%optim_points,2)==d) then
      call greedy_discrete(obj,ctl%optim_points,batchsize,.true.,res%par,bv,vals); res%allvalues=vals
    else if(trim(ctl%method)=='discrete' .and. allocated(ctl%optim_points) .and. size(ctl%optim_points,2)==d*batchsize) then
      call discrete_optimum(obj,ctl%optim_points,.true.,best,bv,vals); res%allvalues=vals; res%par=flat_to_points(best,d,batchsize)
    else
      call bounded_de(obj,lo,hi,.true.,ctl,best,bv); res%par=flat_to_points(best,d,batchsize)
    end if
    res%value=bv; res%alpha=alpha
  contains
    real(dp) function obj(v) result(f)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: xx(:,:),nv(:)
      xx=flat_to_points(v,d,size(v)/d); allocate(nv(size(v)/d)); nv=nn
      if(integration%has_weights) then
        f=vorobvol_optim_parallel(xx,integration%points,integration%weights,pr%mean,pr%sd,model,threshold,alpha, &
                                  merge('>','<',above),nv,current)
      else
        f=vorobvol_optim_parallel(xx,integration%points,oldmean=pr%mean,oldsd=pr%sd,model=model,threshold=threshold, &
                                  alpha=alpha,type_ex=merge('>','<',above),newnoise=nv,current_crit=current)
      end if
    end function obj
  end function max_futurevol_parallel

  subroutine conservative_alpha(model,lower,upper,threshold,above,alpha)
    type(krig_model), intent(in) :: model
    real(dp), intent(in) :: lower(:),upper(:),threshold
    logical, intent(in) :: above
    real(dp), intent(out) :: alpha
    type(krig_prediction) :: prc
    type(conservative_result) :: ce
    real(dp), allocatable :: u(:,:),design(:,:)
    integer :: i,j,n
    n=max(20,500*model%d); u=sobol_points(n,model%d); design=u
    do j=1,model%d; design(:,j)=lower(j)+u(:,j)*(upper(j)-lower(j)); end do
    prc=predict_nobias_km(model,design,'UK',.true.)
    do i=1,n; prc%covariance(i,i)=prc%covariance(i,i)+1.0e-7_dp; end do
    ce=conservative_estimate(alpha=0.95_dp,mean=prc%mean,covariance=prc%covariance,design=design, &
                             threshold=threshold,excursion_type=merge('>','<',above),algo='GANMC',verb=0)
    if(ce%ok) then; alpha=ce%level; else; alpha=0.5_dp; end if
  end subroutine conservative_alpha

  subroutine greedy_discrete(fun,candidates,batch,maximize,selected,best_v,allvalues)
    procedure(local_objective) :: fun
    real(dp), intent(in) :: candidates(:,:)
    integer, intent(in) :: batch
    logical, intent(in) :: maximize
    real(dp), allocatable, intent(out) :: selected(:,:),allvalues(:)
    real(dp), intent(out) :: best_v
    real(dp), allocatable :: trial(:,:),flat(:)
    integer :: i,j,ibest,d,n
    d=size(candidates,2); n=size(candidates,1)
    allocate(selected(batch,d),allvalues(n))
    do j=1,batch
      allocate(trial(j,d))
      if(j>1) trial(1:j-1,:)=selected(1:j-1,:)
      do i=1,n
        trial(j,:)=candidates(i,:); flat=points_to_flat(trial); allvalues(i)=fun(flat)
      end do
      ibest=1
      do i=2,n
        if((maximize .and. allvalues(i)>allvalues(ibest)) .or. &
           (.not.maximize .and. allvalues(i)<allvalues(ibest))) ibest=i
      end do
      selected(j,:)=candidates(ibest,:); best_v=allvalues(ibest)
      deallocate(trial)
    end do
  end subroutine greedy_discrete

  function points_to_flat(x) result(v)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: v(:)
    integer :: i,d
    d=size(x,2); allocate(v(size(x,1)*d))
    do i=1,size(x,1); v((i-1)*d+1:i*d)=x(i,:); end do
  end function points_to_flat

  function flat_to_points(v,d,batch) result(x)
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: d,batch
    real(dp), allocatable :: x(:,:)
    integer :: j
    allocate(x(batch,d))
    do j=1,batch; x(j,:)=v((j-1)*d+1:j*d); end do
  end function flat_to_points

  function scale_points(u,lower,upper) result(x)
    real(dp), intent(in) :: u(:,:),lower(:),upper(:)
    real(dp), allocatable :: x(:,:)
    integer :: j
    x=u
    do j=1,size(u,2); x(:,j)=lower(j)+u(:,j)*(upper(j)-lower(j)); end do
  end function scale_points
end module kriginv_maximize
