module kriginv_integration
  use kriginv_kinds, only : dp
  use kriginv_sobol, only : sobol_points
  use kriginv_model, only : krig_model, krig_prediction, predict_nobias_km
  use kriginv_criteria, only : excursion_probability, vorob_threshold
  use kriginv_math, only : normal_pdf
  implicit none
  private
  type, public :: integration_control
    integer :: n_points=0
    integer :: n_candidates=0
    character(len=16) :: distrib='sobol'
    character(len=16) :: init_distrib='sobol'
    real(dp) :: min_prob=0.001_dp
    integer :: seed=12345
    real(dp), allocatable :: specified_points(:,:)
  end type integration_control
  type, public :: integration_result
    real(dp), allocatable :: points(:,:), weights(:)
    logical :: has_weights=.false.
    real(dp) :: alpha=0.5_dp
    logical :: has_alpha=.false.
    logical :: ok=.true.
  end type integration_result
  public :: integration_design
contains
  function integration_design(lower,upper,model,thresholds,control) result(res)
    real(dp), intent(in) :: lower(:),upper(:)
    type(krig_model), intent(in), optional :: model
    real(dp), intent(in), optional :: thresholds(:)
    type(integration_control), intent(in), optional :: control
    type(integration_result) :: res
    type(integration_control) :: ctl
    type(krig_prediction) :: pr
    real(dp), allocatable :: cand(:,:),pn(:),tau(:),prob(:),wn(:),wfun(:)
    integer, allocatable :: idx(:),idx2(:)
    integer :: d,i
    real(dp) :: s
    if(size(lower)/=size(upper) .or. any(upper<=lower)) then; res%ok=.false.; return; end if
    d=size(lower); ctl%n_points=100*d; ctl%n_candidates=1000*d
    if(present(control)) ctl=control
    if(ctl%n_points<=0) ctl%n_points=100*d
    if(ctl%n_candidates<=0) ctl%n_candidates=10*ctl%n_points
    call set_seed(ctl%seed)
    if(allocated(ctl%specified_points)) then
      res%points=ctl%specified_points; res%has_weights=.false.; return
    end if
    select case(trim(ctl%distrib))
    case('sobol')
      res%points=scale_unit(sobol_points(ctl%n_points,d),lower,upper)
    case('MC','mc')
      allocate(res%points(ctl%n_points,d)); call random_number(res%points); res%points=scale_unit(res%points,lower,upper)
    case('sur','vorob','jn','timse','imse')
      if(.not.present(model) .or. .not.present(thresholds)) then; res%ok=.false.; return; end if
      if(trim(ctl%init_distrib)=='MC' .or. trim(ctl%init_distrib)=='mc') then
        allocate(cand(ctl%n_candidates,d)); call random_number(cand)
      else
        cand=sobol_points(ctl%n_candidates,d)
      end if
      cand=scale_unit(cand,lower,upper); pr=predict_nobias_km(model,cand,'UK',.false.)
      pn=excursion_probability(pr%mean,pr%sd,thresholds)
      allocate(tau(size(pn))); tau=0.0_dp
      select case(trim(ctl%distrib))
      case('sur')
        tau=pn*(1.0_dp-pn)
      case('vorob')
        res%alpha=vorob_threshold(pn); res%has_alpha=.true.
        do i=1,size(pn)
          if(pn(i)>res%alpha) then; tau(i)=1.0_dp-pn(i); else; tau(i)=pn(i); end if
        end do
      case('jn')
        tau=pn
      case('timse','imse')
        if(trim(ctl%distrib)=='imse') then
          tau=pr%sd*pr%sd
        else
          allocate(wfun(size(pn))); wfun=0.0_dp
          do i=1,size(thresholds)
            where(pr%sd>0.0_dp)
              wfun=wfun+normal_pdf((pr%mean-thresholds(i))/max(pr%sd,tiny(1.0_dp)))/max(pr%sd,tiny(1.0_dp))
            end where
          end do
          tau=wfun*pr%sd*pr%sd
        end if
      end select
      s=sum(tau); if(s<=0.0_dp) then; tau=1.0_dp; s=real(size(tau),dp); end if
      allocate(prob(size(tau))); prob=max(tau/s,ctl%min_prob/real(size(tau),dp)); prob=prob/sum(prob)
      allocate(wn(size(prob))); wn=1.0_dp/(prob*real(size(prob)*ctl%n_points,dp))
      if(trim(ctl%distrib)=='jn') then
        idx=sample_indices(prob,ctl%n_points); idx2=sample_indices(prob,ctl%n_points)
        allocate(res%points(2*ctl%n_points,d),res%weights(ctl%n_points))
        res%points(1:ctl%n_points,:)=cand(idx,:); res%points(ctl%n_points+1:,:)=cand(idx2,:)
        res%weights=wn(idx)*wn(idx2)*real(ctl%n_points,dp); res%has_weights=.true.
      else
        idx=sample_indices(prob,ctl%n_points); res%points=cand(idx,:); res%weights=wn(idx); res%has_weights=.true.
      end if
    case default
      res%ok=.false.
    end select
  end function integration_design

  function scale_unit(u,lower,upper) result(x)
    real(dp), intent(in) :: u(:,:),lower(:),upper(:)
    real(dp), allocatable :: x(:,:)
    integer :: j
    x=u
    do j=1,size(u,2); x(:,j)=lower(j)+u(:,j)*(upper(j)-lower(j)); end do
  end function scale_unit

  function sample_indices(prob,n) result(idx)
    real(dp), intent(in) :: prob(:)
    integer, intent(in) :: n
    integer, allocatable :: idx(:)
    real(dp), allocatable :: cdf(:),u(:)
    integer :: i,j
    allocate(cdf(size(prob)),u(n),idx(n)); cdf(1)=prob(1)
    do i=2,size(prob); cdf(i)=cdf(i-1)+prob(i); end do
    cdf(size(cdf))=1.0_dp; call random_number(u)
    do i=1,n
      idx(i)=size(prob)
      do j=1,size(prob)
        if(u(i)<=cdf(j)) then; idx(i)=j; exit; end if
      end do
    end do
  end function sample_indices

  subroutine set_seed(seed)
    integer, intent(in) :: seed
    integer, allocatable :: put(:)
    integer :: n,i
    call random_seed(size=n); allocate(put(n))
    do i=1,n; put(i)=modulo(seed+104729*i,2147483646)+1; end do
    call random_seed(put=put)
  end subroutine set_seed
end module kriginv_integration
