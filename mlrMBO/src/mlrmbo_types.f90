module mlrmbo_types
  use mlrmbo_kinds, only : dp, i8
  implicit none
  private
  integer, parameter, public :: mbo_real=1, mbo_integer=2, mbo_categorical=3
  integer, parameter, public :: crit_mean=1, crit_se=2, crit_ei=3, crit_cb=4, crit_aei=5, crit_eqi=6, crit_adacb=7
  integer, parameter, public :: mo_dib=1, mo_parego=2, mo_mspot=3
  integer, parameter, public :: batch_none=0, batch_cb=1, batch_cl=2, batch_moimbo=3
  integer, parameter, public :: lie_min=1, lie_max=2, lie_mean=3
  integer, parameter, public :: ref_all=1, ref_front=2, ref_const=3
  integer, parameter, public :: dib_sms=1, dib_eps=2

  type, public :: mbo_space
    integer :: d=0
    integer, allocatable :: kind(:)
    real(dp), allocatable :: lower(:),upper(:)
    integer, allocatable :: nlevels(:)
    integer, allocatable :: condition_parent(:),condition_level(:)
  contains
    procedure :: repair => repair_point
    procedure :: active => point_active
  end type mbo_space

  type, public :: mbo_control
    integer :: n_objectives=1
    logical, allocatable :: minimize(:)
    integer :: n_init=0
    integer :: max_iter=20
    integer :: max_evals=0
    integer :: propose_points=1
    integer :: interleave_random_points=0
    integer :: infill_criterion=crit_cb
    real(dp) :: se_threshold=1.0e-6_dp
    real(dp) :: cb_lambda=1.0_dp
    real(dp) :: cb_lambda_start=2.0_dp
    real(dp) :: cb_lambda_end=0.1_dp
    real(dp) :: eqi_beta=0.75_dp
    logical :: aei_use_nugget=.false.
    integer :: focus_restarts=3
    integer :: focus_iterations=5
    integer :: focus_points=250
    logical :: filter_proposed=.true.
    real(dp) :: filter_tol=1.0e-6_dp
    integer :: batch_method=batch_none
    integer :: liar=lie_min
    integer :: multiobj_method=mo_dib
    integer :: dib_indicator=dib_sms
    real(dp) :: sms_eps=-1.0_dp
    integer :: ref_point_method=ref_all
    real(dp) :: ref_point_offset=1.0_dp
    real(dp), allocatable :: ref_point(:)
    real(dp) :: parego_rho=0.05_dp
    integer :: parego_s=100
    logical :: parego_normalize_front=.false.
    integer :: mspot_pool=1000
    character(len=16) :: covariance='matern5_2'
    logical :: cov_reestimate=.true.
    integer :: km_multistart=2
    integer :: km_pop_size=20
    integer :: km_max_iter=200
    real(dp) :: km_tol=1.0e-7_dp
    real(dp) :: target_value=0.0_dp
    logical :: use_target_value=.false.
    integer(i8) :: seed=104729_i8
  end type mbo_control

  type, public :: mbo_path
    real(dp), allocatable :: x(:,:),y(:,:)
    integer, allocatable :: dob(:)
    real(dp), allocatable :: crit(:)
    integer :: n=0
  end type mbo_path

  type, public :: mbo_result
    type(mbo_path) :: path
    real(dp), allocatable :: best_x(:),best_y(:)
    real(dp), allocatable :: pareto_x(:,:),pareto_y(:,:)
    integer :: iterations=0
    integer :: evaluations=0
    logical :: terminated=.false.
    character(len=24) :: termination=''
  end type mbo_result

  abstract interface
    subroutine mbo_objective(x,y)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(:)
    end subroutine mbo_objective
  end interface
  public :: mbo_objective
  public :: init_space, init_control, path_append
contains
  subroutine init_space(space,kind,lower,upper,nlevels,condition_parent,condition_level)
    type(mbo_space), intent(out) :: space
    integer, intent(in) :: kind(:)
    real(dp), intent(in) :: lower(:),upper(:)
    integer, intent(in), optional :: nlevels(:),condition_parent(:),condition_level(:)
    integer :: d
    d=size(kind)
    if(size(lower)/=d .or. size(upper)/=d) error stop 'init_space: dimension mismatch'
    if(any(upper<lower)) error stop 'init_space: upper < lower'
    space%d=d; space%kind=kind; space%lower=lower; space%upper=upper
    allocate(space%nlevels(d)); space%nlevels=0
    if(present(nlevels)) then
      if(size(nlevels)/=d) error stop 'init_space: nlevels length mismatch'
      space%nlevels=nlevels
    end if
    allocate(space%condition_parent(d),space%condition_level(d))
    space%condition_parent=0; space%condition_level=0
    if(present(condition_parent)) then
      if(size(condition_parent)/=d) error stop 'init_space: condition_parent length mismatch'
      space%condition_parent=condition_parent
    end if
    if(present(condition_level)) then
      if(size(condition_level)/=d) error stop 'init_space: condition_level length mismatch'
      space%condition_level=condition_level
    end if
  end subroutine init_space

  subroutine init_control(control,n_objectives,minimize)
    type(mbo_control), intent(out) :: control
    integer, intent(in) :: n_objectives
    logical, intent(in), optional :: minimize(:)
    control%n_objectives=n_objectives
    allocate(control%minimize(n_objectives)); control%minimize=.true.
    if(present(minimize)) then
      if(size(minimize)/=n_objectives) error stop 'init_control: minimize length mismatch'
      control%minimize=minimize
    end if
  end subroutine init_control

  subroutine repair_point(self,x)
    class(mbo_space), intent(in) :: self
    real(dp), intent(inout) :: x(:)
    integer :: j,nl
    if(size(x)/=self%d) error stop 'repair_point: dimension mismatch'
    do j=1,self%d
      x(j)=min(max(x(j),self%lower(j)),self%upper(j))
      select case(self%kind(j))
      case(mbo_integer)
        x(j)=real(nint(x(j)),dp)
        x(j)=min(max(x(j),self%lower(j)),self%upper(j))
      case(mbo_categorical)
        nl=self%nlevels(j)
        if(nl<=0) nl=max(1,nint(self%upper(j)-self%lower(j))+1)
        x(j)=real(min(nl,max(1,nint(x(j)))),dp)
      case default
      end select
      if(self%condition_parent(j)>0) then
        if(nint(x(self%condition_parent(j)))/=self%condition_level(j)) x(j)=self%lower(j)
      end if
    end do
  end subroutine repair_point

  logical function point_active(self,x,j) result(a)
    class(mbo_space), intent(in) :: self
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: j
    if(j<1 .or. j>self%d) error stop 'point_active: invalid index'
    if(self%condition_parent(j)==0) then
      a=.true.
    else
      a=(nint(x(self%condition_parent(j)))==self%condition_level(j))
    end if
  end function point_active

  subroutine path_append(path,xnew,ynew,dob,crit)
    type(mbo_path), intent(inout) :: path
    real(dp), intent(in) :: xnew(:,:),ynew(:,:)
    integer, intent(in) :: dob
    real(dp), intent(in), optional :: crit(:)
    real(dp), allocatable :: x(:,:),y(:,:),c(:)
    integer, allocatable :: d(:)
    integer :: old,k
    if(size(xnew,1)/=size(ynew,1)) error stop 'path_append: row mismatch'
    k=size(xnew,1); old=path%n
    if(old==0) then
      path%x=xnew; path%y=ynew
      allocate(path%dob(k),path%crit(k)); path%dob=dob; path%crit=0.0_dp
      if(present(crit)) path%crit=crit
    else
      allocate(x(old+k,size(path%x,2)),y(old+k,size(path%y,2)),d(old+k),c(old+k))
      x(1:old,:)=path%x; x(old+1:,:)=xnew
      y(1:old,:)=path%y; y(old+1:,:)=ynew
      d(1:old)=path%dob; d(old+1:)=dob
      c(1:old)=path%crit; c(old+1:)=0.0_dp
      if(present(crit)) c(old+1:)=crit
      call move_alloc(x,path%x); call move_alloc(y,path%y)
      call move_alloc(d,path%dob); call move_alloc(c,path%crit)
    end if
    path%n=old+k
  end subroutine path_append
end module mlrmbo_types
