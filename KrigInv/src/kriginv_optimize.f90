module kriginv_optimize
  use kriginv_kinds, only : dp
  implicit none
  private
  abstract interface
    function objective_function(x) result(v)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: v
    end function objective_function
  end interface
  type, public :: optimizer_control
    character(len=16) :: method='de'
    integer :: pop_size=0
    integer :: max_generations=0
    integer :: seed=12345
    real(dp) :: differential_weight=0.8_dp
    real(dp) :: crossover=0.9_dp
    real(dp), allocatable :: optim_points(:,:)
  end type optimizer_control
  public :: bounded_de, discrete_optimum
contains
  subroutine bounded_de(fun,lower,upper,maximize,ctl,best_x,best_v)
    procedure(objective_function) :: fun
    real(dp), intent(in) :: lower(:),upper(:)
    logical, intent(in) :: maximize
    type(optimizer_control), intent(in) :: ctl
    real(dp), allocatable, intent(out) :: best_x(:)
    real(dp), intent(out) :: best_v
    real(dp), allocatable :: pop(:,:),val(:),trial(:),r(:)
    real(dp) :: vt,f,cr,u
    integer :: d,np,ng,i,j,a,b,c,jrand,g
    d=size(lower); np=ctl%pop_size; if(np<=0) np=max(20,20*d)
    ng=ctl%max_generations; if(ng<=0) ng=max(20,15*d)
    f=ctl%differential_weight; cr=ctl%crossover
    call set_seed(ctl%seed)
    allocate(pop(np,d),val(np),trial(d),r(d)); call random_number(pop)
    do j=1,d; pop(:,j)=lower(j)+pop(:,j)*(upper(j)-lower(j)); end do
    do i=1,np; val(i)=fun(pop(i,:)); end do
    do g=1,ng
      do i=1,np
        call choose_three(np,i,a,b,c)
        call random_number(r); call random_number(u); jrand=1+int(u*real(d,dp)); jrand=min(d,jrand)
        trial=pop(i,:)
        do j=1,d
          if(r(j)<cr .or. j==jrand) trial(j)=pop(a,j)+f*(pop(b,j)-pop(c,j))
          trial(j)=max(lower(j),min(upper(j),trial(j)))
        end do
        vt=fun(trial)
        if(better(vt,val(i),maximize)) then; pop(i,:)=trial; val(i)=vt; end if
      end do
    end do
    i=1
    do j=2,np; if(better(val(j),val(i),maximize)) i=j; end do
    best_x=pop(i,:); best_v=val(i)
  end subroutine bounded_de

  subroutine discrete_optimum(fun,points,maximize,best_x,best_v,all_values)
    procedure(objective_function) :: fun
    real(dp), intent(in) :: points(:,:)
    logical, intent(in) :: maximize
    real(dp), allocatable, intent(out) :: best_x(:),all_values(:)
    real(dp), intent(out) :: best_v
    integer :: i,ibest
    allocate(all_values(size(points,1)))
    do i=1,size(points,1); all_values(i)=fun(points(i,:)); end do
    ibest=1
    do i=2,size(points,1); if(better(all_values(i),all_values(ibest),maximize)) ibest=i; end do
    best_x=points(ibest,:); best_v=all_values(ibest)
  end subroutine discrete_optimum

  logical function better(a,b,maximize) result(v)
    real(dp), intent(in) :: a,b
    logical, intent(in) :: maximize
    if(maximize) then; v=a>b; else; v=a<b; end if
  end function better

  subroutine choose_three(n,skip,a,b,c)
    integer, intent(in) :: n,skip
    integer, intent(out) :: a,b,c
    real(dp) :: u
    do
      call random_number(u); a=1+int(u*real(n,dp)); a=min(n,a)
      if(a/=skip) exit
    end do
    do
      call random_number(u); b=1+int(u*real(n,dp)); b=min(n,b)
      if(b/=skip .and. b/=a) exit
    end do
    do
      call random_number(u); c=1+int(u*real(n,dp)); c=min(n,c)
      if(c/=skip .and. c/=a .and. c/=b) exit
    end do
  end subroutine choose_three

  subroutine set_seed(seed)
    integer, intent(in) :: seed
    integer, allocatable :: put(:)
    integer :: n,i
    call random_seed(size=n); allocate(put(n))
    do i=1,n; put(i)=modulo(seed+7919*i,2147483646)+1; end do
    call random_seed(put=put)
  end subroutine set_seed
end module kriginv_optimize
