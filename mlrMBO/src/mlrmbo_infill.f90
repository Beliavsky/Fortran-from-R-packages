module mlrmbo_infill
  use mlrmbo_kinds, only : dp
  use mlrmbo_types, only : mbo_space, mbo_control, mbo_path, mbo_real, mbo_integer
  use mlrmbo_rng, only : mbo_rng
  use mlrmbo_math, only : lhs_design
  implicit none
  private
  abstract interface
    subroutine infill_evaluator(x,v)
      import dp
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: v(:)
    end subroutine infill_evaluator
  end interface
  public :: sample_space_lhs, sample_space_random, focus_search, points_too_close
contains
  subroutine sample_space_lhs(space,rng,n,x,lo,hi)
    type(mbo_space), intent(in) :: space
    type(mbo_rng), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: x(:,:)
    real(dp), intent(in), optional :: lo(:),hi(:)
    real(dp), allocatable :: l(:),u(:),row(:)
    integer :: i
    allocate(l(space%d),u(space%d))
    if(present(lo)) then
      l=lo
    else
      l=space%lower
    end if
    if(present(hi)) then
      u=hi
    else
      u=space%upper
    end if
    call lhs_design(rng,n,l,u,x)
    allocate(row(space%d))
    do i=1,n; row=x(i,:); call space%repair(row); x(i,:)=row; end do
  end subroutine sample_space_lhs

  subroutine sample_space_random(space,rng,n,x)
    type(mbo_space), intent(in) :: space
    type(mbo_rng), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: x(:,:)
    real(dp), allocatable :: row(:)
    integer :: i,j
    allocate(x(n,space%d),row(space%d))
    do i=1,n
      do j=1,space%d; row(j)=space%lower(j)+(space%upper(j)-space%lower(j))*rng%uniform(); end do
      call space%repair(row); x(i,:)=row
    end do
  end subroutine sample_space_random

  subroutine focus_search(space,rng,control,evaluator,xbest,bestval,avoid)
    type(mbo_space), intent(in) :: space
    type(mbo_rng), intent(inout) :: rng
    type(mbo_control), intent(in) :: control
    procedure(infill_evaluator) :: evaluator
    real(dp), allocatable, intent(out) :: xbest(:)
    real(dp), intent(out) :: bestval
    real(dp), intent(in), optional :: avoid(:,:)
    real(dp), allocatable :: lo(:),hi(:),x(:,:),v(:),localx(:),row(:)
    real(dp) :: localv,range
    integer :: r,it,j,idx
    bestval=huge(1.0_dp); allocate(xbest(space%d)); xbest=space%lower
    do r=1,max(1,control%focus_restarts)
      lo=space%lower; hi=space%upper
      do it=1,max(1,control%focus_iterations)
        call sample_space_lhs(space,rng,max(2,control%focus_points),x,lo,hi)
        if(present(avoid)) call replace_close_points(space,rng,x,avoid,control%filter_tol)
        call evaluator(x,v)
        idx=minloc(v,dim=1); localv=v(idx); localx=x(idx,:)
        if(localv<bestval) then; bestval=localv; xbest=localx; end if
        do j=1,space%d
          if(space%kind(j)==mbo_real .or. space%kind(j)==mbo_integer) then
            range=hi(j)-lo(j)
            lo(j)=max(space%lower(j),localx(j)-0.25_dp*range)
            hi(j)=min(space%upper(j),localx(j)+0.25_dp*range)
            if(space%kind(j)==mbo_integer) then
              lo(j)=real(floor(lo(j)),dp); hi(j)=real(ceiling(hi(j)),dp)
            end if
          end if
        end do
      end do
    end do
    if(present(avoid)) then
      if(points_too_close(xbest,avoid,control%filter_tol)) then
        allocate(row(space%d)); call random_nonclose(space,rng,avoid,control%filter_tol,row); xbest=row
        block
          real(dp), allocatable :: xx(:,:), vv(:)
          allocate(xx(1,space%d)); xx(1,:)=xbest; call evaluator(xx,vv); bestval=vv(1)
        end block
      end if
    end if
  end subroutine focus_search

  logical function points_too_close(x,points,tol) result(close)
    real(dp), intent(in) :: x(:),points(:,:),tol
    integer :: i
    close=.false.; if(size(points,1)==0) return
    do i=1,size(points,1)
      if(maxval(abs(x-points(i,:)))<tol) then; close=.true.; return; end if
    end do
  end function points_too_close

  subroutine replace_close_points(space,rng,x,avoid,tol)
    type(mbo_space), intent(in) :: space
    type(mbo_rng), intent(inout) :: rng
    real(dp), intent(inout) :: x(:,:)
    real(dp), intent(in) :: avoid(:,:),tol
    real(dp), allocatable :: row(:)
    integer :: i
    allocate(row(space%d))
    do i=1,size(x,1)
      if(points_too_close(x(i,:),avoid,tol)) then
        call random_nonclose(space,rng,avoid,tol,row); x(i,:)=row
      end if
    end do
  end subroutine replace_close_points

  subroutine random_nonclose(space,rng,avoid,tol,row)
    type(mbo_space), intent(in) :: space
    type(mbo_rng), intent(inout) :: rng
    real(dp), intent(in) :: avoid(:,:),tol
    real(dp), intent(out) :: row(:)
    integer :: j,tries
    do tries=1,1000
      do j=1,space%d; row(j)=space%lower(j)+(space%upper(j)-space%lower(j))*rng%uniform(); end do
      call space%repair(row)
      if(.not.points_too_close(row,avoid,tol)) return
    end do
  end subroutine random_nonclose
end module mlrmbo_infill
