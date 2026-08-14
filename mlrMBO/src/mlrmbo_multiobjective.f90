module mlrmbo_multiobjective
  use mlrmbo_kinds, only : dp
  use mlrmbo_types, only : mbo_control, mbo_path, ref_all, ref_front, ref_const
  use mlrmbo_rng, only : mbo_rng
  implicit none
  private
  public :: dominated_mask, nondominated_points, dominated_hypervolume, hypervolume_contributions
  public :: eps_indicator_values, sms_indicator_values, reference_point
  public :: parego_weights, parego_scalarize
contains
  function dominated_mask(points) result(dom)
    real(dp), intent(in) :: points(:,:)
    logical, allocatable :: dom(:)
    integer :: n,i,j
    logical :: weak,strict
    n=size(points,1); allocate(dom(n)); dom=.false.
    do i=1,n
      do j=1,n
        if(j==i) cycle
        weak=all(points(j,:)<=points(i,:))
        strict=any(points(j,:)<points(i,:))
        if(weak .and. strict) then; dom(i)=.true.; exit; end if
      end do
    end do
  end function dominated_mask

  function nondominated_points(points) result(front)
    real(dp), intent(in) :: points(:,:)
    real(dp), allocatable :: front(:,:)
    logical, allocatable :: dom(:)
    integer :: i,k
    dom=dominated_mask(points); allocate(front(count(.not.dom),size(points,2)))
    k=0; do i=1,size(points,1); if(.not.dom(i)) then; k=k+1; front(k,:)=points(i,:); end if; end do
  end function nondominated_points

  recursive real(dp) function dominated_hypervolume(points,ref) result(hv)
    real(dp), intent(in) :: points(:,:),ref(:)
    real(dp), allocatable :: p(:,:),front(:,:),z(:),active(:,:)
    integer :: i,j,k,m,nz,na
    real(dp) :: width
    m=size(points,2)
    if(size(ref)/=m) error stop 'dominated_hypervolume: dimension mismatch'
    nfilter: block
      logical, allocatable :: keep(:)
      integer :: n
      allocate(keep(size(points,1))); keep=.true.
      do i=1,size(points,1)
        if(any(points(i,:)>=ref)) keep(i)=.false.
      end do
      n=count(keep); allocate(p(n,m)); k=0
      do i=1,size(points,1); if(keep(i)) then; k=k+1; p(k,:)=points(i,:); end if; end do
    end block nfilter
    if(size(p,1)==0) then; hv=0.0_dp; return; end if
    front=nondominated_points(p)
    if(m==1) then
      hv=max(0.0_dp,ref(1)-minval(front(:,1))); return
    end if
    allocate(z(size(front,1))); z=front(:,m); call sort_real(z)
    nz=unique_in_place(z)
    hv=0.0_dp
    do i=1,nz
      if(z(i)>=ref(m)) cycle
      if(i<nz) then; width=min(z(i+1),ref(m))-z(i); else; width=ref(m)-z(i); end if
      if(width<=0.0_dp) cycle
      na=count(front(:,m)<=z(i))
      allocate(active(na,m-1)); k=0
      do j=1,size(front,1)
        if(front(j,m)<=z(i)) then; k=k+1; active(k,:)=front(j,1:m-1); end if
      end do
      hv=hv+width*dominated_hypervolume(active,ref(1:m-1))
      deallocate(active)
    end do
  end function dominated_hypervolume

  function hypervolume_contributions(candidates,front,ref) result(c)
    real(dp), intent(in) :: candidates(:,:),front(:,:),ref(:)
    real(dp), allocatable :: c(:),aug(:,:)
    real(dp) :: base
    integer :: i,n
    n=size(front,1); allocate(c(size(candidates,1)))
    base=dominated_hypervolume(front,ref)
    allocate(aug(n+1,size(front,2))); if(n>0) aug(1:n,:)=front
    do i=1,size(candidates,1)
      aug(n+1,:)=candidates(i,:); c(i)=dominated_hypervolume(aug,ref)-base
    end do
  end function hypervolume_contributions

  function eps_indicator_values(candidates,front) result(v)
    real(dp), intent(in) :: candidates(:,:),front(:,:)
    real(dp), allocatable :: v(:)
    real(dp) :: d,dist
    integer :: i,j
    allocate(v(size(candidates,1)))
    do i=1,size(candidates,1)
      dist=huge(1.0_dp)
      do j=1,size(front,1)
        d=maxval(front(j,:)-candidates(i,:)); dist=min(dist,d)
      end do
      v(i)=-dist
    end do
  end function eps_indicator_values

  function sms_indicator_values(candidates,front,eps,ref) result(v)
    real(dp), intent(in) :: candidates(:,:),front(:,:),eps(:),ref(:)
    real(dp), allocatable :: v(:),aug(:,:)
    real(dp) :: penalty,pf,d,base
    integer :: i,j,k,n
    logical :: real_greater,valid
    allocate(v(size(candidates,1))); n=size(front,1)
    base=dominated_hypervolume(front,ref)
    allocate(aug(n+1,size(front,2))); if(n>0) aug(1:n,:)=front
    do i=1,size(candidates,1)
      penalty=0.0_dp
      do j=1,n
        pf=1.0_dp; real_greater=.false.; valid=.true.
        do k=1,size(front,2)
          d=candidates(i,k)-front(j,k)
          if(d < -eps(k)) then; valid=.false.; exit; end if
          if(d > -eps(k)) real_greater=.true.
          pf=pf*(1.0_dp+max(d,0.0_dp))
        end do
        if(valid .and. real_greater) then; pf=pf-1.0_dp; else; pf=0.0_dp; end if
        penalty=max(penalty,pf)
      end do
      if(penalty<=epsilon(1.0_dp)) then
        aug(n+1,:)=candidates(i,:)
        v(i)=-(dominated_hypervolume(aug,ref)-base)
      else
        v(i)=penalty
      end if
    end do
  end function sms_indicator_values

  function reference_point(y,control) result(ref)
    real(dp), intent(in) :: y(:,:)
    type(mbo_control), intent(in) :: control
    real(dp), allocatable :: ref(:),front(:,:)
    integer :: j
    allocate(ref(size(y,2)))
    select case(control%ref_point_method)
    case(ref_const)
      if(.not.allocated(control%ref_point)) error stop 'reference_point: missing constant reference point'
      ref=control%ref_point
    case(ref_front)
      front=nondominated_points(y)
      do j=1,size(y,2); ref(j)=maxval(front(:,j))+control%ref_point_offset; end do
    case default
      do j=1,size(y,2); ref(j)=maxval(y(:,j))+control%ref_point_offset; end do
    end select
  end function reference_point

  subroutine parego_weights(rng,q,m,w,s_grid)
    type(mbo_rng), intent(inout) :: rng
    integer, intent(in) :: q,m
    real(dp), allocatable, intent(out) :: w(:,:)
    integer, intent(in), optional :: s_grid
    integer, allocatable :: sep(:),cuts(:),counts(:)
    integer :: i,j,k,s,npos,tries,it
    logical :: duplicate
    if(m<1 .or. q<1) error stop 'parego_weights: invalid dimensions'
    s=100; if(present(s_grid)) s=max(1,s_grid)
    allocate(w(q,m),counts(m))
    if(m==1) then; w=1.0_dp; return; end if
    npos=s+m-1; allocate(sep(npos),cuts(m-1))
    do i=1,q
      tries=0
      do
        tries=tries+1; sep=[(j,j=1,npos)]
        do j=npos,2,-1
          k=rng%randint(1,j); it=sep(j); sep(j)=sep(k); sep(k)=it
        end do
        cuts=sep(1:m-1); call sort_int(cuts)
        counts(1)=cuts(1)-1
        do j=2,m-1; counts(j)=cuts(j)-cuts(j-1)-1; end do
        counts(m)=npos-cuts(m-1)
        w(i,:)=real(counts,dp)/real(s,dp)
        duplicate=.false.
        if(i>1) then
          do j=1,i-1
            if(maxval(abs(w(i,:)-w(j,:)))<0.5_dp/real(s,dp)) then; duplicate=.true.; exit; end if
          end do
        end if
        if(.not.duplicate .or. tries>200) exit
      end do
    end do
  end subroutine parego_weights

  subroutine sort_int(x)
    integer, intent(inout) :: x(:)
    integer :: i,j,t
    do i=2,size(x)
      t=x(i); j=i-1
      do while(j>=1)
        if(x(j)<=t) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=t
    end do
  end subroutine sort_int

  subroutine parego_scalarize(y,minimize,w,rho,ys)
    real(dp), intent(in) :: y(:,:),w(:),rho
    logical, intent(in) :: minimize(:)
    real(dp), allocatable, intent(out) :: ys(:)
    real(dp), allocatable :: z(:,:),lo(:),hi(:),row(:)
    integer :: i,j,n,m
    n=size(y,1); m=size(y,2)
    if(size(w)/=m .or. size(minimize)/=m) error stop 'parego_scalarize: dimension mismatch'
    allocate(z(n,m),lo(m),hi(m),ys(n),row(m)); z=y
    do j=1,m
      if(.not.minimize(j)) z(:,j)=-z(:,j)
      lo(j)=minval(z(:,j)); hi(j)=maxval(z(:,j))
      if(hi(j)>lo(j)) then; z(:,j)=(z(:,j)-lo(j))/(hi(j)-lo(j)); else; z(:,j)=0.0_dp; end if
    end do
    do i=1,n
      row=w*z(i,:); ys(i)=maxval(row)+rho*sum(row)
    end do
  end subroutine parego_scalarize

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: t
    do i=2,size(x); t=x(i); j=i-1; do while(j>=1); if(x(j)<=t) exit; x(j+1)=x(j); j=j-1; end do; x(j+1)=t; end do
  end subroutine sort_real

  integer function unique_in_place(x) result(n)
    real(dp), intent(inout) :: x(:)
    integer :: i
    if(size(x)==0) then; n=0; return; end if
    n=1
    do i=2,size(x)
      if(x(i)>x(n)) then; n=n+1; x(n)=x(i); end if
    end do
  end function unique_in_place
end module mlrmbo_multiobjective
