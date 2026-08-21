module relsurv_utils
  use relsurv_kinds, only : dp
  use relsurv_ratetable, only : ratetable_type, make_ratetable
  implicit none
  private
  public :: transrate, join_ratetables, epanechnikov_smooth, years_auc
  public :: epa_smooth, epanechnikov_boundary_matrix
contains
  function transrate(men,women,year_cut,days_per_year) result(tab)
    real(dp),intent(in)::men(:,:),women(:,:),year_cut(:)
    real(dp),intent(in),optional::days_per_year
    type(ratetable_type)::tab
    real(dp)::dpy
    integer::na,ny,a,s,y,idx,maxcut
    integer,allocatable::dims(:),factor(:),ncuts(:)
    real(dp),allocatable::cuts(:,:),rate(:)
    dpy=365.241_dp; if(present(days_per_year))dpy=days_per_year
    na=size(men,1); ny=size(men,2)
    if(size(women,1)/=na .or. size(women,2)/=ny .or. size(year_cut)/=ny) error stop 'transrate: dimension mismatch'
    allocate(dims(3),factor(3),ncuts(3)); dims=[na,2,ny]; factor=[0,1,0]; ncuts=[na,0,ny]
    maxcut=max(na,ny); allocate(cuts(maxcut,3)); cuts=0.0_dp
    do a=1,na; cuts(a,1)=real(a-1,dp)*dpy; end do
    cuts(1:ny,3)=year_cut
    allocate(rate(na*2*ny)); rate=0.0_dp
    do y=1,ny
      do s=1,2
        do a=1,na
          idx=a+(s-1)*na+(y-1)*na*2
          if(s==1) rate(idx)=-log(max(men(a,y),tiny(1.0_dp)))/dpy
          if(s==2) rate(idx)=-log(max(women(a,y),tiny(1.0_dp)))/dpy
        end do
      end do
    end do
    tab=make_ratetable(dims,factor,cuts,ncuts,rate)
  end function transrate

  function join_ratetables(tables) result(out)
    type(ratetable_type),intent(in)::tables(:)
    type(ratetable_type)::out
    integer::nt,nd,nbase,i,maxc
    integer,allocatable::dims(:),factor(:),ncuts(:)
    real(dp),allocatable::cuts(:,:),rates(:)
    nt=size(tables); if(nt<1) error stop 'join_ratetables: empty list'
    nd=tables(1)%ndim; nbase=size(tables(1)%rate)
    do i=2,nt
      if(tables(i)%ndim/=nd .or. size(tables(i)%rate)/=nbase) error stop 'join_ratetables: incompatible tables'
      if(any(tables(i)%dims/=tables(1)%dims) .or. &
         any(tables(i)%factor/=tables(1)%factor)) then
        error stop 'join_ratetables: incompatible dimensions'
      end if
    end do
    allocate(dims(nd+1),factor(nd+1),ncuts(nd+1)); dims(1:nd)=tables(1)%dims; dims(nd+1)=nt
    factor(1:nd)=tables(1)%factor; factor(nd+1)=1; ncuts(1:nd)=tables(1)%ncuts; ncuts(nd+1)=0
    maxc=size(tables(1)%cuts,1); allocate(cuts(maxc,nd+1)); cuts=0.0_dp; cuts(:,1:nd)=tables(1)%cuts
    allocate(rates(nbase*nt)); do i=1,nt; rates((i-1)*nbase+1:i*nbase)=tables(i)%rate; end do
    out=make_ratetable(dims,factor,cuts,ncuts,rates)
  end function join_ratetables

  subroutine epanechnikov_smooth(source_time,source_value,target_time,bandwidth,value,left_only)
    real(dp),intent(in)::source_time(:),source_value(:),target_time(:),bandwidth
    real(dp),intent(out)::value(:)
    logical,intent(in),optional::left_only
    logical::left
    real(dp)::u,w,sw
    integer::i,j
    left=.false.; if(present(left_only))left=left_only
    do i=1,size(target_time)
      value(i)=0.0_dp; sw=0.0_dp
      do j=1,size(source_time)
        if(left .and. source_time(j)>target_time(i))cycle
        u=(target_time(i)-source_time(j))/bandwidth
        if(abs(u)<=1.0_dp) then
          w=0.75_dp*(1.0_dp-u*u)/bandwidth
          value(i)=value(i)+w*source_value(j); sw=sw+w
        end if
      end do
      if(sw>0.0_dp)value(i)=max(0.0_dp,value(i)/sw)
    end do
  end subroutine epanechnikov_smooth

  subroutine years_auc(time,survival,auc_lost,auc_lost_cumulative,scale)
    real(dp),intent(in)::time(:),survival(:)
    real(dp),intent(out)::auc_lost
    real(dp),intent(out),optional::auc_lost_cumulative(:)
    real(dp),intent(in),optional::scale
    real(dp)::sc,prev,inc,cum
    integer::i
    sc=365.241_dp; if(present(scale))sc=scale; prev=0.0_dp; cum=0.0_dp
    do i=1,size(time)
      inc=(time(i)-prev)*(1.0_dp-survival(i)); cum=cum+inc
      if(present(auc_lost_cumulative))auc_lost_cumulative(i)=cum/sc
      prev=time(i)
    end do
    auc_lost=cum/sc
  end subroutine years_auc


  subroutine epa_smooth(event_time, lambda_ns, target_time, value, bwin, n_bwin, left)
    real(dp), intent(in) :: event_time(:), lambda_ns(:), target_time(:)
    real(dp), intent(out) :: value(:)
    real(dp), intent(in), optional :: bwin
    integer, intent(in), optional :: n_bwin
    logical, intent(in), optional :: left
    integer :: n, nb, i, j, k, i1, i2, nseg
    real(dp) :: bwrel, gap
    logical :: lft
    integer, allocatable :: cuts(:), cuts2(:)
    real(dp), allocatable :: bw(:), kmat(:,:), eval_time(:)

    n = size(event_time)
    if (size(lambda_ns) /= n) error stop 'epa_smooth: event/lambda shape'
    if (size(value) /= size(target_time)) error stop 'epa_smooth: target/value shape'
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    nb = 16
    if (present(n_bwin)) nb = max(1,n_bwin)
    bwrel = real(n,dp)/100.0_dp
    if (present(bwin)) bwrel = bwin*real(n,dp)/100.0_dp
    lft = .false.
    if (present(left)) lft = left

    allocate(cuts(nb+1))
    cuts(1)=1
    do i=1,nb
      cuts(i+1)=max(1,min(n,ceiling(real(n*i,dp)/real(nb,dp))))
    end do
    call compress_cuts(cuts,cuts2)
    call move_alloc(cuts2,cuts)
    nseg=max(1,size(cuts)-1)
    allocate(bw(nseg))
    do i=1,nseg
      i1=cuts(i);i2=cuts(i+1)
      gap=0.0_dp
      if(i2>i1)then
        do j=i1+1,i2
          gap=max(gap,event_time(j)-event_time(j-1))
        end do
      end if
      if(gap<=0.0_dp)gap=max(event_time(max(i2,2))-event_time(max(1,i2-1)),epsilon(1.0_dp))
      bw(i)=max(bwrel*gap,epsilon(1.0_dp))
    end do

    do while(size(cuts)>2 .and. event_time(cuts(2))<bw(1))
      allocate(cuts2(size(cuts)-1))
      cuts2(1)=cuts(1);cuts2(2:)=cuts(3:)
      call move_alloc(cuts2,cuts)
      nseg=size(cuts)-1
      if(size(bw)/=nseg)then
        deallocate(bw);allocate(bw(nseg))
        do i=1,nseg
          i1=cuts(i);i2=cuts(i+1);gap=0.0_dp
          do j=i1+1,i2;gap=max(gap,event_time(j)-event_time(j-1));end do
          if(gap<=0.0_dp)gap=epsilon(1.0_dp)
          bw(i)=max(bwrel*gap,epsilon(1.0_dp))
        end do
      end if
    end do

    if(lft)then
      ! Upstream epa(left=TRUE) smooths only at event times.  For arbitrary
      ! target_time, evaluate the same predictable kernel at the requested points.
      allocate(kmat(size(target_time),n))
      call left_kernel_matrix(target_time,event_time,bw,cuts,kmat)
    else
      allocate(kmat(size(target_time),n))
      call epanechnikov_boundary_matrix(target_time,event_time,bw,cuts,kmat)
    end if
    value=max(matmul(kmat,lambda_ns),0.0_dp)
  end subroutine epa_smooth

  subroutine epanechnikov_boundary_matrix(target, source, bandwidth, cuts, kernel)
    real(dp),intent(in)::target(:),source(:),bandwidth(:)
    integer,intent(in)::cuts(:)
    real(dp),intent(out)::kernel(:,:)
    integer::i,j,s
    real(dp)::b,q,u,rb,uu
    rb=maxval(source);kernel=0.0_dp
    do i=1,size(target)
      s=segment_for_time(target(i),source,cuts)
      if(s<1)cycle
      b=bandwidth(min(s,size(bandwidth)))
      q=min(target(i)/b,1.0_dp,(rb-target(i))/b)
      q=max(0.0_dp,q)
      do j=1,size(source)
        u=(target(i)-source(j))/b
        if(q<1.0_dp)then
          if(target(i)>b)then
            uu=-u
          else
            uu=u
          end if
          if(uu>=-1.0_dp .and. uu<=q)then
            kernel(i,j)=boundary_kernel(q,uu)/b
          end if
        else if(abs(u)<=1.0_dp)then
          kernel(i,j)=0.75_dp*(1.0_dp-u*u)/b
        end if
      end do
    end do
  end subroutine epanechnikov_boundary_matrix

  subroutine left_kernel_matrix(target,source,bandwidth,cuts,kernel)
    real(dp),intent(in)::target(:),source(:),bandwidth(:)
    integer,intent(in)::cuts(:)
    real(dp),intent(out)::kernel(:,:)
    integer::i,j,s
    real(dp)::b,u
    kernel=0.0_dp
    do i=1,size(target)
      s=segment_for_time(target(i),source,cuts)
      if(s<1)cycle
      b=bandwidth(min(s,size(bandwidth)))
      do j=1,size(source)
        u=(target(i)-source(j))/b
        if(u>=0.0_dp .and. u<=1.0_dp)kernel(i,j)=1.5_dp*(1.0_dp-u*u)/b
      end do
    end do
  end subroutine left_kernel_matrix

  pure function boundary_kernel(q,t) result(v)
    real(dp),intent(in)::q,t
    real(dp)::v
    v=12.0_dp*(t+1.0_dp)/(1.0_dp+q)**4 * &
      ((1.0_dp-2.0_dp*q)*t + (3.0_dp*q*q-2.0_dp*q+1.0_dp)/2.0_dp)
  end function boundary_kernel

  integer function segment_for_time(t,source,cuts) result(seg)
    real(dp),intent(in)::t,source(:)
    integer,intent(in)::cuts(:)
    integer::i
    seg=0
    do i=1,size(cuts)-1
      if(t>source(cuts(i)) .and. t<=source(cuts(i+1)))then;seg=i;return;end if
    end do
    if(t<=source(1))seg=1
    if(t>source(size(source)))seg=size(cuts)-1
  end function segment_for_time

  subroutine compress_cuts(cuts,out)
    integer,intent(in)::cuts(:)
    integer,allocatable,intent(out)::out(:)
    integer::i,n
    n=1
    do i=2,size(cuts);if(cuts(i)/=cuts(i-1))n=n+1;end do
    allocate(out(n));out(1)=cuts(1);n=1
    do i=2,size(cuts)
      if(cuts(i)/=cuts(i-1))then;n=n+1;out(n)=cuts(i);end if
    end do
    if(size(out)==1)then
      deallocate(out);allocate(out(2));out=[1,max(1,cuts(size(cuts)))]
    end if
  end subroutine compress_cuts

end module relsurv_utils
