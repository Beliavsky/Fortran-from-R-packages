! SPDX-License-Identifier: LGPL-2.0-or-later
module survival_utils
  use survival_kinds, only : dp
  implicit none
  private
  public :: surv_split, finegray_expand, pseudo_survival
contains

  subroutine surv_split(start,stop,status,cut,out_start,out_stop,out_status,row_id)
    real(dp), intent(in) :: start(:), stop(:), cut(:)
    integer, intent(in) :: status(:)
    real(dp), allocatable, intent(out) :: out_start(:), out_stop(:)
    integer, allocatable, intent(out) :: out_status(:), row_id(:)
    integer :: n, i, j, m, k
    real(dp) :: left

    n = size(stop)
    m = 0
    do i=1,n
      m = m + 1 + count(cut>start(i) .and. cut<stop(i))
    end do
    allocate(out_start(m),out_stop(m),out_status(m),row_id(m))
    k=0
    do i=1,n
      left=start(i)
      do j=1,size(cut)
        if(cut(j)>left .and. cut(j)<stop(i)) then
          k=k+1
          out_start(k)=left
          out_stop(k)=cut(j)
          out_status(k)=0
          row_id(k)=i
          left=cut(j)
        end if
      end do
      k=k+1
      out_start(k)=left
      out_stop(k)=stop(i)
      out_status(k)=status(i)
      row_id(k)=i
    end do
  end subroutine surv_split

  subroutine finegray_expand(tstart,tstop,ctime,cprob,extend,keep,row_id, &
                             out_start,out_stop,out_weight,add)
    real(dp), intent(in) :: tstart(:),tstop(:),ctime(:),cprob(:)
    logical, intent(in) :: extend(:),keep(:)
    integer, allocatable, intent(out) :: row_id(:),add(:)
    real(dp), allocatable, intent(out) :: out_start(:),out_stop(:),out_weight(:)
    integer :: n,ncut,extra,nout,i,j,k,iadd
    real(dp) :: tempwt

    n=size(tstart)
    ncut=size(cprob)
    extra=0
    do i=1,n
      if(.not.extend(i)) cycle
      j=1
      do while(j<=ncut .and. ctime(j)<tstop(i))
        j=j+1
      end do
      j=j+1
      do while(j<=ncut)
        if(keep(j)) extra=extra+1
        j=j+1
      end do
    end do
    nout=n+extra
    allocate(row_id(nout),add(nout),out_start(nout),out_stop(nout),out_weight(nout))
    k=0
    do i=1,n
      k=k+1
      row_id(k)=i
      out_start(k)=tstart(i)
      out_stop(k)=tstop(i)
      out_weight(k)=1.0_dp
      add(k)=0
      if(.not.extend(i)) cycle
      j=1
      do while(j<=ncut .and. ctime(j)<tstop(i))
        j=j+1
      end do
      if(j>ncut) cycle
      out_stop(k)=ctime(j)
      tempwt=cprob(j)
      j=j+1
      iadd=0
      do while(j<=ncut)
        if(keep(j)) then
          k=k+1
          iadd=iadd+1
          row_id(k)=i
          out_start(k)=ctime(j-1)
          out_stop(k)=ctime(j)
          out_weight(k)=cprob(j)/tempwt
          add(k)=iadd
        end if
        j=j+1
      end do
    end do
  end subroutine finegray_expand

  subroutine pseudo_survival(time,status,eval_time,pseudo)
    real(dp), intent(in) :: time(:), eval_time
    integer, intent(in) :: status(:)
    real(dp), intent(out) :: pseudo(:)
    integer :: n, i, j, k
    real(dp) :: full, loo
    real(dp), allocatable :: t(:)
    integer, allocatable :: s(:)

    n=size(time)
    full=km_at_time(time,status,eval_time)
    allocate(t(n-1),s(n-1))
    do i=1,n
      k=0
      do j=1,n
        if(j/=i) then
          k=k+1
          t(k)=time(j)
          s(k)=status(j)
        end if
      end do
      loo=km_at_time(t,s,eval_time)
      pseudo(i)=real(n,dp)*full-real(n-1,dp)*loo
    end do
  end subroutine pseudo_survival

  real(dp) function km_at_time(time,status,t0) result(surv)
    real(dp), intent(in) :: time(:), t0
    integer, intent(in) :: status(:)
    real(dp), allocatable :: u(:)
    integer :: i,k,m
    real(dp) :: nrisk,deaths

    allocate(u(size(time)))
    m=0
    do i=1,size(time)
      if(time(i)<=t0 .and. .not.contains_time(u,m,time(i))) then
        m=m+1
        u(m)=time(i)
      end if
    end do
    call sort_real(u(1:m))
    surv=1.0_dp
    do k=1,m
      nrisk=0.0_dp
      deaths=0.0_dp
      do i=1,size(time)
        if(time(i)>=u(k)) nrisk=nrisk+1.0_dp
        if(same_time(time(i),u(k)) .and. status(i)/=0) deaths=deaths+1.0_dp
      end do
      if(nrisk>0.0_dp) surv=surv*(nrisk-deaths)/nrisk
    end do
  end function km_at_time

  logical function contains_time(x,nused,value) result(found)
    real(dp), intent(in) :: x(:),value
    integer, intent(in) :: nused
    integer :: i
    found=.false.
    do i=1,nused
      if(same_time(x(i),value)) then
        found=.true.
        return
      end if
    end do
  end function contains_time

  pure logical function same_time(a,b) result(equal)
    real(dp), intent(in) :: a,b
    real(dp) :: scale
    scale=max(1.0_dp,abs(a),abs(b))
    equal=abs(a-b)<=8.0_dp*epsilon(1.0_dp)*scale
  end function same_time

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: value
    do i=2,size(x)
      value=x(i)
      j=i-1
      do while(j>=1)
        if(x(j)<=value) exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=value
    end do
  end subroutine sort_real
end module survival_utils
