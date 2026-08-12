! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_geometry
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use adagio_kinds, only : dp
  use adagio_types, only : maxsub_result, maxsub2d_result, maxempty_result, count_result
  use adagio_utils, only : argsort_real, sort_real
  implicit none
  private
  public :: maxsub, maxsub2d, maxempty, count_values, occurs
  public :: transfinite_forward, transfinite_inverse

contains

  function maxsub(x) result(res)
    real(dp), intent(in) :: x(:)
    type(maxsub_result) :: res
    real(dp) :: best, cur
    integer :: i, start
    best = 0.0_dp; cur = 0.0_dp; start = 1
    res%first = 0; res%last = 0
    do i = 1, size(x)
       if (cur > -x(i)) then
          cur = cur + x(i)
          if (cur > best) then
             best = cur; res%first = start; res%last = i
          end if
       else
          cur = 0.0_dp; start = i + 1
       end if
    end do
    res%sum = best
  end function maxsub

  function maxsub2d(a) result(res)
    real(dp), intent(in) :: a(:,:)
    type(maxsub2d_result) :: res
    real(dp), allocatable :: colsum(:)
    type(maxsub_result) :: one
    real(dp) :: best
    integer :: top, bot, nr, nc
    nr = size(a,1); nc = size(a,2)
    allocate(colsum(nc)); best = 0.0_dp
    res%inds = 0
    if (all(a >= 0.0_dp)) then
       res%sum = sum(a); res%inds = [1,nr,1,nc]
       allocate(res%submat(nr,nc)); res%submat = a; return
    end if
    do top = 1, nr
       colsum = 0.0_dp
       do bot = top, nr
          colsum = colsum + a(bot,:)
          one = maxsub(colsum)
          if (one%sum > best) then
             best = one%sum
             res%inds = [top, bot, one%first, one%last]
          end if
       end do
    end do
    res%sum = best
    if (best > 0.0_dp) then
       allocate(res%submat(res%inds(2)-res%inds(1)+1, res%inds(4)-res%inds(3)+1))
       res%submat = a(res%inds(1):res%inds(2),res%inds(3):res%inds(4))
    else
       allocate(res%submat(0,0))
    end if
  end function maxsub2d

  function maxempty(x, y, ax, ay) result(res)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(in), optional :: ax(2), ay(2)
    type(maxempty_result) :: res
    real(dp) :: bx(2), by(2), area, tl, tr, ri, li
    real(dp), allocatable :: d(:), xs(:), ys(:)
    integer, allocatable :: ord(:)
    integer :: i, j, n
    bx = [0.0_dp,1.0_dp]; by = [0.0_dp,1.0_dp]
    if (present(ax)) bx=ax; if (present(ay)) by=ay
    n = size(x)
    allocate(d(n+2)); d=[bx(1),x,bx(2)]; call sort_real(d)
    res%area=-1.0_dp
    do i=1,n+1
       area=(d(i+1)-d(i))*(by(2)-by(1))
       if (area>res%area) then
          res%area=area; res%rect=[d(i),by(1),d(i+1),by(2)]
       end if
    end do
    allocate(ord(n),xs(n),ys(n)); call argsort_real(y,ord)
    xs=x(ord); ys=y(ord)
    do i=1,n
       tl=bx(1); tr=bx(2)
       do j=i+1,n
          if (xs(j)>tl .and. xs(j)<tr) then
             area=(tr-tl)*(ys(j)-ys(i))
             if (area>res%area) then
                res%area=area; res%rect=[tl,ys(i),tr,ys(j)]
             end if
             if (xs(j)>xs(i)) then; tr=xs(j); else; tl=xs(j); end if
          end if
       end do
       area=(tr-tl)*(by(2)-ys(i))
       if(area>res%area) then
          res%area=area; res%rect=[tl,ys(i),tr,by(2)]
       end if
       ri=bx(2); li=bx(1)
       do j=1,n
          if (ys(j)<ys(i) .and. xs(j)>xs(i)) ri=min(ri,xs(j))
          if (ys(j)<ys(i) .and. xs(j)<xs(i)) li=max(li,xs(j))
       end do
       area=(ri-li)*(ys(i)-by(1))
       if(area>res%area) then
          res%area=area; res%rect=[li,by(1),ri,ys(i)]
       end if
    end do
  end function maxempty

  function count_values(x, sorted) result(res)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: sorted
    type(count_result) :: res
    real(dp), allocatable :: z(:), vals(:)
    integer, allocatable :: cnt(:)
    logical :: dosort
    integer :: i, k, n
    n=size(x); dosort=.true.; if(present(sorted)) dosort=sorted
    allocate(z(n)); z=x
    if(dosort) call sort_real(z)
    allocate(vals(n),cnt(n)); k=1; vals(1)=z(1); cnt(1)=1
    do i=2,n
       if(z(i) <= vals(k) .and. z(i) >= vals(k)) then
          cnt(k)=cnt(k)+1
       else
          k=k+1; vals(k)=z(i); cnt(k)=1
       end if
    end do
    allocate(res%values(k),res%counts(k)); res%values=vals(1:k); res%counts=cnt(1:k)
  end function count_values

  function occurs(subseq, series) result(inds)
    real(dp), intent(in) :: subseq(:), series(:)
    integer, allocatable :: inds(:)
    integer, allocatable :: tmp(:)
    integer :: i, j, k, m, n
    m=size(subseq); n=size(series)
    if(m>n) then; allocate(inds(0)); return; end if
    allocate(tmp(n-m+1)); k=0
    do i=1,n-m+1
       do j=1,m
          if (.not. (series(i+j-1) <= subseq(j) .and. series(i+j-1) >= subseq(j))) exit
       end do
       if(j>m) then; k=k+1; tmp(k)=i; end if
    end do
    allocate(inds(k)); if(k>0) inds=tmp(1:k)
  end function occurs

  function transfinite_forward(x, lower, upper) result(hx)
    real(dp), intent(in) :: x(:), lower(:), upper(:)
    real(dp) :: hx(size(x)), qnan
    integer :: i
    qnan=ieee_value(0.0_dp,ieee_quiet_nan)
    do i=1,size(x)
       if(x(i)<lower(i) .or. x(i)>upper(i)) then; hx=qnan; return; end if
       if(ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i))) then
          hx(i)=atanh(2.0_dp*(x(i)-lower(i))/(upper(i)-lower(i))-1.0_dp)
       else if(.not.ieee_is_finite(lower(i)) .and. .not.ieee_is_finite(upper(i))) then
          hx(i)=x(i)
       else if(ieee_is_finite(lower(i))) then
          hx(i)=log(x(i)-lower(i))
       else
          hx(i)=log(upper(i)-x(i))
       end if
    end do
  end function transfinite_forward

  function transfinite_inverse(x, lower, upper) result(hx)
    real(dp), intent(in) :: x(:), lower(:), upper(:)
    real(dp) :: hx(size(x))
    integer :: i
    do i=1,size(x)
       if(ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i))) then
          hx(i)=lower(i)+(upper(i)-lower(i))/2.0_dp*(1.0_dp+tanh(x(i)))
       else if(.not.ieee_is_finite(lower(i)) .and. .not.ieee_is_finite(upper(i))) then
          hx(i)=x(i)
       else if(ieee_is_finite(lower(i))) then
          hx(i)=lower(i)+exp(x(i))
       else
          hx(i)=upper(i)-exp(x(i))
       end if
    end do
  end function transfinite_inverse

end module adagio_geometry
