! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_data
  use acdm_kinds, only : dp, tiny_pos, ACDM_SUCCESS, ACDM_BAD_INPUT, &
                         ACDM_NUMERIC_FAILURE
  use acdm_math, only : natural_spline_second, natural_spline_eval, &
                        least_squares, solve_linear
  implicit none
  private

  integer, parameter, public :: DURATION_TRADE = 1
  integer, parameter, public :: DURATION_PRICE = 2
  integer, parameter, public :: DURATION_VOLUME = 3
  integer, parameter, public :: DIURNAL_CUBIC_SPLINE = 1
  integer, parameter, public :: DIURNAL_SMOOTH_SPLINE = 2
  integer, parameter, public :: DIURNAL_SUPER_SMOOTHER = 3
  integer, parameter, public :: DIURNAL_FFF = 4

  type, public :: duration_result
    integer, allocatable :: year(:), month(:), day(:)
    real(dp), allocatable :: second_of_day(:)
    real(dp), allocatable :: durations(:), price(:), volume(:)
    integer, allocatable :: transaction_count(:)
    integer :: status = ACDM_BAD_INPUT
  end type duration_result

  type, public :: diurnal_result
    real(dp), allocatable :: adjusted(:)
    real(dp), allocatable :: fitted(:)
    real(dp), allocatable :: grid_time(:)
    real(dp), allocatable :: grid_fit(:, :)
    integer, allocatable :: group_labels(:)
    integer :: status = ACDM_BAD_INPUT
  end type diurnal_result

  public :: compute_durations, diurnal_adjust, flexible_fourier_fit
  public :: super_smoother_fit

contains

  subroutine compute_durations(year, month, day, second_of_day, result, &
                               duration_type, open_second, close_second, &
                               remove_zero, price, volume, price_difference, &
                               cumulative_volume)
    integer, intent(in) :: year(:), month(:), day(:)
    real(dp), intent(in) :: second_of_day(:)
    type(duration_result), intent(out) :: result
    integer, intent(in), optional :: duration_type
    real(dp), intent(in), optional :: open_second, close_second
    logical, intent(in), optional :: remove_zero
    real(dp), intent(in), optional :: price(:), volume(:)
    real(dp), intent(in), optional :: price_difference, cumulative_volume

    integer :: n, typ, i, j, k, first, last, count, nout
    real(dp) :: opn, cls, pdiff, cvol, last_time, last_price
    real(dp) :: this_time, sum_volume
    logical :: rmzero, emit
    integer, allocatable :: oy(:), om(:), od(:), on(:)
    real(dp), allocatable :: ot(:), ox(:), op(:), ov(:)

    result%status = ACDM_BAD_INPUT
    n = size(second_of_day)
    if (n < 1 .or. size(year) /= n .or. size(month) /= n .or. size(day) /= n) return
    if (any(second_of_day < 0.0_dp) .or. any(second_of_day >= 86400.0_dp)) return
    typ = DURATION_TRADE
    if (present(duration_type)) typ = duration_type
    if (typ < DURATION_TRADE .or. typ > DURATION_VOLUME) return
    if (typ == DURATION_PRICE .and. .not. present(price)) return
    if (typ == DURATION_VOLUME .and. .not. present(volume)) return
    if (present(price)) then
      if (size(price) /= n) return
    end if
    if (present(volume)) then
      if (size(volume) /= n) return
    end if
    opn = 10.0_dp * 3600.0_dp
    cls = 18.0_dp * 3600.0_dp + 25.0_dp * 60.0_dp
    if (present(open_second)) opn = open_second
    if (present(close_second)) cls = close_second
    if (opn < 0.0_dp .or. cls <= opn .or. cls >= 86400.0_dp) return
    rmzero = .true.
    if (present(remove_zero)) rmzero = remove_zero
    pdiff = 0.1_dp
    if (present(price_difference)) pdiff = price_difference
    cvol = 10000.0_dp
    if (present(cumulative_volume)) cvol = cumulative_volume
    if (pdiff < 0.0_dp .or. cvol <= 0.0_dp) return

    ! Reject unordered input. Cross-day ordering is checked lexicographically.
    do i = 2, n
      if (date_key(year(i), month(i), day(i)) < date_key(year(i-1), month(i-1), day(i-1))) return
      if (date_key(year(i), month(i), day(i)) == date_key(year(i-1), month(i-1), day(i-1))) then
        if (second_of_day(i) < second_of_day(i-1)) return
      end if
    end do

    allocate(oy(n), om(n), od(n), on(n), ot(n), ox(n), op(n), ov(n))
    oy = 0; om = 0; od = 0; on = 0
    ot = 0.0_dp; ox = 0.0_dp; op = 0.0_dp; ov = 0.0_dp
    nout = 0
    last_time = opn
    last_price = 0.0_dp
    this_time = 0.0_dp
    sum_volume = 0.0_dp
    first = 1
    do while (first <= n)
      last = first
      do while (last < n)
        if (year(last+1) /= year(first) .or. month(last+1) /= month(first) .or. &
            day(last+1) /= day(first)) exit
        last = last + 1
      end do

      i = first
      do while (i <= last)
        if (second_of_day(i) >= opn) exit
        i = i + 1
      end do
      if (i <= last) then
        if (second_of_day(i) > cls) then
          first = last + 1
          cycle
        end if
        last_time = opn
        if (typ == DURATION_PRICE) then
          last_price = price(i)
          ! A quote exactly at the opening establishes the baseline.
          do while (i <= last)
            if (second_of_day(i) > opn) exit
            last_price = price(i)
            i = i + 1
          end do
          if (i <= last) then
            if (last_time < second_of_day(i) .and. second_of_day(i) > opn) then
              last_time = second_of_day(i)
              last_price = price(i)
              i = i + 1
            end if
          end if
        end if
        sum_volume = 0.0_dp
        count = 0
        do while (i <= last)
          if (second_of_day(i) > cls) exit
          this_time = second_of_day(i)
          j = i
          do while (j < last)
            if (abs(second_of_day(j+1) - this_time) > epsilon(1.0_dp)) exit
            j = j + 1
          end do
          emit = .false.
          select case (typ)
          case (DURATION_TRADE)
            count = j - i + 1
            emit = (.not. rmzero) .or. (this_time > last_time)
            if (.not. rmzero) then
              do k = i, j
                call append_duration(k, second_of_day(k) - last_time, 1)
                last_time = second_of_day(k)
              end do
            else if (emit) then
              call append_duration(j, this_time - last_time, count)
              last_time = this_time
            end if
          case (DURATION_PRICE)
            do k = i, j
              count = count + 1
              if (present(volume)) sum_volume = sum_volume + volume(k)
              if (abs(price(k) - last_price) >= pdiff) then
                if ((.not. rmzero) .or. this_time > last_time) then
                  call append_duration(k, this_time - last_time, count, sum_volume)
                  last_time = this_time
                  last_price = price(k)
                  count = 0
                  sum_volume = 0.0_dp
                end if
              end if
            end do
          case (DURATION_VOLUME)
            do k = i, j
              count = count + 1
              sum_volume = sum_volume + volume(k)
              if (sum_volume >= cvol) then
                if ((.not. rmzero) .or. this_time > last_time) then
                  call append_duration(k, this_time - last_time, count, sum_volume)
                  last_time = this_time
                  count = 0
                  sum_volume = 0.0_dp
                end if
              end if
            end do
          end select
          i = j + 1
        end do
      end if
      first = last + 1
    end do

    allocate(result%year(nout), result%month(nout), result%day(nout))
    allocate(result%second_of_day(nout), result%durations(nout))
    allocate(result%price(nout), result%volume(nout), result%transaction_count(nout))
    if (nout > 0) then
      result%year = oy(1:nout); result%month = om(1:nout); result%day = od(1:nout)
      result%second_of_day = ot(1:nout); result%durations = ox(1:nout)
      result%price = op(1:nout); result%volume = ov(1:nout)
      result%transaction_count = on(1:nout)
    end if
    result%status = ACDM_SUCCESS

  contains

    subroutine append_duration(idx, duration, ntrans, aggregate_volume)
      integer, intent(in) :: idx, ntrans
      real(dp), intent(in) :: duration
      real(dp), intent(in), optional :: aggregate_volume
      nout = nout + 1
      oy(nout) = year(idx); om(nout) = month(idx); od(nout) = day(idx)
      ot(nout) = second_of_day(idx); ox(nout) = duration; on(nout) = ntrans
      if (present(price)) op(nout) = price(idx)
      if (present(aggregate_volume)) then
        ov(nout) = aggregate_volume
      else if (present(volume)) then
        ov(nout) = sum(volume(max(1, idx-ntrans+1):idx))
      end if
    end subroutine append_duration

  end subroutine compute_durations

  pure integer(kind=8) function date_key(y, m, d) result(key)
    integer, intent(in) :: y, m, d
    key = int(y, 8) * 10000_8 + int(m, 8) * 100_8 + int(d, 8)
  end function date_key

  subroutine diurnal_adjust(time_seconds, durations, result, method, nodes, &
                            group, fourier_order, smooth_penalty)
    real(dp), intent(in) :: time_seconds(:), durations(:)
    type(diurnal_result), intent(out) :: result
    integer, intent(in), optional :: method, group(:), fourier_order
    real(dp), intent(in), optional :: nodes(:), smooth_penalty

    integer :: n, meth, gmax, g, i, ng, q, st, grid_n
    integer, allocatable :: labels(:), idx(:)
    real(dp), allocatable :: x(:), y(:), fit(:), grid(:), gfit(:)
    real(dp) :: penalty

    result%status = ACDM_BAD_INPUT
    n = size(durations)
    if (n < 4 .or. size(time_seconds) /= n) return
    if (any(durations <= 0.0_dp)) return
    meth = DIURNAL_CUBIC_SPLINE
    if (present(method)) meth = method
    if (meth < DIURNAL_CUBIC_SPLINE .or. meth > DIURNAL_FFF) return
    q = 4
    if (present(fourier_order)) q = max(0, fourier_order)
    penalty = 0.0_dp
    if (meth == DIURNAL_SMOOTH_SPLINE) penalty = 1.0_dp
    if (present(smooth_penalty)) penalty = max(0.0_dp, smooth_penalty)

    allocate(labels(n))
    labels = 1
    if (present(group)) then
      if (size(group) /= n .or. any(group < 1)) return
      labels = group
    end if
    gmax = maxval(labels)
    allocate(result%fitted(n), result%adjusted(n), result%group_labels(gmax))
    result%fitted = 0.0_dp
    result%group_labels = [(g, g=1,gmax)]
    grid_n = 201
    allocate(result%grid_time(grid_n), result%grid_fit(grid_n, gmax))
    do i = 1, grid_n
      result%grid_time(i) = minval(time_seconds) + real(i-1,dp) * &
        (maxval(time_seconds)-minval(time_seconds))/real(grid_n-1,dp)
    end do
    result%grid_fit = 0.0_dp

    do g = 1, gmax
      ng = count(labels == g)
      if (ng < 4) return
      allocate(idx(ng), x(ng), y(ng), fit(ng), grid(grid_n), gfit(grid_n))
      idx = pack([(i, i=1,n)], labels == g)
      x = time_seconds(idx)
      y = durations(idx)
      grid = result%grid_time
      select case (meth)
      case (DIURNAL_CUBIC_SPLINE, DIURNAL_SMOOTH_SPLINE)
        if (present(nodes)) then
          call binned_spline_fit(x, y, nodes, fit, grid, gfit, st, penalty)
        else
          call automatic_spline_fit(x, y, fit, grid, gfit, st, penalty)
        end if
      case (DIURNAL_SUPER_SMOOTHER)
        call super_smoother_fit(x, y, fit, grid, gfit, st)
      case (DIURNAL_FFF)
        call flexible_fourier_fit(x, y, q, fit, grid, gfit, st)
      end select
      if (st /= ACDM_SUCCESS .or. any(fit <= tiny_pos)) return
      result%fitted(idx) = fit
      result%grid_fit(:, g) = gfit
      deallocate(idx, x, y, fit, grid, gfit)
    end do
    result%adjusted = durations / result%fitted
    result%status = ACDM_SUCCESS
  end subroutine diurnal_adjust

  subroutine automatic_spline_fit(x, y, fitted, grid, grid_fit, status, penalty)
    real(dp), intent(in) :: x(:), y(:), grid(:), penalty
    real(dp), intent(out) :: fitted(:), grid_fit(:)
    integer, intent(out) :: status
    real(dp), allocatable :: nodes(:)
    integer :: i, nn

    nn = min(12, max(4, int(sqrt(real(size(x), dp)))))
    allocate(nodes(nn))
    do i = 1, nn
      nodes(i) = minval(x) + real(i-1,dp) * (maxval(x)-minval(x))/real(nn-1,dp)
    end do
    call binned_spline_fit(x, y, nodes, fitted, grid, grid_fit, status, penalty)
  end subroutine automatic_spline_fit

  subroutine binned_spline_fit(x, y, nodes, fitted, grid, grid_fit, status, penalty)
    real(dp), intent(in) :: x(:), y(:), nodes(:), grid(:), penalty
    real(dp), intent(out) :: fitted(:), grid_fit(:)
    integer, intent(out) :: status
    integer :: nb, i, j, k, used, st
    real(dp), allocatable :: bx(:), by(:), y2(:), a(:, :), b(:), sol(:)
    integer, allocatable :: cnt(:)
    real(dp) :: mid

    status = ACDM_BAD_INPUT
    nb = size(nodes)-1
    if (nb < 3 .or. any(nodes(2:) <= nodes(:size(nodes)-1))) return
    if (minval(x) < nodes(1) .or. maxval(x) > nodes(size(nodes))) return
    allocate(bx(nb), by(nb), cnt(nb))
    bx = 0.0_dp; by = 0.0_dp; cnt = 0
    do i = 1, size(x)
      j = nb
      do k = 1, nb
        if (x(i) >= nodes(k) .and. x(i) < nodes(k+1)) then
          j = k; exit
        end if
      end do
      mid = 0.5_dp * (nodes(j)+nodes(j+1))
      bx(j) = mid
      by(j) = by(j) + y(i)
      cnt(j) = cnt(j) + 1
    end do
    used = count(cnt > 0)
    if (used < 4) return
    bx = pack(bx, cnt > 0)
    by = pack(by, cnt > 0) / real(pack(cnt, cnt > 0), dp)
    if (penalty > 0.0_dp .and. used > 2) then
      allocate(a(used,used), b(used), sol(used))
      a = 0.0_dp; b = by
      do i = 1, used
        a(i,i) = 1.0_dp
      end do
      do i = 1, used-2
        a(i,i) = a(i,i) + penalty
        a(i,i+1) = a(i,i+1) - 2.0_dp*penalty
        a(i,i+2) = a(i,i+2) + penalty
        a(i+1,i) = a(i+1,i) - 2.0_dp*penalty
        a(i+1,i+1) = a(i+1,i+1) + 4.0_dp*penalty
        a(i+1,i+2) = a(i+1,i+2) - 2.0_dp*penalty
        a(i+2,i) = a(i+2,i) + penalty
        a(i+2,i+1) = a(i+2,i+1) - 2.0_dp*penalty
        a(i+2,i+2) = a(i+2,i+2) + penalty
      end do
      call solve_linear(a,b,sol,st)
      if (st /= ACDM_SUCCESS) return
      by = sol
    end if
    allocate(y2(used))
    call natural_spline_second(bx, by, y2, st)
    if (st /= ACDM_SUCCESS) return
    do i = 1, size(x)
      fitted(i) = max(tiny_pos, natural_spline_eval(bx,by,y2,x(i)))
    end do
    do i = 1, size(grid)
      grid_fit(i) = max(tiny_pos, natural_spline_eval(bx,by,y2,grid(i)))
    end do
    status = ACDM_SUCCESS
  end subroutine binned_spline_fit

  subroutine flexible_fourier_fit(x, y, q, fitted, grid, grid_fit, status)
    real(dp), intent(in) :: x(:), y(:), grid(:)
    integer, intent(in) :: q
    real(dp), intent(out) :: fitted(:), grid_fit(:)
    integer, intent(out) :: status
    real(dp), allocatable :: design(:, :), beta(:), res(:)
    real(dp) :: xmin, xrange, t
    integer :: i, p

    status = ACDM_BAD_INPUT
    if (size(x) /= size(y) .or. size(fitted) /= size(y) .or. q < 0) return
    xmin = minval(x); xrange = maxval(x)-xmin
    if (xrange <= 0.0_dp) return
    p = 2 + 2*q
    if (size(x) <= p) return
    allocate(design(size(x),p), beta(p), res(size(x)))
    do i = 1, size(x)
      t = (x(i)-xmin)/xrange
      call fourier_row(t,q,design(i,:))
    end do
    call least_squares(design,y,beta,res,status,ridge=1.0e-12_dp)
    if (status /= ACDM_SUCCESS) return
    fitted = matmul(design,beta)
    do i = 1, size(grid)
      t = (grid(i)-xmin)/xrange
      call fourier_value(t,q,beta,grid_fit(i))
    end do
    fitted = max(fitted,tiny_pos); grid_fit = max(grid_fit,tiny_pos)
  end subroutine flexible_fourier_fit

  subroutine fourier_row(t,q,row)
    real(dp), intent(in) :: t
    integer, intent(in) :: q
    real(dp), intent(out) :: row(:)
    integer :: j
    row(1)=1.0_dp; row(2)=t
    do j=1,q
      row(2+j)=cos(2.0_dp*acos(-1.0_dp)*real(j,dp)*t)
      row(2+q+j)=sin(2.0_dp*acos(-1.0_dp)*real(j,dp)*t)
    end do
  end subroutine fourier_row

  subroutine fourier_value(t,q,beta,value)
    real(dp), intent(in) :: t, beta(:)
    integer, intent(in) :: q
    real(dp), intent(out) :: value
    real(dp), allocatable :: row(:)
    allocate(row(size(beta)))
    call fourier_row(t,q,row)
    value=dot_product(row,beta)
  end subroutine fourier_value

  subroutine super_smoother_fit(x,y,fitted,grid,grid_fit,status)
    real(dp), intent(in) :: x(:),y(:),grid(:)
    real(dp), intent(out) :: fitted(:),grid_fit(:)
    integer, intent(out) :: status
    real(dp), parameter :: spans(3)=[0.05_dp,0.20_dp,0.50_dp]
    real(dp), allocatable :: cand(:, :), cv(:, :), local_span(:), smooth_span(:)
    integer :: i,j,k,best,n

    status=ACDM_BAD_INPUT; n=size(x)
    if (size(y)/=n .or. size(fitted)/=n .or. n<8) return
    allocate(cand(n,3),cv(n,3),local_span(n),smooth_span(n))
    do j=1,3
      call running_local_linear(x,y,spans(j),cand(:,j),cv(:,j))
    end do
    do i=1,n
      best=1
      do j=2,3
        if (cv(i,j)<cv(i,best)) best=j
      end do
      local_span(i)=spans(best)
    end do
    call running_mean_x(x,local_span,0.20_dp,smooth_span)
    do i=1,n
      k=1
      if (smooth_span(i)>0.125_dp) k=2
      if (smooth_span(i)>0.35_dp) k=3
      fitted(i)=cand(i,k)
    end do
    do i=1,size(grid)
      call local_linear_at(x,y,grid(i),median_span(smooth_span),grid_fit(i))
    end do
    fitted=max(fitted,tiny_pos); grid_fit=max(grid_fit,tiny_pos)
    status=ACDM_SUCCESS
  end subroutine super_smoother_fit

  subroutine running_local_linear(x,y,span,fit,cv)
    real(dp),intent(in)::x(:),y(:),span
    real(dp),intent(out)::fit(:),cv(:)
    integer::i,n,k
    real(dp)::lev
    n=size(x); k=max(3,ceiling(span*real(n,dp)))
    do i=1,n
      call local_linear_at(x,y,x(i),real(k,dp)/real(n,dp),fit(i),lev)
      cv(i)=((y(i)-fit(i))/max(0.05_dp,1.0_dp-lev))**2
    end do
  end subroutine running_local_linear

  subroutine local_linear_at(x,y,x0,span,value,leverage)
    real(dp),intent(in)::x(:),y(:),x0,span
    real(dp),intent(out)::value
    real(dp),intent(out),optional::leverage
    real(dp),allocatable::dist(:)
    real(dp)::h,u,w,s0,s1,s2,t0,t1,den,lev
    integer::i,n,k
    n=size(x); k=max(3,min(n,ceiling(span*real(n,dp))))
    allocate(dist(n)); dist=abs(x-x0); call sort_real(dist)
    h=max(dist(k),epsilon(1.0_dp)*max(1.0_dp,abs(x0)))
    s0=0.0_dp;s1=0.0_dp;s2=0.0_dp;t0=0.0_dp;t1=0.0_dp
    do i=1,n
      u=abs(x(i)-x0)/h
      if (u<1.0_dp) then
        w=(1.0_dp-u**3)**3
        s0=s0+w; s1=s1+w*(x(i)-x0); s2=s2+w*(x(i)-x0)**2
        t0=t0+w*y(i); t1=t1+w*(x(i)-x0)*y(i)
      end if
    end do
    den=s0*s2-s1*s1
    if (abs(den)<=tiny_pos) then
      value=t0/max(s0,tiny_pos); lev=1.0_dp/max(1.0_dp,s0)
    else
      value=(s2*t0-s1*t1)/den
      lev=s2/den
    end if
    if (present(leverage)) leverage=min(0.99_dp,max(0.0_dp,lev))
  end subroutine local_linear_at

  subroutine running_mean_x(x,y,span,out)
    real(dp),intent(in)::x(:),y(:),span
    real(dp),intent(out)::out(:)
    integer::i
    do i=1,size(x)
      call local_linear_at(x,y,x(i),span,out(i))
    end do
  end subroutine running_mean_x

  function median_span(x) result(m)
    real(dp),intent(in)::x(:)
    real(dp)::m
    real(dp),allocatable::z(:)
    integer::n
    z=x; call sort_real(z); n=size(z)
    if(mod(n,2)==1) then;m=z((n+1)/2);else;m=0.5_dp*(z(n/2)+z(n/2+1));end if
  end function median_span

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::v
    do i=2,size(x)
      v=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=v) exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=v
    end do
  end subroutine sort_real

end module acdm_data
