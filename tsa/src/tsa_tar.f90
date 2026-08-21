! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_tar
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use tsa_kinds, only : dp
  use tsa_types, only : tar_result, tar_multi_result, tsa_test_result
  use tsa_utils, only : sort_with_index, quantile_sorted
  use tsa_simulation, only : tar_sim
  use tseries_special, only : normal_cdf
  use tseries_linalg, only : least_squares
  implicit none
  private

  public :: tar_fit, tar_fit_multi, tar_skeleton, tar_predict, tlrt_test, tlrt_p_value

contains

  subroutine tar_fit(y, p1, p2, d, result, constant1, constant2, &
      estimate_threshold, threshold, a, b, order_select, method, status, &
      transform, center, standard)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: p1, p2, d
    type(tar_result), intent(out) :: result
    logical, intent(in), optional :: constant1, constant2
    logical, intent(in), optional :: estimate_threshold, order_select
    real(dp), intent(in), optional :: threshold, a, b
    character(len=*), intent(in), optional :: method, transform
    logical, intent(in), optional :: center, standard
    integer, intent(out), optional :: status

    logical :: c1, c2, est, os, use_maic
    integer :: p, start, m, maxm, nvalid, i, j, k, lo, hi, bestidx, st
    integer :: idxsplit, p1best, p2best, q1, q2
    real(dp) :: aa, bb, best, criterion, rss1, rss2, aic1, aic2, th
    real(dp), allocatable :: response(:), thvar(:), sorted(:), xfull(:,:), yw(:)
    real(dp) :: ymean, ysd, nanv
    logical :: ctr, stdz
    character(len=5) :: trans
    integer, allocatable :: ord(:), tid(:)

    c1 = .true.
    c2 = .true.
    est = .true.
    os = .true.
    aa = 0.05_dp
    bb = 0.95_dp
    use_maic = .true.
    if (present(constant1)) c1 = constant1
    if (present(constant2)) c2 = constant2
    if (present(estimate_threshold)) est = estimate_threshold
    if (present(order_select)) os = order_select
    if (present(a)) aa = a
    if (present(b)) bb = b
    if (present(method)) use_maic = method(1:1) == 'M' .or. method(1:1) == 'm'

    ctr = .false.
    stdz = .false.
    trans = 'no'
    if (present(center)) ctr = center
    if (present(standard)) stdz = standard
    if (present(transform)) trans = adjustl(transform)
    allocate(yw(size(y)))
    yw = y
    select case (trim(trans))
    case ('no')
      continue
    case ('log')
      do i = 1, size(yw)
        if (.not. ieee_is_finite(yw(i))) cycle
        if (yw(i) <= 0.0_dp) then
          call set_status(6)
          return
        end if
        yw(i) = log(yw(i))
      end do
    case ('log10')
      do i = 1, size(yw)
        if (.not. ieee_is_finite(yw(i))) cycle
        if (yw(i) <= 0.0_dp) then
          call set_status(6)
          return
        end if
        yw(i) = log10(yw(i))
      end do
    case ('sqrt')
      do i = 1, size(yw)
        if (.not. ieee_is_finite(yw(i))) cycle
        if (yw(i) < 0.0_dp) then
          call set_status(6)
          return
        end if
        yw(i) = sqrt(yw(i))
      end do
    case default
      call set_status(6)
      return
    end select
    nvalid = count(ieee_is_finite(yw))
    if (nvalid == 0) then
      call set_status(6)
      return
    end if
    ymean = sum(pack(yw,ieee_is_finite(yw)))/real(nvalid,dp)
    if (nvalid > 1) then
      ysd = sqrt(sum(pack((yw-ymean)**2,ieee_is_finite(yw)))/real(nvalid-1,dp))
    else
      ysd = 1.0_dp
    end if
    if (ysd <= epsilon(1.0_dp)) ysd = 1.0_dp
    do i = 1, size(yw)
      if (.not. ieee_is_finite(yw(i))) cycle
      if (stdz) then
        yw(i) = (yw(i)-ymean)/ysd
      else if (ctr) then
        yw(i) = yw(i)-ymean
      end if
    end do
    result%transform = trans
    result%transform_mean = ymean
    result%transform_sd = ysd
    result%centered = ctr
    result%standardized = stdz

    p = max(p1, p2)
    start = max(p, d) + 1
    result%p1 = p1
    result%p2 = p2
    result%d = d
    result%constant1 = c1
    result%constant2 = c2

    if (p1 < 0 .or. p2 < 0 .or. d < 1 .or. size(yw) < start+5) then
      call set_status(1)
      return
    end if
    if (aa < 0.0_dp .or. bb > 1.0_dp .or. aa >= bb) then
      call set_status(2)
      return
    end if

    maxm = size(yw)-start+1
    allocate(response(maxm), thvar(maxm), xfull(maxm,p), tid(maxm))
    m = 0
    do i = start, size(yw)
      if (.not. ieee_is_finite(yw(i))) cycle
      if (.not. ieee_is_finite(yw(i-d))) cycle
      if (p > 0) then
        if (.not. all(ieee_is_finite(yw(i-p:i-1)))) cycle
      end if
      m = m+1
      tid(m) = i
      response(m) = yw(i)
      thvar(m) = yw(i-d)
      do j = 1, p
        xfull(m,j) = yw(i-j)
      end do
    end do
    if (m < p+6) then
      call set_status(1)
      return
    end if
    response=response(:m)
    thvar=thvar(:m)
    xfull=xfull(:m,:)
    tid=tid(:m)

    allocate(sorted(m), ord(m))
    call sort_with_index(thvar, sorted, ord)
    lo = max(1, ceiling(aa*real(m,dp)))
    hi = min(m-1, floor(bb*real(m,dp)))
    lo = max(lo, max(min_regime_size(p1,c1), min_regime_size(p2,c2)))
    hi = min(hi, m-max(min_regime_size(p1,c1), min_regime_size(p2,c2)))
    if (lo >= hi) then
      call set_status(3)
      return
    end if

    if (est) then
      best = huge(1.0_dp)
      bestidx = lo
      do k = lo, hi
        call regime_score(ord(:k), p1, c1, use_maic .and. os, &
          rss1, aic1, q1, st)
        if (st /= 0) cycle
        call regime_score(ord(k+1:), p2, c2, use_maic .and. os, &
          rss2, aic2, q2, st)
        if (st /= 0) cycle
        if (use_maic) then
          criterion = aic1 + aic2
        else
          criterion = rss1 + rss2
        end if
        if (criterion < best) then
          best = criterion
          bestidx = k
        end if
      end do
      if (best >= huge(1.0_dp)/2.0_dp) then
        call set_status(4)
        return
      end if
      idxsplit = bestidx
      th = sorted(idxsplit)
    else
      if (.not. present(threshold)) then
        call set_status(5)
        return
      end if
      th = threshold
      idxsplit = count(sorted <= th)
      idxsplit = max(lo, min(hi, idxsplit))
    end if

    result%threshold = th
    result%threshold_index = idxsplit
    call final_fit(idxsplit, p1best, p2best, st)
    result%p1 = p1best
    result%p2 = p2best
    call set_status(st)

  contains

    pure integer function min_regime_size(pmax, inc) result(nmin)
      integer, intent(in) :: pmax
      logical, intent(in) :: inc
      nmin = max(3, pmax + merge(1,0,inc) + 2)
    end function min_regime_size

    subroutine make_design(indices, pord, inc, xm, ym)
      integer, intent(in) :: indices(:), pord
      logical, intent(in) :: inc
      real(dp), allocatable, intent(out) :: xm(:,:), ym(:)
      integer :: ii, jj, cols, row, off

      cols = pord + merge(1,0,inc)
      allocate(xm(size(indices),cols), ym(size(indices)))
      off = merge(1,0,inc)
      do ii = 1, size(indices)
        row = indices(ii)
        ym(ii) = response(row)
        if (inc) xm(ii,1) = 1.0_dp
        do jj = 1, pord
          xm(ii,jj+off) = xfull(row,jj)
        end do
      end do
    end subroutine make_design

    subroutine fit_one_order(indices, pord, inc, rss, aic, beta, resid, istat)
      integer, intent(in) :: indices(:), pord
      logical, intent(in) :: inc
      real(dp), intent(out) :: rss, aic
      real(dp), allocatable, intent(out) :: beta(:), resid(:)
      integer, intent(out) :: istat
      real(dp), allocatable :: xm(:,:), ym(:)
      integer :: nobs, npar

      call make_design(indices, pord, inc, xm, ym)
      nobs = size(ym)
      npar = size(xm,2)
      allocate(beta(npar), resid(nobs))
      if (npar == 0) then
        resid = ym
        istat = 0
      else
        call least_squares(xm, ym, beta, residuals=resid, status=istat)
      end if
      if (istat /= 0) then
        rss = huge(1.0_dp)
        aic = huge(1.0_dp)
        return
      end if
      rss = sum(resid**2)
      aic = real(nobs,dp)*log(max(rss/real(nobs,dp),tiny(1.0_dp))) + &
        2.0_dp*real(npar,dp)
    end subroutine fit_one_order

    subroutine regime_score(indices, pmax, inc, select_order, &
        rss, aic, pbest, istat)
      integer, intent(in) :: indices(:), pmax
      logical, intent(in) :: inc, select_order
      real(dp), intent(out) :: rss, aic
      integer, intent(out) :: pbest, istat
      real(dp), allocatable :: beta(:), resid(:)
      real(dp) :: rr, aaic, best_aic
      integer :: po, first, last, st1

      if (select_order) then
        first = 0
        last = pmax
      else
        first = pmax
        last = pmax
      end if
      best_aic = huge(1.0_dp)
      rss = huge(1.0_dp)
      aic = huge(1.0_dp)
      pbest = first
      istat = 1
      do po = first, last
        if (size(indices) <= po + merge(1,0,inc)) cycle
        call fit_one_order(indices, po, inc, rr, aaic, beta, resid, st1)
        if (st1 == 0 .and. aaic < best_aic) then
          best_aic = aaic
          rss = rr
          aic = aaic
          pbest = po
          istat = 0
        end if
      end do
    end subroutine regime_score

    subroutine final_fit(ksplit, po1, po2, status_local)
      integer, intent(in) :: ksplit
      integer, intent(out) :: po1, po2, status_local
      real(dp), allocatable :: b1(:), b2(:), r1(:), r2(:)
      real(dp) :: rsslo, rsshi, aiclo, aichi
      integer :: s1, s2, ii, row, off1, off2

      call regime_score(ord(:ksplit), p1, c1, use_maic .and. os, &
        rsslo, aiclo, po1, s1)
      call regime_score(ord(ksplit+1:), p2, c2, use_maic .and. os, &
        rsshi, aichi, po2, s2)
      status_local = max(s1,s2)
      if (status_local /= 0) return

      call fit_one_order(ord(:ksplit), po1, c1, rsslo, aiclo, b1, r1, s1)
      call fit_one_order(ord(ksplit+1:), po2, c2, rsshi, aichi, b2, r2, s2)
      status_local = max(s1,s2)
      if (status_local /= 0) return

      allocate(result%phi1(po1+1), result%phi2(po2+1))
      result%phi1 = 0.0_dp
      result%phi2 = 0.0_dp
      off1 = merge(1,0,c1)
      off2 = merge(1,0,c2)
      if (c1) result%phi1(1) = b1(1)
      if (c2) result%phi2(1) = b2(1)
      if (po1 > 0) result%phi1(2:) = b1(off1+1:)
      if (po2 > 0) result%phi2(2:) = b2(off2+1:)

      result%n1 = ksplit
      result%n2 = m-ksplit
      result%rms1 = rsslo / &
        real(max(1,ksplit-size(b1)),dp)
      result%rms2 = rsshi / &
        real(max(1,m-ksplit-size(b2)),dp)

      allocate(result%residuals(size(yw)), &
        result%standardized_residuals(size(yw)), result%fitted(size(yw)))
      nanv = ieee_value(0.0_dp,ieee_quiet_nan)
      result%residuals = nanv
      result%standardized_residuals = nanv
      result%fitted = yw
      do ii = 1, ksplit
        row = ord(ii)
        result%residuals(tid(row)) = r1(ii)
        result%standardized_residuals(tid(row)) = &
          r1(ii)/sqrt(max(result%rms1,tiny(1.0_dp)))
        result%fitted(tid(row)) = response(row)-r1(ii)
      end do
      do ii = ksplit+1, m
        row = ord(ii)
        result%residuals(tid(row)) = r2(ii-ksplit)
        result%standardized_residuals(tid(row)) = &
          r2(ii-ksplit)/sqrt(max(result%rms2,tiny(1.0_dp)))
        result%fitted(tid(row)) = response(row)-r2(ii-ksplit)
      end do

      result%aic = real(ksplit,dp)* &
        log(max(rsslo/real(ksplit,dp),tiny(1.0_dp))) + &
        real(m-ksplit,dp)* &
        log(max(rsshi/real(m-ksplit,dp),tiny(1.0_dp))) + &
        real(m,dp)*(1.0_dp+log(2.0_dp*acos(-1.0_dp))) + &
        2.0_dp*real(size(b1)+size(b2)+1,dp)
    end subroutine final_fit

    subroutine set_status(istat)
      integer, intent(in) :: istat
      result%status = istat
      if (present(status)) status = istat
    end subroutine set_status
  end subroutine tar_fit


  subroutine tar_fit_multi(y, p1, p2, d, result, constant1, constant2, &
      estimate_threshold, threshold, a, b, order_select, method, status, &
      transform, center, standard)
    real(dp), intent(in) :: y(:,:)
    integer, intent(in) :: p1, p2, d
    type(tar_multi_result), intent(out) :: result
    logical, intent(in), optional :: constant1, constant2
    logical, intent(in), optional :: estimate_threshold, order_select
    real(dp), intent(in), optional :: threshold, a, b
    character(len=*), intent(in), optional :: method, transform
    logical, intent(in), optional :: center, standard
    integer, intent(out), optional :: status

    real(dp), allocatable :: work(:,:), resp(:), thv(:), lags(:,:), sorted(:)
    real(dp), allocatable :: b1(:), b2(:), r1(:), r2(:)
    integer, allocatable :: sid(:), tid(:), ord(:)
    real(dp) :: aa, bb, crit, best, rr1, rr2, ac1, ac2, th, mu, sdv, nanv
    integer :: n, ns, p, start, maxrows, m, i, j, sidx, k, lo, hi, split
    integer :: st1, st2, k1, k2, npar1, npar2, colcount
    logical :: c1, c2, est, os, maic, ctr, stdz
    character(len=5) :: trans

    c1 = .true.; c2 = .true.; est = .true.; os = .true.; maic = .true.
    ctr = .false.; stdz = .false.; aa = 0.05_dp; bb = 0.95_dp; trans = 'no'
    if (present(constant1)) c1 = constant1
    if (present(constant2)) c2 = constant2
    if (present(estimate_threshold)) est = estimate_threshold
    if (present(order_select)) os = order_select
    if (present(a)) aa = a
    if (present(b)) bb = b
    if (present(center)) ctr = center
    if (present(standard)) stdz = standard
    if (present(method)) maic = method(1:1) == 'M' .or. method(1:1) == 'm'
    if (present(transform)) trans = adjustl(transform)

    n = size(y,1); ns = size(y,2); p = max(p1,p2); start = max(p,d)+1
    result%nseries = ns; result%p1 = p1; result%p2 = p2; result%d = d
    result%constant1 = c1; result%constant2 = c2
    result%centered = ctr; result%standardized = stdz; result%transform = trans
    if (n < start+4 .or. ns < 1 .or. p1 < 0 .or. p2 < 0 .or. d < 1) then
      call set_status_multi(1); return
    end if
    if (aa < 0.0_dp .or. bb > 1.0_dp .or. aa >= bb) then
      call set_status_multi(2); return
    end if

    allocate(work(n,ns), result%transform_mean(ns), result%transform_sd(ns))
    work = y
    do sidx = 1, ns
      select case (trim(trans))
      case ('no')
        continue
      case ('log')
        do i=1,n
          if (ieee_is_finite(work(i,sidx))) then
            if (work(i,sidx) <= 0.0_dp) then
              call set_status_multi(3); return
            end if
            work(i,sidx) = log(work(i,sidx))
          end if
        end do
      case ('log10')
        do i=1,n
          if (ieee_is_finite(work(i,sidx))) then
            if (work(i,sidx) <= 0.0_dp) then
              call set_status_multi(3); return
            end if
            work(i,sidx) = log10(work(i,sidx))
          end if
        end do
      case ('sqrt')
        do i=1,n
          if (ieee_is_finite(work(i,sidx))) then
            if (work(i,sidx) < 0.0_dp) then
              call set_status_multi(3); return
            end if
            work(i,sidx) = sqrt(work(i,sidx))
          end if
        end do
      case default
        call set_status_multi(3); return
      end select
      colcount = count(ieee_is_finite(work(:,sidx)))
      if (colcount == 0) then
        call set_status_multi(4); return
      end if
      mu = sum(pack(work(:,sidx),ieee_is_finite(work(:,sidx))))/real(colcount,dp)
      if (colcount > 1) then
        sdv = sqrt(sum(pack((work(:,sidx)-mu)**2,ieee_is_finite(work(:,sidx))))/ &
          real(colcount-1,dp))
      else
        sdv = 1.0_dp
      end if
      if (sdv <= epsilon(1.0_dp)) sdv = 1.0_dp
      result%transform_mean(sidx) = mu; result%transform_sd(sidx) = sdv
      do i=1,n
        if (.not. ieee_is_finite(work(i,sidx))) cycle
        if (stdz) then
          work(i,sidx) = (work(i,sidx)-mu)/sdv
        else if (ctr) then
          work(i,sidx) = work(i,sidx)-mu
        end if
      end do
    end do

    maxrows = (n-start+1)*ns
    allocate(resp(maxrows),thv(maxrows),lags(maxrows,p),sid(maxrows),tid(maxrows))
    m = 0
    do sidx = 1, ns
      do i = start, n
        if (.not. ieee_is_finite(work(i,sidx))) cycle
        if (.not. ieee_is_finite(work(i-d,sidx))) cycle
        if (p > 0) then
          if (.not. all(ieee_is_finite(work(i-p:i-1,sidx)))) cycle
        end if
        m = m+1; resp(m)=work(i,sidx); thv(m)=work(i-d,sidx)
        sid(m)=sidx; tid(m)=i
        do j=1,p
          lags(m,j)=work(i-j,sidx)
        end do
      end do
    end do
    if (m < 10) then
      call set_status_multi(5); return
    end if
    resp=resp(:m); thv=thv(:m); lags=lags(:m,:); sid=sid(:m); tid=tid(:m)
    allocate(sorted(m),ord(m)); call sort_with_index(thv,sorted,ord)

    npar1 = ns*(p1+merge(1,0,c1)); npar2 = ns*(p2+merge(1,0,c2))
    lo=max(1,ceiling(aa*real(m,dp))); hi=min(m-1,floor(bb*real(m,dp)))
    lo=max(lo,max(npar1+1,2*p1+1)); hi=min(hi,m-max(npar2+1,2*p2+1))
    do while (lo < hi .and. .not. each_series_enough(ord(:lo),p1+2))
      lo=lo+1
    end do
    do while (hi > lo .and. .not. each_series_enough(ord(hi+1:),p2+2))
      hi=hi-1
    end do
    if (lo >= hi) then
      call set_status_multi(6); return
    end if

    if (est) then
      best=huge(1.0_dp); split=lo
      do k=lo,hi
        call regime_fit(ord(:k),p1,c1,maic .and. os,rr1,ac1,k1,b1,r1,st1)
        if (st1/=0) cycle
        call regime_fit(ord(k+1:),p2,c2,maic .and. os,rr2,ac2,k2,b2,r2,st2)
        if (st2/=0) cycle
        if (maic) then; crit=ac1+ac2; else; crit=rr1+rr2; end if
        if (crit < best) then; best=crit; split=k; end if
      end do
      if (best >= huge(1.0_dp)/2.0_dp) then
        call set_status_multi(7); return
      end if
      th=sorted(split)
    else
      if (.not. present(threshold)) then
        call set_status_multi(8); return
      end if
      th=threshold; split=count(sorted<=th); split=max(lo,min(hi,split))
    end if

    call regime_fit(ord(:split),p1,c1,maic .and. os,rr1,ac1,k1,b1,r1,st1)
    call regime_fit(ord(split+1:),p2,c2,maic .and. os,rr2,ac2,k2,b2,r2,st2)
    if (st1/=0 .or. st2/=0) then
      call set_status_multi(9); return
    end if
    result%threshold=th; result%threshold_index=split
    if (maic) then
      result%p1=max(0,k1-merge(1,0,c1)); result%p2=max(0,k2-merge(1,0,c2))
    end if
    result%coefficients1=b1; result%coefficients2=b2
    allocate(result%rms1(ns),result%rms2(ns),result%n1(ns),result%n2(ns), &
      result%aic_series(ns))
    result%rms1=0.0_dp; result%rms2=0.0_dp; result%n1=0; result%n2=0
    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    allocate(result%residuals(n,ns),result%standardized_residuals(n,ns),result%fitted(n,ns))
    result%residuals=nanv; result%standardized_residuals=nanv; result%fitted=work
    do i=1,split
      j=ord(i); sidx=sid(j); result%n1(sidx)=result%n1(sidx)+1
      result%rms1(sidx)=result%rms1(sidx)+r1(i)**2
      result%residuals(tid(j),sidx)=r1(i); result%fitted(tid(j),sidx)=resp(j)-r1(i)
    end do
    do i=split+1,m
      j=ord(i); sidx=sid(j); result%n2(sidx)=result%n2(sidx)+1
      result%rms2(sidx)=result%rms2(sidx)+r2(i-split)**2
      result%residuals(tid(j),sidx)=r2(i-split); result%fitted(tid(j),sidx)=resp(j)-r2(i-split)
    end do
    do sidx=1,ns
      result%rms1(sidx)=result%rms1(sidx)/real(max(1,result%n1(sidx)-result%p1-merge(1,0,c1)),dp)
      result%rms2(sidx)=result%rms2(sidx)/real(max(1,result%n2(sidx)-result%p2-merge(1,0,c2)),dp)
      do i=1,n
        if (ieee_is_finite(result%residuals(i,sidx))) then
          if (work(i-d,sidx) <= th) then
            result%standardized_residuals(i,sidx)=result%residuals(i,sidx)/sqrt(max(result%rms1(sidx),tiny(1.0_dp)))
          else
            result%standardized_residuals(i,sidx)=result%residuals(i,sidx)/sqrt(max(result%rms2(sidx),tiny(1.0_dp)))
          end if
        end if
      end do
      result%aic_series(sidx)=real(result%n1(sidx),dp)*log(max( &
        result%rms1(sidx)*real(max(1,result%n1(sidx)-result%p1-merge(1,0,c1)),dp)/ &
        real(max(1,result%n1(sidx)),dp),tiny(1.0_dp))) + &
        real(result%n2(sidx),dp)*log(max( &
        result%rms2(sidx)*real(max(1,result%n2(sidx)-result%p2-merge(1,0,c2)),dp)/ &
        real(max(1,result%n2(sidx)),dp),tiny(1.0_dp))) + &
        real(result%n1(sidx)+result%n2(sidx),dp)*(1.0_dp+log(2.0_dp*acos(-1.0_dp))) + &
        2.0_dp*real(result%p1+result%p2+merge(1,0,c1)+merge(1,0,c2)+1,dp)
    end do
    result%aic=sum(result%aic_series); call set_status_multi(0)

  contains
    logical function each_series_enough(indices,need) result(ok)
      integer,intent(in)::indices(:),need
      integer::ss
      ok=.true.
      do ss=1,ns
        if(count(sid(indices)==ss)<need) then; ok=.false.; return; end if
      end do
    end function each_series_enough

    subroutine regime_fit(indices,pmax,inc,select_order,rss,aic,kbest,beta,resid,istat)
      integer,intent(in)::indices(:),pmax
      logical,intent(in)::inc,select_order
      real(dp),intent(out)::rss,aic
      integer,intent(out)::kbest,istat
      real(dp),allocatable,intent(out)::beta(:),resid(:)
      real(dp),allocatable::xf(:,:),xx(:,:),bbeta(:),rres(:),yy(:)
      real(dp)::rr,aa,besta
      integer::ncols,kk,first,last,st
      call design_full(indices,pmax,inc,xf,yy)
      ncols=size(xf,2)
      if(ncols==0) then
        allocate(beta(0),resid(size(yy))); resid=yy; rss=sum(yy**2)
        aic=real(size(yy),dp)*log(max(rss/real(size(yy),dp),tiny(1.0_dp))); kbest=0; istat=0; return
      end if
      if(select_order) then; first=1; last=ncols; else; first=ncols; last=ncols; end if
      besta=huge(1.0_dp); istat=1; kbest=first
      do kk=first,last
        allocate(xx(size(yy),kk),bbeta(kk),rres(size(yy))); xx=xf(:,:kk)
        call least_squares(xx,yy,bbeta,residuals=rres,status=st)
        if(st==0) then
          rr=sum(rres**2); aa=real(size(yy),dp)*log(max(rr/real(size(yy),dp),tiny(1.0_dp)))+2.0_dp*real(kk,dp)
          if(aa<besta) then; besta=aa; kbest=kk; rss=rr; aic=aa; beta=bbeta; resid=rres; istat=0; end if
        end if
        deallocate(xx,bbeta,rres)
      end do
    end subroutine regime_fit

    subroutine design_full(indices,pord,inc,xm,ym)
      integer,intent(in)::indices(:),pord
      logical,intent(in)::inc
      real(dp),allocatable,intent(out)::xm(:,:),ym(:)
      integer::block,ncols,ii,jj,row,ss,off
      block=pord+merge(1,0,inc); ncols=block*ns
      allocate(xm(size(indices),ncols),ym(size(indices))); xm=0.0_dp
      do ii=1,size(indices)
        row=indices(ii); ss=sid(row); ym(ii)=resp(row); off=0
        if(inc) then; xm(ii,1)=1.0_dp; off=1; end if
        do jj=1,pord; xm(ii,off+jj)=lags(row,jj); end do
        if(ss>1) then
          off=(ss-1)*block
          if(inc) then; xm(ii,off+1)=1.0_dp; off=off+1; end if
          do jj=1,pord; xm(ii,off+jj)=lags(row,jj); end do
        end if
      end do
    end subroutine design_full

    subroutine set_status_multi(istat)
      integer,intent(in)::istat
      result%status=istat; if(present(status)) status=istat
    end subroutine set_status_multi
  end subroutine tar_fit_multi

  subroutine tar_skeleton(phi1, phi2, threshold, d, n, ntrans, xstart, &
      tail, status, cycle_length)
    real(dp), intent(in) :: phi1(:), phi2(:), threshold
    integer, intent(in) :: d, n, ntrans
    real(dp), intent(in), optional :: xstart(:)
    real(dp), allocatable, intent(out) :: tail(:)
    integer, intent(out) :: status, cycle_length
    real(dp), allocatable :: x(:), startv(:)
    integer :: p, mp, i, j, cl, m

    p = max(size(phi1),size(phi2))-1
    mp = max(p,d)
    allocate(x(n+ntrans), startv(mp))
    startv = 0.0_dp
    if (present(xstart)) then
      startv(1:min(mp,size(xstart))) = xstart(1:min(mp,size(xstart)))
    end if
    x = 0.0_dp
    x(1:mp) = startv
    do i = mp+1, size(x)
      if (x(i-d) <= threshold) then
        x(i) = phi1(1)
        do j = 1, min(p,size(phi1)-1)
          x(i) = x(i) + phi1(j+1)*x(i-j)
        end do
      else
        x(i) = phi2(1)
        do j = 1, min(p,size(phi2)-1)
          x(i) = x(i) + phi2(j+1)*x(i-j)
        end do
      end if
    end do
    m = min(50,size(x))
    allocate(tail(m))
    tail = x(size(x)-m+1:)
    status = 0
    cycle_length = 0
    if (abs(tail(m)) > 1.0e7_dp) then
      status = 2
      return
    end if
    do cl = 1, m/2
      if (all(abs(tail(cl+1:)-tail(:m-cl)) <= 1.0e-4_dp)) then
        cycle_length = cl
        exit
      end if
    end do
  end subroutine tar_skeleton

  subroutine tar_predict(model, n_ahead, n_sim, median, lower, upper, status)
    type(tar_result), intent(in) :: model
    integer, intent(in) :: n_ahead, n_sim
    real(dp), allocatable, intent(out) :: median(:), lower(:), upper(:)
    integer, intent(out) :: status
    real(dp), allocatable :: sim(:,:), path(:), startv(:)
    integer :: i, j, mp, st

    mp = max(model%d,max(model%p1,model%p2))
    allocate(startv(mp))
    startv = 0.0_dp
    if (.not. allocated(model%fitted)) then
      status = 1
      allocate(median(0),lower(0),upper(0))
      return
    end if
    startv = model%fitted(size(model%fitted)-mp+1:) + &
      model%residuals(size(model%residuals)-mp+1:)
    allocate(sim(n_sim,n_ahead))
    do i = 1, n_sim
      call tar_sim(model%phi1, model%phi2, model%threshold, model%d, &
        sqrt(model%rms1), sqrt(model%rms2), n_ahead, 0, path, &
        startv, status=st)
      if (st /= 0) then
        status = st
        return
      end if
      sim(i,:) = path
    end do
    allocate(median(n_ahead), lower(n_ahead), upper(n_ahead))
    do j = 1, n_ahead
      median(j) = quantile_sorted(sim(:,j),0.5_dp)
      lower(j) = quantile_sorted(sim(:,j),0.025_dp)
      upper(j) = quantile_sorted(sim(:,j),0.975_dp)
    end do
    status = 0
  end subroutine tar_predict

  pure real(dp) function normal_pdf(x) result(v)
    real(dp), intent(in) :: x
    v = exp(-0.5_dp*x*x)/sqrt(2.0_dp*acos(-1.0_dp))
  end function normal_pdf

  pure real(dp) function tlrt_p_value(y, a, b, p) result(pv)
    real(dp), intent(in) :: y, a, b
    integer, intent(in) :: p
    real(dp) :: lower, upper, temp, z, df, dens

    lower = invnorm(a)
    upper = invnorm(b)
    if (p == 0) then
      temp = t1(upper)-t1(lower)
      z = sqrt(max(y,tiny(1.0_dp)))
      pv = sqrt(2.0_dp/acos(-1.0_dp))*exp(-y/2.0_dp) * &
        (temp*(z-1.0_dp/z)+1.0_dp/z)
    else
      temp = real(p-1,dp)*(t1(upper)-t1(lower)) + &
        tp1(upper)-tp1(lower) + tp2(upper)-tp2(lower)
      df = real(p+1,dp)
      dens = chi_pdf(y,df)
      pv = 1.0_dp-exp(-2.0_dp*dens*(y/df-1.0_dp)*temp)
    end if
    pv = max(0.0_dp,min(1.0_dp,pv))

  contains
    pure real(dp) function t1(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: f
      f = normal_cdf(x)
      v = 0.5_dp*log(f/(1.0_dp-f))
    end function t1

    pure real(dp) function tp1(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: f, ff, bb, cc, r
      f = normal_cdf(x)
      ff = normal_pdf(x)
      bb = 2.0_dp*f-x*ff
      cc = f*(f-x*ff)-ff*ff
      r = 0.5_dp*(bb+sqrt(max(0.0_dp,bb*bb-4.0_dp*cc)))
      v = 0.5_dp*log(r/(1.0_dp-r))
    end function tp1

    pure real(dp) function tp2(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: f, ff, bb, cc, r
      f = normal_cdf(x)
      ff = normal_pdf(x)
      bb = 2.0_dp*f-x*ff
      cc = f*(f-x*ff)-ff*ff
      r = 0.5_dp*(bb-sqrt(max(0.0_dp,bb*bb-4.0_dp*cc)))
      v = 0.5_dp*log(r/(1.0_dp-r))
    end function tp2

    pure real(dp) function chi_pdf(x,df) result(v)
      real(dp), intent(in) :: x, df
      v = exp((df/2.0_dp-1.0_dp)*log(max(x,tiny(1.0_dp))) - x/2.0_dp - &
        (df/2.0_dp)*log(2.0_dp) - log_gamma(df/2.0_dp))
    end function chi_pdf

    pure real(dp) function invnorm(pp) result(v)
      use tsa_utils, only : normal_quantile
      real(dp), intent(in) :: pp
      v = normal_quantile(pp)
    end function invnorm
  end function tlrt_p_value

  function tlrt_test(y, p, d, a, b) result(res)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: p, d
    real(dp), intent(in), optional :: a, b
    type(tsa_test_result) :: res
    real(dp) :: aa, bb, rss0, rss1
    type(tar_result) :: tr
    real(dp), allocatable :: xmat(:,:), yy(:), beta(:), rr(:)
    integer :: start, n, i, j, st

    aa = 0.25_dp
    bb = 0.75_dp
    if (present(a)) aa = a
    if (present(b)) bb = b
    call tar_fit(y, p, p, d, tr, a=aa, b=bb, order_select=.false., &
      method='CLS', status=st)
    if (st /= 0) then
      res%status = st
      return
    end if
    start = max(p,d)+1
    n = size(y)-start+1
    allocate(xmat(n,p+1), yy(n), beta(p+1), rr(n))
    xmat(:,1) = 1.0_dp
    yy = y(start:)
    do i = 1, n
      do j = 1, p
        xmat(i,j+1) = y(start+i-1-j)
      end do
    end do
    call least_squares(xmat, yy, beta, residuals=rr, status=st)
    if (st /= 0) then
      res%status = st
      return
    end if
    rss0 = sum(rr**2)
    rss1 = tr%rms1*real(max(1,tr%n1-(p+1)),dp) + &
      tr%rms2*real(max(1,tr%n2-(p+1)),dp)
    res%statistic = real(n,dp)*(rss0-rss1)/max(rss1,tiny(1.0_dp))
    res%p_value = tlrt_p_value(res%statistic,aa,bb,p)
    res%order = p
    res%status = 0
  end function tlrt_test
end module tsa_tar
