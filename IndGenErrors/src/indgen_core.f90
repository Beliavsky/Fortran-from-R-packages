! SPDX-License-Identifier: GPL-3.0-only
module indgen_core
  use indgen_kinds, only : dp
  use indgen_types
  use indgen_special, only : chi_square_survival
  use indgen_cvm_tables, only : empirical_cvm_cdf
  use indgen_moebius, only : moebius_full_stat
  implicit none
  private

  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp

  public :: cvm_2series_core, cvm_3series_core
  public :: crosscor_2series_core, crosscor_3series_core
  public :: crossdep_2series_core, crossdep_3series_core

contains

  pure integer function circular_index(i,lag,n) result(j)
    integer, intent(in) :: i, lag, n
    j = modulo(i-1+lag,n)+1
  end function circular_index

  subroutine set_pair_lags(lagmax,lags)
    integer, intent(in) :: lagmax
    integer, allocatable, intent(out) :: lags(:,:)
    integer :: k
    allocate(lags(2,2*lagmax+1))
    do k = -lagmax, lagmax
      lags(1,k+lagmax+1) = k
      lags(2,k+lagmax+1) = 0
    end do
  end subroutine set_pair_lags

  subroutine set_triple_lags(lagmax,lags)
    integer, intent(in) :: lagmax
    integer, allocatable, intent(out) :: lags(:,:)
    integer :: i, j, k, m
    m = (2*lagmax+1)**2
    allocate(lags(2,m))
    k = 0
    do i = -lagmax, lagmax
      do j = -lagmax, lagmax
        k = k+1
        lags(:,k) = [i,j]
      end do
    end do
  end subroutine set_triple_lags

  subroutine rank_maximum(x,r)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: r(:)
    integer :: i
    do i = 1, size(x)
      r(i) = real(count(x <= x(i)),dp)
    end do
  end subroutine rank_maximum

  function cvm_pair_raw(r1,r2,lag) result(value)
    real(dp), intent(in) :: r1(:), r2(:)
    integer, intent(in) :: lag
    real(dp) :: value
    real(dp) :: de, cte, a, b
    integer :: i, j, n, ii, jj

    n = size(r1)
    de = 2.0_dp*real(n,dp)*real(n+1,dp)
    cte = real(2*n+1,dp)/(6.0_dp*real(n,dp))
    value = 0.0_dp
    do i = 1, n
      ii = circular_index(i,lag,n)
      do j = 1, n
        jj = circular_index(j,lag,n)
        a = r1(i)*(r1(i)-1.0_dp)/de+r1(j)*(r1(j)-1.0_dp)/de+cte- &
          max(r1(i),r1(j))/real(n+1,dp)
        b = r2(ii)*(r2(ii)-1.0_dp)/de+r2(jj)*(r2(jj)-1.0_dp)/de+cte- &
          max(r2(ii),r2(jj))/real(n+1,dp)
        value = value+a*b
      end do
    end do
    value = value/real(n,dp)
  end function cvm_pair_raw

  function cvm_triple_raw(r1,r2,r3,lag2,lag3) result(value)
    real(dp), intent(in) :: r1(:), r2(:), r3(:)
    integer, intent(in) :: lag2, lag3
    real(dp) :: value
    real(dp) :: de, cte, a, b, c
    integer :: i, j, n, i2, j2, i3, j3

    n = size(r1)
    de = 2.0_dp*real(n,dp)*real(n+1,dp)
    cte = real(2*n+1,dp)/(6.0_dp*real(n,dp))
    value = 0.0_dp
    do i = 1, n
      i2 = circular_index(i,lag2,n)
      i3 = circular_index(i,lag3,n)
      do j = 1, n
        j2 = circular_index(j,lag2,n)
        j3 = circular_index(j,lag3,n)
        a = r1(i)*(r1(i)-1.0_dp)/de+r1(j)*(r1(j)-1.0_dp)/de+cte- &
          max(r1(i),r1(j))/real(n+1,dp)
        b = r2(i2)*(r2(i2)-1.0_dp)/de+r2(j2)*(r2(j2)-1.0_dp)/de+cte- &
          max(r2(i2),r2(j2))/real(n+1,dp)
        c = r3(i3)*(r3(i3)-1.0_dp)/de+r3(j3)*(r3(j3)-1.0_dp)/de+cte- &
          max(r3(i3),r3(j3))/real(n+1,dp)
        value = value+a*b*c
      end do
    end do
    value = value/real(n,dp)
  end function cvm_triple_raw

  pure function bias_tdn(d,n) result(value)
    integer, intent(in) :: d, n
    real(dp) :: value
    value = (real(n-1,dp)/real(6*n,dp))**d+ &
      real(n-1,dp)/(-real(6*n,dp))**d-(1.0_dp/6.0_dp)**d
  end function bias_tdn

  pure function kappa_xik(k,j) result(value)
    integer, intent(in) :: k, j
    real(dp) :: value
    real(dp), parameter :: a(6) = [1.0_dp,2.0_dp,8.0_dp,48.0_dp,384.0_dp,3840.0_dp]
    real(dp), parameter :: zeta(6) = [ &
      1.0_dp/6.0_dp,1.0_dp/90.0_dp,1.0_dp/945.0_dp, &
      1.0_dp/9450.0_dp,1.0_dp/93555.0_dp,691.0_dp/638512875.0_dp ]
    value = a(j)*zeta(j)**k
  end function kappa_xik

  pure function hermite_prob(n,x) result(value)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: value, h0, h1, h2
    integer :: k

    if (n == 0) then
      value = 1.0_dp
    else if (n == 1) then
      value = x
    else
      h0 = 1.0_dp
      h1 = x
      do k = 2, n
        h2 = x*h1-real(k-1,dp)*h0
        h0 = h1
        h1 = h2
      end do
      value = h1
    end if
  end function hermite_prob

  pure function phi_approx(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value, ax, t, y
    integer :: sgn
    real(dp), parameter :: p = 0.3275911_dp
    real(dp), parameter :: a1 = 0.254829592_dp
    real(dp), parameter :: a2 = -0.284496736_dp
    real(dp), parameter :: a3 = 1.421413741_dp
    real(dp), parameter :: a4 = -1.453152027_dp
    real(dp), parameter :: a5 = 1.061405429_dp

    sgn = 1
    if (x < 0.0_dp) sgn = -1
    ax = abs(x)/sqrt(2.0_dp)
    t = 1.0_dp/(1.0_dp+p*ax)
    y = 1.0_dp-(((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t)*exp(-ax*ax)
    value = 0.5_dp*(1.0_dp+real(sgn,dp)*y)
  end function phi_approx

  pure function edgeworth_cdf(y,kappa) result(value)
    real(dp), intent(in) :: y, kappa(6)
    real(dp) :: value, mu, sigma, z, gam(4), p1, p11, p2, p3, p4
    integer :: j

    mu = kappa(1)
    sigma = sqrt(kappa(2))
    z = (y-mu)/sigma
    do j = 1, 4
      gam(j) = kappa(j+2)/sigma**(j+2)
    end do
    p1 = gam(1)*hermite_prob(2,z)/6.0_dp+gam(2)*hermite_prob(3,z)/24.0_dp+ &
      gam(1)**2*hermite_prob(5,z)/72.0_dp
    p11 = gam(3)*hermite_prob(4,z)/120.0_dp+ &
      gam(1)*gam(2)*hermite_prob(6,z)/144.0_dp
    p2 = gam(1)**3*hermite_prob(8,z)/1296.0_dp+ &
      gam(4)*hermite_prob(5,z)/720.0_dp
    p3 = gam(2)**2*hermite_prob(7,z)/1152.0_dp+ &
      gam(1)*gam(3)*hermite_prob(7,z)/720.0_dp
    p4 = gam(1)**2*gam(2)*hermite_prob(9,z)/1728.0_dp+ &
      gam(1)**4*hermite_prob(11,z)/31104.0_dp
    value = phi_approx(z)-(p1+p11+p2+p3+p4)*exp(-0.5_dp*z*z)/sqrt(2.0_dp*pi)
  end function edgeworth_cdf

  function cvm_2series_core(x,y,lagmax) result(out)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: lagmax
    type(cvm_test_result) :: out
    real(dp), allocatable :: r1(:), r2(:)
    real(dp) :: raw, p
    real(dp) :: cum_w(6), cum_f(6)
    integer :: n, lag, k, j, m

    n = size(x)
    if (n < 2 .or. size(y) /= n .or. lagmax < 0) then
      out%status = indgen_invalid_argument
      return
    end if
    m = 2*lagmax+1
    allocate(out%cvm(m),out%p_cvm(m),r1(n),r2(n))
    call set_pair_lags(lagmax,out%lags)
    call rank_maximum(x,r1)
    call rank_maximum(y,r2)
    out%n = n
    out%wstat = 0.0_dp
    out%fstat = 0.0_dp
    k = 0
    do lag = -lagmax, lagmax
      k = k+1
      raw = cvm_pair_raw(r1,r2,lag)
      out%cvm(k) = 90.0_dp*(raw-1.0_dp/36.0_dp)
      p = 1.0_dp-empirical_cvm_cdf(n,out%cvm(k),2)
      out%p_cvm(k) = p
      out%wstat = out%wstat+raw
      out%fstat = out%fstat-2.0_dp*log(p)
    end do
    out%wstat = out%wstat-real(m,dp)*bias_tdn(2,n)
    do j = 1, 6
      cum_w(j) = real(m,dp)*kappa_xik(2,j)
      cum_f(j) = real(2*m,dp)*kappa_xik(0,j)
    end do
    out%p_wstat = 1.0_dp-edgeworth_cdf(out%wstat,cum_w)
    out%p_fstat = 1.0_dp-edgeworth_cdf(out%fstat,cum_f)
  end function cvm_2series_core

  function cvm_3series_core(x,y,z,lag2,lag3) result(out)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: lag2, lag3
    type(cvm_three_result) :: out
    real(dp), allocatable :: r1(:), r2(:), r3(:)
    real(dp) :: raw, p
    real(dp) :: cum_w12(6), cum_w123(6), cum_w(6)
    real(dp) :: cum_f12(6), cum_f123(6), cum_f(6)
    integer :: n, i, j, k, m1, m2

    n = size(x)
    if (n < 2 .or. size(y) /= n .or. size(z) /= n .or. lag2 < 0 .or. lag3 < 0) then
      out%status = indgen_invalid_argument
      return
    end if
    out%xy = cvm_2series_core(x,y,lag2)
    out%xz = cvm_2series_core(x,z,lag2)
    out%yz = cvm_2series_core(y,z,lag2)
    m1 = 2*lag2+1
    m2 = (2*lag3+1)**2
    allocate(out%xyz%cvm(m2),out%xyz%p_cvm(m2),r1(n),r2(n),r3(n))
    call set_triple_lags(lag3,out%xyz%lags)
    call rank_maximum(x,r1)
    call rank_maximum(y,r2)
    call rank_maximum(z,r3)
    out%xyz%n = n
    out%xyz%wstat = 0.0_dp
    out%xyz%fstat = 0.0_dp
    k = 0
    do i = -lag3, lag3
      do j = -lag3, lag3
        k = k+1
        raw = cvm_triple_raw(r1,r2,r3,i,j)
        out%xyz%cvm(k) = 90.0_dp*sqrt(90.0_dp)*(raw-1.0_dp/216.0_dp)
        p = 1.0_dp-empirical_cvm_cdf(n,out%xyz%cvm(k),3)
        out%xyz%p_cvm(k) = p
        out%xyz%wstat = out%xyz%wstat+raw
        out%xyz%fstat = out%xyz%fstat-2.0_dp*log(p)
      end do
    end do
    out%xyz%wstat = pi*pi*(out%xyz%wstat-real(m2,dp)*bias_tdn(3,n))
    do j = 1, 6
      cum_w12(j) = real(m1,dp)*kappa_xik(2,j)
      cum_w123(j) = real(m2,dp)*kappa_xik(3,j)*(pi*pi)**j
      cum_w(j) = 3.0_dp*cum_w12(j)+cum_w123(j)
      cum_f12(j) = real(2*m1,dp)*kappa_xik(0,j)
      cum_f123(j) = real(2*m2,dp)*kappa_xik(0,j)
      cum_f(j) = 3.0_dp*cum_f12(j)+cum_f123(j)
    end do
    out%xyz%p_wstat = 1.0_dp-edgeworth_cdf(out%xyz%wstat,cum_w123)
    out%xyz%p_fstat = 1.0_dp-edgeworth_cdf(out%xyz%fstat,cum_f123)
    out%wstat = out%xy%wstat+out%xz%wstat+out%yz%wstat+out%xyz%wstat
    out%fstat = out%xy%fstat+out%xz%fstat+out%yz%fstat+out%xyz%fstat
    out%p_wstat = 1.0_dp-edgeworth_cdf(out%wstat,cum_w)
    out%p_fstat = 1.0_dp-edgeworth_cdf(out%fstat,cum_f)
  end function cvm_3series_core

  function crosscor_pair_value(x,y,lag,status) result(value)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: lag
    integer, intent(out) :: status
    real(dp) :: value, mx, my, sx, sy, den
    integer :: i, n

    n = size(x)
    mx = sum(x)/real(n,dp)
    my = sum(y)/real(n,dp)
    sx = sum((x-mx)**2)
    sy = sum((y-my)**2)
    den = sqrt(sx*sy)
    if (den <= tiny(1.0_dp)) then
      value = 0.0_dp
      status = indgen_numerical_error
      return
    end if
    value = 0.0_dp
    do i = 1, n
      value = value+(x(i)-mx)*(y(circular_index(i,lag,n))-my)
    end do
    value = value/den
    status = indgen_success
  end function crosscor_pair_value

  function crosscor_triple_value(x,y,z,lag2,lag3,status) result(value)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: lag2, lag3
    integer, intent(out) :: status
    real(dp) :: value, mx, my, mz, vx, vy, vz, den
    integer :: i, n

    n = size(x)
    mx = sum(x)/real(n,dp)
    my = sum(y)/real(n,dp)
    mz = sum(z)/real(n,dp)
    vx = sum((x-mx)**2)/real(n,dp)
    vy = sum((y-my)**2)/real(n,dp)
    vz = sum((z-mz)**2)/real(n,dp)
    den = sqrt(vx*vy*vz)
    if (den <= tiny(1.0_dp)) then
      value = 0.0_dp
      status = indgen_numerical_error
      return
    end if
    value = 0.0_dp
    do i = 1, n
      value = value+(x(i)-mx)*(y(circular_index(i,lag2,n))-my)* &
        (z(circular_index(i,lag3,n))-mz)
    end do
    value = value/real(n,dp)/den
    status = indgen_success
  end function crosscor_triple_value

  function crosscor_2series_core(x,y,lagmax) result(out)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: lagmax
    type(lag_test_result) :: out
    integer :: n, lag, k, status

    n = size(x)
    if (n < 2 .or. size(y) /= n .or. lagmax < 0) then
      out%status = indgen_invalid_argument
      return
    end if
    allocate(out%stat(2*lagmax+1))
    call set_pair_lags(lagmax,out%lags)
    out%n = n
    out%aggregate = 0.0_dp
    k = 0
    do lag = -lagmax, lagmax
      k = k+1
      out%stat(k) = crosscor_pair_value(x,y,lag,status)
      if (status /= indgen_success) out%status = status
      out%aggregate = out%aggregate+out%stat(k)**2
    end do
    out%aggregate = real(n,dp)*out%aggregate
    out%p_aggregate = chi_square_survival(out%aggregate,real(2*lagmax+1,dp))
  end function crosscor_2series_core

  function crosscor_3series_core(x,y,z,lag2,lag3) result(out)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: lag2, lag3
    type(four_lag_test_result) :: out
    integer :: n, i, j, k, status, m2

    n = size(x)
    if (n < 2 .or. size(y) /= n .or. size(z) /= n .or. lag2 < 0 .or. lag3 < 0) then
      out%status = indgen_invalid_argument
      return
    end if
    out%xy = crosscor_2series_core(x,y,lag2)
    out%xz = crosscor_2series_core(x,z,lag2)
    out%yz = crosscor_2series_core(y,z,lag2)
    m2 = (2*lag3+1)**2
    allocate(out%xyz%stat(m2))
    call set_triple_lags(lag3,out%xyz%lags)
    out%xyz%n = n
    out%xyz%aggregate = 0.0_dp
    k = 0
    do i = -lag3, lag3
      do j = -lag3, lag3
        k = k+1
        out%xyz%stat(k) = crosscor_triple_value(x,y,z,i,j,status)
        if (status /= indgen_success) out%xyz%status = status
        out%xyz%aggregate = out%xyz%aggregate+out%xyz%stat(k)**2
      end do
    end do
    out%xyz%aggregate = real(n,dp)*out%xyz%aggregate
    out%xyz%p_aggregate = chi_square_survival(out%xyz%aggregate,real(m2,dp))
    out%aggregate = out%xy%aggregate+out%xz%aggregate+out%yz%aggregate+out%xyz%aggregate
    out%p_aggregate = chi_square_survival(out%aggregate,real(3*(2*lag2+1)+m2,dp))
  end function crosscor_3series_core

  subroutine dependence_at_lag(x,y,lag,s,g,e,status)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: lag
    real(dp), intent(out) :: s, g, e
    integer, intent(out) :: status
    real(dp), allocatable :: work(:,:)
    integer :: i, n

    n = size(x)
    allocate(work(n,2))
    work(:,1) = x
    do i = 1, n
      work(i,2) = y(circular_index(i,lag,n))
    end do
    call moebius_full_stat(work,s,g,e,status)
  end subroutine dependence_at_lag

  subroutine dependence_at_two_lags(x,y,z,lag2,lag3,s,g,e,status)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: lag2, lag3
    real(dp), intent(out) :: s, g, e
    integer, intent(out) :: status
    real(dp), allocatable :: work(:,:)
    integer :: i, n

    n = size(x)
    allocate(work(n,3))
    work(:,1) = x
    do i = 1, n
      work(i,2) = y(circular_index(i,lag2,n))
      work(i,3) = z(circular_index(i,lag3,n))
    end do
    call moebius_full_stat(work,s,g,e,status)
  end subroutine dependence_at_two_lags

  subroutine finish_lag_test(out,n,df)
    type(lag_test_result), intent(inout) :: out
    integer, intent(in) :: n, df
    out%n = n
    out%aggregate = real(n,dp)*sum(out%stat**2)
    out%p_aggregate = chi_square_survival(out%aggregate,real(df,dp))
  end subroutine finish_lag_test

  function crossdep_2series_core(x,y,lagmax) result(out)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in) :: lagmax
    type(dependence_two_result) :: out
    integer :: n, lag, k, status, m

    n = size(x)
    if (n < 2 .or. size(y) /= n .or. lagmax < 0) then
      out%status = indgen_invalid_argument
      return
    end if
    m = 2*lagmax+1
    allocate(out%spearman%stat(m),out%vdw%stat(m),out%savage%stat(m))
    call set_pair_lags(lagmax,out%spearman%lags)
    call set_pair_lags(lagmax,out%vdw%lags)
    call set_pair_lags(lagmax,out%savage%lags)
    k = 0
    do lag = -lagmax, lagmax
      k = k+1
      call dependence_at_lag(x,y,lag,out%spearman%stat(k),out%vdw%stat(k), &
        out%savage%stat(k),status)
      if (status /= indgen_success) out%status = status
    end do
    call finish_lag_test(out%spearman,n,m)
    call finish_lag_test(out%vdw,n,m)
    call finish_lag_test(out%savage,n,m)
  end function crossdep_2series_core

  subroutine assign_pair_measure(target,source,which)
    type(lag_test_result), intent(out) :: target
    type(dependence_two_result), intent(in) :: source
    integer, intent(in) :: which
    select case (which)
    case (1)
      target = source%spearman
    case (2)
      target = source%vdw
    case (3)
      target = source%savage
    end select
  end subroutine assign_pair_measure

  subroutine finish_four_result(out,n,df)
    type(four_lag_test_result), intent(inout) :: out
    integer, intent(in) :: n, df
    out%aggregate = out%xy%aggregate+out%xz%aggregate+out%yz%aggregate+out%xyz%aggregate
    out%p_aggregate = chi_square_survival(out%aggregate,real(df,dp))
    out%xyz%n = n
  end subroutine finish_four_result

  function crossdep_3series_core(x,y,z,lag2,lag3) result(out)
    real(dp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: lag2, lag3
    type(dependence_three_result) :: out
    type(dependence_two_result) :: xy, xz, yz
    integer :: n, i, j, k, status, m2, df

    n = size(x)
    if (n < 2 .or. size(y) /= n .or. size(z) /= n .or. lag2 < 0 .or. lag3 < 0) then
      out%status = indgen_invalid_argument
      return
    end if
    xy = crossdep_2series_core(x,y,lag2)
    xz = crossdep_2series_core(x,z,lag2)
    yz = crossdep_2series_core(y,z,lag2)
    call assign_pair_measure(out%spearman%xy,xy,1)
    call assign_pair_measure(out%spearman%xz,xz,1)
    call assign_pair_measure(out%spearman%yz,yz,1)
    call assign_pair_measure(out%vdw%xy,xy,2)
    call assign_pair_measure(out%vdw%xz,xz,2)
    call assign_pair_measure(out%vdw%yz,yz,2)
    call assign_pair_measure(out%savage%xy,xy,3)
    call assign_pair_measure(out%savage%xz,xz,3)
    call assign_pair_measure(out%savage%yz,yz,3)

    m2 = (2*lag3+1)**2
    allocate(out%spearman%xyz%stat(m2),out%vdw%xyz%stat(m2), &
      out%savage%xyz%stat(m2))
    call set_triple_lags(lag3,out%spearman%xyz%lags)
    call set_triple_lags(lag3,out%vdw%xyz%lags)
    call set_triple_lags(lag3,out%savage%xyz%lags)
    k = 0
    do i = -lag3, lag3
      do j = -lag3, lag3
        k = k+1
        call dependence_at_two_lags(x,y,z,i,j,out%spearman%xyz%stat(k), &
          out%vdw%xyz%stat(k),out%savage%xyz%stat(k),status)
        if (status /= indgen_success) out%status = status
      end do
    end do
    call finish_lag_test(out%spearman%xyz,n,m2)
    call finish_lag_test(out%vdw%xyz,n,m2)
    call finish_lag_test(out%savage%xyz,n,m2)
    df = 3*(2*lag2+1)+m2
    call finish_four_result(out%spearman,n,df)
    call finish_four_result(out%vdw,n,df)
    call finish_four_result(out%savage,n,df)
  end function crossdep_3series_core

end module indgen_core
