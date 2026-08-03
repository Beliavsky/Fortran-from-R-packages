! SPDX-License-Identifier: MIT
module jumptest_statistics
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use jumptest_kinds, only : dp, pi
  use jumptest_probability, only : normal_cdf, normal_quantile, chi_square_cdf, &
    beta_cdf, sample_quantile_type7
  use jumptest_status, only : JT_SUCCESS, JT_INVALID_ARGUMENT, JT_INVALID_DIMENSION, &
    JT_NONFINITE_INPUT, JT_DEGENERATE_SAMPLE
  implicit none
  private

  integer, parameter, public :: METHOD_BNS = 1
  integer, parameter, public :: METHOD_AMED = 2
  integer, parameter, public :: METHOD_AMIN = 3

  type, public :: statp_result
    real(dp) :: stat = 0.0_dp
    real(dp) :: pvalue = 1.0_dp
    integer :: status = JT_SUCCESS
  end type statp_result

  type, public :: adjp_result
    real(dp), allocatable :: stat(:)
    real(dp), allocatable :: pvalue(:)
    real(dp), allocatable :: adjp(:)
    integer :: status = JT_SUCCESS
  end type adjp_result

  type, public :: pcombine_result
    real(dp), allocatable :: pvalue(:,:)
    character(len=8), allocatable :: methods(:)
    integer :: status = JT_SUCCESS
  end type pcombine_result

  public :: bns_statistic, amin_statistic, amed_statistic
  public :: jumptestday, jumptestperiod, pcombine, ppool
  public :: bh_adjust

contains

  function bns_statistic(ret, status) result(stat)
    real(dp), intent(in) :: ret(:)
    integer, intent(out), optional :: status
    real(dp) :: stat
    integer :: n, i, stat_local
    real(dp) :: rv, bv, tp, b, a, mu43, triproduct

    stat_local = JT_SUCCESS
    stat = ieee_value(0.0_dp, ieee_quiet_nan)
    n = size(ret)
    if (n < 3) then
      stat_local = JT_INVALID_DIMENSION
    else if (.not. all(ieee_is_finite(ret))) then
      stat_local = JT_NONFINITE_INPUT
    else
      rv = sum(ret*ret)
      bv = 0.0_dp
      do i = 1, n - 1
        bv = bv + abs(ret(i))*abs(ret(i + 1))
      end do
      bv = 0.5_dp*pi*real(n, dp)/real(n - 1, dp)*bv
      mu43 = 2.0_dp**(2.0_dp/3.0_dp)*gamma(7.0_dp/6.0_dp)/sqrt(pi)
      tp = 0.0_dp
      do i = 1, n - 2
        triproduct = abs(ret(i))*abs(ret(i + 1))*abs(ret(i + 2))
        tp = tp + triproduct**(4.0_dp/3.0_dp)
      end do
      tp = mu43**(-3.0_dp)*real(n, dp)**2/real(n - 2, dp)*tp
      if (rv <= tiny(1.0_dp) .or. bv <= tiny(1.0_dp)) then
        stat_local = JT_DEGENERATE_SAMPLE
      else
        b = tp/(bv*bv)
        a = (0.5_dp*pi)**2 + pi - 5.0_dp
        stat = sqrt(real(n, dp))*(rv - bv)/rv/sqrt(a*max(1.0_dp, b))
      end if
    end if
    if (present(status)) status = stat_local
  end function bns_statistic

  function amin_statistic(ret, status) result(stat)
    real(dp), intent(in) :: ret(:)
    integer, intent(out), optional :: status
    real(dp) :: stat
    integer :: n, i, stat_local
    real(dp) :: rv, mrv, mrq, value

    stat_local = JT_SUCCESS
    stat = ieee_value(0.0_dp, ieee_quiet_nan)
    n = size(ret)
    if (n < 2) then
      stat_local = JT_INVALID_DIMENSION
    else if (.not. all(ieee_is_finite(ret))) then
      stat_local = JT_NONFINITE_INPUT
    else
      rv = sum(ret*ret)
      mrv = 0.0_dp
      mrq = 0.0_dp
      do i = 1, n - 1
        value = min(abs(ret(i)), abs(ret(i + 1)))
        mrv = mrv + value*value
        mrq = mrq + value**4
      end do
      mrv = pi/(pi - 2.0_dp)*real(n, dp)/real(n - 1, dp)*mrv
      mrq = pi/(3.0_dp*pi - 8.0_dp)*real(n, dp)**2/real(n - 1, dp)*mrq
      if (rv <= tiny(1.0_dp) .or. mrv <= tiny(1.0_dp)) then
        stat_local = JT_DEGENERATE_SAMPLE
      else
        stat = (1.0_dp - mrv/rv)/sqrt(1.81_dp/real(n, dp)*max(1.0_dp, mrq/(mrv*mrv)))
      end if
    end if
    if (present(status)) status = stat_local
  end function amin_statistic

  function amed_statistic(ret, status) result(stat)
    real(dp), intent(in) :: ret(:)
    integer, intent(out), optional :: status
    real(dp) :: stat
    integer :: n, i, stat_local
    real(dp) :: rv, merv, merq, value

    stat_local = JT_SUCCESS
    stat = ieee_value(0.0_dp, ieee_quiet_nan)
    n = size(ret)
    if (n < 3) then
      stat_local = JT_INVALID_DIMENSION
    else if (.not. all(ieee_is_finite(ret))) then
      stat_local = JT_NONFINITE_INPUT
    else
      rv = sum(ret*ret)
      merv = 0.0_dp
      merq = 0.0_dp
      do i = 1, n - 2
        value = median3(abs(ret(i)), abs(ret(i + 1)), abs(ret(i + 2)))
        merv = merv + value*value
        merq = merq + value**4
      end do
      merv = pi/(6.0_dp - 4.0_dp*sqrt(3.0_dp) + pi)* &
        real(n, dp)/real(n - 2, dp)*merv
      merq = 3.0_dp*pi/(9.0_dp*pi + 72.0_dp - 52.0_dp*sqrt(3.0_dp))* &
        real(n, dp)**2/real(n - 2, dp)*merq
      if (rv <= tiny(1.0_dp) .or. merv <= tiny(1.0_dp)) then
        stat_local = JT_DEGENERATE_SAMPLE
      else
        stat = (1.0_dp - merv/rv)/sqrt(0.96_dp/real(n, dp)* &
          max(1.0_dp, merq/(merv*merv)))
      end if
    end if
    if (present(status)) status = stat_local
  end function amed_statistic

  subroutine jumptestday(ret, result, method)
    real(dp), intent(in) :: ret(:)
    type(statp_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    integer :: code

    code = parse_test_method(method)
    if (code < 1) then
      result%stat = ieee_value(0.0_dp, ieee_quiet_nan)
      result%pvalue = result%stat
      result%status = JT_INVALID_ARGUMENT
      return
    end if
    select case (code)
    case (METHOD_BNS)
      result%stat = bns_statistic(ret, result%status)
    case (METHOD_AMED)
      result%stat = amed_statistic(ret, result%status)
    case default
      result%stat = amin_statistic(ret, result%status)
    end select
    if (result%status == JT_SUCCESS) then
      result%pvalue = 1.0_dp - normal_cdf(result%stat)
    else
      result%pvalue = ieee_value(0.0_dp, ieee_quiet_nan)
    end if
  end subroutine jumptestday

  subroutine jumptestperiod(retmat, result, method)
    real(dp), intent(in) :: retmat(:,:)
    type(adjp_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    type(statp_result) :: one
    integer :: j

    allocate(result%stat(size(retmat, 2)), result%pvalue(size(retmat, 2)), &
      result%adjp(size(retmat, 2)))
    if (size(retmat, 1) < 1 .or. size(retmat, 2) < 1) then
      result%status = JT_INVALID_DIMENSION
      return
    end if
    result%status = JT_SUCCESS
    do j = 1, size(retmat, 2)
      call jumptestday(retmat(:, j), one, method)
      result%stat(j) = one%stat
      result%pvalue(j) = one%pvalue
      if (one%status /= JT_SUCCESS .and. result%status == JT_SUCCESS) result%status = one%status
    end do
    if (result%status == JT_SUCCESS) then
      call bh_adjust(result%pvalue, result%adjp)
    else
      result%adjp = ieee_value(0.0_dp, ieee_quiet_nan)
    end if
  end subroutine jumptestperiod

  subroutine pcombine(retmat, methods, result)
    real(dp), intent(in) :: retmat(:,:)
    character(len=*), intent(in) :: methods(:)
    type(pcombine_result), intent(out) :: result
    type(adjp_result) :: one
    integer :: j

    if (size(methods) < 2 .or. size(retmat, 1) < 1 .or. size(retmat, 2) < 1) then
      result%status = JT_INVALID_DIMENSION
      allocate(result%pvalue(0, 0), result%methods(0))
      return
    end if
    allocate(result%pvalue(size(retmat, 2), size(methods)))
    allocate(result%methods(size(methods)))
    result%status = JT_SUCCESS
    do j = 1, size(methods)
      result%methods(j) = adjustl(methods(j))
      call jumptestperiod(retmat, one, methods(j))
      if (one%status /= JT_SUCCESS) then
        result%status = one%status
        result%pvalue = ieee_value(0.0_dp, ieee_quiet_nan)
        return
      end if
      result%pvalue(:, j) = one%pvalue
    end do
  end subroutine pcombine

  subroutine ppool(pmat, result, method)
    real(dp), intent(in) :: pmat(:,:)
    type(adjp_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    character(len=2) :: code
    real(dp), allocatable :: p(:,:), transformed(:,:), rowmax(:), work(:,:)
    real(dp), allocatable :: corr(:,:)
    real(dp) :: threshold, upper_sum, cc, ff
    integer :: nrow, nmet, i, j, nkeep

    nrow = size(pmat, 1)
    nmet = size(pmat, 2)
    allocate(result%stat(nrow), result%pvalue(nrow), result%adjp(nrow))
    if (nrow < 1 .or. nmet < 1) then
      result%status = JT_INVALID_DIMENSION
      return
    end if
    if (.not. all(ieee_is_finite(pmat))) then
      result%status = JT_NONFINITE_INPUT
      result%stat = ieee_value(0.0_dp, ieee_quiet_nan)
      result%pvalue = result%stat
      result%adjp = result%stat
      return
    end if
    code = parse_pool_method(method)
    if (code == '??') then
      result%status = JT_INVALID_ARGUMENT
      result%stat = ieee_value(0.0_dp, ieee_quiet_nan)
      result%pvalue = result%stat
      result%adjp = result%stat
      return
    end if

    allocate(p(nrow, nmet))
    p = min(1.0_dp - 1.0e-5_dp, max(1.0e-5_dp, pmat))
    select case (code)
    case ('SD')
      allocate(transformed(nrow, nmet), rowmax(nrow))
      transformed = normal_quantile(1.0_dp - p)
      rowmax = maxval(p, dim=2)
      threshold = sample_quantile_type7(rowmax, 0.2_dp)
      nkeep = count(rowmax >= threshold)
      if (nkeep < 2) nkeep = nrow
      allocate(work(nkeep, nmet))
      call selected_rows(transformed, rowmax, threshold, work, nkeep)
      allocate(corr(nmet, nmet))
      call correlation_matrix(work, corr)
      upper_sum = sum_upper(corr)
      do i = 1, nrow
        result%stat(i) = sum(transformed(i, :))/sqrt(max(tiny(1.0_dp), &
          real(nmet, dp) + 2.0_dp*upper_sum))
      end do
      result%pvalue = 1.0_dp - normal_cdf(result%stat)
    case ('FD')
      allocate(transformed(nrow, nmet), rowmax(nrow))
      transformed = -2.0_dp*log(p)
      rowmax = maxval(p, dim=2)
      threshold = sample_quantile_type7(rowmax, 0.2_dp)
      nkeep = count(rowmax >= threshold)
      if (nkeep < 2) nkeep = nrow
      allocate(work(nkeep, nmet))
      call selected_rows(transformed, rowmax, threshold, work, nkeep)
      allocate(corr(nmet, nmet))
      call correlation_matrix(work, corr)
      do j = 1, nmet
        do i = 1, nmet
          corr(i, j) = 3.263_dp*corr(i, j) + 0.71_dp*corr(i, j)**2 + &
            0.027_dp*corr(i, j)**3
        end do
      end do
      upper_sum = sum_upper(corr)
      cc = (2.0_dp*real(nmet, dp) + upper_sum)/(2.0_dp*real(nmet, dp))
      ff = (2.0_dp*real(nmet, dp))**2/ &
        max(tiny(1.0_dp), 2.0_dp*real(nmet, dp) + upper_sum)
      do i = 1, nrow
        result%stat(i) = sum(transformed(i, :))/cc
        result%pvalue(i) = 1.0_dp - chi_square_cdf(result%stat(i), ff)
      end do
    case ('SI')
      allocate(transformed(nrow, nmet))
      transformed = normal_quantile(1.0_dp - p)
      do i = 1, nrow
        result%stat(i) = sum(transformed(i, :))/sqrt(real(nmet, dp))
      end do
      result%pvalue = 1.0_dp - normal_cdf(result%stat)
    case ('FI')
      do i = 1, nrow
        result%stat(i) = -2.0_dp*sum(log(p(i, :)))
        result%pvalue(i) = 1.0_dp - chi_square_cdf(result%stat(i), 2.0_dp*real(nmet, dp))
      end do
    case ('MI')
      do i = 1, nrow
        result%stat(i) = minval(p(i, :))
        result%pvalue(i) = beta_cdf(result%stat(i), 1.0_dp, real(nmet, dp))
      end do
    case default
      do i = 1, nrow
        result%stat(i) = maxval(p(i, :))
        result%pvalue(i) = beta_cdf(result%stat(i), real(nmet, dp), 1.0_dp)
      end do
    end select
    call bh_adjust(result%pvalue, result%adjp)
    result%status = JT_SUCCESS
  end subroutine ppool

  subroutine bh_adjust(pvalue, adjusted)
    real(dp), intent(in) :: pvalue(:)
    real(dp), intent(out) :: adjusted(:)
    real(dp), allocatable :: sorted(:), sorted_adjusted(:)
    integer, allocatable :: order(:)
    integer :: n, i

    n = size(pvalue)
    if (size(adjusted) /= n .or. n < 1) return
    allocate(sorted(n), sorted_adjusted(n), order(n))
    sorted = pvalue
    order = [(i, i = 1, n)]
    call sort_with_order(sorted, order)
    sorted_adjusted(n) = min(1.0_dp, sorted(n))
    do i = n - 1, 1, -1
      sorted_adjusted(i) = min(sorted_adjusted(i + 1), &
        min(1.0_dp, real(n, dp)*sorted(i)/real(i, dp)))
    end do
    do i = 1, n
      adjusted(order(i)) = sorted_adjusted(i)
    end do
  end subroutine bh_adjust

  pure function median3(a, b, c) result(median)
    real(dp), intent(in) :: a, b, c
    real(dp) :: median

    median = a + b + c - min(a, min(b, c)) - max(a, max(b, c))
  end function median3

  function parse_test_method(method) result(code)
    character(len=*), intent(in), optional :: method
    integer :: code
    character(len=:), allocatable :: text

    if (present(method)) then
      text = uppercase(trim(adjustl(method)))
    else
      text = 'BNS'
    end if
    select case (text)
    case ('BNS')
      code = METHOD_BNS
    case ('AMED')
      code = METHOD_AMED
    case ('AMIN')
      code = METHOD_AMIN
    case default
      code = -1
    end select
  end function parse_test_method

  function parse_pool_method(method) result(code)
    character(len=*), intent(in), optional :: method
    character(len=2) :: code
    character(len=:), allocatable :: text

    if (present(method)) then
      text = uppercase(trim(adjustl(method)))
    else
      text = 'SD'
    end if
    select case (text)
    case ('SD', 'FD', 'SI', 'FI', 'MI', 'MA')
      code = text
    case default
      code = '??'
    end select
  end function parse_pool_method

  pure function uppercase(text) result(result)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: result
    integer :: i, value

    result = text
    do i = 1, len(text)
      value = iachar(result(i:i))
      if (value >= iachar('a') .and. value <= iachar('z')) then
        result(i:i) = achar(value - iachar('a') + iachar('A'))
      end if
    end do
  end function uppercase

  subroutine selected_rows(values, selector, threshold, selected, desired)
    real(dp), intent(in) :: values(:,:), selector(:), threshold
    real(dp), intent(out) :: selected(:,:)
    integer, intent(in) :: desired
    integer :: i, k

    if (desired == size(values, 1)) then
      selected = values
      return
    end if
    k = 0
    do i = 1, size(values, 1)
      if (selector(i) >= threshold) then
        k = k + 1
        selected(k, :) = values(i, :)
      end if
    end do
  end subroutine selected_rows

  subroutine correlation_matrix(x, corr)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: corr(:,:)
    real(dp), allocatable :: mean(:), centered(:,:), ss(:)
    real(dp) :: denominator
    integer :: i, j, n

    n = size(x, 1)
    allocate(mean(size(x, 2)), centered(size(x, 1), size(x, 2)), ss(size(x, 2)))
    mean = sum(x, dim=1)/real(max(1, n), dp)
    do j = 1, size(x, 2)
      centered(:, j) = x(:, j) - mean(j)
      ss(j) = sum(centered(:, j)**2)
    end do
    corr = 0.0_dp
    do j = 1, size(x, 2)
      corr(j, j) = 1.0_dp
      do i = 1, j - 1
        denominator = sqrt(ss(i)*ss(j))
        if (denominator > tiny(1.0_dp)) then
          corr(i, j) = sum(centered(:, i)*centered(:, j))/denominator
          corr(j, i) = corr(i, j)
        end if
      end do
    end do
  end subroutine correlation_matrix

  pure function sum_upper(matrix) result(value)
    real(dp), intent(in) :: matrix(:,:)
    real(dp) :: value
    integer :: i, j

    value = 0.0_dp
    do j = 2, size(matrix, 2)
      do i = 1, j - 1
        value = value + matrix(i, j)
      end do
    end do
  end function sum_upper

  subroutine sort_with_order(values, order)
    real(dp), intent(inout) :: values(:)
    integer, intent(inout) :: order(:)
    real(dp) :: key
    integer :: key_order, i, j

    do i = 2, size(values)
      key = values(i)
      key_order = order(i)
      j = i - 1
      do while (j >= 1)
        if (values(j) <= key) exit
        values(j + 1) = values(j)
        order(j + 1) = order(j)
        j = j - 1
      end do
      values(j + 1) = key
      order(j + 1) = key_order
    end do
  end subroutine sort_with_order

end module jumptest_statistics
