! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
module creditr_curve
  use creditr_kinds, only: dp, creditr_ok, creditr_invalid_input, creditr_no_bracket, creditr_io_error
  use creditr_dates, only: date_t, conventions_t, add_months, modified_following, year_fraction, parse_date
  implicit none
  private

  type, public :: rate_quote_t
    integer :: expiry_months = 0
    character(len=1) :: instrument_type = 'M'
    real(kind=dp) :: rate = 0.0_dp
  end type rate_quote_t

  type, public :: zero_curve_t
    type(date_t) :: base_date
    integer, allocatable :: node_serial(:)
    real(kind=dp), allocatable :: log_discount(:)
  contains
    procedure :: discount_serial
    procedure :: discount_date
    generic :: discount => discount_serial, discount_date
    procedure :: forward_discount
  end type zero_curve_t

  public :: build_zero_curve, read_rate_quotes_csv

contains

  pure real(kind=dp) function discount_serial(self, serial) result(df)
    class(zero_curve_t), intent(in) :: self
    integer, intent(in) :: serial
    integer :: n, lo, hi, mid
    real(kind=dp) :: w, slope
    n = size(self%node_serial)
    if (serial <= self%base_date%serial()) then
      df = 1.0_dp
    else if (serial <= self%node_serial(1)) then
      w = real(serial - self%base_date%serial(), dp) / &
        real(self%node_serial(1) - self%base_date%serial(), dp)
      df = exp(w * self%log_discount(1))
    else if (serial >= self%node_serial(n)) then
      if (n == 1) then
        slope = self%log_discount(1) / real(self%node_serial(1) - self%base_date%serial(), dp)
      else
        slope = (self%log_discount(n) - self%log_discount(n - 1)) / &
          real(self%node_serial(n) - self%node_serial(n - 1), dp)
      end if
      df = exp(self%log_discount(n) + slope * real(serial - self%node_serial(n), dp))
    else
      lo = 1
      hi = n
      do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (self%node_serial(mid) <= serial) then
          lo = mid
        else
          hi = mid
        end if
      end do
      w = real(serial - self%node_serial(lo), dp) / real(self%node_serial(hi) - self%node_serial(lo), dp)
      df = exp((1.0_dp - w) * self%log_discount(lo) + w * self%log_discount(hi))
    end if
  end function discount_serial

  pure real(kind=dp) function discount_date(self, dt) result(df)
    class(zero_curve_t), intent(in) :: self
    type(date_t), intent(in) :: dt
    df = self%discount_serial(dt%serial())
  end function discount_date

  pure real(kind=dp) function forward_discount(self, start_date, end_date) result(df)
    class(zero_curve_t), intent(in) :: self
    type(date_t), intent(in) :: start_date, end_date
    df = self%discount(end_date) / self%discount(start_date)
  end function forward_discount

  subroutine build_zero_curve(base_date, quotes, conventions, curve, status)
    type(date_t), intent(in) :: base_date
    type(rate_quote_t), intent(in) :: quotes(:)
    type(conventions_t), intent(in) :: conventions
    type(zero_curve_t), intent(out) :: curve
    integer, intent(out), optional :: status
    integer, allocatable :: nodes(:)
    real(kind=dp), allocatable :: logs(:)
    integer :: i, n, max_nodes, mat_serial, local_status
    type(date_t) :: maturity
    real(kind=dp) :: yf, df

    if (size(quotes) < 1) then
      if (present(status)) status = creditr_invalid_input
      return
    end if
    max_nodes = 4 + sum(max(1, quotes%expiry_months / max(1, conventions%fixed_frequency_months)))
    allocate(nodes(max_nodes), logs(max_nodes))
    n = 0
    curve%base_date = base_date

    do i = 1, size(quotes)
      if (quotes(i)%expiry_months <= 0 .or. quotes(i)%rate <= -1.0_dp) then
        if (present(status)) status = creditr_invalid_input
        return
      end if
      if (quotes(i)%instrument_type == 'M' .or. quotes(i)%instrument_type == 'm') then
        maturity = add_months(base_date, quotes(i)%expiry_months)
        yf = year_fraction(base_date, maturity, conventions%mm_dcc)
        df = 1.0_dp / (1.0_dp + quotes(i)%rate * yf)
        call insert_node(nodes, logs, n, maturity%serial(), log(df))
      else
        maturity = modified_following(add_months(base_date, quotes(i)%expiry_months))
        mat_serial = maturity%serial()
        call bootstrap_swap(base_date, maturity, quotes(i)%rate, conventions, nodes, logs, n, local_status)
        if (local_status /= creditr_ok) then
          if (present(status)) status = local_status
          return
        end if
        if (nodes(n) /= mat_serial) then
          if (present(status)) status = creditr_invalid_input
          return
        end if
      end if
    end do

    allocate(curve%node_serial(n), curve%log_discount(n))
    curve%node_serial = nodes(1:n)
    curve%log_discount = logs(1:n)
    if (present(status)) status = creditr_ok
  end subroutine build_zero_curve

  subroutine bootstrap_swap(base_date, maturity, swap_rate, conventions, nodes, logs, n, status)
    type(date_t), intent(in) :: base_date, maturity
    real(kind=dp), intent(in) :: swap_rate
    type(conventions_t), intent(in) :: conventions
    integer, intent(inout) :: nodes(:)
    real(kind=dp), intent(inout) :: logs(:)
    integer, intent(inout) :: n
    integer, intent(out) :: status
    integer :: k, npay, i, prev_serial, mat_serial
    integer, allocatable :: pay_serial(:)
    type(date_t) :: pay_date
    real(kind=dp) :: lo, hi, flo, fhi, mid, fmid, logdf_mat, w, previous_log

    mat_serial = maturity%serial()
    npay = ceiling(real(max(1, 12 * (maturity%year - base_date%year) + maturity%month - base_date%month), dp) / &
      real(conventions%fixed_frequency_months, dp))
    allocate(pay_serial(npay))
    k = 0
    do i = 1, npay
      pay_date = modified_following(add_months(base_date, i * conventions%fixed_frequency_months))
      if (pay_date%serial() >= mat_serial) then
        k = k + 1
        pay_serial(k) = mat_serial
        exit
      else
        k = k + 1
        pay_serial(k) = pay_date%serial()
      end if
    end do
    npay = k

    lo = -100.0_dp
    if (n > 0) then
      hi = min(-1.0e-14_dp, logs(n))
    else
      hi = -1.0e-8_dp
    end if
    flo = swap_equation(lo, pay_serial(1:npay), base_date, swap_rate, conventions, nodes, logs, n)
    fhi = swap_equation(hi, pay_serial(1:npay), base_date, swap_rate, conventions, nodes, logs, n)
    if (flo * fhi > 0.0_dp) then
      status = creditr_no_bracket
      return
    end if
    do i = 1, 200
      mid = 0.5_dp * (lo + hi)
      fmid = swap_equation(mid, pay_serial(1:npay), base_date, swap_rate, conventions, nodes, logs, n)
      if (abs(fmid) < 1.0e-14_dp .or. abs(hi - lo) < 1.0e-14_dp) exit
      if (flo * fmid <= 0.0_dp) then
        hi = mid
        fhi = fmid
      else
        lo = mid
        flo = fmid
      end if
    end do
    logdf_mat = 0.5_dp * (lo + hi)

    prev_serial = base_date%serial()
    previous_log = 0.0_dp
    if (n > 0) then
      prev_serial = nodes(n)
      previous_log = logs(n)
    end if
    do i = 1, npay
      if (pay_serial(i) > prev_serial) then
        w = real(pay_serial(i) - prev_serial, dp) / real(mat_serial - prev_serial, dp)
        call insert_node(nodes, logs, n, pay_serial(i), &
          (1.0_dp - w) * previous_log + w * logdf_mat)
      end if
    end do
    status = creditr_ok
  end subroutine bootstrap_swap

  real(kind=dp) function swap_equation(logdf_mat, pay_serial, base_date, swap_rate, conventions, &
      nodes, logs, n) result(value)
    real(kind=dp), intent(in) :: logdf_mat, swap_rate
    integer, intent(in) :: pay_serial(:), nodes(:), n
    real(kind=dp), intent(in) :: logs(:)
    type(date_t), intent(in) :: base_date
    type(conventions_t), intent(in) :: conventions
    integer :: i, prev
    real(kind=dp) :: annuity, df, alpha
    annuity = 0.0_dp
    prev = base_date%serial()
    do i = 1, size(pay_serial)
      df = trial_discount(pay_serial(i), pay_serial(size(pay_serial)), logdf_mat, base_date%serial(), nodes, logs, n)
      alpha = year_fraction_serial(prev, pay_serial(i), conventions%fixed_dcc)
      annuity = annuity + alpha * df
      prev = pay_serial(i)
    end do
    value = swap_rate * annuity + exp(logdf_mat) - 1.0_dp
  end function swap_equation

  pure real(kind=dp) function year_fraction_serial(s0, s1, convention) result(yf)
    use creditr_dates, only: date_from_serial
    integer, intent(in) :: s0, s1
    character(len=*), intent(in) :: convention
    yf = year_fraction(date_from_serial(s0), date_from_serial(s1), convention)
  end function year_fraction_serial

  pure real(kind=dp) function trial_discount(serial, maturity_serial, logdf_mat, base_serial, nodes, logs, n) result(df)
    integer, intent(in) :: serial, maturity_serial, base_serial, nodes(:), n
    real(kind=dp), intent(in) :: logdf_mat, logs(:)
    integer :: lo, hi, mid, left_serial
    real(kind=dp) :: left_log, w
    if (n == 0) then
      w = real(serial - base_serial, dp) / real(maturity_serial - base_serial, dp)
      df = exp(w * logdf_mat)
    else if (serial > nodes(n)) then
      left_serial = nodes(n)
      left_log = logs(n)
      w = real(serial - left_serial, dp) / real(maturity_serial - left_serial, dp)
      df = exp((1.0_dp - w) * left_log + w * logdf_mat)
    else if (serial <= nodes(1)) then
      w = real(serial - base_serial, dp) / real(nodes(1) - base_serial, dp)
      df = exp(w * logs(1))
    else
      lo = 1
      hi = n
      do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (nodes(mid) <= serial) then
          lo = mid
        else
          hi = mid
        end if
      end do
      if (nodes(lo) == serial) then
        df = exp(logs(lo))
      else
        w = real(serial - nodes(lo), dp) / real(nodes(hi) - nodes(lo), dp)
        df = exp((1.0_dp - w) * logs(lo) + w * logs(hi))
      end if
    end if
  end function trial_discount

  subroutine insert_node(nodes, logs, n, serial, logdf)
    integer, intent(inout) :: nodes(:)
    real(kind=dp), intent(inout) :: logs(:)
    integer, intent(inout) :: n
    integer, intent(in) :: serial
    real(kind=dp), intent(in) :: logdf
    if (n > 0) then
      if (serial == nodes(n)) then
        logs(n) = logdf
        return
      end if
    end if
    n = n + 1
    nodes(n) = serial
    logs(n) = logdf
  end subroutine insert_node

  subroutine read_rate_quotes_csv(filename, quotes, status)
    character(len=*), intent(in) :: filename
    type(rate_quote_t), allocatable, intent(out) :: quotes(:)
    integer, intent(out), optional :: status
    integer :: unit, ios, n, i
    character(len=256) :: line
    character(len=16) :: expiry
    character(len=1) :: kind
    real(kind=dp) :: rate
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      if (present(status)) status = creditr_io_error
      return
    end if
    read(unit, '(a)', iostat=ios) line
    n = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) n = n + 1
    end do
    rewind(unit)
    read(unit, '(a)') line
    allocate(quotes(n))
    do i = 1, n
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      call parse_quote_line(line, expiry, kind, rate, ios)
      if (ios /= 0) exit
      quotes(i)%expiry_months = expiry_to_months(expiry)
      quotes(i)%instrument_type = kind
      quotes(i)%rate = rate
    end do
    close(unit)
    if (present(status)) status = merge(creditr_ok, creditr_io_error, ios == 0 .or. ios < 0)
  end subroutine read_rate_quotes_csv

  subroutine parse_quote_line(line, expiry, kind, rate, ios)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: expiry
    character(len=1), intent(out) :: kind
    real(kind=dp), intent(out) :: rate
    integer, intent(out) :: ios
    integer :: p1, p2
    p1 = index(line, ',')
    p2 = index(line(p1 + 1:), ',') + p1
    if (p1 <= 0 .or. p2 <= p1) then
      ios = 1
      return
    end if
    expiry = adjustl(line(1:p1 - 1))
    kind = adjustl(line(p1 + 1:p2 - 1))
    read(line(p2 + 1:), *, iostat=ios) rate
  end subroutine parse_quote_line

  pure integer function expiry_to_months(expiry) result(months)
    character(len=*), intent(in) :: expiry
    integer :: n, ios
    character(len=1) :: suffix
    n = 0
    suffix = expiry(len_trim(expiry):len_trim(expiry))
    read(expiry(1:len_trim(expiry) - 1), *, iostat=ios) n
    if (ios /= 0) then
      months = 0
    else if (suffix == 'Y' .or. suffix == 'y') then
      months = 12 * n
    else
      months = n
    end if
  end function expiry_to_months

end module creditr_curve
