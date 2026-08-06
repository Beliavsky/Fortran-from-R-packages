module arfima_polynomial
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok, arfima_invalid_input
  use arfima_types, only : arfima_error, set_error
  implicit none
  private
  public :: ar_to_pacf, pacf_to_ar, is_stationary_polynomial
  public :: polynomial_convolution, shift_polynomial, fractional_coefficients
  public :: difference_series, integrate_series, differencing_polynomial
  public :: psi_weights_ar, psi_weights, exact_integration_weights

contains

  function ar_to_pacf(phi) result(pacf)
    real(dp), intent(in) :: phi(:)
    real(dp), allocatable :: pacf(:)
    real(dp), allocatable :: work(:), next(:)
    real(dp) :: a
    integer :: l, k, ll, i

    l = size(phi)
    allocate(pacf(l))
    if (l == 0) return
    allocate(work(l))
    work = phi
    do k = 1, l
      ll = l + 1 - k
      a = work(ll)
      pacf(ll) = a
      if (ll == 1) exit
      if (abs(abs(a) - 1.0_dp) <= 10.0_dp * epsilon(1.0_dp)) exit
      allocate(next(ll-1))
      do i = 1, ll - 1
        next(i) = (work(i) + a * work(ll-i)) / (1.0_dp - a*a)
      end do
      work(1:ll-1) = next
      deallocate(next)
    end do
  end function ar_to_pacf

  function pacf_to_ar(pacf) result(phi)
    real(dp), intent(in) :: pacf(:)
    real(dp), allocatable :: phi(:)
    real(dp), allocatable :: old(:)
    integer :: l, k, i

    l = size(pacf)
    allocate(phi(l))
    if (l == 0) return
    phi(1) = pacf(1)
    do k = 2, l
      allocate(old(k-1))
      old = phi(1:k-1)
      do i = 1, k - 1
        phi(i) = old(i) - pacf(k) * old(k-i)
      end do
      phi(k) = pacf(k)
      deallocate(old)
    end do
  end function pacf_to_ar

  logical function is_stationary_polynomial(phi, margin)
    real(dp), intent(in) :: phi(:)
    real(dp), intent(in), optional :: margin
    real(dp), allocatable :: pacf(:)
    real(dp) :: m
    m = 0.0_dp
    if (present(margin)) m = margin
    if (size(phi) == 0) then
      is_stationary_polynomial = .true.
      return
    end if
    pacf = ar_to_pacf(phi)
    is_stationary_polynomial = all(abs(pacf) < 1.0_dp - m)
  end function is_stationary_polynomial

  function polynomial_convolution(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp), allocatable :: c(:)
    integer :: i, j
    if (size(a) == 0 .or. size(b) == 0) then
      allocate(c(0))
      return
    end if
    allocate(c(size(a)+size(b)-1))
    c = 0.0_dp
    do i = 1, size(a)
      do j = 1, size(b)
        c(i+j-1) = c(i+j-1) + a(i)*b(j)
      end do
    end do
  end function polynomial_convolution

  function shift_polynomial(a, period) result(b)
    real(dp), intent(in) :: a(:)
    integer, intent(in) :: period
    real(dp), allocatable :: b(:)
    integer :: i
    if (size(a) == 0) then
      allocate(b(0))
      return
    end if
    if (period == 0) then
      allocate(b(size(a)))
      b = a
    else if (period > 0) then
      allocate(b((size(a)-1)*period+1))
      b = 0.0_dp
      do i = 1, size(a)
        b(1+(i-1)*period) = a(i)
      end do
    else
      allocate(b((size(a)-1)/abs(period)+1))
      do i = 1, size(b)
        b(i) = a(1+(i-1)*abs(period))
      end do
    end if
  end function shift_polynomial

  function fractional_coefficients(d, n, bterm) result(c)
    real(dp), intent(in) :: d
    integer, intent(in) :: n
    real(dp), intent(in), optional :: bterm
    real(dp), allocatable :: c(:)
    real(dp) :: bt
    integer :: k
    bt = -1.0_dp
    if (present(bterm)) bt = bterm
    allocate(c(n+1))
    c(1) = 1.0_dp
    do k = 1, n
      c(k+1) = c(k) * (d-real(k-1,dp))/real(k,dp) * bt
    end do
  end function fractional_coefficients

  function differencing_polynomial(dint, dseas, period) result(a)
    integer, intent(in) :: dint, dseas, period
    real(dp), allocatable :: a(:), one(:), seas(:), tmp(:)
    integer :: i
    allocate(a(1))
    a = 1.0_dp
    if (dint > 0) then
      allocate(one(2))
      one = [1.0_dp, -1.0_dp]
      do i = 1, dint
        tmp = polynomial_convolution(a, one)
        call move_alloc(tmp, a)
      end do
    end if
    if (dseas > 0) then
      allocate(seas(period+1))
      seas = 0.0_dp
      seas(1) = 1.0_dp
      seas(period+1) = -1.0_dp
      do i = 1, dseas
        tmp = polynomial_convolution(a, seas)
        call move_alloc(tmp, a)
      end do
    end if
  end function differencing_polynomial

  subroutine difference_series(z, dint, dseas, period, y, error)
    real(dp), intent(in) :: z(:)
    integer, intent(in) :: dint, dseas, period
    real(dp), allocatable, intent(out) :: y(:)
    type(arfima_error), intent(out) :: error
    real(dp), allocatable :: a(:)
    integer :: lag, t, j

    call set_error(error, arfima_ok, '')
    if (dint < 0 .or. dseas < 0 .or. (dseas > 0 .and. period < 2)) then
      allocate(y(0))
      call set_error(error, arfima_invalid_input, 'invalid differencing order or period')
      return
    end if
    a = differencing_polynomial(dint, dseas, period)
    lag = size(a)-1
    if (size(z) <= lag) then
      allocate(y(0))
      call set_error(error, arfima_invalid_input, 'series too short for requested differencing')
      return
    end if
    allocate(y(size(z)-lag))
    do t = lag+1, size(z)
      y(t-lag) = 0.0_dp
      do j = 0, lag
        y(t-lag) = y(t-lag) + a(j+1)*z(t-j)
      end do
    end do
  end subroutine difference_series

  subroutine integrate_series(y, zinit, dint, dseas, period, z, error)
    real(dp), intent(in) :: y(:), zinit(:)
    integer, intent(in) :: dint, dseas, period
    real(dp), allocatable, intent(out) :: z(:)
    type(arfima_error), intent(out) :: error
    real(dp), allocatable :: a(:)
    integer :: lag, t, j

    call set_error(error, arfima_ok, '')
    a = differencing_polynomial(dint, dseas, period)
    lag = size(a)-1
    if (size(zinit) /= lag) then
      allocate(z(0))
      call set_error(error, arfima_invalid_input, 'zinit length must equal dint + period*dseas')
      return
    end if
    allocate(z(lag+size(y)))
    if (lag > 0) z(1:lag) = zinit
    do t = lag+1, size(z)
      z(t) = y(t-lag)
      do j = 1, lag
        z(t) = z(t) - a(j+1)*z(t-j)
      end do
    end do
  end subroutine integrate_series

  function psi_weights_ar(phi, maxlag) result(x)
    real(dp), intent(in) :: phi(:)
    integer, intent(in) :: maxlag
    real(dp), allocatable :: x(:)
    integer :: i, j, p
    p = size(phi)
    allocate(x(maxlag+1))
    x = 0.0_dp
    x(1) = 1.0_dp
    do i = 1, maxlag
      do j = 1, min(i,p)
        x(i+1) = x(i+1) + phi(j)*x(i-j+1)
      end do
    end do
  end function psi_weights_ar

  function psi_weights(phi, theta, phiseas, thetaseas, dfrac, dfs, dint, dseas, period, n) result(w)
    real(dp), intent(in) :: phi(:), theta(:), phiseas(:), thetaseas(:)
    real(dp), intent(in) :: dfrac, dfs
    integer, intent(in) :: dint, dseas, period, n
    real(dp), allocatable :: w(:), tmp(:), q(:), ps(:)

    allocate(w(1)); w = 1.0_dp
    if (size(theta)>0) then
      allocate(q(size(theta)+1)); q = [1.0_dp, -theta]
      tmp = polynomial_convolution(w,q); call move_alloc(tmp,w)
    end if
    if (size(phi)>0) then
      ps = psi_weights_ar(phi,n)
      tmp = polynomial_convolution(w,ps); call move_alloc(tmp,w)
    end if
    if (dint /= 0 .or. abs(dfrac)>0.0_dp) then
      q = fractional_coefficients(-real(dint,dp)-dfrac,n)
      tmp = polynomial_convolution(w,q); call move_alloc(tmp,w)
    end if
    if (period>0) then
      if (size(phiseas)>0) then
        ps = shift_polynomial(psi_weights_ar(phiseas,n),period)
        tmp = polynomial_convolution(w,ps); call move_alloc(tmp,w)
      end if
      if (size(thetaseas)>0) then
        q = shift_polynomial([1.0_dp,-thetaseas],period)
        tmp = polynomial_convolution(w,q); call move_alloc(tmp,w)
      end if
      if (dseas /= 0 .or. abs(dfs)>0.0_dp) then
        q = shift_polynomial(fractional_coefficients(-real(dseas,dp)-dfs,n),period)
        tmp = polynomial_convolution(w,q); call move_alloc(tmp,w)
      end if
    end if
    if (size(w)<n) then
      allocate(tmp(n)); tmp=0.0_dp; tmp(1:size(w))=w; call move_alloc(tmp,w)
    else
      w=w(1:n)
    end if
  end function psi_weights

  function exact_integration_weights(dint,dseas,period,len) result(w)
    integer, intent(in) :: dint,dseas,period,len
    real(dp), allocatable :: w(:), a(:)
    a = differencing_polynomial(dint,dseas,period)
    allocate(w(len)); w=0.0_dp; w(1)=1.0_dp
    call inverse_filter_weights(a,w)
  contains
    subroutine inverse_filter_weights(poly,out)
      real(dp), intent(in) :: poly(:)
      real(dp), intent(inout) :: out(:)
      integer :: i,j
      do i=2,size(out)
        do j=1,min(i-1,size(poly)-1)
          out(i)=out(i)-poly(j+1)*out(i-j)
        end do
      end do
    end subroutine inverse_filter_weights
  end function exact_integration_weights

end module arfima_polynomial
