module arfima_autocov
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok, arfima_invalid_input, arfima_not_stationary, arfima_singular
  use arfima_types, only : arfima_error, arfima_spec, arfima_parameters, set_error, &
    long_memory_none, long_memory_fd, long_memory_fgn, long_memory_pla
  use arfima_polynomial, only : is_stationary_polynomial, shift_polynomial
  use arfima_linalg, only : solve_linear
  implicit none
  private
  public :: tacvf_arma, tacvf_fdwn, tacvf_fgn, tacvf_pla
  public :: mix_autocov, tacvf_farma, tacvf_arfima, riemann_zeta_borwein

contains

  subroutine tacvf_arma(phi, theta, maxlag, sigma2, gamma_out, error)
    real(dp), intent(in) :: phi(:), theta(:)
    integer, intent(in) :: maxlag
    real(dp), intent(in), optional :: sigma2
    real(dp), allocatable, intent(out) :: gamma_out(:)
    type(arfima_error), intent(out) :: error
    integer :: p, q, rr, i, j, k, info
    real(dp) :: s2
    real(dp), allocatable :: c(:), b(:), theta2(:), phi2(:), a(:,:), g(:), sol(:)

    call set_error(error, arfima_ok, '')
    s2 = 1.0_dp
    if (present(sigma2)) s2 = sigma2
    if (maxlag < 0 .or. s2 < 0.0_dp) then
      allocate(gamma_out(0))
      call set_error(error, arfima_invalid_input, 'maxlag and sigma2 must be nonnegative')
      return
    end if
    if (.not. is_stationary_polynomial(phi) .or. .not. is_stationary_polynomial(theta)) then
      allocate(gamma_out(0))
      call set_error(error, arfima_not_stationary, 'AR or MA polynomial is outside the PACF stability region')
      return
    end if
    p = size(phi)
    q = size(theta)
    allocate(gamma_out(maxlag+1))
    gamma_out = 0.0_dp
    if (max(p,q) == 0) then
      gamma_out(1) = s2
      return
    end if
    rr = max(p,q)+1
    allocate(c(q+1), b(rr), theta2(q+1), phi2(3*rr))
    c = 0.0_dp
    b = 0.0_dp
    theta2 = 0.0_dp
    phi2 = 0.0_dp
    c(1) = 1.0_dp
    theta2(1) = -1.0_dp
    if (q>0) theta2(2:q+1) = theta
    phi2(rr) = -1.0_dp
    if (p>0) phi2(rr+1:rr+p) = phi
    do k = 1, q
      c(k+1) = -theta(k)
      do i = 1, min(p,k)
        c(k+1) = c(k+1) + phi(i)*c(k+1-i)
      end do
    end do
    do k = 0, q
      do i = k, q
        b(k+1) = b(k+1) - theta2(i+1)*c(i-k+1)
      end do
    end do
    if (p == 0) then
      gamma_out(1:min(maxlag+1,q+1)) = s2*b(1:min(maxlag+1,q+1))
      return
    end if
    allocate(a(rr,rr))
    do i = 1, rr
      do j = 1, rr
        if (j == 1) then
          a(i,j) = phi2(rr+i-1)
        else
          a(i,j) = phi2(rr+i-j) + phi2(rr+i+j-2)
        end if
      end do
    end do
    call solve_linear(a, -b, sol, info)
    if (info /= 0) then
      gamma_out = 0.0_dp
      call set_error(error, arfima_singular, 'singular ARMA autocovariance equations')
      return
    end if
    allocate(g(max(maxlag+1,rr)))
    g = 0.0_dp
    g(1:rr) = sol
    do i = rr+1, maxlag+1
      do j = 1, p
        g(i) = g(i) + phi(j)*g(i-j)
      end do
    end do
    gamma_out = s2*g(1:maxlag+1)
  end subroutine tacvf_arma

  subroutine tacvf_fdwn(dfrac, maxlag, gamma_out, error)
    real(dp), intent(in) :: dfrac
    integer, intent(in) :: maxlag
    real(dp), allocatable, intent(out) :: gamma_out(:)
    type(arfima_error), intent(out) :: error
    integer :: i
    call set_error(error, arfima_ok, '')
    if (dfrac <= -1.0_dp .or. dfrac >= 0.5_dp .or. maxlag < 0) then
      allocate(gamma_out(0))
      call set_error(error, arfima_not_stationary, 'fractional d must be in (-1,0.5)')
      return
    end if
    allocate(gamma_out(maxlag+1))
    gamma_out(1) = gamma(1.0_dp-2.0_dp*dfrac)/(gamma(1.0_dp-dfrac)**2)
    do i = 1, maxlag
      gamma_out(i+1) = (real(i-1,dp)+dfrac)/(real(i,dp)-dfrac)*gamma_out(i)
    end do
  end subroutine tacvf_fdwn

  subroutine tacvf_fgn(hurst, maxlag, gamma_out, error)
    real(dp), intent(in) :: hurst
    integer, intent(in) :: maxlag
    real(dp), allocatable, intent(out) :: gamma_out(:)
    type(arfima_error), intent(out) :: error
    integer :: k
    real(dp) :: h2
    call set_error(error, arfima_ok, '')
    if (hurst <= 0.0_dp .or. hurst >= 1.0_dp .or. maxlag < 0) then
      allocate(gamma_out(0))
      call set_error(error, arfima_not_stationary, 'Hurst parameter must be in (0,1)')
      return
    end if
    h2 = 2.0_dp*hurst
    allocate(gamma_out(maxlag+1))
    do k = 0, maxlag
      gamma_out(k+1) = 0.5_dp*(abs(real(k+1,dp))**h2 - 2.0_dp*abs(real(k,dp))**h2 + &
        abs(real(k-1,dp))**h2)
    end do
  end subroutine tacvf_fgn

  real(dp) function riemann_zeta_borwein(s, n) result(zeta)
    real(dp), intent(in) :: s
    integer, intent(in), optional :: n
    integer :: nn, k
    real(dp), allocatable :: d(:)
    real(dp) :: temp
    nn = 20
    if (present(n)) nn = n
    allocate(d(0:nn))
    d(0) = 1.0_dp
    do k = 1, nn
      temp = gamma(real(nn+k,dp))/(gamma(real(nn-k+1,dp))*gamma(real(2*k+1,dp)))
      d(k) = d(k-1) + real(nn,dp)*temp*4.0_dp**k
    end do
    zeta = 0.0_dp
    do k = 0, nn-1
      zeta = zeta + (-1.0_dp)**k*(d(k)-d(nn))/real(k+1,dp)**s
    end do
    zeta = -zeta/(d(nn)*(1.0_dp-2.0_dp**(1.0_dp-s)))
  end function riemann_zeta_borwein

  subroutine tacvf_pla(alpha, maxlag, gamma_out, error)
    real(dp), intent(in) :: alpha
    integer, intent(in) :: maxlag
    real(dp), allocatable, intent(out) :: gamma_out(:)
    type(arfima_error), intent(out) :: error
    integer :: k
    real(dp) :: c
    call set_error(error, arfima_ok, '')
    if (alpha <= 0.0_dp .or. alpha >= 3.0_dp .or. maxlag < 0) then
      allocate(gamma_out(0))
      call set_error(error, arfima_not_stationary, 'PLA alpha must be in (0,3)')
      return
    end if
    allocate(gamma_out(maxlag+1))
    gamma_out(1) = 1.0_dp
    if (maxlag > 0) then
      c = -1.0_dp/(2.0_dp*riemann_zeta_borwein(alpha))
      do k = 1, maxlag
        gamma_out(k+1) = c*real(k,dp)**(-alpha)
      end do
    end if
  end subroutine tacvf_pla

  function mix_autocov(x, y, maxlag) result(z)
    real(dp), intent(in) :: x(:), y(:)
    integer, intent(in), optional :: maxlag
    real(dp), allocatable :: z(:)
    integer :: m, k, j, nx, ny
    nx = size(x)-1
    ny = size(y)-1
    m = min(nx,ny)
    if (present(maxlag)) m = min(maxlag,min(nx,ny))
    allocate(z(m+1))
    z = 0.0_dp
    do k = 0, m
      do j = -nx, nx
        if (abs(k-j) <= ny) z(k+1) = z(k+1) + x(abs(j)+1)*y(abs(k-j)+1)
      end do
    end do
  end function mix_autocov

  subroutine tacvf_farma(phi, theta, lmodel, long_parameter, maxlag, sigma2, gamma_out, error, no_lag_truncation)
    real(dp), intent(in) :: phi(:), theta(:)
    integer, intent(in) :: lmodel, maxlag
    real(dp), intent(in) :: long_parameter, sigma2
    real(dp), allocatable, intent(out) :: gamma_out(:)
    type(arfima_error), intent(out) :: error
    logical, intent(in), optional :: no_lag_truncation
    logical :: no_trunc, only_long
    integer :: lag_trunc
    real(dp), allocatable :: x(:), y(:), z(:)
    type(arfima_error) :: e

    call set_error(error, arfima_ok, '')
    no_trunc = .false.
    if (present(no_lag_truncation)) no_trunc = no_lag_truncation
    only_long = size(phi)==0 .and. size(theta)==0
    if (lmodel == long_memory_none) then
      call tacvf_arma(phi,theta,maxlag,sigma2,gamma_out,error)
      return
    end if
    if (no_trunc .or. only_long) then
      lag_trunc = maxlag
    else
      lag_trunc = max(maxlag,256)
    end if
    select case(lmodel)
    case(long_memory_fd)
      call tacvf_fdwn(long_parameter,lag_trunc,x,e)
    case(long_memory_fgn)
      call tacvf_fgn(long_parameter,lag_trunc,x,e)
    case(long_memory_pla)
      call tacvf_pla(long_parameter,lag_trunc,x,e)
    case default
      allocate(gamma_out(0))
      call set_error(error, arfima_invalid_input, 'unknown long-memory model')
      return
    end select
    if (e%code /= arfima_ok) then
      allocate(gamma_out(0)); error=e; return
    end if
    if (only_long) then
      allocate(gamma_out(maxlag+1)); gamma_out=sigma2*x(1:maxlag+1); return
    end if
    call tacvf_arma(phi,theta,lag_trunc,1.0_dp,y,e)
    if (e%code /= arfima_ok) then
      allocate(gamma_out(0)); error=e; return
    end if
    z = mix_autocov(x,y,maxlag)
    allocate(gamma_out(maxlag+1))
    gamma_out = sigma2*z(1:maxlag+1)
  end subroutine tacvf_farma

  subroutine tacvf_arfima(spec, params, maxlag, sigma2, gamma_out, error)
    type(arfima_spec), intent(in) :: spec
    type(arfima_parameters), intent(in) :: params
    integer, intent(in) :: maxlag
    real(dp), intent(in), optional :: sigma2
    real(dp), allocatable, intent(out) :: gamma_out(:)
    type(arfima_error), intent(out) :: error
    real(dp) :: s2, lp, slp
    real(dp), allocatable :: phi(:),theta(:),phs(:),ths(:),model(:),seas(:),shifted(:),mixed(:)
    type(arfima_error) :: e
    integer :: lag_trunc

    call set_error(error, arfima_ok, '')
    s2=1.0_dp; if(present(sigma2)) s2=sigma2
    call copy_or_empty(params%phi,phi)
    call copy_or_empty(params%theta,theta)
    call copy_or_empty(params%phiseas,phs)
    call copy_or_empty(params%thetaseas,ths)
    lp = get_long_parameter(spec%lmodel,params,.false.)
    if (spec%period <= 0) then
      if (size(phs)>0 .or. size(ths)>0 .or. spec%slmodel/=long_memory_none) then
        allocate(gamma_out(0)); call set_error(error,arfima_invalid_input,'seasonal terms require period >= 2'); return
      end if
      call tacvf_farma(phi,theta,spec%lmodel,lp,maxlag,s2,gamma_out,error)
      return
    end if
    if (spec%period < 2) then
      allocate(gamma_out(0)); call set_error(error,arfima_invalid_input,'seasonal period must be >= 2'); return
    end if
    lag_trunc=max(maxlag,256)*max(2,spec%period)
    call tacvf_farma(phi,theta,spec%lmodel,lp,lag_trunc,1.0_dp,model,e,.true.)
    if(e%code/=arfima_ok) then; allocate(gamma_out(0)); error=e; return; end if
    slp=get_long_parameter(spec%slmodel,params,.true.)
    call tacvf_farma(phs,ths,spec%slmodel,slp,2*max(maxlag,256),1.0_dp,seas,e,.true.)
    if(e%code/=arfima_ok) then; allocate(gamma_out(0)); error=e; return; end if
    shifted=shift_polynomial(seas,spec%period)
    if(size(shifted)<size(model)) then
      allocate(gamma_out(0)); call set_error(error,arfima_invalid_input,'seasonal autocovariance truncation is too short'); return
    end if
    mixed=mix_autocov(model,shifted,maxlag)
    allocate(gamma_out(maxlag+1)); gamma_out=s2*mixed(1:maxlag+1)
  contains
    subroutine copy_or_empty(source,target)
      real(dp), allocatable, intent(in) :: source(:)
      real(dp), allocatable, intent(out) :: target(:)
      if(allocated(source)) then
        allocate(target(size(source))); target=source
      else
        allocate(target(0))
      end if
    end subroutine copy_or_empty
    real(dp) function get_long_parameter(model,p,seasonal) result(value)
      integer,intent(in)::model
      type(arfima_parameters),intent(in)::p
      logical,intent(in)::seasonal
      value=0.0_dp
      if(.not.seasonal) then
        select case(model)
        case(long_memory_fd); value=p%dfrac
        case(long_memory_fgn); value=p%hurst
        case(long_memory_pla); value=p%alpha
        end select
      else
        select case(model)
        case(long_memory_fd); value=p%dfs
        case(long_memory_fgn); value=p%hurst_seasonal
        case(long_memory_pla); value=p%alpha_seasonal
        end select
      end if
    end function get_long_parameter
  end subroutine tacvf_arfima

end module arfima_autocov
