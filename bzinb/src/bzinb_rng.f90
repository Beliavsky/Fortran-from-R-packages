module bzinb_rng
  use bzinb_kinds, only : dp, pi_dp
  implicit none
  private
  public :: set_bzinb_seed, random_normal, random_gamma, random_poisson, categorical4
contains
  subroutine set_bzinb_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i + 7919*i*i, huge(1)-1)
      if (put(i) <= 0) put(i) = i + 17
    end do
    call random_seed(put=put)
  end subroutine set_bzinb_seed

  real(dp) function random_normal() result(z)
    real(dp) :: u1, u2
    call random_number(u1); call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi_dp*u2)
  end function random_normal

  recursive real(dp) function random_gamma(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d, c, z, v, u
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      x = random_gamma(shape + 1.0_dp, scale)*u**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = random_normal()
        v = 1.0_dp + c*z
        if (v > 0.0_dp) exit
      end do
      v = v**3
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(max(u,tiny(1.0_dp))) < 0.5_dp*z*z + d*(1.0_dp - v + log(v))) exit
    end do
    x = scale*d*v
  end function random_gamma

  integer function random_poisson(lambda) result(k)
    real(dp), intent(in) :: lambda
    real(dp) :: l, p, u, c, beta, alpha, x, lhs, rhs
    integer :: n
    if (lambda <= 0.0_dp) then
      k = 0
    else if (lambda < 30.0_dp) then
      l = exp(-lambda); p = 1.0_dp; n = 0
      do
        n = n + 1
        call random_number(u); p = p*u
        if (p <= l) exit
      end do
      k = n - 1
    else
      c = 0.767_dp - 3.36_dp/lambda
      beta = pi_dp/sqrt(3.0_dp*lambda)
      alpha = beta*lambda
      l = log(c) - lambda - log(beta)
      do
        call random_number(u)
        x = (alpha - log((1.0_dp-u)/u))/beta
        n = floor(x + 0.5_dp)
        if (n < 0) cycle
        call random_number(u)
        lhs = alpha - beta*x + log(u/((1.0_dp + exp(alpha-beta*x))**2))
        rhs = l + n*log(lambda) - log_gamma(real(n+1,dp))
        if (lhs <= rhs) exit
      end do
      k = n
    end if
  end function random_poisson

  integer function categorical4(p) result(k)
    real(dp), intent(in) :: p(4)
    real(dp) :: u, s
    integer :: i
    call random_number(u)
    s = 0.0_dp
    do i = 1, 4
      s = s + p(i)
      if (u <= s) then
        k = i
        return
      end if
    end do
    k = 4
  end function categorical4
end module bzinb_rng
