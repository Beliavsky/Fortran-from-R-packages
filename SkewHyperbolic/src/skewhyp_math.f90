module skewhyp_math
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  integer, parameter, public :: dp = real64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  public :: bessel_k, normal_rand, gamma_rand, adaptive_simpson, binom_coeff, double_factorial_odd

  abstract interface
    function scalar_fun(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fun
  end interface
contains
  function bessel_k(nu, x) result(v)
    real(dp), intent(in) :: nu, x
    real(dp) :: v, h, t, s, term, mu
    integer :: i, n
    if (x <= 0.0_dp) then
      v = huge(1.0_dp)
      return
    end if
    if (x > 50.0_dp) then
      mu = 4.0_dp*nu*nu
      v = sqrt(pi/(2.0_dp*x))*exp(-x) * &
          (1.0_dp + (mu-1.0_dp)/(8.0_dp*x) + &
          (mu-1.0_dp)*(mu-9.0_dp)/(2.0_dp*(8.0_dp*x)**2))
      return
    end if
    n = 320
    h = 18.0_dp/real(n,dp)
    s = 0.0_dp
    do i = 0, n
      t = h*real(i,dp)
      term = exp(-x*cosh(t))*cosh(nu*t)
      if (i == 0 .or. i == n) then
        s = s + term
      else if (mod(i,2) == 0) then
        s = s + 2.0_dp*term
      else
        s = s + 4.0_dp*term
      end if
    end do
    v = h*s/3.0_dp
  end function bessel_k

  function normal_rand() result(z)
    real(dp) :: z, u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function normal_rand

  recursive function gamma_rand(shape) result(x)
    real(dp), intent(in) :: shape
    real(dp) :: x, d, c, z, u
    if (shape <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      x = gamma_rand(shape + 1.0_dp)*max(u,tiny(1.0_dp))**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = normal_rand()
        if (1.0_dp + c*z > 0.0_dp) exit
      end do
      x = (1.0_dp + c*z)**3
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(max(u,tiny(1.0_dp))) < 0.5_dp*z*z + d*(1.0_dp-x+log(x))) exit
    end do
    x = d*x
  end function gamma_rand

  function adaptive_simpson(f,a,b,tol,maxdepth) result(res)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a,b,tol
    integer, intent(in), optional :: maxdepth
    real(dp) :: res, fa, fb, fm, whole
    integer :: md
    md = 18
    if (present(maxdepth)) md = maxdepth
    fa = f(a); fb = f(b); fm = f(0.5_dp*(a+b))
    whole = (b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
    res = as_rec(f,a,b,fa,fm,fb,whole,tol,md)
  end function adaptive_simpson

  recursive function as_rec(f,a,b,fa,fm,fb,whole,tol,depth) result(res)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a,b,fa,fm,fb,whole,tol
    integer, intent(in) :: depth
    real(dp) :: res,c,lm,rm,flm,frm,left,right
    c = 0.5_dp*(a+b); lm = 0.5_dp*(a+c); rm = 0.5_dp*(c+b)
    flm = f(lm); frm = f(rm)
    left = (c-a)*(fa+4.0_dp*flm+fm)/6.0_dp
    right = (b-c)*(fm+4.0_dp*frm+fb)/6.0_dp
    if (depth <= 0 .or. abs(left+right-whole) <= 15.0_dp*tol) then
      res = left + right + (left+right-whole)/15.0_dp
    else
      res = as_rec(f,a,c,fa,flm,fm,left,tol/2.0_dp,depth-1) + &
            as_rec(f,c,b,fm,frm,fb,right,tol/2.0_dp,depth-1)
    end if
  end function as_rec

  pure function binom_coeff(n,k) result(v)
    integer, intent(in) :: n,k
    real(dp) :: v
    integer :: i, kk
    if (k < 0 .or. k > n) then
      v = 0.0_dp
      return
    end if
    kk = min(k,n-k)
    v = 1.0_dp
    do i=1,kk
      v = v*real(n-kk+i,dp)/real(i,dp)
    end do
  end function binom_coeff

  pure function double_factorial_odd(r) result(v)
    integer, intent(in) :: r
    real(dp) :: v
    integer :: j
    if (r == 0) then
      v = 1.0_dp
      return
    end if
    v = 1.0_dp
    do j=1,r
      v = v*real(2*j-1,dp)
    end do
  end function double_factorial_odd
end module skewhyp_math
