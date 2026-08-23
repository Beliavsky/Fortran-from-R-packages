module gh_math
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  integer, parameter, public :: dp = real64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  public :: bessel_k, normal_rand, normal_cdf, adaptive_simpson
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
    real(dp) :: v, h, t, s, term
    integer :: i, n
    if (x <= 0.0_dp) then
      v = huge(1.0_dp); return
    end if
    if (x > 40.0_dp) then
      v = sqrt(pi/(2.0_dp*x))*exp(-x)*(1.0_dp + (4.0_dp*nu*nu-1.0_dp)/(8.0_dp*x))
      return
    end if
    n = 240
    h = 16.0_dp/real(n,dp)
    s = 0.0_dp
    do i=0,n
      t = h*real(i,dp)
      term = exp(-x*cosh(t))*cosh(nu*t)
      if (i==0 .or. i==n) then
        s = s + term
      else if (mod(i,2)==0) then
        s = s + 2.0_dp*term
      else
        s = s + 4.0_dp*term
      end if
    end do
    v = h*s/3.0_dp
  end function bessel_k

  function normal_rand() result(z)
    real(dp) :: z, u1, u2
    call random_number(u1); call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function normal_rand

  pure function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  function adaptive_simpson(f,a,b,tol,maxdepth) result(res)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a,b,tol
    integer, intent(in), optional :: maxdepth
    real(dp) :: res, fa, fb, fm, whole
    integer :: md
    md=18; if (present(maxdepth)) md=maxdepth
    fa=f(a); fb=f(b); fm=f(0.5_dp*(a+b))
    whole=(b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
    res=as_rec(f,a,b,fa,fm,fb,whole,tol,md)
  end function adaptive_simpson

  recursive function as_rec(f,a,b,fa,fm,fb,whole,tol,depth) result(res)
    procedure(scalar_fun) :: f
    real(dp),intent(in)::a,b,fa,fm,fb,whole,tol
    integer,intent(in)::depth
    real(dp)::res,c,lm,rm,flm,frm,left,right
    c=0.5_dp*(a+b); lm=0.5_dp*(a+c); rm=0.5_dp*(c+b)
    flm=f(lm); frm=f(rm)
    left=(c-a)*(fa+4.0_dp*flm+fm)/6.0_dp
    right=(b-c)*(fm+4.0_dp*frm+fb)/6.0_dp
    if (depth<=0 .or. abs(left+right-whole)<=15.0_dp*tol) then
      res=left+right+(left+right-whole)/15.0_dp
    else
      res=as_rec(f,a,c,fa,flm,fm,left,tol/2.0_dp,depth-1)+ &
          as_rec(f,c,b,fm,frm,fb,right,tol/2.0_dp,depth-1)
    end if
  end function as_rec
end module gh_math
