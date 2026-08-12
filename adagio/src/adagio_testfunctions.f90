! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_testfunctions
  use adagio_kinds, only : dp, pi
  implicit none
  private
  public :: fn_rosenbrock, gr_rosenbrock, fn_rastrigin, gr_rastrigin
  public :: fn_nesterov, gr_nesterov, fn_nesterov1, fn_nesterov2
  public :: fn_hald, gr_hald, fn_shor, gr_shor, fn_trefethen, fn_wagon

contains

  pure function fn_rosenbrock(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    integer :: i
    f = 0.0_dp
    do i = 1, size(x)-1
       f = f + 100.0_dp*(x(i+1)-x(i)**2)**2 + (1.0_dp-x(i))**2
    end do
  end function

  pure function gr_rosenbrock(x) result(g)
    real(dp), intent(in) :: x(:)
    real(dp) :: g(size(x))
    integer :: i, n
    n = size(x); g = 0.0_dp
    if (n == 1) then
       g(1) = 2.0_dp*(x(1)-1.0_dp); return
    end if
    g(1) = 2.0_dp*(x(1)-1.0_dp) + 400.0_dp*x(1)*(x(1)**2-x(2))
    do i = 2, n-1
       g(i) = 2.0_dp*(x(i)-1.0_dp) + 400.0_dp*x(i)*(x(i)**2-x(i+1)) &
            + 200.0_dp*(x(i)-x(i-1)**2)
    end do
    g(n) = 200.0_dp*(x(n)-x(n-1)**2)
  end function

  pure function fn_rastrigin(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = 10.0_dp*real(size(x),dp) + sum(x*x - 10.0_dp*cos(2.0_dp*pi*x))
  end function

  pure function gr_rastrigin(x) result(g)
    real(dp), intent(in) :: x(:)
    real(dp) :: g(size(x))
    g = 2.0_dp*x + 20.0_dp*pi*sin(2.0_dp*pi*x)
  end function

  pure function fn_nesterov(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    integer :: i
    f = (1.0_dp-x(1))**2/4.0_dp
    do i = 1, size(x)-1
       f = f + (1.0_dp+x(i+1)-2.0_dp*x(i)**2)**2
    end do
  end function

  pure function gr_nesterov(x) result(g)
    real(dp), intent(in) :: x(:)
    real(dp) :: g(size(x)), r
    integer :: i
    g = 0.0_dp
    g(1) = (x(1)-1.0_dp)/2.0_dp
    do i = 1, size(x)-1
       r = 1.0_dp+x(i+1)-2.0_dp*x(i)**2
       g(i+1) = g(i+1) + 2.0_dp*r
       g(i) = g(i) - 8.0_dp*x(i)*r
    end do
  end function

  pure function fn_nesterov1(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    integer :: i
    f = (1.0_dp-x(1))**2/4.0_dp
    do i = 1, size(x)-1
       f = f + abs(1.0_dp+x(i+1)-2.0_dp*x(i)**2)
    end do
  end function

  pure function fn_nesterov2(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    integer :: i
    f = abs(1.0_dp-x(1))/4.0_dp
    do i = 1, size(x)-1
       f = f + abs(1.0_dp+x(i+1)-2.0_dp*abs(x(i)))
    end do
  end function

  pure function fn_hald(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f, t, r
    integer :: k
    f = 0.0_dp
    do k = 1, 21
       t = -1.0_dp + real(k-1,dp)/10.0_dp
       r = (x(1)+x(2)*t)/(1.0_dp+x(3)*t+x(4)*t*t+x(5)*t**3) - exp(t)
       f = max(f, abs(r))
    end do
  end function

  function gr_hald(x) result(g)
    real(dp), intent(in) :: x(:)
    real(dp) :: g(5), t, a, b, r, best, s1, s2
    integer :: k, kb
    g = 0.0_dp; best = -1.0_dp; kb = 1
    do k = 1, 21
       t = -1.0_dp + real(k-1,dp)/10.0_dp
       r = (x(1)+x(2)*t)/(1.0_dp+x(3)*t+x(4)*t*t+x(5)*t**3) - exp(t)
       if (abs(r) > best) then; best = abs(r); kb = k; end if
    end do
    t = -1.0_dp + real(kb-1,dp)*0.1_dp
    a = 1.0_dp+x(3)*t+x(4)*t*t+x(5)*t**3
    b = x(1)+x(2)*t
    if (abs(a) <= 0.0_dp) then; g = ieee_nan(); return; end if
    r = b/a-exp(t); s1 = sign(1.0_dp,r); s2 = sign(1.0_dp,t*a*b)
    g(1:2) = s1*[1.0_dp,t]/a
    g(3:5) = s2*[t,t*t,t**3]*b/(a*a)
  contains
    function ieee_nan() result(v)
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      real(dp) :: v(5)
      v = ieee_value(0.0_dp, ieee_quiet_nan)
    end function
  end function

  pure function fn_shor(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f, d
    real(dp), parameter :: a(5,10) = reshape([ &
      0._dp,0._dp,0._dp,0._dp,0._dp, 2._dp,1._dp,1._dp,1._dp,3._dp, &
      1._dp,2._dp,1._dp,1._dp,2._dp, 1._dp,4._dp,1._dp,2._dp,2._dp, &
      3._dp,2._dp,1._dp,0._dp,1._dp, 0._dp,2._dp,1._dp,0._dp,1._dp, &
      1._dp,1._dp,1._dp,1._dp,1._dp, 1._dp,0._dp,1._dp,2._dp,1._dp, &
      0._dp,0._dp,2._dp,1._dp,0._dp, 1._dp,1._dp,2._dp,0._dp,0._dp ], [5,10])
    real(dp), parameter :: b(10) = [1._dp,5._dp,10._dp,2._dp,4._dp,3._dp,1.7_dp,2.5_dp,6._dp,4.5_dp]
    integer :: i
    f = 0.0_dp
    do i = 1, 10
       d = b(i)*sum((x-a(:,i))**2)
       if (d > f) f = d
    end do
  end function

  pure function gr_shor(x) result(g)
    real(dp), intent(in) :: x(:)
    real(dp) :: g(5), f, d
    real(dp), parameter :: a(5,10) = reshape([ &
      0._dp,0._dp,0._dp,0._dp,0._dp, 2._dp,1._dp,1._dp,1._dp,3._dp, &
      1._dp,2._dp,1._dp,1._dp,2._dp, 1._dp,4._dp,1._dp,2._dp,2._dp, &
      3._dp,2._dp,1._dp,0._dp,1._dp, 0._dp,2._dp,1._dp,0._dp,1._dp, &
      1._dp,1._dp,1._dp,1._dp,1._dp, 1._dp,0._dp,1._dp,2._dp,1._dp, &
      0._dp,0._dp,2._dp,1._dp,0._dp, 1._dp,1._dp,2._dp,0._dp,0._dp ], [5,10])
    real(dp), parameter :: b(10) = [1._dp,5._dp,10._dp,2._dp,4._dp,3._dp,1.7_dp,2.5_dp,6._dp,4.5_dp]
    integer :: i, k
    f = 0.0_dp; k = 1
    do i = 1, 10
       d = b(i)*sum((x-a(:,i))**2)
       if (d > f) then; f = d; k = i; end if
    end do
    g = 2.0_dp*b(k)*(x-a(:,k))
  end function

  pure function fn_trefethen(p) result(f)
    real(dp), intent(in) :: p(:)
    real(dp) :: f, x, y
    x=p(1); y=p(2)
    f=exp(sin(50*x))+sin(60*exp(y))+sin(70*sin(x))+sin(sin(80*y)) &
      -sin(10*(x+y))+(x*x+y*y)/4
  end function

  pure function fn_wagon(p) result(f)
    real(dp), intent(in) :: p(:)
    real(dp) :: f, x, y, z
    x=p(1); y=p(2); z=p(3)
    f=exp(sin(50*x))+sin(60*exp(y))*sin(60*z)+sin(70*sin(x))*cos(10*z) &
      +sin(sin(80*y))-sin(10*(x+z))+(x*x+y*y+z*z)/4
  end function
end module adagio_testfunctions
