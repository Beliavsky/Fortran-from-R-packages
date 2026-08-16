! Modern Fortran translation of the computational core of R package expint.
!
! Upstream expint code:
!   Copyright (C) 2016-2026 Vincent Goulet
! GSL-derived algorithms:
!   Copyright (C) 2007 Brian Gough
!   Copyright (C) 1996-2002 Gerard Jungman
!
! This translation is distributed under the GNU GPL, version 3 or later.
!
! The exponential-integral Chebyshev expansions and incomplete-gamma
! algorithms below are translated from the corresponding upstream C files,
! which are derived from GSL/SLATEC as documented by the original package.
module expint_mod
   use, intrinsic :: iso_fortran_env, only : real64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, &
                                              ieee_is_nan
   implicit none
   private

   integer, parameter, public :: dp = real64
   real(dp), parameter :: eps = epsilon(1.0_dp)
   real(dp), parameter :: log_tiny = log(tiny(1.0_dp))

   public :: expint, expint_e1, expint_e2, expint_en, expint_ei
   public :: expint_recycle
   public :: gammainc, gamma_inc, gammainc_recycle

   interface expint
      module procedure expint_value
   end interface

   interface gammainc
      module procedure gamma_inc
   end interface

   real(dp), parameter :: ae11(39) = [ &
      0.121503239716065790_dp, -0.065088778513550150_dp, 0.004897651357459670_dp, &
     -0.000649237843027216_dp, 0.000093840434587471_dp, 0.000000420236380882_dp, &
     -0.000008113374735904_dp, 0.000002804247688663_dp, 0.000000056487164441_dp, &
     -0.000000344809174450_dp, 0.000000058209273578_dp, 0.000000038711426349_dp, &
     -0.000000012453235014_dp, -0.000000005118504888_dp, 0.000000002148771527_dp, &
      0.000000000868459898_dp, -0.000000000343650105_dp, -0.000000000179796603_dp, &
      0.000000000047442060_dp, 0.000000000040423282_dp, -0.000000000003543928_dp, &
     -0.000000000008853444_dp, -0.000000000000960151_dp, 0.000000000001692921_dp, &
      0.000000000000607990_dp, -0.000000000000224338_dp, -0.000000000000200327_dp, &
     -0.000000000000006246_dp, 0.000000000000045571_dp, 0.000000000000016383_dp, &
     -0.000000000000005561_dp, -0.000000000000006074_dp, -0.000000000000000862_dp, &
      0.000000000000001223_dp, 0.000000000000000716_dp, -0.000000000000000024_dp, &
     -0.000000000000000201_dp, -0.000000000000000082_dp, 0.000000000000000017_dp ]

   real(dp), parameter :: ae12(25) = [ &
      0.582417495134726740_dp, -0.158348850905782750_dp, -0.006764275590323141_dp, &
      0.005125843950185725_dp, 0.000435232492169391_dp, -0.000143613366305483_dp, &
     -0.000041801320556301_dp, -0.000002713395758640_dp, 0.000001151381913647_dp, &
      0.000000420650022012_dp, 0.000000066581901391_dp, 0.000000000662143777_dp, &
     -0.000000002844104870_dp, -0.000000000940724197_dp, -0.000000000177476602_dp, &
     -0.000000000015830222_dp, 0.000000000002905732_dp, 0.000000000001769356_dp, &
      0.000000000000492735_dp, 0.000000000000093709_dp, 0.000000000000010707_dp, &
     -0.000000000000000537_dp, -0.000000000000000716_dp, -0.000000000000000244_dp, &
     -0.000000000000000058_dp ]

   real(dp), parameter :: e11(19) = [ &
     -16.11346165557149402600_dp, 7.79407277874268027690_dp, -1.95540581886314195070_dp, &
      0.37337293866277945612_dp, -0.05692503191092901938_dp, 0.00721107776966009185_dp, &
     -0.00078104901449841593_dp, 0.00007388093356262168_dp, -0.00000620286187580820_dp, &
      0.00000046816002303176_dp, -0.00000003209288853329_dp, 0.00000000201519974874_dp, &
     -0.00000000011673686816_dp, 0.00000000000627627066_dp, -0.00000000000031481541_dp, &
      0.00000000000001479904_dp, -0.00000000000000065457_dp, 0.00000000000000002733_dp, &
     -0.00000000000000000108_dp ]

   real(dp), parameter :: e12(16) = [ &
     -0.03739021479220279500_dp, 0.04272398606220957700_dp, -0.13031820798497005440_dp, &
      0.01441912402469889073_dp, -0.00134617078051068022_dp, 0.00010731029253063780_dp, &
     -0.00000742999951611943_dp, 0.00000045377325690753_dp, -0.00000002476417211390_dp, &
      0.00000000122076581374_dp, -0.00000000005485141480_dp, 0.00000000000226362142_dp, &
     -0.00000000000008635897_dp, 0.00000000000000306291_dp, -0.00000000000000010148_dp, &
      0.00000000000000000315_dp ]

   real(dp), parameter :: ae13(25) = [ &
     -0.605773246640603460_dp, -0.112535243483660900_dp, 0.013432266247902779_dp, &
     -0.001926845187381145_dp, 0.000309118337720603_dp, -0.000053564132129618_dp, &
      0.000009827812880247_dp, -0.000001885368984916_dp, 0.000000374943193568_dp, &
     -0.000000076823455870_dp, 0.000000016143270567_dp, -0.000000003466802211_dp, &
      0.000000000758754209_dp, -0.000000000168864333_dp, 0.000000000038145706_dp, &
     -0.000000000008733026_dp, 0.000000000002023672_dp, -0.000000000000474132_dp, &
      0.000000000000112211_dp, -0.000000000000026804_dp, 0.000000000000006457_dp, &
     -0.000000000000001568_dp, 0.000000000000000383_dp, -0.000000000000000094_dp, &
      0.000000000000000023_dp ]

   real(dp), parameter :: ae14(26) = [ &
     -0.18929180007530170_dp, -0.08648117855259871_dp, 0.00722410154374659_dp, &
     -0.00080975594575573_dp, 0.00010999134432661_dp, -0.00001717332998937_dp, &
      0.00000298562751447_dp, -0.00000056596491457_dp, 0.00000011526808397_dp, &
     -0.00000002495030440_dp, 0.00000000569232420_dp, -0.00000000135995766_dp, &
      0.00000000033846628_dp, -0.00000000008737853_dp, 0.00000000002331588_dp, &
     -0.00000000000641148_dp, 0.00000000000181224_dp, -0.00000000000052538_dp, &
      0.00000000000015592_dp, -0.00000000000004729_dp, 0.00000000000001463_dp, &
     -0.00000000000000461_dp, 0.00000000000000148_dp, -0.00000000000000048_dp, &
      0.00000000000000016_dp, -0.00000000000000005_dp ]

contains

   pure real(dp) function quiet_nan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure real(dp) function positive_inf() result(x)
      x = ieee_value(0.0_dp, ieee_positive_inf)
   end function positive_inf

   pure real(dp) function cheb_eval(c, x) result(value)
      real(dp), intent(in) :: c(:), x
      real(dp) :: d, dd, tmp, y, y2
      integer :: j

      y = x
      y2 = 2.0_dp * y
      d = 0.0_dp
      dd = 0.0_dp
      do j = size(c), 2, -1
         tmp = d
         d = y2 * d - dd + c(j)
         dd = tmp
      end do
      value = y * d - dd + 0.5_dp * c(1)
   end function cheb_eval

   pure elemental real(dp) function expint_e1(x, scale) result(res)
      real(dp), intent(in) :: x
      logical, intent(in), optional :: scale
      logical :: scaled
      real(dp) :: xmaxt, xmax, s, ch, ln_term

      if (ieee_is_nan(x)) then
         res = x
         return
      end if
      scaled = .false.
      if (present(scale)) scaled = scale

      xmaxt = -log_tiny
      xmax = xmaxt - log(xmaxt)

      if (x < -xmax .and. .not. scaled) then
         res = positive_inf()
      else if (x <= -10.0_dp) then
         s = merge(1.0_dp, exp(-x), scaled) / x
         ch = cheb_eval(ae11, 20.0_dp / x + 1.0_dp)
         res = s * (1.0_dp + ch)
      else if (x <= -4.0_dp) then
         s = merge(1.0_dp, exp(-x), scaled) / x
         ch = cheb_eval(ae12, (40.0_dp / x + 7.0_dp) / 3.0_dp)
         res = s * (1.0_dp + ch)
      else if (x <= -1.0_dp) then
         s = merge(exp(x), 1.0_dp, scaled)
         ln_term = -log(abs(x))
         ch = cheb_eval(e11, (2.0_dp * x + 5.0_dp) / 3.0_dp)
         res = s * (ln_term + ch)
      else if (x == 0.0_dp) then
         res = quiet_nan()
      else if (x <= 1.0_dp) then
         s = merge(exp(x), 1.0_dp, scaled)
         ln_term = -log(abs(x))
         ch = cheb_eval(e12, x)
         res = s * (ln_term - 0.6875_dp + x + ch)
      else if (x <= 4.0_dp) then
         s = merge(1.0_dp, exp(-x), scaled) / x
         ch = cheb_eval(ae13, (8.0_dp / x - 5.0_dp) / 3.0_dp)
         res = s * (1.0_dp + ch)
      else if (x <= xmax .or. scaled) then
         s = merge(1.0_dp, exp(-x), scaled) / x
         ch = cheb_eval(ae14, 8.0_dp / x - 1.0_dp)
         res = s * (1.0_dp + ch)
         if (abs(res) < tiny(1.0_dp)) res = 0.0_dp
      else
         res = 0.0_dp
      end if
   end function expint_e1

   pure elemental real(dp) function expint_e2(x, scale) result(res)
      real(dp), intent(in) :: x
      logical, intent(in), optional :: scale
      logical :: scaled
      real(dp) :: xmaxt, xmax, s, ex, y, sm, sum6

      if (ieee_is_nan(x)) then
         res = x
         return
      end if
      scaled = .false.
      if (present(scale)) scaled = scale

      xmaxt = -log_tiny
      xmax = xmaxt - log(xmaxt)
      if (x < -xmax .and. .not. scaled) then
         res = positive_inf()
      else if (x == 0.0_dp) then
         res = 1.0_dp
      else if (x < 100.0_dp) then
         ex = merge(1.0_dp, exp(-x), scaled)
         res = ex - x * expint_e1(x, scaled)
      else if (x < xmax .or. scaled) then
         s = merge(1.0_dp, exp(-x), scaled)
         y = 1.0_dp / x
         sum6 = -720.0_dp + y * (5040.0_dp + y * (-40320.0_dp + y * (362880.0_dp + &
                y * (-3628800.0_dp + y * (39916800.0_dp + y * (-479001600.0_dp + &
                y * (6227020800.0_dp + y * (-87178291200.0_dp))))))))
         sm = y * (-2.0_dp + y * (6.0_dp + y * (-24.0_dp + y * (120.0_dp + y * sum6))))
         res = s * (1.0_dp + sm) / x
         if (abs(res) < tiny(1.0_dp)) res = 0.0_dp
      else
         res = 0.0_dp
      end if
   end function expint_e2

   pure real(dp) function gamma_inc_f_cf(a, x) result(hn)
      real(dp), intent(in) :: a, x
      integer, parameter :: nmax = 5000
      real(dp) :: small, cn, dn, an, delta
      integer :: n

      small = eps ** 3
      hn = 1.0_dp
      cn = 1.0_dp / small
      dn = 1.0_dp
      do n = 2, nmax - 1
         if (mod(n, 2) == 1) then
            an = 0.5_dp * real(n - 1, dp) / x
         else
            an = (0.5_dp * real(n, dp) - a) / x
         end if
         dn = 1.0_dp + an * dn
         if (abs(dn) < small) dn = small
         cn = 1.0_dp + an / cn
         if (abs(cn) < small) cn = small
         dn = 1.0_dp / dn
         delta = cn * dn
         hn = hn * delta
         if (abs(delta - 1.0_dp) < eps) exit
      end do
   end function gamma_inc_f_cf

   pure real(dp) function gamma_upper_positive(a, x) result(gax)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax = 10000
      real(dp) :: ap, del, sm, gln, p, b, c, d, h, an, delta
      integer :: n

      if (x == 0.0_dp) then
         gax = gamma(a)
         return
      end if
      gln = log_gamma(a)

      if (x < a + 1.0_dp) then
         ap = a
         sm = 1.0_dp / a
         del = sm
         do n = 1, itmax
            ap = ap + 1.0_dp
            del = del * x / ap
            sm = sm + del
            if (abs(del) <= abs(sm) * eps) exit
         end do
         p = sm * exp(-x + a * log(x) - gln)
         if (p >= 1.0_dp) then
            gax = 0.0_dp
         else
            gax = exp(gln + log1p_safe(-p))
         end if
      else
         b = x + 1.0_dp - a
         c = 1.0_dp / tiny(1.0_dp)
         d = 1.0_dp / b
         h = d
         do n = 1, itmax
            an = -real(n, dp) * (real(n, dp) - a)
            b = b + 2.0_dp
            d = an * d + b
            if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
            c = b + an / c
            if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
            d = 1.0_dp / d
            delta = d * c
            h = h * delta
            if (abs(delta - 1.0_dp) <= eps) exit
         end do
         gax = exp(-x + a * log(x)) * h
      end if
   end function gamma_upper_positive

   pure real(dp) function log1p_safe(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: term, sm
      integer :: k

      if (abs(x) > 0.125_dp) then
         y = log(1.0_dp + x)
      else
         term = x
         sm = term
         do k = 2, 200
            term = -term * x * real(k - 1, dp) / real(k, dp)
            sm = sm + term
            if (abs(term) < eps * max(1.0_dp, abs(sm))) exit
         end do
         y = sm
      end if
   end function log1p_safe

   pure elemental real(dp) function gamma_inc(a, x) result(gax)
      real(dp), intent(in) :: a, x
      real(dp) :: da, fa, shift, alpha

      if (ieee_is_nan(a) .or. ieee_is_nan(x)) then
         gax = a + x
      else if (x < 0.0_dp) then
         gax = quiet_nan()
      else if (x == 0.0_dp) then
         gax = gamma(a)
      else if (a == 0.0_dp) then
         gax = expint_e1(x, .false.)
      else if (a > 0.0_dp) then
         gax = gamma_upper_positive(a, x)
      else if (x > 0.25_dp) then
         gax = exp((a - 1.0_dp) * log(x) - x) * gamma_inc_f_cf(a, x)
      else if (abs(a) < 0.5_dp) then
         da = a + 1.0_dp
         gax = gamma_upper_positive(da, x)
         shift = exp(-x + a * log(x))
         gax = (gax - shift) / a
      else
         fa = floor(a)
         da = a - fa
         if (da > 0.0_dp) then
            gax = gamma_upper_positive(da, x)
         else
            gax = expint_e1(x, .false.)
         end if
         alpha = da
         do
            shift = exp(-x + (alpha - 1.0_dp) * log(x))
            gax = (gax - shift) / (alpha - 1.0_dp)
            alpha = alpha - 1.0_dp
            if (alpha <= a) exit
         end do
      end if
   end function gamma_inc

   pure elemental real(dp) function expint_en(x, order, scale) result(res)
      real(dp), intent(in) :: x
      integer, intent(in) :: order
      logical, intent(in), optional :: scale
      logical :: scaled
      real(dp) :: s

      if (ieee_is_nan(x)) then
         res = x
         return
      end if
      scaled = .false.
      if (present(scale)) scaled = scale

      if (order < 0) then
         res = quiet_nan()
      else if (order == 0) then
         if (x == 0.0_dp) then
            res = quiet_nan()
         else
            res = merge(1.0_dp, exp(-x), scaled) / x
            if (abs(res) < tiny(1.0_dp)) res = 0.0_dp
         end if
      else if (order == 1) then
         res = expint_e1(x, scaled)
      else if (order == 2) then
         res = expint_e2(x, scaled)
      else if (x < 0.0_dp) then
         res = quiet_nan()
      else if (x == 0.0_dp) then
         res = 1.0_dp / real(order - 1, dp)
      else
         s = merge(exp(x), 1.0_dp, scaled)
         res = s * gamma_inc(real(1 - order, dp), x) * x ** (order - 1)
         if (abs(res) < tiny(1.0_dp)) res = 0.0_dp
      end if
   end function expint_en

   pure elemental real(dp) function expint_value(x, order, scale) result(res)
      real(dp), intent(in) :: x
      integer, intent(in), optional :: order
      logical, intent(in), optional :: scale
      integer :: n
      logical :: scaled

      n = 1
      if (present(order)) n = order
      scaled = .false.
      if (present(scale)) scaled = scale
      res = expint_en(x, n, scaled)
   end function expint_value

   pure elemental real(dp) function expint_ei(x, scale) result(res)
      real(dp), intent(in) :: x
      logical, intent(in), optional :: scale
      logical :: scaled

      scaled = .false.
      if (present(scale)) scaled = scale
      res = -expint_e1(-x, scaled)
   end function expint_ei

   pure function expint_recycle(x, order, scale) result(values)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: order(:)
      logical, intent(in), optional :: scale
      real(dp), allocatable :: values(:)
      logical :: scaled
      integer :: i, ix, io, n

      scaled = .false.
      if (present(scale)) scaled = scale
      if (size(x) == 0 .or. size(order) == 0) then
         allocate(values(0))
         return
      end if
      n = max(size(x), size(order))
      allocate(values(n))
      do i = 1, n
         ix = 1 + modulo(i - 1, size(x))
         io = 1 + modulo(i - 1, size(order))
         values(i) = expint_en(x(ix), order(io), scaled)
      end do
   end function expint_recycle

   pure function gammainc_recycle(a, x) result(values)
      real(dp), intent(in) :: a(:), x(:)
      real(dp), allocatable :: values(:)
      integer :: i, ia, ix, n

      if (size(a) == 0 .or. size(x) == 0) then
         allocate(values(0))
         return
      end if
      n = max(size(a), size(x))
      allocate(values(n))
      do i = 1, n
         ia = 1 + modulo(i - 1, size(a))
         ix = 1 + modulo(i - 1, size(x))
         values(i) = gamma_inc(a(ia), x(ix))
      end do
   end function gammainc_recycle

end module expint_mod
