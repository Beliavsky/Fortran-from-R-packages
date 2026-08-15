! Modern Fortran translation of R package skewunit.
! SPDX-License-Identifier: GPL-2.0-or-later
module skewunit_distributions
   use skewunit_kinds, only : dp, pi, sqrt2pi, nan_dp, neg_inf_dp
   use skewunit_special, only : normal_pdf, normal_cdf, regularized_beta, &
      beta_log_pdf, safe_log_prob, safe_log1m_prob, clamp01
   use skewunit_rng, only : rand_uniform, rand_normal, rand_beta
   implicit none
   private

   integer, parameter, public :: family_none = 0
   integer, parameter, public :: family_asin = 1
   integer, parameter, public :: family_uquad = 2
   integer, parameter, public :: family_triang = 3
   integer, parameter, public :: family_jsb = 4
   integer, parameter, public :: family_sbeta = 5
   integer, parameter, public :: n_families = 5

   public :: cuberoot, family_id, family_name, family_has_delta
   public :: dasin, pasin, rasin
   public :: dtriang, ptriang, rtriang
   public :: duquad, puquad, ruquad
   public :: djsb, pjsb, rjsb
   public :: dsbeta, psbeta, rsbeta
   public :: baseline_density, baseline_cdf, baseline_random
   public :: dskewunit, pskewunit, pskewunit_vec, rskewunit
   public :: rasin_vec, rtriang_vec, ruquad_vec, rjsb_vec, rsbeta_vec
   public :: rskewunit_vec

contains

   elemental real(dp) function cuberoot(x) result(y)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         y = x**(1.0_dp/3.0_dp)
      else
         y = -(-x)**(1.0_dp/3.0_dp)
      end if
   end function cuberoot

   pure integer function family_id(name) result(id)
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: s
      s = trim(adjustl(name))
      select case(s)
      case('asin','ArcSin','arcsin','ASIN')
         id = family_asin
      case('Uquad','uquad','UQUAD')
         id = family_uquad
      case('triang','triangular','TRIANG')
         id = family_triang
      case('JSB','jsb')
         id = family_jsb
      case('sbeta','SBETA')
         id = family_sbeta
      case('','none','NONE','null','NULL')
         id = family_none
      case default
         id = -1
      end select
   end function family_id

   pure function family_name(id) result(name)
      integer, intent(in) :: id
      character(len=8) :: name
      select case(id)
      case(family_asin)
         name = 'asin'
      case(family_uquad)
         name = 'Uquad'
      case(family_triang)
         name = 'triang'
      case(family_jsb)
         name = 'JSB'
      case(family_sbeta)
         name = 'sbeta'
      case(family_none)
         name = ''
      case default
         name = 'unknown'
      end select
   end function family_name

   pure logical function family_has_delta(family) result(has_delta)
      integer, intent(in) :: family
      has_delta = family == family_jsb .or. family == family_sbeta
   end function family_has_delta

   elemental real(dp) function dasin(x, log_pdf) result(v)
      real(dp), intent(in) :: x
      logical, intent(in), optional :: log_pdf
      logical :: lp
      real(dp) :: lv
      lp = .false.
      if (present(log_pdf)) lp = log_pdf
      if (x <= 0.0_dp .or. x >= 1.0_dp) then
         v = 0.0_dp
         return
      end if
      lv = -log(pi)-0.5_dp*log(x)-0.5_dp*log(1.0_dp-x)
      if (lp) then
         v = lv
      else
         v = exp(lv)
      end if
   end function dasin

   elemental real(dp) function pasin(q, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      real(dp) :: p
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      p = regularized_beta(q,0.5_dp,0.5_dp)
      if (.not.lower) p = 1.0_dp-p
      if (lp) then
         v = safe_log_prob(p)
      else
         v = p
      end if
   end function pasin

   real(dp) function rasin() result(x)
      x = rand_beta(0.5_dp,0.5_dp)
   end function rasin

   elemental real(dp) function dtriang(x, log_pdf) result(v)
      real(dp), intent(in) :: x
      logical, intent(in), optional :: log_pdf
      logical :: lp
      real(dp) :: lv
      lp = .false.
      if (present(log_pdf)) lp = log_pdf
      if (x <= 0.0_dp .or. x >= 1.0_dp) then
         v = 0.0_dp
         return
      end if
      if (x <= 0.5_dp) then
         lv = log(4.0_dp)+log(x)
      else
         lv = log(4.0_dp)+log(1.0_dp-x)
      end if
      if (lp) then
         v = lv
      else
         v = exp(lv)
      end if
   end function dtriang

   elemental real(dp) function ptriang(q, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      real(dp) :: p
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (q <= 0.0_dp) then
         p = 0.0_dp
      else if (q >= 1.0_dp) then
         p = 1.0_dp
      else if (q <= 0.5_dp) then
         p = 2.0_dp*q*q
      else
         p = 1.0_dp-2.0_dp*(1.0_dp-q)**2
      end if
      if (.not.lower) p = 1.0_dp-p
      if (lp) then
         v = safe_log_prob(p)
      else
         v = p
      end if
   end function ptriang

   real(dp) function rtriang() result(x)
      real(dp) :: u
      u = rand_uniform()
      if (u <= 0.5_dp) then
         x = sqrt(u/2.0_dp)
      else
         x = 1.0_dp-sqrt(0.5_dp*(1.0_dp-u))
      end if
   end function rtriang

   elemental real(dp) function duquad(x, a, b, log_pdf) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: a, b
      logical, intent(in), optional :: log_pdf
      real(dp) :: aa, bb, alpha, beta, lv
      logical :: lp
      aa = 0.0_dp
      bb = 1.0_dp
      lp = .false.
      if (present(a)) aa = a
      if (present(b)) bb = b
      if (present(log_pdf)) lp = log_pdf
      if (aa >= bb) then
         v = nan_dp()
         return
      end if
      if (x <= aa .or. x >= bb) then
         v = 0.0_dp
         return
      end if
      alpha = 12.0_dp/(bb-aa)**3
      beta = 0.5_dp*(aa+bb)
      if (x == beta) then
         if (lp) then
            v = neg_inf_dp()
         else
            v = 0.0_dp
         end if
         return
      end if
      lv = log(alpha)+2.0_dp*log(abs(x-beta))
      if (lp) then
         v = lv
      else
         v = exp(lv)
      end if
   end function duquad

   elemental real(dp) function puquad(q, a, b, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q
      real(dp), intent(in), optional :: a, b
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: aa, bb, alpha, beta, p
      logical :: lower, lp
      aa = 0.0_dp
      bb = 1.0_dp
      lower = .true.
      lp = .false.
      if (present(a)) aa = a
      if (present(b)) bb = b
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (aa >= bb) then
         v = nan_dp()
         return
      end if
      if (q <= aa) then
         p = 0.0_dp
      else if (q >= bb) then
         p = 1.0_dp
      else
         alpha = 12.0_dp/(bb-aa)**3
         beta = 0.5_dp*(aa+bb)
         p = alpha*((q-beta)**3+(beta-aa)**3)/3.0_dp
         p = clamp01(p)
      end if
      if (.not.lower) p = 1.0_dp-p
      if (lp) then
         v = safe_log_prob(p)
      else
         v = p
      end if
   end function puquad

   real(dp) function ruquad(a, b) result(x)
      real(dp), intent(in), optional :: a, b
      real(dp) :: aa, bb, alpha, beta
      aa = 0.0_dp
      bb = 1.0_dp
      if (present(a)) aa = a
      if (present(b)) bb = b
      if (aa >= bb) then
         x = nan_dp()
         return
      end if
      alpha = 12.0_dp/(bb-aa)**3
      beta = 0.5_dp*(aa+bb)
      x = beta+cuberoot(3.0_dp*rand_uniform()/alpha-(beta-aa)**3)
   end function ruquad

   elemental real(dp) function djsb(x, delta, log_pdf) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: delta
      logical, intent(in), optional :: log_pdf
      real(dp) :: d, eta, lv
      logical :: lp
      d = 1.0_dp
      lp = .false.
      if (present(delta)) d = delta
      if (present(log_pdf)) lp = log_pdf
      if (d <= 0.0_dp) then
         v = nan_dp()
         return
      end if
      if (x <= 0.0_dp .or. x >= 1.0_dp) then
         v = 0.0_dp
         return
      end if
      eta = log(x)-log(1.0_dp-x)
      lv = log(d)-log(x)-log(1.0_dp-x) &
           -0.5_dp*(d*eta)**2-log(sqrt2pi)
      if (lp) then
         v = lv
      else
         v = exp(lv)
      end if
   end function djsb

   elemental real(dp) function pjsb(q, delta, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q
      real(dp), intent(in), optional :: delta
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: d, eta, p
      logical :: lower, lp
      d = 1.0_dp
      lower = .true.
      lp = .false.
      if (present(delta)) d = delta
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (d <= 0.0_dp) then
         v = nan_dp()
         return
      end if
      if (q <= 0.0_dp) then
         p = 0.0_dp
      else if (q >= 1.0_dp) then
         p = 1.0_dp
      else
         eta = log(q)-log(1.0_dp-q)
         p = normal_cdf(d*eta)
      end if
      if (.not.lower) p = 1.0_dp-p
      if (lp) then
         v = safe_log_prob(p)
      else
         v = p
      end if
   end function pjsb

   real(dp) function rjsb(delta) result(x)
      real(dp), intent(in), optional :: delta
      real(dp) :: d, z
      d = 1.0_dp
      if (present(delta)) d = delta
      if (d <= 0.0_dp) then
         x = nan_dp()
         return
      end if
      z = rand_normal()/d
      if (z >= 0.0_dp) then
         x = 1.0_dp/(1.0_dp+exp(-z))
      else
         x = exp(z)/(1.0_dp+exp(z))
      end if
   end function rjsb

   elemental real(dp) function dsbeta(x, delta, log_pdf) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: delta
      logical, intent(in), optional :: log_pdf
      real(dp) :: d, lv
      logical :: lp
      d = 1.0_dp
      lp = .false.
      if (present(delta)) d = delta
      if (present(log_pdf)) lp = log_pdf
      if (d <= 0.0_dp) then
         v = nan_dp()
         return
      end if
      if (x <= 0.0_dp .or. x >= 1.0_dp) then
         v = 0.0_dp
         return
      end if
      lv = beta_log_pdf(x,d,d)
      if (lp) then
         v = lv
      else
         v = exp(lv)
      end if
   end function dsbeta

   elemental real(dp) function psbeta(q, delta, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q
      real(dp), intent(in), optional :: delta
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: d, p
      logical :: lower, lp
      d = 1.0_dp
      lower = .true.
      lp = .false.
      if (present(delta)) d = delta
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (d <= 0.0_dp) then
         v = nan_dp()
         return
      end if
      p = regularized_beta(q,d,d)
      if (.not.lower) p = 1.0_dp-p
      if (lp) then
         v = safe_log_prob(p)
      else
         v = p
      end if
   end function psbeta

   real(dp) function rsbeta(delta) result(x)
      real(dp), intent(in), optional :: delta
      real(dp) :: d
      d = 1.0_dp
      if (present(delta)) d = delta
      x = rand_beta(d,d)
   end function rsbeta


   subroutine rasin_vec(n, x)
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n)
      integer :: i
      do i = 1, n
         x(i) = rasin()
      end do
   end subroutine rasin_vec

   subroutine rtriang_vec(n, x)
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n)
      integer :: i
      do i = 1, n
         x(i) = rtriang()
      end do
   end subroutine rtriang_vec

   subroutine ruquad_vec(n, x, a, b)
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n)
      real(dp), intent(in), optional :: a, b
      integer :: i
      do i = 1, n
         x(i) = ruquad(a,b)
      end do
   end subroutine ruquad_vec

   subroutine rjsb_vec(n, x, delta)
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n)
      real(dp), intent(in), optional :: delta
      integer :: i
      do i = 1, n
         x(i) = rjsb(delta)
      end do
   end subroutine rjsb_vec

   subroutine rsbeta_vec(n, x, delta)
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n)
      real(dp), intent(in), optional :: delta
      integer :: i
      do i = 1, n
         x(i) = rsbeta(delta)
      end do
   end subroutine rsbeta_vec

   elemental real(dp) function baseline_density(x, family, delta, log_pdf) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: family
      real(dp), intent(in), optional :: delta
      logical, intent(in), optional :: log_pdf
      real(dp) :: d
      logical :: lp
      d = 1.0_dp
      lp = .false.
      if (present(delta)) d = delta
      if (present(log_pdf)) lp = log_pdf
      select case(family)
      case(family_asin)
         v = dasin(x,lp)
      case(family_uquad)
         v = duquad(x,log_pdf=lp)
      case(family_triang)
         v = dtriang(x,lp)
      case(family_jsb)
         v = djsb(x,d,lp)
      case(family_sbeta)
         v = dsbeta(x,d,lp)
      case default
         v = nan_dp()
      end select
   end function baseline_density

   elemental real(dp) function baseline_cdf(x, family, delta, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: family
      real(dp), intent(in), optional :: delta
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: d
      logical :: lower, lp
      d = 1.0_dp
      lower = .true.
      lp = .false.
      if (present(delta)) d = delta
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      select case(family)
      case(family_asin)
         v = pasin(x,lower,lp)
      case(family_uquad)
         v = puquad(x,lower_tail=lower,log_p=lp)
      case(family_triang)
         v = ptriang(x,lower,lp)
      case(family_jsb)
         v = pjsb(x,d,lower,lp)
      case(family_sbeta)
         v = psbeta(x,d,lower,lp)
      case default
         v = nan_dp()
      end select
   end function baseline_cdf

   real(dp) function baseline_random(family, delta) result(x)
      integer, intent(in) :: family
      real(dp), intent(in), optional :: delta
      real(dp) :: d
      d = 1.0_dp
      if (present(delta)) d = delta
      select case(family)
      case(family_asin)
         x = rasin()
      case(family_uquad)
         x = ruquad()
      case(family_triang)
         x = rtriang()
      case(family_jsb)
         x = rjsb(d)
      case(family_sbeta)
         x = rsbeta(d)
      case default
         x = nan_dp()
      end select
   end function baseline_random

   elemental real(dp) function dskewunit(x, lambda, delta, delta2, &
      family1, family2, log_pdf) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: lambda, delta, delta2
      integer, intent(in), optional :: family1, family2
      logical, intent(in), optional :: log_pdf
      real(dp) :: lam, d1, d2, df, dg, z, lf, pg, lv
      integer :: f1, f2
      logical :: lp

      lam = 0.0_dp
      d1 = 1.0_dp
      d2 = 1.0_dp
      f1 = family_asin
      f2 = family_asin
      lp = .false.
      if (present(lambda)) lam = lambda
      if (present(delta)) d1 = delta
      if (present(delta2)) d2 = delta2
      if (present(family1)) f1 = family1
      if (present(family2)) f2 = family2
      if (present(log_pdf)) lp = log_pdf

      if (f1 < family_asin .or. f1 > family_sbeta .or. &
          f2 < family_none .or. f2 > family_sbeta .or. &
          abs(lam) > 1.0_dp .or. d1 <= 0.0_dp .or. d2 <= 0.0_dp) then
         v = nan_dp()
         return
      end if
      if (x <= 0.0_dp .or. x >= 1.0_dp) then
         v = 0.0_dp
         return
      end if

      if (f2 == family_none) then
         v = baseline_density(x,f1,d1,lp)
         return
      end if

      df = d1
      if (.not.family_has_delta(f1)) df = 1.0_dp
      if (family_has_delta(f2)) then
         if (family_has_delta(f1)) then
            dg = d2
         else
            dg = d1
         end if
      else
         dg = 1.0_dp
      end if

      lf = baseline_density(x,f1,df,.true.)
      z = lam*(x-0.5_dp)+0.5_dp
      pg = baseline_cdf(z,f2,dg)
      if (pg <= 0.0_dp) then
         lv = neg_inf_dp()
      else
         lv = log(2.0_dp)+log(pg)+lf
      end if
      if (lp) then
         v = lv
      else
         v = exp(lv)
      end if
   end function dskewunit

   real(dp) function pskewunit(x, lambda, delta, delta2, family1, family2, &
      log_p, rel_tol) result(v)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: lambda, delta, delta2
      integer, intent(in), optional :: family1, family2
      logical, intent(in), optional :: log_p
      real(dp), intent(in), optional :: rel_tol
      real(dp) :: lam, d1, d2, tol, p
      integer :: f1, f2
      logical :: lp

      lam = 0.0_dp
      d1 = 1.0_dp
      d2 = 1.0_dp
      f1 = family_asin
      f2 = family_asin
      lp = .false.
      tol = 1.0e-10_dp
      if (present(lambda)) lam = lambda
      if (present(delta)) d1 = delta
      if (present(delta2)) d2 = delta2
      if (present(family1)) f1 = family1
      if (present(family2)) f2 = family2
      if (present(log_p)) lp = log_p
      if (present(rel_tol)) tol = rel_tol

      if (x <= 0.0_dp) then
         p = 0.0_dp
      else if (x >= 1.0_dp) then
         p = 1.0_dp
      else if (f2 == family_none) then
         p = baseline_cdf(x,f1,d1)
      else if (x <= 0.5_dp) then
         p = skew_integral_theta(0.0_dp,asin(sqrt(x)),lam,d1,d2,f1,f2,tol)
         p = clamp01(p)
      else
         p = 1.0_dp-skew_integral_theta(asin(sqrt(x)),0.5_dp*pi, &
            lam,d1,d2,f1,f2,tol)
         p = clamp01(p)
      end if

      if (lp) then
         v = safe_log_prob(p)
      else
         v = p
      end if

   end function pskewunit


   real(dp) function skew_integral_theta(a, b, lambda, delta, delta2, &
      family1, family2, rel_tol) result(ans)
      real(dp), intent(in) :: a, b, lambda, delta, delta2, rel_tol
      integer, intent(in) :: family1, family2
      real(dp) :: err
      call skew_adapt(a,b,lambda,delta,delta2,family1,family2, &
         max(rel_tol,1.0e-14_dp),1.0e-13_dp,20,ans,err)
   end function skew_integral_theta

   recursive subroutine skew_adapt(a, b, lambda, delta, delta2, family1, &
      family2, rtol, atol, depth, ans, err)
      real(dp), intent(in) :: a, b, lambda, delta, delta2, rtol, atol
      integer, intent(in) :: family1, family2, depth
      real(dp), intent(out) :: ans, err
      real(dp) :: whole, ewhole, left, right, el, er, mid, tolerance

      call skew_gk15(a,b,lambda,delta,delta2,family1,family2,whole,ewhole)
      tolerance = max(atol,rtol*abs(whole))
      if (ewhole <= tolerance .or. depth <= 1) then
         ans = whole
         err = ewhole
         return
      end if
      mid = 0.5_dp*(a+b)
      call skew_adapt(a,mid,lambda,delta,delta2,family1,family2,rtol, &
         0.5_dp*atol,depth-1,left,el)
      call skew_adapt(mid,b,lambda,delta,delta2,family1,family2,rtol, &
         0.5_dp*atol,depth-1,right,er)
      ans = left+right
      err = el+er
   end subroutine skew_adapt

   subroutine skew_gk15(a, b, lambda, delta, delta2, family1, family2, &
      resk, abserr)
      real(dp), intent(in) :: a, b, lambda, delta, delta2
      integer, intent(in) :: family1, family2
      real(dp), intent(out) :: resk, abserr
      real(dp), parameter :: xgk(8) = [ &
         0.991455371120812639206854697526329_dp, &
         0.949107912342758524526189684047851_dp, &
         0.864864423359769072789712788640926_dp, &
         0.741531185599394439863864773280788_dp, &
         0.586087235467691130294144838258730_dp, &
         0.405845151377397166906606412076961_dp, &
         0.207784955007898467600689403773245_dp,0.0_dp ]
      real(dp), parameter :: wgk(8) = [ &
         0.022935322010529224963732008058970_dp, &
         0.063092092629978553290700663189204_dp, &
         0.104790010322250183839876322541518_dp, &
         0.140653259715525918745189590510238_dp, &
         0.169004726639267902826583426598550_dp, &
         0.190350578064785409913256402421014_dp, &
         0.204432940075298892414161999234649_dp, &
         0.209482141084727828012999174891714_dp ]
      real(dp), parameter :: wg(4) = [ &
         0.129484966168869693270611432679082_dp, &
         0.279705391489276667901467771423780_dp, &
         0.381830050505118944950369775488975_dp, &
         0.417959183673469387755102040816327_dp ]
      real(dp) :: center, half, fc, f1v, f2v, resg
      integer :: j

      center = 0.5_dp*(a+b)
      half = 0.5_dp*(b-a)
      fc = skew_theta_value(center,lambda,delta,delta2,family1,family2)
      resk = wgk(8)*fc
      resg = wg(4)*fc
      do j = 1, 7
         f1v = skew_theta_value(center-half*xgk(j),lambda,delta,delta2, &
            family1,family2)
         f2v = skew_theta_value(center+half*xgk(j),lambda,delta,delta2, &
            family1,family2)
         resk = resk+wgk(j)*(f1v+f2v)
      end do
      resg = resg+wg(1)*( &
         skew_theta_value(center-half*xgk(2),lambda,delta,delta2,family1,family2) &
         + skew_theta_value(center+half*xgk(2),lambda,delta,delta2,family1,family2))
      resg = resg+wg(2)*( &
         skew_theta_value(center-half*xgk(4),lambda,delta,delta2,family1,family2) &
         + skew_theta_value(center+half*xgk(4),lambda,delta,delta2,family1,family2))
      resg = resg+wg(3)*( &
         skew_theta_value(center-half*xgk(6),lambda,delta,delta2,family1,family2) &
         + skew_theta_value(center+half*xgk(6),lambda,delta,delta2,family1,family2))
      resk = resk*half
      resg = resg*half
      abserr = abs(resk-resg)
   end subroutine skew_gk15

   real(dp) function skew_theta_value(theta, lambda, delta, delta2, &
      family1, family2) result(y)
      real(dp), intent(in) :: theta, lambda, delta, delta2
      integer, intent(in) :: family1, family2
      real(dp) :: t
      t = sin(theta)**2
      y = dskewunit(t,lambda,delta,delta2,family1,family2) &
          *2.0_dp*sin(theta)*cos(theta)
   end function skew_theta_value

   subroutine pskewunit_vec(x, p, lambda, delta, delta2, family1, family2, log_p, rel_tol)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: p(size(x))
      real(dp), intent(in), optional :: lambda, delta, delta2
      integer, intent(in), optional :: family1, family2
      logical, intent(in), optional :: log_p
      real(dp), intent(in), optional :: rel_tol
      integer :: i
      do i = 1, size(x)
         p(i) = pskewunit(x(i),lambda,delta,delta2,family1,family2,log_p,rel_tol)
      end do
   end subroutine pskewunit_vec

   real(dp) function rskewunit(lambda, delta, delta2, family1, family2) result(x)
      real(dp), intent(in), optional :: lambda, delta, delta2
      integer, intent(in), optional :: family1, family2
      real(dp) :: lam, d1, d2, df, dg, y, prob
      integer :: f1, f2

      lam = 0.0_dp
      d1 = 1.0_dp
      d2 = 1.0_dp
      f1 = family_asin
      f2 = family_asin
      if (present(lambda)) lam = lambda
      if (present(delta)) d1 = delta
      if (present(delta2)) d2 = delta2
      if (present(family1)) f1 = family1
      if (present(family2)) f2 = family2

      if (f2 == family_none) then
         x = baseline_random(f1,d1)
         return
      end if

      df = 1.0_dp
      if (family_has_delta(f1)) df = d1
      dg = 1.0_dp
      if (family_has_delta(f2)) then
         if (family_has_delta(f1)) then
            dg = d2
         else
            dg = d1
         end if
      end if

      do
         y = baseline_random(f1,df)
         prob = baseline_cdf(lam*(y-0.5_dp)+0.5_dp,f2,dg)
         ! Upstream accepts with prob/2; using prob gives the same accepted
         ! distribution with twice the acceptance rate.
         if (rand_uniform() <= prob) exit
      end do
      x = y
   end function rskewunit

   subroutine rskewunit_vec(n, x, lambda, delta, delta2, family1, family2)
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n)
      real(dp), intent(in), optional :: lambda, delta, delta2
      integer, intent(in), optional :: family1, family2
      integer :: i
      do i = 1, n
         x(i) = rskewunit(lambda,delta,delta2,family1,family2)
      end do
   end subroutine rskewunit_vec

end module skewunit_distributions
