! SPDX-License-Identifier: MIT
module zero_one_distributions
   use zero_one_kinds, only : dp, pi
   use zero_one_special, only : beta_pdf, beta_cdf, normal_pdf, normal_cdf, normal_quantile, quiet_nan
   implicit none
   private
   public :: dber, pber, qber, rber
   public :: dber2, pber2, qber2, rber2, ber2_transform
   public :: duhlg, puhlg, quhlg, ruhlg
   public :: dumb, pumb, qumb, rumb
   public :: duphn, puphn, quphn, ruphn
   interface rber
      module procedure rber_scalar, rber_vector
   end interface rber
   interface rber2
      module procedure rber2_scalar, rber2_vector
   end interface rber2
   interface ruhlg
      module procedure ruhlg_scalar, ruhlg_vector
   end interface ruhlg
   interface rumb
      module procedure rumb_scalar, rumb_vector
   end interface rumb
   interface ruphn
      module procedure ruphn_scalar, ruphn_vector
   end interface ruphn
contains
   elemental real(dp) function dber(x, mu, sigma, nu, log_value) result(v)
      real(dp), intent(in) :: x, mu, sigma, nu
      logical, intent(in), optional :: log_value
      logical :: want_log
      real(dp) :: den
      want_log = .false.
      if (present(log_value)) want_log = log_value
      if (mu <= 0.0_dp .or. mu >= 1.0_dp .or. sigma <= 0.0_dp .or. nu < 0.0_dp .or. nu > 1.0_dp) then
         v = quiet_nan()
         return
      end if
      if (x < 0.0_dp .or. x > 1.0_dp) then
         den = 0.0_dp
      else
         den = nu+(1.0_dp-nu)*beta_pdf(x,mu*sigma,(1.0_dp-mu)*sigma)
      end if
      if (want_log) then
         if (den > 0.0_dp) then
            v = log(den)
         else
            v = -huge(1.0_dp)
         end if
      else
         v = den
      end if
   end function dber

   elemental real(dp) function pber(q, mu, sigma, nu, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q, mu, sigma, nu
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (mu <= 0.0_dp .or. mu >= 1.0_dp .or. sigma <= 0.0_dp .or. nu < 0.0_dp .or. nu > 1.0_dp) then
         v = quiet_nan()
         return
      end if
      if (q <= 0.0_dp) then
         v = 0.0_dp
      else if (q >= 1.0_dp) then
         v = 1.0_dp
      else
         v = nu*q+(1.0_dp-nu)*beta_cdf(q,mu*sigma,(1.0_dp-mu)*sigma)
      end if
      if (.not.lower) v = 1.0_dp-v
      if (lp) then
         if (v > 0.0_dp) then
            v = log(v)
         else
            v = -huge(1.0_dp)
         end if
      end if
   end function pber

   elemental real(dp) function qber(p, mu, sigma, nu, lower_tail, log_p) result(x)
      real(dp), intent(in) :: p, mu, sigma, nu
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      real(dp) :: prob, lo, hi, f, den, cand
      integer :: it
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      prob = p
      if (lp) prob = exp(prob)
      if (.not.lower) prob = 1.0_dp-prob
      if (prob < 0.0_dp .or. prob > 1.0_dp .or. mu <= 0.0_dp .or. mu >= 1.0_dp .or. &
          sigma <= 0.0_dp .or. nu < 0.0_dp .or. nu > 1.0_dp) then
         x = quiet_nan()
         return
      end if
      if (prob <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (prob >= 1.0_dp) then
         x = 1.0_dp
         return
      end if
      lo = 0.0_dp
      hi = 1.0_dp
      x = min(1.0_dp,max(0.0_dp,mu))
      do it = 1, 120
         f = pber(x,mu,sigma,nu)-prob
         if (f < 0.0_dp) then
            lo = x
         else
            hi = x
         end if
         den = dber(x,mu,sigma,nu)
         if (den > 1.0e-300_dp) then
            cand = x-f/den
         else
            cand = 0.5_dp*(lo+hi)
         end if
         if (cand <= lo .or. cand >= hi) cand = 0.5_dp*(lo+hi)
         x = cand
         if (hi-lo < 2.0e-14_dp) exit
      end do
   end function qber

   subroutine rber_scalar(x, mu, sigma, nu)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu, sigma, nu
      real(dp) :: u(size(x))
      integer :: i
      call random_number(u)
      do i = 1, size(x)
         x(i) = qber(u(i),mu,sigma,nu)
      end do
   end subroutine rber_scalar

   subroutine rber_vector(x, mu, sigma, nu)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu(:), sigma(:), nu(:)
      real(dp) :: u(size(x))
      integer :: i
      if (size(mu)/=size(x) .or. size(sigma)/=size(x) .or. size(nu)/=size(x)) then
         x = quiet_nan()
         return
      end if
      call random_number(u)
      do i = 1, size(x)
         x(i) = qber(u(i),mu(i),sigma(i),nu(i))
      end do
   end subroutine rber_vector

   pure subroutine ber2_transform(gamma, alpha, base_mu, theta)
      real(dp), intent(in) :: gamma, alpha
      real(dp), intent(out) :: base_mu, theta
      theta = alpha*(1.0_dp-abs(2.0_dp*gamma-1.0_dp))
      if (abs(1.0_dp-theta) <= 16.0_dp*epsilon(1.0_dp)) then
         base_mu = 0.5_dp
      else
         base_mu = (gamma-0.5_dp*theta)/(1.0_dp-theta)
      end if
   end subroutine ber2_transform

   real(dp) function dber2(x, mu, sigma, nu, log_value) result(v)
      real(dp), intent(in) :: x, mu, sigma, nu
      logical, intent(in), optional :: log_value
      real(dp) :: base_mu, theta
      logical :: want_log
      want_log = .false.
      if (present(log_value)) want_log = log_value
      if (mu < 0.0_dp .or. mu > 1.0_dp .or. sigma <= 0.0_dp .or. nu < 0.0_dp .or. nu > 1.0_dp) then
         v = quiet_nan()
         return
      end if
      call ber2_transform(mu,nu,base_mu,theta)
      if (theta >= 1.0_dp-16.0_dp*epsilon(1.0_dp)) then
         if (x >= 0.0_dp .and. x <= 1.0_dp) then
            v = 1.0_dp
         else
            v = 0.0_dp
         end if
         if (want_log) then
            if (v > 0.0_dp) then
               v = 0.0_dp
            else
               v = -huge(1.0_dp)
            end if
         end if
      else
         v = dber(x,base_mu,sigma,theta,want_log)
      end if
   end function dber2

   real(dp) function pber2(q, mu, sigma, nu, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q, mu, sigma, nu
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: base_mu, theta
      logical :: lower, lp
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (mu < 0.0_dp .or. mu > 1.0_dp .or. sigma <= 0.0_dp .or. nu < 0.0_dp .or. nu > 1.0_dp) then
         v = quiet_nan()
         return
      end if
      call ber2_transform(mu,nu,base_mu,theta)
      if (theta >= 1.0_dp-16.0_dp*epsilon(1.0_dp)) then
         v = min(1.0_dp,max(0.0_dp,q))
         if (.not.lower) v = 1.0_dp-v
         if (lp) then
            if (v > 0.0_dp) then
               v = log(v)
            else
               v = -huge(1.0_dp)
            end if
         end if
      else
         v = pber(q,base_mu,sigma,theta,lower,lp)
      end if
   end function pber2

   elemental real(dp) function qber2(p, mu, sigma, nu, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p, mu, sigma, nu
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: base_mu, theta, prob
      logical :: lower, lp
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (mu < 0.0_dp .or. mu > 1.0_dp .or. sigma <= 0.0_dp .or. nu < 0.0_dp .or. nu > 1.0_dp) then
         v = quiet_nan()
         return
      end if
      call ber2_transform(mu,nu,base_mu,theta)
      if (theta >= 1.0_dp-16.0_dp*epsilon(1.0_dp)) then
         prob = p
         if (lp) prob = exp(prob)
         if (.not.lower) prob = 1.0_dp-prob
         v = prob
      else
         v = qber(p,base_mu,sigma,theta,lower,lp)
      end if
   end function qber2

   subroutine rber2_scalar(x, mu, sigma, nu)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu, sigma, nu
      real(dp) :: u(size(x))
      integer :: i
      call random_number(u)
      do i = 1, size(x)
         x(i) = qber2(u(i),mu,sigma,nu)
      end do
   end subroutine rber2_scalar

   subroutine rber2_vector(x, mu, sigma, nu)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu(:), sigma(:), nu(:)
      real(dp) :: u(size(x))
      integer :: i
      if (size(mu)/=size(x) .or. size(sigma)/=size(x) .or. size(nu)/=size(x)) then
         x = quiet_nan()
         return
      end if
      call random_number(u)
      do i = 1, size(x)
         x(i) = qber2(u(i),mu(i),sigma(i),nu(i))
      end do
   end subroutine rber2_vector

   elemental real(dp) function duhlg(x, mu, log_value) result(v)
      real(dp), intent(in) :: x, mu
      logical, intent(in), optional :: log_value
      logical :: want_log
      real(dp) :: lp
      want_log = .false.
      if (present(log_value)) want_log = log_value
      if (mu <= 0.0_dp .or. x <= 0.0_dp .or. x >= 1.0_dp) then
         v = quiet_nan()
         return
      end if
      lp = log(2.0_dp*mu)-2.0_dp*log(mu+(2.0_dp-mu)*x)
      if (want_log) then
         v = lp
      else
         v = exp(lp)
      end if
   end function duhlg

   elemental real(dp) function puhlg(q, mu, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q, mu
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (mu <= 0.0_dp) then
         v = quiet_nan()
         return
      end if
      if (q <= 0.0_dp) then
         v = 0.0_dp
      else if (q >= 1.0_dp) then
         v = 1.0_dp
      else
         v = 1.0_dp-mu*(1.0_dp-q)/(mu+(2.0_dp-mu)*q)
      end if
      if (.not.lower) v = 1.0_dp-v
      if (lp) then
         if (v > 0.0_dp) then
            v = log(v)
         else
            v = -huge(1.0_dp)
         end if
      end if
   end function puhlg

   elemental real(dp) function quhlg(p, mu, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p, mu
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      real(dp) :: prob
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      prob = p
      if (lp) prob = exp(prob)
      if (.not.lower) prob = 1.0_dp-prob
      if (mu <= 0.0_dp .or. prob < 0.0_dp .or. prob > 1.0_dp) then
         v = quiet_nan()
      else
         v = prob*mu/(2.0_dp-2.0_dp*prob+prob*mu)
      end if
   end function quhlg

   subroutine ruhlg_scalar(x, mu)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu
      real(dp) :: u(size(x))
      call random_number(u)
      x = quhlg(u,mu)
   end subroutine ruhlg_scalar

   subroutine ruhlg_vector(x, mu)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu(:)
      real(dp) :: u(size(x))
      if (size(mu)/=size(x)) then
         x = quiet_nan()
         return
      end if
      call random_number(u)
      x = quhlg(u,mu)
   end subroutine ruhlg_vector

   elemental real(dp) function dumb(x, mu, log_value) result(v)
      real(dp), intent(in) :: x, mu
      logical, intent(in), optional :: log_value
      logical :: want_log
      real(dp) :: l, lp
      want_log = .false.
      if (present(log_value)) want_log = log_value
      if (mu <= 0.0_dp) then
         v = quiet_nan()
         return
      end if
      if (x <= 0.0_dp .or. x >= 1.0_dp) then
         if (want_log) then
            v = -huge(1.0_dp)
         else
            v = 0.0_dp
         end if
         return
      end if
      l = log(1.0_dp/x)
      lp = 0.5_dp*log(2.0_dp/pi)+2.0_dp*log(l)-3.0_dp*log(mu)-log(x)-l*l/(2.0_dp*mu*mu)
      if (want_log) then
         v = lp
      else
         v = exp(lp)
      end if
   end function dumb

   elemental real(dp) function pumb(q, mu, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q, mu
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      real(dp) :: l, z
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (mu <= 0.0_dp) then
         v = quiet_nan()
         return
      end if
      if (q <= 0.0_dp) then
         v = 0.0_dp
      else if (q >= 1.0_dp) then
         v = 1.0_dp
      else
         l = log(1.0_dp/q)
         z = l/(sqrt(2.0_dp)*mu)
         v = erfc(z)+sqrt(2.0_dp/pi)*(l/mu)*exp(-0.5_dp*(l/mu)**2)
         v = min(1.0_dp,max(0.0_dp,v))
      end if
      if (.not.lower) v = 1.0_dp-v
      if (lp) then
         if (v > 0.0_dp) then
            v = log(v)
         else
            v = -huge(1.0_dp)
         end if
      end if
   end function pumb

   elemental real(dp) function qumb(p, mu, lower_tail, log_p) result(x)
      real(dp), intent(in) :: p, mu
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      real(dp) :: prob, lo, hi, mid
      integer :: it
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      prob = p
      if (lp) prob = exp(prob)
      if (.not.lower) prob = 1.0_dp-prob
      if (mu <= 0.0_dp .or. prob < 0.0_dp .or. prob > 1.0_dp) then
         x = quiet_nan()
         return
      end if
      if (prob <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (prob >= 1.0_dp) then
         x = 1.0_dp
         return
      end if
      lo = 0.0_dp
      hi = 1.0_dp
      do it = 1, 120
         mid = 0.5_dp*(lo+hi)
         if (pumb(mid,mu) < prob) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function qumb

   subroutine rumb_scalar(x, mu)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu
      real(dp) :: u(size(x))
      integer :: i
      call random_number(u)
      do i = 1, size(x)
         x(i) = qumb(u(i),mu)
      end do
   end subroutine rumb_scalar

   subroutine rumb_vector(x, mu)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu(:)
      real(dp) :: u(size(x))
      integer :: i
      if (size(mu)/=size(x)) then
         x = quiet_nan()
         return
      end if
      call random_number(u)
      do i = 1, size(x)
         x(i) = qumb(u(i),mu(i))
      end do
   end subroutine rumb_vector

   elemental real(dp) function duphn(x, mu, sigma, log_value) result(v)
      real(dp), intent(in) :: x, mu, sigma
      logical, intent(in), optional :: log_value
      logical :: want_log
      real(dp) :: z, h, lp
      want_log = .false.
      if (present(log_value)) want_log = log_value
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. x <= 0.0_dp .or. x >= 1.0_dp) then
         v = quiet_nan()
         return
      end if
      z = (1.0_dp-x)/(sigma*x)
      h = erf(z/sqrt(2.0_dp))
      if (h <= 0.0_dp) then
         lp = -huge(1.0_dp)
      else
         lp = log(2.0_dp*mu)-log(sigma)-2.0_dp*log(x)+log(normal_pdf(z))+(mu-1.0_dp)*log(h)
      end if
      if (want_log) then
         v = lp
      else
         v = exp(lp)
      end if
   end function duphn

   elemental real(dp) function puphn(q, mu, sigma, lower_tail, log_p) result(v)
      real(dp), intent(in) :: q, mu, sigma
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      real(dp) :: z, h
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp) then
         v = quiet_nan()
         return
      end if
      if (q <= 0.0_dp) then
         v = 0.0_dp
      else if (q >= 1.0_dp) then
         v = 1.0_dp
      else
         z = (1.0_dp-q)/(sigma*q)
         h = erf(z/sqrt(2.0_dp))
         v = 1.0_dp-h**mu
      end if
      if (.not.lower) v = 1.0_dp-v
      if (lp) then
         if (v > 0.0_dp) then
            v = log(v)
         else
            v = -huge(1.0_dp)
         end if
      end if
   end function puphn

   elemental real(dp) function quphn(p, mu, sigma, lower_tail, log_p) result(v)
      real(dp), intent(in) :: p, mu, sigma
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lp
      real(dp) :: prob, inside, z
      lower = .true.
      lp = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) lp = log_p
      prob = p
      if (lp) prob = exp(prob)
      if (.not.lower) prob = 1.0_dp-prob
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. prob <= 0.0_dp .or. prob >= 1.0_dp) then
         v = quiet_nan()
         return
      end if
      inside = 0.5_dp*(1.0_dp+(1.0_dp-prob)**(1.0_dp/mu))
      z = normal_quantile(inside)
      v = 1.0_dp/(sigma*z+1.0_dp)
   end function quphn

   subroutine ruphn_scalar(x, mu, sigma)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu, sigma
      real(dp) :: u(size(x))
      integer :: i
      call random_number(u)
      do i = 1, size(x)
         x(i) = quphn(u(i),mu,sigma)
      end do
   end subroutine ruphn_scalar
   subroutine ruphn_vector(x, mu, sigma)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: mu(:), sigma(:)
      real(dp) :: u(size(x))
      if (size(mu)/=size(x) .or. size(sigma)/=size(x)) then
         x = quiet_nan()
         return
      end if
      call random_number(u)
      x = quphn(u,mu,sigma)
   end subroutine ruphn_vector
end module zero_one_distributions
