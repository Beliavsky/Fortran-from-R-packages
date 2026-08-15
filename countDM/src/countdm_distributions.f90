module countdm_distributions
   use countdm_kinds, only: dp
   use countdm_math, only: logsumexp2
   implicit none
   private
   public :: touchard_polynomial, tp, bell_number
   public :: dbell, dborel, dpoisson_count
   public :: dbellt, pbellt, qbellt, rbellt
   public :: dzibellt, pzibellt, qzibellt, rzibellt
   public :: dzip, dzibell, dzoip, dzoibell

contains

   pure real(dp) function touchard_polynomial(x, theta) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: theta
      real(dp), allocatable :: s(:), next(:)
      integer :: i, k
      if (x < 0 .or. theta < 0.0_dp) then
         v = 0.0_dp
         return
      end if
      if (x == 0) then
         v = 1.0_dp
         return
      end if
      allocate(s(0:x), next(0:x))
      s = 0.0_dp
      s(0) = 1.0_dp
      do i = 1, x
         next = 0.0_dp
         do k = 1, i
            next(k) = real(k, dp) * s(k) + s(k-1)
         end do
         s = next
      end do
      v = 0.0_dp
      do k = 1, x
         v = v + s(k) * theta**k
      end do
   end function touchard_polynomial


   pure real(dp) function tp(x, theta) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: theta
      v = touchard_polynomial(x, theta)
   end function tp

   pure real(dp) function bell_number(x) result(v)
      integer, intent(in) :: x
      v = touchard_polynomial(x, 1.0_dp)
   end function bell_number

   pure real(dp) function dbell(x, theta, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: theta
      logical, intent(in), optional :: log_p
      real(dp) :: lv, b
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      if (x < 0 .or. theta <= 0.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      b = bell_number(x)
      lv = real(x, dp) * log(theta) + 1.0_dp - exp(theta) + log(b) - log_gamma(real(x + 1, dp))
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dbell

   pure real(dp) function dborel(x, alpha, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: alpha
      logical, intent(in), optional :: log_p
      real(dp) :: lv
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      if (x < 1 .or. alpha <= 0.0_dp .or. alpha > 1.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      lv = -alpha * real(x, dp) + real(x - 1, dp) * log(alpha * real(x, dp)) &
         - log_gamma(real(x + 1, dp))
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dborel

   pure real(dp) function dpoisson_count(x, theta, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: theta
      logical, intent(in), optional :: log_p
      real(dp) :: lv
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      if (x < 0 .or. theta <= 0.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      lv = real(x, dp) * log(theta) - theta - log_gamma(real(x + 1, dp))
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dpoisson_count

   pure real(dp) function dbellt(x, lambda, theta, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda, theta
      logical, intent(in), optional :: log_p
      real(dp) :: lv, tp
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      if (x < 0 .or. lambda <= 0.0_dp .or. theta <= 0.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      tp = touchard_polynomial(x, theta)
      lv = real(x, dp) * log(lambda) + theta * (1.0_dp - exp(lambda)) &
         + log(tp) - log_gamma(real(x + 1, dp))
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dbellt

   pure real(dp) function pbellt(q, lambda, theta, lower_tail, log_p) result(v)
      integer, intent(in) :: q
      real(dp), intent(in) :: lambda, theta
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, lg
      real(dp) :: c
      integer :: k
      lower = .true.; if (present(lower_tail)) lower = lower_tail
      lg = .false.; if (present(log_p)) lg = log_p
      if (q < 0) then
         c = 0.0_dp
      else
         c = 0.0_dp
         do k = 0, q
            c = c + dbellt(k, lambda, theta)
         end do
         c = min(1.0_dp, max(0.0_dp, c))
      end if
      if (.not. lower) c = 1.0_dp - c
      if (lg) then
         if (c <= 0.0_dp) then; v = -huge(1.0_dp); else; v = log(c); end if
      else
         v = c
      end if
   end function pbellt

   integer function qbellt(p, lambda, theta, lower_tail, log_p) result(q)
      real(dp), intent(in) :: p, lambda, theta
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: pp, c
      logical :: lower, lg
      integer :: k
      lower = .true.; if (present(lower_tail)) lower = lower_tail
      lg = .false.; if (present(log_p)) lg = log_p
      pp = p
      if (lg) pp = exp(pp)
      if (.not. lower) pp = 1.0_dp - pp
      pp = min(1.0_dp, max(0.0_dp, pp))
      if (pp <= 0.0_dp) then
         q = 0
         return
      end if
      c = 0.0_dp
      do k = 0, 100000
         c = c + dbellt(k, lambda, theta)
         if (c >= pp .or. 1.0_dp - c <= 8.0_dp * epsilon(1.0_dp)) then
            q = k
            return
         end if
      end do
      q = 100000
   end function qbellt

   subroutine rbellt(n, lambda, theta, x)
      integer, intent(in) :: n
      real(dp), intent(in) :: lambda, theta
      integer, intent(out) :: x(n)
      real(dp) :: u
      integer :: i
      do i = 1, n
         call random_number(u)
         x(i) = qbellt(u, lambda, theta)
      end do
   end subroutine rbellt

   pure real(dp) function dzibellt(x, lambda, theta, pi0, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: lambda, theta, pi0
      logical, intent(in), optional :: log_p
      real(dp) :: lv
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      if (pi0 < 0.0_dp .or. pi0 > 1.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      lv = log(max(1.0_dp - pi0, tiny(1.0_dp))) + dbellt(x, lambda, theta, .true.)
      if (x == 0) lv = logsumexp2(log(max(pi0, tiny(1.0_dp))), lv)
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dzibellt

   pure real(dp) function pzibellt(q, lambda, theta, pi0, lower_tail, log_p) result(v)
      integer, intent(in) :: q
      real(dp), intent(in) :: lambda, theta, pi0
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: c
      logical :: lower, lg
      lower = .true.; if (present(lower_tail)) lower = lower_tail
      lg = .false.; if (present(log_p)) lg = log_p
      if (q < 0) then
         c = 0.0_dp
      else
         c = pi0 + (1.0_dp - pi0) * pbellt(q, lambda, theta)
      end if
      if (.not. lower) c = 1.0_dp - c
      if (lg) then
         if (c <= 0.0_dp) then; v = -huge(1.0_dp); else; v = log(c); end if
      else
         v = c
      end if
   end function pzibellt

   integer function qzibellt(p, lambda, theta, pi0, lower_tail, log_p) result(q)
      real(dp), intent(in) :: p, lambda, theta, pi0
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: pp
      logical :: lower, lg
      lower = .true.; if (present(lower_tail)) lower = lower_tail
      lg = .false.; if (present(log_p)) lg = log_p
      pp = p
      if (lg) pp = exp(pp)
      if (.not. lower) pp = 1.0_dp - pp
      pp = max(0.0_dp, (pp - pi0) / max(1.0_dp - pi0, tiny(1.0_dp)))
      q = qbellt(pp, lambda, theta)
   end function qzibellt

   subroutine rzibellt(n, lambda, theta, pi0, x)
      integer, intent(in) :: n
      real(dp), intent(in) :: lambda, theta, pi0
      integer, intent(out) :: x(n)
      real(dp) :: u
      integer :: i
      do i = 1, n
         call random_number(u)
         x(i) = qzibellt(u, lambda, theta, pi0)
      end do
   end subroutine rzibellt

   pure real(dp) function dzip(x, alpha, theta, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: alpha, theta
      logical, intent(in), optional :: log_p
      real(dp) :: lv
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      if (x < 0 .or. alpha < 0.0_dp .or. alpha > 1.0_dp .or. theta <= 0.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      lv = log(max(1.0_dp - alpha, tiny(1.0_dp))) + dpoisson_count(x, theta, .true.)
      if (x == 0) lv = logsumexp2(log(max(alpha, tiny(1.0_dp))), lv)
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dzip

   pure real(dp) function dzibell(x, alpha, lambda, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: alpha, lambda
      logical, intent(in), optional :: log_p
      real(dp) :: lv
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      if (x < 0 .or. alpha < 0.0_dp .or. alpha > 1.0_dp .or. lambda <= 0.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      lv = log(max(1.0_dp - alpha, tiny(1.0_dp))) + dbell(x, lambda, .true.)
      if (x == 0) lv = logsumexp2(log(max(alpha, tiny(1.0_dp))), lv)
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dzibell

   pure real(dp) function dzoip(x, alpha, beta, theta, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: alpha, beta, theta
      logical, intent(in), optional :: log_p
      real(dp) :: lv, base
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      base = 1.0_dp - alpha - beta
      if (x < 0 .or. alpha < 0.0_dp .or. beta < 0.0_dp .or. base < 0.0_dp .or. theta <= 0.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      lv = log(max(base, tiny(1.0_dp))) + dpoisson_count(x, theta, .true.)
      if (x == 0) lv = logsumexp2(log(max(alpha, tiny(1.0_dp))), lv)
      if (x == 1) lv = logsumexp2(log(max(beta, tiny(1.0_dp))), lv)
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dzoip

   pure real(dp) function dzoibell(x, alpha, beta, theta, log_p) result(v)
      integer, intent(in) :: x
      real(dp), intent(in) :: alpha, beta, theta
      logical, intent(in), optional :: log_p
      real(dp) :: lv, base, norm
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      base = 1.0_dp - alpha - beta
      if (x < 0 .or. alpha <= 0.0_dp .or. beta <= 0.0_dp .or. base <= 0.0_dp .or. theta <= 0.0_dp) then
         if (lg) then; v = -huge(1.0_dp); else; v = 0.0_dp; end if
         return
      end if
      if (x == 0) then
         lv = log(alpha)
      else if (x == 1) then
         lv = log(beta)
      else
         norm = 1.0_dp - dbell(0, theta) - dbell(1, theta)
         lv = log(base) + dbell(x, theta, .true.) - log(max(norm, tiny(1.0_dp)))
      end if
      if (lg) then; v = lv; else; v = exp(lv); end if
   end function dzoibell

end module countdm_distributions
