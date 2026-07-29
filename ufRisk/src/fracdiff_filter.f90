! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_filter
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fracdiff_kinds, only : dp
   use fracdiff_status, only : fd_ok, fd_invalid_input, fd_gamma_error
   implicit none
   private

   public :: haslett_raftery_filter
   public :: arma_residuals, arma_residual_jacobian
   public :: conditional_arma_series_residuals

contains

   subroutine haslett_raftery_filter(x, d, m_terms, has_arma, y, sum_log_v, estimated_mean, status)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: d
      integer, intent(in) :: m_terms
      logical, intent(in) :: has_arma
      real(dp), intent(out) :: y(:)
      real(dp), intent(out) :: sum_log_v
      real(dp), intent(out) :: estimated_mean
      integer, intent(out) :: status

      real(dp), allocatable :: amk(:), ak(:), vk(:), phi(:), pi_weights(:)
      real(dp) :: r, s, t, u, v, z, g0, denominator
      integer :: n, mcap, j, k, km

      n = size(x)
      status = fd_ok
      sum_log_v = 0.0_dp
      estimated_mean = 0.0_dp
      if (n < 2 .or. size(y) /= n .or. m_terms < 1 .or. d >= 0.5_dp .or. d <= -0.5_dp) then
         status = fd_invalid_input
         y = 0.0_dp
         return
      end if

      mcap = min(m_terms, n)
      allocate(amk(n), ak(n), vk(mcap), phi(max(1,mcap)), pi_weights(max(1,mcap)))
      amk = 0.0_dp
      ak = 0.0_dp
      vk = 0.0_dp
      phi = 0.0_dp
      pi_weights = 0.0_dp

      amk(1) = 0.0_dp
      ak(1) = 1.0_dp
      g0 = gamma(1.0_dp - 2.0_dp*d)/(gamma(1.0_dp - d)**2)
      if (.not. ieee_is_finite(g0) .or. g0 <= 0.0_dp) then
         status = fd_gamma_error
         y = 0.0_dp
         return
      end if
      vk(1) = g0

      if (n >= 2) then
         z = d/(1.0_dp - d)
         amk(2) = z*x(1)
         ak(2) = 1.0_dp - z
         phi(1) = z
         if (mcap >= 2) vk(2) = g0*(1.0_dp - z*z)
      end if

      do k = 3, mcap
         km = k - 1
         t = real(km, dp)
         u = t - d
         do j = 1, km - 1
            s = t - real(j, dp)
            phi(j) = phi(j)*t*(s - d)/(u*s)
         end do
         v = d/u
         phi(km) = v
         vk(k) = vk(km)*(1.0_dp - v*v)
         if (vk(k) <= 0.0_dp .or. .not. ieee_is_finite(vk(k))) then
            status = fd_gamma_error
            y = 0.0_dp
            return
         end if
         u = 0.0_dp
         v = 1.0_dp
         do j = 1, km
            t = phi(j)
            u = u + t*x(k - j)
            v = v - t
         end do
         amk(k) = u
         ak(k) = v
      end do

      if (m_terms < n) then
         pi_weights(1) = d
         s = d
         do j = 2, mcap
            u = real(j, dp)
            t = pi_weights(j - 1)*(u - 1.0_dp - d)/u
            s = s + t
            pi_weights(j) = t
         end do
         s = 1.0_dp - s
         r = 0.0_dp
         u = real(mcap, dp)
         t = u*pi_weights(mcap)
         do k = mcap + 1, n
            km = k - mcap
            z = 0.0_dp
            do j = 1, mcap
               z = z + pi_weights(j)*x(k - j)
            end do
            if (abs(r) <= tiny(1.0_dp)) then
               amk(k) = z
               ak(k) = s
            else
               if (abs(d) <= sqrt(epsilon(1.0_dp))) then
                  v = -t*log(u/real(k,dp))
               else
                  v = t*(1.0_dp - (u/real(k,dp))**d)/d
               end if
               denominator = real(km, dp) - 1.0_dp
               amk(k) = z + v*r/denominator
               ak(k) = s - v
            end if
            r = r + x(km)
         end do
      end if

      r = 0.0_dp
      s = 0.0_dp
      z = vk(mcap)
      do k = 1, n
         t = ak(k)
         u = (x(k) - amk(k))*t
         v = t*t
         if (k <= mcap) then
            z = vk(k)
            u = u/z
            v = v/z
         end if
         r = r + u
         s = s + v
      end do
      if (s <= 0.0_dp) then
         status = fd_gamma_error
         y = 0.0_dp
         return
      end if
      estimated_mean = r/s

      sum_log_v = sum(log(vk))
      do k = 1, n
         t = x(k) - amk(k) - estimated_mean*ak(k)
         if (k <= mcap) t = t/sqrt(vk(k))
         y(k) = t
      end do

      if (has_arma) y = y - z/real(n, dp)
   end subroutine haslett_raftery_filter

   subroutine arma_residuals(y, ar, ma, residuals, status)
      real(dp), intent(in) :: y(:), ar(:), ma(:)
      real(dp), intent(out) :: residuals(:)
      integer, intent(out) :: status

      integer :: n, p, q, maxpq, nm, k, km, l
      real(dp) :: value

      n = size(y)
      p = size(ar)
      q = size(ma)
      maxpq = max(p, q)
      nm = n - maxpq
      status = fd_ok
      if (nm < 1 .or. size(residuals) /= nm) then
         status = fd_invalid_input
         residuals = 0.0_dp
         return
      end if

      residuals = 0.0_dp
      do k = maxpq + 1, n
         km = k - maxpq
         value = y(k)
         do l = 1, p
            value = value - ar(l)*y(k - l)
         end do
         do l = 1, q
            if (km <= l) exit
            value = value + ma(l)*residuals(km - l)
         end do
         residuals(km) = value
      end do
   end subroutine arma_residuals

   subroutine arma_residual_jacobian(y, ar, ma, residuals, jacobian, status)
      real(dp), intent(in) :: y(:), ar(:), ma(:)
      real(dp), intent(out) :: residuals(:)
      real(dp), intent(out) :: jacobian(:,:)
      integer, intent(out) :: status

      integer :: n, p, q, pq, maxpq, nm, i, k, km, l
      real(dp) :: value, propagated

      n = size(y)
      p = size(ar)
      q = size(ma)
      pq = p + q
      maxpq = max(p, q)
      nm = n - maxpq
      status = fd_ok
      if (nm < 1 .or. size(residuals) /= nm .or. size(jacobian,1) /= nm .or. &
          size(jacobian,2) /= pq) then
         status = fd_invalid_input
         residuals = 0.0_dp
         jacobian = 0.0_dp
         return
      end if

      residuals = 0.0_dp
      jacobian = 0.0_dp
      do k = maxpq + 1, n
         km = k - maxpq
         value = y(k)
         do l = 1, p
            value = value - ar(l)*y(k - l)
         end do
         do l = 1, q
            if (km <= l) exit
            value = value + ma(l)*residuals(km - l)
         end do
         residuals(km) = value
      end do

      do i = 1, pq
         do k = maxpq + 1, n
            km = k - maxpq
            propagated = 0.0_dp
            do l = 1, q
               if (km <= l) exit
               propagated = propagated + ma(l)*jacobian(km - l, i)
            end do
            if (i <= q) then
               if (km > i) propagated = propagated + residuals(km - i)
            else
               propagated = propagated - y(k - (i - q))
            end if
            jacobian(km,i) = propagated
         end do
      end do
   end subroutine arma_residual_jacobian

   subroutine conditional_arma_series_residuals(x, ar, ma, residuals, fitted, status)
      real(dp), intent(in) :: x(:), ar(:), ma(:)
      real(dp), intent(out) :: residuals(:), fitted(:)
      integer, intent(out) :: status

      integer :: n, p, q, t, lag
      real(dp) :: prediction

      n = size(x)
      p = size(ar)
      q = size(ma)
      status = fd_ok
      if (size(residuals) /= n .or. size(fitted) /= n) then
         status = fd_invalid_input
         return
      end if

      residuals = 0.0_dp
      fitted = 0.0_dp
      do t = 1, n
         prediction = 0.0_dp
         do lag = 1, min(p, t - 1)
            prediction = prediction + ar(lag)*x(t - lag)
         end do
         do lag = 1, min(q, t - 1)
            prediction = prediction - ma(lag)*residuals(t - lag)
         end do
         fitted(t) = prediction
         residuals(t) = x(t) - prediction
      end do
   end subroutine conditional_arma_series_residuals

end module fracdiff_filter
