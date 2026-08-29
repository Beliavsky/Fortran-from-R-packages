module lmtest_pan
   use lmtest_kinds, only : dp
   implicit none
   private
   public :: pan_probability

contains

   function pan_probability(aeig, x, c, nterms) result(prob)
      ! Modern free-format translation of the amended Applied Statistics
      ! Algorithm AS 153 (AS R52), Farebrother (1984), as distributed in
      ! lmtest src/pan.f. The original comments and source are retained under
      ! upstream/src/pan.f.
      real(dp), intent(in) :: aeig(:), x, c
      integer, intent(in), optional :: nterms
      real(dp) :: prob
      real(dp), allocatable :: a(:)
      integer :: m, n, d, h, i, j1, j2, j3, j4, k, l1, l2, nu, n2
      real(dp) :: num, pin, prod, sgn, sumv, sum1, u, v, yv

      m = size(aeig)
      n = 15
      if (present(nterms)) n = nterms
      if (m < 1 .or. n < 1) then
         prob = 0.0_dp
         return
      end if
      allocate(a(0:m))
      a(0) = x
      a(1:m) = aeig

      if (a(1) > a(m)) then
         h = m
         k = -1
         i = 1
      else
         h = 1
         k = 1
         i = m
      end if
      nu = h
      do
         if (a(nu) >= x) exit
         if (nu == i) then
            if (c >= 0.0_dp) then
               prob = 1.0_dp
            else
               prob = 0.0_dp
            end if
            return
         end if
         nu = nu + k
      end do

      if (nu == h .and. c <= 0.0_dp) then
         prob = 0.0_dp
         return
      end if
      if (k == 1) nu = nu - 1
      h = m - nu
      if (abs(c) <= tiny(1.0_dp)) then
         yv = real(h-nu,dp)
      else
         yv = c * (a(1)-a(m))
      end if

      if (yv >= 0.0_dp) then
         d = 2
         h = nu
         k = -k
         j1 = 0
         j2 = 2
         j3 = 3
         j4 = 1
      else
         d = -2
         nu = nu + 1
         j1 = m - 2
         j2 = m - 1
         j3 = m + 1
         j4 = m
      end if
      pin = 2.0_dp * atan(1.0_dp) / real(n,dp)
      sumv = 0.5_dp * real(k+1,dp)
      sgn = real(k,dp) / real(n,dp)
      n2 = 2*n - 1

      do l1 = h-2*(h/2), 0, -1
         l2 = j2
         do
            if ((d > 0 .and. l2 > nu) .or. (d < 0 .and. l2 < nu)) exit
            sum1 = a(j4)
            prod = a(l2)
            u = 0.5_dp * (sum1+prod)
            v = 0.5_dp * (sum1-prod)
            sum1 = 0.0_dp
            do i = 1, n2, 2
               yv = u - v*cos(real(i,dp)*pin)
               num = yv - x
               if (abs(num) <= tiny(1.0_dp)) cycle
               prod = exp(-c/num)
               do k = 1, j1
                  prod = prod * num/(yv-a(k))
               end do
               do k = j3, m
                  prod = prod * num/(yv-a(k))
               end do
               sum1 = sum1 + sqrt(abs(prod))
            end do
            sgn = -sgn
            sumv = sumv + sgn*sum1
            j1 = j1 + d
            j3 = j3 + d
            j4 = j4 + d
            l2 = l2 + d
         end do
         if (d == 2) then
            j3 = j3 - 1
         else
            j1 = j1 + 1
         end if
         j2 = 0
         nu = 0
      end do
      prob = sumv
   end function pan_probability

end module lmtest_pan
