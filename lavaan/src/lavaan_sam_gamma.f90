module lavaan_sam_gamma
   use lavaan_kinds, only : dp
   implicit none
   private
   public :: sam_continuous_gamma, sam_browne_unbiased_gamma
contains

   subroutine sam_continuous_gamma(data, gamma, unbiased, meanstructure, status)
      real(dp), intent(in) :: data(:, :)
      real(dp), allocatable, intent(out) :: gamma(:, :)
      logical, intent(in), optional :: unbiased, meanstructure
      integer, intent(out), optional :: status
      real(dp), allocatable :: mu(:), s(:, :), res(:, :), z(:, :), zm(:), gcov(:, :), cross(:, :)
      real(dp), allocatable :: gnt(:, :), sv(:), gu(:, :)
      logical :: ub, means
      integer :: n, p, q, r, a, b, c, d, k, l
      real(dp) :: scale1, scale2

      n = size(data, 1)
      p = size(data, 2)
      ub = .false.
      if (present(unbiased)) ub = unbiased
      means = .true.
      if (present(meanstructure)) means = meanstructure
      if (present(status)) status = 0
      if (n < 4 .or. p < 1) then
         allocate(gamma(0, 0))
         if (present(status)) status = -1
         return
      end if
      allocate(mu(p), s(p, p), res(n, p))
      mu = sum(data, dim=1) / real(n, dp)
      do r = 1, n
      res(r, :) = data(r, :) - mu
      end do
      s = matmul(transpose(res), res) / real(n - 1, dp)
      q = p * (p + 1) / 2
      allocate(z(n, q))
      k = 0
      do b = 1, p
         do a = b, p
            k = k + 1
            do r = 1, n
            z(r, k) = res(r, a) * res(r, b)
            end do
         end do
      end do
      zm = sum(z, dim=1) / real(n, dp)
      allocate(gcov(q, q))
      gcov = 0.0_dp
      do r = 1, n
         gcov = gcov + spread(z(r, :) - zm, 2, q) * spread(z(r, :) - zm, 1, q)
      end do
      gcov = gcov / real(n, dp)
      allocate(cross(p, q))
      cross = 0.0_dp
      do r = 1, n
         cross = cross + spread(res(r, :), 2, q) * spread(z(r, :) - zm, 1, p)
      end do
      cross = cross / real(n, dp)

      if (ub) then
         allocate(gnt(q, q), sv(q))
         gnt = 0.0_dp
         k = 0
         do b = 1, p
            do a = b, p
               k = k + 1
               sv(k) = s(a, b)
               l = 0
               do d = 1, p
                  do c = d, p
                     l = l + 1
                     gnt(k, l) = s(a, c) * s(b, d) + s(a, d) * s(b, c)
                  end do
               end do
            end do
         end do
         scale1 = real(n * (n - 1), dp) / real((n - 2) * (n - 3), dp)
         scale2 = real(n, dp) / real((n - 2) * (n - 3), dp)
         gu = scale1 * gcov - scale2 * (gnt - 2.0_dp / real(n - 1, dp) * &
              spread(sv, 2, q) * spread(sv, 1, q))
         gcov = 0.5_dp * (gu + transpose(gu))
         cross = cross * real(n, dp) / real(n - 2, dp)
      end if

      if (means) then
         allocate(gamma(p + q, p + q))
         gamma = 0.0_dp
         gamma(1:p, 1:p) = s
         gamma(1:p, p + 1:p + q) = cross
         gamma(p + 1:p + q, 1:p) = transpose(cross)
         gamma(p + 1:p + q, p + 1:p + q) = gcov
      else
         gamma = gcov
      end if
   end subroutine sam_continuous_gamma

   subroutine sam_browne_unbiased_gamma(data, gamma, status, meanstructure)
      real(dp), intent(in) :: data(:, :)
      real(dp), allocatable, intent(out) :: gamma(:, :)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: meanstructure
      logical :: ms
      ms = .true.
      if (present(meanstructure)) ms = meanstructure
      call sam_continuous_gamma(data, gamma, unbiased=.true., meanstructure=ms, status=status)
   end subroutine sam_browne_unbiased_gamma

end module lavaan_sam_gamma
