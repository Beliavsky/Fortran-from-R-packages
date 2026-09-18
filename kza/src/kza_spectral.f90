! SPDX-License-Identifier: GPL-3.0-only
module kza_spectral
   use kza_kinds, only : dp, pi
   use kza_utils, only : complex_is_nan, quiet_nan, r_round_even
   implicit none
   private

   public :: kzft
   public :: kzs
   public :: kztp
   public :: periodogram
   public :: transfer_function

contains

   pure function kzft(x, f, m, k) result(zout)
      real(dp), intent(in) :: x(:) !! Raw one-dimensional data; NaNs are treated as missing observations in each local transform.
      real(dp), intent(in), optional :: f !! Target frequency in cycles per time unit; defaults to zero.
      real(dp), intent(in), optional :: m !! KZFT window size; defaults to 1 and may be non-integer as in the R implementation.
      integer, intent(in), optional :: k !! Number of KZFT iterations; defaults to 1.
      complex(dp), allocatable :: zout(:)
      complex(dp), allocatable :: next(:)
      complex(dp), allocatable :: work(:)
      complex(dp) :: phase
      complex(dp) :: total
      integer :: count_valid
      integer :: h
      integer :: i
      integer :: iter
      integer :: j
      integer :: left_pad
      integer :: niter
      integer :: right_pad
      integer :: source
      integer :: total_n
      integer :: t
      real(dp) :: freq
      real(dp) :: mval

      freq = 0.0_dp
      if (present(f)) freq = f
      mval = 1.0_dp
      if (present(m)) mval = m
      niter = 1
      if (present(k)) niter = max(k, 0)
      h = max(floor((mval - 1.0_dp) / 2.0_dp), 0)
      t = max(ceiling((mval - 1.0_dp) / 2.0_dp), 0)
      left_pad = h * niter
      right_pad = t * niter
      total_n = size(x) + left_pad + right_pad

      allocate(work(total_n), next(total_n))
      work = cmplx(quiet_nan(), 0.0_dp, kind=dp)
      if (size(x) > 0) work(left_pad + 1:left_pad + size(x)) = cmplx(x, 0.0_dp, kind=dp)

      do iter = 1, niter
         do i = 1, total_n
            total = cmplx(0.0_dp, 0.0_dp, kind=dp)
            count_valid = 0
            if (.not. complex_is_nan(work(i))) then
               total = total + work(i)
               count_valid = count_valid + 1
            end if
            do j = 1, h
               source = i - j
               if (source >= 1) then
                  if (.not. complex_is_nan(work(source))) then
                     phase = exp(cmplx(0.0_dp, 2.0_dp * pi * freq * real(j, dp), kind=dp))
                     total = total + work(source) * phase
                     count_valid = count_valid + 1
                  end if
               end if
            end do
            do j = 1, t
               source = i + j
               if (source <= total_n) then
                  if (.not. complex_is_nan(work(source))) then
                     phase = exp(cmplx(0.0_dp, -2.0_dp * pi * freq * real(j, dp), kind=dp))
                     total = total + work(source) * phase
                     count_valid = count_valid + 1
                  end if
               end if
            end do
            if (count_valid == 0) then
               next(i) = cmplx(quiet_nan(), 0.0_dp, kind=dp)
            else
               next(i) = total / real(count_valid, dp)
            end if
         end do
         work = next
      end do

      allocate(zout(size(x)))
      if (size(x) > 0) then
         zout = work(t * niter + 1:total_n - h * niter)
      end if
   end function kzft

   pure function kzs(y, m, k, t) result(smoothed)
      real(dp), intent(in) :: y(:) !! One-dimensional data smoothed by KZFT at zero frequency.
      real(dp), intent(in), optional :: m !! KZFT window size; when absent, upstream's data-driven default is used.
      integer, intent(in), optional :: k !! Number of zero-frequency KZFT iterations; defaults to 3.
      real(dp), intent(in), optional :: t(:) !! Optional indexing set retained for R API parity; upstream currently does not use it.
      real(dp), allocatable :: smoothed(:)
      complex(dp), allocatable :: z(:)
      integer :: i
      integer :: niter
      real(dp) :: mean_diff2
      real(dp) :: mval

      niter = 3
      if (present(k)) niter = k
      if (present(m)) then
         mval = m
      else if (size(y) <= 1) then
         mval = 2.0_dp
      else
         mean_diff2 = 0.0_dp
         do i = 1, size(y) - 1
            mean_diff2 = mean_diff2 + (y(i + 1) - y(i)) ** 2
         end do
         mean_diff2 = mean_diff2 / real(size(y) - 1, dp)
         mval = 100.0_dp * sqrt(mean_diff2)
         if (mval > real(size(y), dp)) mval = 2.0_dp
      end if
      z = kzft(y, f=0.0_dp, m=mval, k=niter)
      allocate(smoothed(size(y)))
      smoothed = real(z, dp)
   end function kzs

   pure function transfer_function(m, k, lamda, omega) result(tf)
      integer, intent(in) :: m !! Number of points in the moving window; must be positive for the standard transfer function.
      integer, intent(in) :: k !! Number of KZFT iterations, used as the power on the single-pass magnitude response.
      real(dp), intent(in), optional :: lamda(:) !! Frequencies in cycles per time unit; defaults to -0.5:0.01:0.5.
      real(dp), intent(in), optional :: omega !! Center frequency in cycles per time unit; defaults to zero.
      real(dp), allocatable :: tf(:)
      real(dp), allocatable :: frequencies(:)
      complex(dp) :: total
      integer :: i
      integer :: j
      real(dp) :: center
      real(dp) :: angle

      center = 0.0_dp
      if (present(omega)) center = omega
      if (present(lamda)) then
         frequencies = lamda
      else
         allocate(frequencies(101))
         do i = 1, 101
            frequencies(i) = -0.5_dp + 0.01_dp * real(i - 1, dp)
         end do
      end if
      allocate(tf(size(frequencies)))
      if (m <= 0) then
         tf = quiet_nan()
         return
      end if
      do i = 1, size(frequencies)
         total = cmplx(0.0_dp, 0.0_dp, kind=dp)
         angle = 2.0_dp * pi * (frequencies(i) - center)
         do j = 1, m
            total = total + exp(cmplx(0.0_dp, angle * real(j, dp), kind=dp))
         end do
         tf(i) = (abs(total) / real(m, dp)) ** k
      end do
   end function transfer_function

   pure function periodogram(y) result(p)
      real(dp), intent(in) :: y(:) !! Raw data transformed with R's unnormalized forward-FFT sign convention.
      real(dp), allocatable :: p(:, :)
      complex(dp) :: fourier
      integer :: half
      integer :: j
      integer :: k
      integer :: n
      real(dp) :: angle

      n = size(y)
      half = n / 2
      allocate(p(half, 2))
      do k = 0, half - 1
         fourier = cmplx(0.0_dp, 0.0_dp, kind=dp)
         do j = 0, n - 1
            angle = -2.0_dp * pi * real(j * k, dp) / real(max(n, 1), dp)
            fourier = fourier + cmplx(y(j + 1), 0.0_dp, kind=dp) * &
               exp(cmplx(0.0_dp, angle, kind=dp))
         end do
         p(k + 1, 1) = real(k + 1, dp) / real(max(n, 1), dp)
         p(k + 1, 2) = abs(fourier)
      end do
   end function periodogram

   function kztp(x, m, k, box) result(tp)
      real(dp), intent(in) :: x(:) !! One-dimensional signal used to construct the third-order KZ periodogram.
      integer, intent(in) :: m !! KZFT window and number of frequency channels.
      integer, intent(in) :: k !! Number of KZFT iterations.
      real(dp), intent(in), optional :: box(:) !! Four fractional limits (row low/high, column low/high); defaults to [0,0.5,0,0.5].
      complex(dp), allocatable :: tp(:, :)
      complex(dp), allocatable :: z(:, :)
      complex(dp), allocatable :: zcol(:)
      real(dp) :: limits(4)
      integer :: cm1
      integer :: cm2
      integer :: delta_cm
      integer :: delta_rm
      integer :: freq3
      integer :: i
      integer :: j
      integer :: rm1
      integer :: rm2
      integer :: time_index

      if ((m - 1) * k + 1 > size(x)) then
         error stop "kztp: (m-1)*k+1 must not exceed the input length"
      end if
      if (m <= 0) error stop "kztp: m must be positive"
      limits = [0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp]
      if (present(box)) then
         if (size(box) /= 4) error stop "kztp: box must contain four fractional limits"
         limits = box
      end if

      allocate(z(size(x), m))
      do i = 1, m
         zcol = kzft(x, f=real(i - 1, dp) / real(m, dp), m=real(m, dp), k=k)
         z(:, i) = zcol
      end do

      rm1 = max(r_round_even(real(m, dp) * limits(1)), 1)
      rm2 = r_round_even(real(m, dp) * limits(2))
      cm1 = max(r_round_even(real(m, dp) * limits(3)), 1)
      cm2 = r_round_even(real(m, dp) * limits(4))
      delta_rm = rm2 - rm1 + 1
      delta_cm = cm2 - cm1 + 1
      if (delta_rm <= 0 .or. delta_cm <= 0) error stop "kztp: box selects an empty frequency region"
      allocate(tp(delta_rm, delta_cm))

      do j = 1, delta_cm
         do i = 1, delta_rm
            freq3 = i + j + rm1 + cm1 - 2
            if (freq3 < 1 .or. freq3 > m) error stop "kztp: box requires a frequency index outside 1:m"
            tp(i, j) = cmplx(0.0_dp, 0.0_dp, kind=dp)
            do time_index = 1, size(x)
               tp(i, j) = tp(i, j) + z(time_index, i + rm1 - 1) * &
                  z(time_index, j + cm1 - 1) * conjg(z(time_index, freq3)) * real(m, dp) ** 2
            end do
            tp(i, j) = tp(i, j) / real(size(x), dp)
         end do
      end do
   end function kztp

end module kza_spectral
