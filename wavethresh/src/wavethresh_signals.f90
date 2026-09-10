! SPDX-License-Identifier: GPL-2.0-or-later
! Test signals and random Haar constructions translated from wavethresh 4.7.3.
module wavethresh_signals
   use wavethresh_types, only : dp, wd_t, signal_set_t, chirp_t, test_signal_t
   use wavethresh_transform_1d, only : wst, wr_wst
   use wavethresh_spectrum, only : ewspec
   implicit none
   private
   public :: doppler, simchirp, dj_ex, test_data_ct, haar_ma, haar_concat, lsw_sim, check_my_ews
contains

   elemental function doppler(t) result(value)
      real(dp), intent(in) :: t !! Position in [0,1] at which the wavethresh Doppler test function is evaluated.
      real(dp) :: value
      value = sqrt(max(t * (1.0_dp - t), 0.0_dp)) * &
         sin((2.0_dp * acos(-1.0_dp) * 1.05_dp) / (t + 0.05_dp))
   end function doppler

   function simchirp(n) result(chirp)
      integer, intent(in), optional :: n !! Number of samples in the chirp; default is 1024.
      type(chirp_t) :: chirp
      integer :: count
      integer :: i
      count = 1024
      if (present(n)) count = n
      if (count < 1) then
         allocate(chirp%x(0), chirp%y(0))
         return
      end if
      allocate(chirp%x(count), chirp%y(count))
      do i = 1, count
         chirp%x(i) = 1.0e-5_dp - 1.0_dp + 2.0_dp * real(i - 1, dp) / real(count, dp)
         chirp%y(i) = sin(acos(-1.0_dp) / chirp%x(i))
      end do
   end function simchirp

   function dj_ex(n, signal_sd, rsnr, noisy, seed) result(signals)
      integer, intent(in), optional :: n !! Number of samples in each Donoho-Johnstone test signal; default is 1024.
      real(dp), intent(in), optional :: signal_sd !! Target sample standard deviation; default is 7.
      real(dp), intent(in), optional :: rsnr !! Signal-to-noise standard-deviation ratio; default is 7.
      logical, intent(in), optional :: noisy !! Add Gaussian noise to every signal when true; default is false.
      integer, intent(in), optional :: seed !! Optional deterministic seed for the Fortran intrinsic RNG.
      type(signal_set_t) :: signals
      real(dp), parameter :: knots(11) = [0.10_dp, 0.13_dp, 0.15_dp, 0.23_dp, 0.25_dp, &
         0.40_dp, 0.44_dp, 0.65_dp, 0.76_dp, 0.78_dp, 0.81_dp]
      real(dp), parameter :: block_height(11) = [4.0_dp, -5.0_dp, 3.0_dp, -4.0_dp, 5.0_dp, &
         -4.2_dp, 2.1_dp, 4.3_dp, -3.1_dp, 2.1_dp, -4.2_dp]
      real(dp), parameter :: bump_height(11) = [4.0_dp, 5.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
         4.2_dp, 2.1_dp, 4.3_dp, 3.1_dp, 5.1_dp, 4.2_dp]
      real(dp), parameter :: bump_width(11) = [0.005_dp, 0.005_dp, 0.006_dp, 0.01_dp, 0.01_dp, &
         0.03_dp, 0.01_dp, 0.01_dp, 0.005_dp, 0.008_dp, 0.005_dp]
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: noise(:)
      real(dp) :: target
      real(dp) :: ratio
      real(dp) :: term
      integer :: count
      integer :: i
      integer :: j
      logical :: add_noise

      count = 1024
      if (present(n)) count = n
      target = 7.0_dp
      if (present(signal_sd)) target = signal_sd
      ratio = 7.0_dp
      if (present(rsnr)) ratio = rsnr
      add_noise = .false.
      if (present(noisy)) add_noise = noisy
      if (count < 2) then
         allocate(signals%blocks(0), signals%bumps(0), signals%heavi_sine(0), signals%doppler(0))
         return
      end if
      allocate(x(count), signals%blocks(count), signals%bumps(count))
      allocate(signals%heavi_sine(count), signals%doppler(count))
      do i = 1, count
         x(i) = real(i, dp) / real(count, dp)
      end do
      signals%blocks = 0.0_dp
      signals%bumps = 0.0_dp
      do j = 1, size(knots)
         do i = 1, count
            if (x(i) > knots(j)) signals%blocks(i) = signals%blocks(i) + block_height(j)
            term = max(0.0_dp, 1.0_dp - abs((x(i) - knots(j)) / bump_width(j)))
            signals%bumps(i) = signals%bumps(i) + bump_height(j) * term**4
         end do
      end do
      signals%heavi_sine = 4.0_dp * sin(4.0_dp * acos(-1.0_dp) * x) - &
         sign(1.0_dp, x - 0.30_dp) - sign(1.0_dp, 0.72_dp - x)
      signals%doppler = sqrt(x * (1.0_dp - x)) * &
         sin((2.0_dp * acos(-1.0_dp) * 0.95_dp) / (x + 0.05_dp))
      call normalize_sample_sd(signals%blocks, target)
      call normalize_sample_sd(signals%bumps, target)
      call normalize_sample_sd(signals%heavi_sine, target)
      call normalize_sample_sd(signals%doppler, target)
      if (add_noise .and. ratio > 0.0_dp) then
         if (present(seed)) call seed_rng(seed)
         noise = normal_vector(count) * target / ratio
         signals%blocks = signals%blocks + noise
         noise = normal_vector(count) * target / ratio
         signals%bumps = signals%bumps + noise
         noise = normal_vector(count) * target / ratio
         signals%heavi_sine = signals%heavi_sine + noise
         noise = normal_vector(count) * target / ratio
         signals%doppler = signals%doppler + noise
      end if
   end function dj_ex

   function test_data_ct(signal_type, n, signal_sd, rsnr, seed) result(result)
      character(len=*), intent(in), optional :: signal_type !! ppoly, blocks, bumps, heavi, or doppler; default is ppoly.
      integer, intent(in), optional :: n !! Number of test-signal samples; default is 512.
      real(dp), intent(in), optional :: signal_sd !! Target sample standard deviation; default is 1.
      real(dp), intent(in), optional :: rsnr !! Signal-to-noise standard-deviation ratio; default is 7.
      integer, intent(in), optional :: seed !! Optional deterministic seed for the added Gaussian noise.
      type(test_signal_t) :: result
      real(dp), parameter :: knots(11) = [0.10_dp, 0.13_dp, 0.15_dp, 0.23_dp, 0.25_dp, &
         0.40_dp, 0.44_dp, 0.65_dp, 0.76_dp, 0.78_dp, 0.81_dp]
      real(dp), parameter :: block_height(11) = [4.0_dp, -5.0_dp, 3.0_dp, -4.0_dp, 5.0_dp, &
         -4.2_dp, 2.1_dp, 4.3_dp, -3.1_dp, 2.1_dp, -4.2_dp]
      real(dp), parameter :: bump_height(11) = [4.0_dp, 5.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
         4.2_dp, 2.1_dp, 4.3_dp, 3.1_dp, 5.1_dp, 4.2_dp]
      real(dp), parameter :: bump_width(11) = [0.005_dp, 0.005_dp, 0.006_dp, 0.01_dp, 0.01_dp, &
         0.03_dp, 0.01_dp, 0.01_dp, 0.005_dp, 0.008_dp, 0.005_dp]
      character(len=24) :: kind
      real(dp) :: target
      real(dp) :: ratio
      integer :: count
      integer :: i
      integer :: j

      count = 512
      if (present(n)) count = n
      target = 1.0_dp
      if (present(signal_sd)) target = signal_sd
      ratio = 7.0_dp
      if (present(rsnr)) ratio = rsnr
      kind = "ppoly"
      if (present(signal_type)) kind = signal_type
      if (count < 2 .or. ratio <= 0.0_dp) then
         result%message = "n must exceed one and rsnr must be positive"
         return
      end if
      allocate(result%x(count), result%signal(count), result%noisy(count))
      do i = 1, count
         result%x(i) = real(i - 1, dp) / real(count, dp)
      end do
      result%signal = 0.0_dp
      select case (trim(kind))
      case ("ppoly")
         do i = 1, count
            if (result%x(i) <= 0.5_dp) then
               result%signal(i) = -16.0_dp * result%x(i)**3 + 12.0_dp * result%x(i)**2
            else if (result%x(i) <= 0.75_dp) then
               result%signal(i) = result%x(i) * &
                  (16.0_dp * result%x(i)**2 - 40.0_dp * result%x(i) + 28.0_dp) / 3.0_dp - 1.5_dp
            else
               result%signal(i) = result%x(i) * &
                  (16.0_dp * result%x(i)**2 - 32.0_dp * result%x(i) + 16.0_dp) / 3.0_dp
            end if
         end do
      case ("blocks")
         do j = 1, size(knots)
            where (result%x > knots(j)) result%signal = result%signal + block_height(j)
         end do
      case ("bumps")
         do j = 1, size(knots)
            result%signal = result%signal + bump_height(j) / &
               (1.0_dp + abs((result%x - knots(j)) / bump_width(j)))**4
         end do
      case ("heavi")
         result%signal = 4.0_dp * sin(4.0_dp * acos(-1.0_dp) * result%x) - &
            sign(1.0_dp, result%x - 0.30_dp) - sign(1.0_dp, 0.72_dp - result%x)
      case ("doppler")
         result%signal = doppler(result%x)
      case default
         result%message = "unknown wavethresh test.dataCT signal type"
         return
      end select
      call normalize_sample_sd(result%signal, target)
      if (present(seed)) call seed_rng(seed)
      result%noisy = result%signal + normal_vector(count) * target / ratio
      result%signal_type = kind
      result%ok = .true.
      result%message = "ok"
   end function test_data_ct

   function haar_ma(n, standard_deviation, order, seed) result(x)
      integer, intent(in) :: n !! Number of output observations in the Haar moving-average realization.
      real(dp), intent(in), optional :: standard_deviation !! Standard deviation of Gaussian increments; default is one.
      integer, intent(in), optional :: order !! Haar moving-average order; default is five.
      integer, intent(in), optional :: seed !! Optional deterministic seed for the Fortran intrinsic RNG.
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: z(:)
      real(dp) :: sd
      integer :: jorder
      integer :: i
      sd = 1.0_dp
      if (present(standard_deviation)) sd = standard_deviation
      jorder = 5
      if (present(order)) jorder = order
      if (n < 1 .or. jorder < 1) then
         allocate(x(0))
         return
      end if
      if (present(seed)) call seed_rng(seed)
      z = normal_vector(n + 2**jorder - 1) * sd
      allocate(x(n))
      x = 0.0_dp
      do i = 2**jorder, 2**(jorder - 1) + 1, -1
         x = x + z(i:i + n - 1)
      end do
      do i = 2**(jorder - 1), 1, -1
         x = x - z(i:i + n - 1)
      end do
      x = x * 2.0_dp**(-0.5_dp * real(jorder, dp))
   end function haar_ma

   function haar_concat(seed) result(x)
      integer, intent(in), optional :: seed !! Optional base seed used for the four successive HaarMA segments.
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: one(:)
      real(dp), allocatable :: two(:)
      real(dp), allocatable :: three(:)
      real(dp), allocatable :: four(:)
      integer :: base
      base = 13579
      if (present(seed)) base = seed
      one = haar_ma(128, order=1, seed=base)
      two = haar_ma(128, order=2, seed=base + 1)
      three = haar_ma(128, order=3, seed=base + 2)
      four = haar_ma(128, order=4, seed=base + 3)
      x = [one, two, three, four]
   end function haar_concat

   function check_my_ews(spec, nsim, seed) result(average)
      type(wd_t), intent(in) :: spec !! Nonnegative stationary-wavelet spectrum used to generate repeated LSW realizations.
      integer, intent(in), optional :: nsim !! Number of Monte Carlo realizations averaged; default is 10 and must be positive.
      integer, intent(in), optional :: seed !! Optional deterministic base seed; realization k uses seed+k-1.
      type(wd_t) :: average
      type(wd_t) :: transformed
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: zeros(:)
      integer :: count
      integer :: simulation
      integer :: level
      integer :: base_seed

      count = 10
      if (present(nsim)) count = nsim
      if (.not. spec%ok .or. trim(spec%transform_type) /= "station" .or. count < 1) then
         average%message = "checkmyews requires a valid stationary spectrum and positive nsim"
         return
      end if
      if (spec%n_original < 1 .or. spec%nlevels < 1) then
         average%message = "checkmyews requires a nonempty stationary spectrum"
         return
      end if
      allocate(zeros(spec%n_original), source=0.0_dp)
      average = wst(zeros, filter_number=1.0_dp, family="DaubExPhase")
      if (.not. average%ok) then
         average%message = "checkmyews could not initialize the Haar accumulator"
         return
      end if
      base_seed = 0
      if (present(seed)) base_seed = seed
      do simulation = 1, count
         if (present(seed)) then
            series = lsw_sim(spec, seed=base_seed + simulation - 1)
         else
            series = lsw_sim(spec)
         end if
         if (size(series) /= spec%n_original) then
            average%ok = .false.
            average%message = "checkmyews LSW simulation failed"
            return
         end if
         transformed = wst(series, filter_number=1.0_dp, family="DaubExPhase")
         associate (estimate => ewspec(transformed, filter_number=1.0_dp, family="DaubExPhase"))
            if (.not. estimate%ok) then
               average%ok = .false.
               average%message = estimate%message
               return
            end if
            do level = 0, average%nlevels - 1
               average%detail(level)%values = average%detail(level)%values + estimate%spectrum%detail(level)%values
            end do
            do level = 0, average%nlevels
               average%scaling(level)%values = average%scaling(level)%values + estimate%spectrum%scaling(level)%values
            end do
         end associate
      end do
      do level = 0, average%nlevels - 1
         average%detail(level)%values = average%detail(level)%values / real(count, dp)
      end do
      do level = 0, average%nlevels
         average%scaling(level)%values = average%scaling(level)%values / real(count, dp)
      end do
      average%ok = .true.
      average%message = "ok"
   end function check_my_ews

   function lsw_sim(spec, seed) result(series)
      type(wd_t), intent(in) :: spec !! Nonnegative stationary-wavelet spectrum whose detail rows define local energies.
      integer, intent(in), optional :: seed !! Optional deterministic seed for the intrinsic normal random-number stream.
      real(dp), allocatable :: series(:)
      type(wd_t) :: randomized
      real(dp), allocatable :: noise(:)
      integer :: level
      integer :: n

      if (.not. spec%ok .or. trim(spec%transform_type) /= "station") then
         allocate(series(0))
         return
      end if
      if (spec%nlevels <= 0 .or. spec%n_original /= 2**spec%nlevels) then
         allocate(series(0))
         return
      end if
      n = spec%n_original
      do level = 0, spec%nlevels - 1
         if (.not. allocated(spec%detail(level)%values)) then
            allocate(series(0))
            return
         end if
         if (size(spec%detail(level)%values) /= n .or. any(spec%detail(level)%values < 0.0_dp)) then
            allocate(series(0))
            return
         end if
      end do
      if (.not. allocated(spec%scaling)) then
         allocate(series(0))
         return
      end if
      if (lbound(spec%scaling, 1) > 0 .or. ubound(spec%scaling, 1) < spec%nlevels) then
         allocate(series(0))
         return
      end if
      do level = 0, spec%nlevels
         if (.not. allocated(spec%scaling(level)%values)) then
            allocate(series(0))
            return
         end if
      end do
      if (present(seed)) call seed_rng(seed)
      randomized%n_original = spec%n_original
      randomized%nlevels = spec%nlevels
      randomized%transform_type = spec%transform_type
      randomized%boundary = spec%boundary
      randomized%filter = spec%filter
      randomized%ok = .true.
      randomized%message = "ok"
      allocate(randomized%detail(0:spec%nlevels - 1))
      allocate(randomized%scaling(0:spec%nlevels))
      do level = 0, spec%nlevels
         randomized%scaling(level)%values = spec%scaling(level)%values
      end do
      do level = spec%nlevels - 1, 0, -1
         noise = normal_vector(n)
         randomized%detail(level)%values = sqrt(spec%detail(level)%values) * &
            2.0_dp**(spec%nlevels - level) * noise
      end do
      series = wr_wst(randomized)
   end function lsw_sim

   subroutine normalize_sample_sd(x, target)
      real(dp), intent(inout) :: x(:) !! Signal values rescaled in place to the requested sample standard deviation.
      real(dp), intent(in) :: target !! Desired sample standard deviation after rescaling.
      real(dp) :: mean_value
      real(dp) :: variance
      if (size(x) < 2) return
      mean_value = sum(x) / real(size(x), dp)
      variance = sum((x - mean_value)**2) / real(size(x) - 1, dp)
      if (variance > tiny(1.0_dp)) x = x * target / sqrt(variance)
   end subroutine normalize_sample_sd

   function normal_vector(n) result(z)
      integer, intent(in) :: n !! Number of independent standard-normal variates requested.
      real(dp), allocatable :: z(:)
      real(dp) :: u1
      real(dp) :: u2
      integer :: i
      allocate(z(n))
      i = 1
      do while (i <= n)
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         z(i) = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
         if (i + 1 <= n) z(i + 1) = sqrt(-2.0_dp * log(u1)) * sin(2.0_dp * acos(-1.0_dp) * u2)
         i = i + 2
      end do
   end function normal_vector

   subroutine seed_rng(seed)
      integer, intent(in) :: seed !! User seed expanded deterministically to the compiler RNG seed vector.
      integer, allocatable :: put(:)
      integer :: nseed
      integer :: i
      call random_seed(size=nseed)
      allocate(put(nseed))
      do i = 1, nseed
         put(i) = modulo(abs(seed) + 104729 * i, huge(1) - 1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_rng

end module wavethresh_signals
