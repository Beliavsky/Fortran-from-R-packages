! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_auxiliary
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH
   implicit none
   private

   type, public :: pava_result
      real(dp), allocatable :: values(:)
      integer, allocatable :: blocks(:)
      logical :: increasing = .true.
   end type pava_result

   public :: pava_blocks, pava_fitted, autocorrelation, isoacf

contains

   subroutine pava_blocks(x, result, status, weights, block_sizes, increasing)
      real(dp), intent(in) :: x(:)
      type(pava_result), intent(out) :: result
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: weights(:)
      integer, intent(in), optional :: block_sizes(:)
      logical, intent(in), optional :: increasing
      real(dp), allocatable :: level(:), weight_work(:)
      integer, allocatable :: count(:)
      integer :: n, m, i
      logical :: up, violates

      n = size(x)
      if (n <= 0) then
         allocate(result%values(0), result%blocks(0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            allocate(result%values(0), result%blocks(0))
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
      end if
      if (present(block_sizes)) then
         if (size(block_sizes) /= n .or. any(block_sizes <= 0)) then
            allocate(result%values(0), result%blocks(0))
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
      end if

      up = .true.
      if (present(increasing)) up = increasing
      allocate(level(n), weight_work(n), count(n))
      m = 0
      do i = 1, n
         m = m + 1
         level(m) = x(i)
         weight_work(m) = 1.0_dp
         if (present(weights)) weight_work(m) = weights(i)
         count(m) = 1
         if (present(block_sizes)) count(m) = block_sizes(i)

         do while (m > 1)
            if (up) then
               violates = level(m - 1) > level(m)
            else
               violates = level(m - 1) < level(m)
            end if
            if (.not. violates) exit
            level(m - 1) = (weight_work(m - 1) * level(m - 1) + &
               weight_work(m) * level(m)) / (weight_work(m - 1) + weight_work(m))
            weight_work(m - 1) = weight_work(m - 1) + weight_work(m)
            count(m - 1) = count(m - 1) + count(m)
            m = m - 1
         end do
      end do

      allocate(result%values(m), result%blocks(m))
      result%values = level(1:m)
      result%blocks = count(1:m)
      result%increasing = up
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine pava_blocks

   subroutine pava_fitted(x, fitted, status, weights, increasing)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: fitted(:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: weights(:)
      logical, intent(in), optional :: increasing
      type(pava_result) :: result
      integer :: i, j, position, info

      if (present(weights) .and. present(increasing)) then
         call pava_blocks(x, result, info, weights = weights, increasing = increasing)
      else if (present(weights)) then
         call pava_blocks(x, result, info, weights = weights)
      else if (present(increasing)) then
         call pava_blocks(x, result, info, increasing = increasing)
      else
         call pava_blocks(x, result, info)
      end if
      if (info /= SANDWICH_SUCCESS) then
         allocate(fitted(0))
         if (present(status)) status = info
         return
      end if

      allocate(fitted(size(x)))
      position = 1
      do i = 1, size(result%values)
         do j = 1, result%blocks(i)
            fitted(position) = result%values(i)
            position = position + 1
         end do
      end do
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine pava_fitted

   subroutine autocorrelation(x, lagmax, acf, status, demean, unbiased)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: lagmax
      real(dp), allocatable, intent(out) :: acf(:)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: demean, unbiased
      real(dp), allocatable :: z(:)
      real(dp) :: center, denominator, numerator
      integer :: n, lag
      logical :: use_demean, use_unbiased

      n = size(x)
      if (n <= 1 .or. lagmax < 0 .or. lagmax >= n) then
         allocate(acf(0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      use_demean = .true.
      if (present(demean)) use_demean = demean
      use_unbiased = .false.
      if (present(unbiased)) use_unbiased = unbiased
      center = 0.0_dp
      if (use_demean) center = sum(x) / real(n, dp)
      allocate(z(n))
      z = x - center
      denominator = dot_product(z, z)
      if (denominator <= tiny(1.0_dp)) then
         allocate(acf(0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      allocate(acf(lagmax + 1))
      acf(1) = 1.0_dp
      do lag = 1, lagmax
         numerator = dot_product(z(1:n - lag), z(1 + lag:n))
         if (use_unbiased) numerator = numerator * real(n, dp) / real(n - lag, dp)
         acf(lag + 1) = numerator / denominator
      end do
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine autocorrelation

   subroutine isoacf(x, result, status, lagmax, weave1)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: result(:)
      integer, intent(out), optional :: status
      integer, intent(in), optional :: lagmax
      logical, intent(in), optional :: weave1
      real(dp), allocatable :: raw(:), sequence(:), fitted(:), weights(:)
      real(dp), allocatable :: a(:), b(:)
      real(dp) :: ma, mb
      integer :: n, maximum_lag, lag, m, info
      logical :: use_weave1

      n = size(x)
      if (n <= 2) then
         allocate(result(0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      use_weave1 = .false.
      if (present(weave1)) use_weave1 = weave1

      if (use_weave1) then
         maximum_lag = int(5.0_dp * sqrt(real(n, dp)))
         if (present(lagmax)) maximum_lag = lagmax
         maximum_lag = min(n - 1, max(1, maximum_lag))
         allocate(raw(maximum_lag - 1))
         do lag = 1, maximum_lag - 1
            m = n - lag
            allocate(a(m), b(m))
            a = x(1:m)
            b = x(1 + lag:n)
            ma = sum(a) / real(m, dp)
            mb = sum(b) / real(m, dp)
            raw(lag) = dot_product(a - ma, b - mb) / real(max(m - 1, 1), dp)
            raw(lag) = raw(lag) / (dot_product(x - sum(x) / real(n, dp), &
               x - sum(x) / real(n, dp)) / real(n - 1, dp))
            deallocate(a, b)
         end do
         allocate(sequence(size(raw) + 1), weights(size(raw) + 1))
         sequence(1:size(raw)) = raw
         sequence(size(sequence)) = 0.0_dp
         do lag = 1, size(raw)
            weights(lag) = real(n - lag, dp)
         end do
         weights(size(weights)) = huge(1.0_dp) / 100.0_dp
         call pava_fitted(sequence, fitted, info, weights = weights, increasing = .false.)
      else
         maximum_lag = n - 2
         if (present(lagmax)) maximum_lag = min(maximum_lag, lagmax)
         call autocorrelation(x, maximum_lag, raw, info)
         if (info /= SANDWICH_SUCCESS) then
            allocate(result(0))
            if (present(status)) status = info
            return
         end if
         allocate(sequence(maximum_lag + 1))
         if (maximum_lag > 0) sequence(1:maximum_lag) = raw(2:maximum_lag + 1)
         sequence(maximum_lag + 1) = 0.0_dp
         call pava_fitted(sequence, fitted, info, increasing = .false.)
      end if

      if (info /= SANDWICH_SUCCESS) then
         allocate(result(0))
         if (present(status)) status = info
         return
      end if
      allocate(result(size(fitted) + 1))
      result(1) = 1.0_dp
      result(2:) = fitted
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine isoacf

end module sandwich_auxiliary
