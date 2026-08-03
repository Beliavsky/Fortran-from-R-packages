! SPDX-License-Identifier: GPL-2.0-only
module glmnet_data
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use glmnet_kinds, only : dp, glmnet_eps
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument
   use glmnet_types, only : glmnet_survival_data, glmnet_sparse_csc
   implicit none
   private
   public :: na_replace, na_sparse_fix, prepare_x, make_x, rmult, stratify_surv
contains
   subroutine na_replace(x, result, status, replacement)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: result(:,:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: replacement(:)
      real(dp), allocatable :: fill(:)
      integer :: i, j, count_finite
      status = glmnet_success
      allocate(result(size(x, 1), size(x, 2)), fill(size(x, 2)))
      if (present(replacement)) then
         if (size(replacement) /= size(x, 2)) then
            status = glmnet_invalid_argument
            result = 0.0_dp
            return
         end if
         fill = replacement
      else
         do j = 1, size(x, 2)
            fill(j) = 0.0_dp
            count_finite = 0
            do i = 1, size(x, 1)
               if (ieee_is_finite(x(i, j))) then
                  fill(j) = fill(j) + x(i, j)
                  count_finite = count_finite + 1
               end if
            end do
            if (count_finite > 0) fill(j) = fill(j) / real(count_finite, dp)
         end do
      end if
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            if (ieee_is_finite(x(i, j))) then
               result(i, j) = x(i, j)
            else
               result(i, j) = fill(j)
            end if
         end do
      end do
   end subroutine na_replace

   subroutine na_sparse_fix(x, replacement)
      type(glmnet_sparse_csc), intent(inout) :: x
      real(dp), intent(in), optional :: replacement
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      if (present(replacement)) value = replacement
      if (.not. allocated(x%values)) return
      do i = 1, size(x%values)
         if (.not. ieee_is_finite(x%values(i))) x%values(i) = value
      end do
   end subroutine na_sparse_fix

   subroutine prepare_x(x, result, keep, status, remove_constant)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: result(:,:)
      integer, allocatable, intent(out) :: keep(:)
      integer, intent(out) :: status
      logical, intent(in), optional :: remove_constant
      real(dp), allocatable :: filled(:,:)
      logical :: drop_constant
      integer :: j, count_keep
      call na_replace(x, filled, status)
      if (status /= glmnet_success) then
         allocate(result(0, 0), keep(0))
         return
      end if
      drop_constant = .true.
      if (present(remove_constant)) drop_constant = remove_constant
      count_keep = 0
      do j = 1, size(filled, 2)
         if (.not. drop_constant .or. &
             maxval(filled(:, j)) - minval(filled(:, j)) > 100.0_dp * glmnet_eps) &
            count_keep = count_keep + 1
      end do
      allocate(keep(count_keep), result(size(filled, 1), count_keep))
      count_keep = 0
      do j = 1, size(filled, 2)
         if (.not. drop_constant .or. &
             maxval(filled(:, j)) - minval(filled(:, j)) > 100.0_dp * glmnet_eps) then
            count_keep = count_keep + 1
            keep(count_keep) = j
            result(:, count_keep) = filled(:, j)
         end if
      end do
   end subroutine prepare_x

   subroutine make_x(train, combined, status, test)
      real(dp), intent(in) :: train(:,:)
      real(dp), allocatable, intent(out) :: combined(:,:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: test(:,:)
      status = glmnet_success
      if (present(test)) then
         if (size(test, 2) /= size(train, 2)) then
            status = glmnet_invalid_argument
            allocate(combined(0, 0))
            return
         end if
         allocate(combined(size(train, 1) + size(test, 1), size(train, 2)))
         combined(:size(train, 1), :) = train
         combined(size(train, 1) + 1:, :) = test
      else
         allocate(combined(size(train, 1), size(train, 2)))
         combined = train
      end if
   end subroutine make_x

   subroutine rmult(probabilities, draws, counts, status, seed)
      real(dp), intent(in) :: probabilities(:)
      integer, intent(in) :: draws
      integer, allocatable, intent(out) :: counts(:)
      integer, intent(out) :: status
      integer, intent(in), optional :: seed
      real(dp), allocatable :: cumulative(:)
      real(dp) :: total, u
      integer(kind=8) :: state
      integer :: d, k
      status = glmnet_success
      if (draws < 0 .or. size(probabilities) < 1 .or. any(probabilities < 0.0_dp)) then
         status = glmnet_invalid_argument
         allocate(counts(0))
         return
      end if
      total = sum(probabilities)
      if (total <= glmnet_eps) then
         status = glmnet_invalid_argument
         allocate(counts(0))
         return
      end if
      allocate(counts(size(probabilities)), cumulative(size(probabilities)))
      counts = 0
      cumulative(1) = probabilities(1) / total
      do k = 2, size(probabilities)
         cumulative(k) = cumulative(k - 1) + probabilities(k) / total
      end do
      cumulative(size(cumulative)) = 1.0_dp
      state = 88172645463325252_8
      if (present(seed)) state = int(abs(seed) + 1, kind=8)
      do d = 1, draws
         call next_uniform(state, u)
         do k = 1, size(cumulative)
            if (u <= cumulative(k)) then
               counts(k) = counts(k) + 1
               exit
            end if
         end do
      end do
   end subroutine rmult

   subroutine stratify_surv(start_time, stop_time, event, strata, data, status)
      real(dp), intent(in) :: start_time(:), stop_time(:)
      integer, intent(in) :: event(:), strata(:)
      type(glmnet_survival_data), intent(out) :: data
      integer, intent(out) :: status
      integer :: n
      n = size(stop_time)
      status = glmnet_success
      if (size(start_time) /= n .or. size(event) /= n .or. size(strata) /= n) then
         status = glmnet_invalid_argument
         return
      end if
      if (any(start_time < 0.0_dp) .or. any(stop_time <= start_time) .or. &
          any(event < 0) .or. any(event > 1) .or. any(strata < 1)) then
         status = glmnet_invalid_argument
         return
      end if
      allocate(data%start(n), data%stop(n), data%event(n), data%strata(n))
      data%start = start_time
      data%stop = stop_time
      data%event = event
      data%strata = strata
   end subroutine stratify_surv

   subroutine next_uniform(state, value)
      integer(kind=8), intent(inout) :: state
      real(dp), intent(out) :: value
      integer(kind=8), parameter :: mask53 = int(z'001FFFFFFFFFFFFF', kind=8)
      state = ieor(state, shiftl(state, 13))
      state = ieor(state, shiftr(state, 7))
      state = ieor(state, shiftl(state, 17))
      value = real(iand(state, mask53), dp) / real(mask53 + 1_8, dp)
   end subroutine next_uniform
end module glmnet_data
