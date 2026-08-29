! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_efp
   use r_kinds, only : dp
   use r_linalg, only : inverse_matrix
   use strucchange_recresid, only : recursive_residuals
   use strucchange_regression, only : ols_fit, root_matrix
   use strucchange_utils, only : cumulative_sum, sample_standard_deviation
   implicit none
   private
   public :: moving_estimates_process
   public :: ols_cusum
   public :: ols_mosum
   public :: process_max
   public :: process_max_l2
   public :: process_mean_l2
   public :: process_range
   public :: recursive_cusum
   public :: recursive_estimates_process
   public :: recursive_mosum
   public :: score_cusum
   public :: score_mosum
contains
   subroutine ols_cusum(x, y, process, info)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: process(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: beta(:), e(:), cs(:)
      real(dp) :: rss, sigma
      integer :: i, n, k, rank

      n = size(x, 1)
      k = size(x, 2)
      call ols_fit(x, y, beta, e, rss, rank, info)
      if (info /= 0 .or. rank < k .or. n <= k) then
         allocate(process(0, 0))
         if (info == 0) info = -1
         return
      end if
      sigma = sqrt(rss / real(n - k, dp))
      allocate(process(n + 1, 1), cs(n))
      cs = cumulative_sum(e)
      process(1, 1) = 0.0_dp
      do i = 1, n
         process(i + 1, 1) = cs(i) / (sigma * sqrt(real(n, dp)))
      end do
      info = 0
   end subroutine ols_cusum

   subroutine recursive_cusum(x, y, process, info)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: process(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: w(:), cs(:)
      real(dp) :: sigma
      integer :: i, nw

      call recursive_residuals(x, y, w, info)
      if (info /= 0) then
         allocate(process(0, 0))
         return
      end if
      nw = size(w)
      sigma = sample_standard_deviation(w)
      allocate(process(nw + 1, 1), cs(nw))
      cs = cumulative_sum(w)
      process(1, 1) = 0.0_dp
      do i = 1, nw
         process(i + 1, 1) = cs(i) / (sigma * sqrt(real(nw, dp)))
      end do
      info = 0
   end subroutine recursive_cusum

   subroutine ols_mosum(x, y, h, process, info)
      real(dp), intent(in) :: x(:, :), y(:), h
      real(dp), allocatable, intent(out) :: process(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: beta(:), e(:), cs0(:)
      real(dp) :: rss, sigma
      integer :: i, n, k, nh, rank

      n = size(x, 1)
      k = size(x, 2)
      nh = floor(real(n, dp) * h)
      if (nh < 1 .or. nh > n) then
         allocate(process(0, 0))
         info = -1
         return
      end if
      call ols_fit(x, y, beta, e, rss, rank, info)
      if (info /= 0 .or. rank < k .or. n <= k) then
         allocate(process(0, 0))
         if (info == 0) info = -2
         return
      end if
      sigma = sqrt(rss / real(n - k, dp))
      allocate(cs0(0:n))
      cs0(0) = 0.0_dp
      do i = 1, n
         cs0(i) = cs0(i - 1) + e(i)
      end do
      allocate(process(n - nh + 1, 1))
      do i = 0, n - nh
         process(i + 1, 1) = (cs0(i + nh) - cs0(i)) / &
            (sigma * sqrt(real(n, dp)))
      end do
      info = 0
   end subroutine ols_mosum

   subroutine recursive_mosum(x, y, h, process, info)
      real(dp), intent(in) :: x(:, :), y(:), h
      real(dp), allocatable, intent(out) :: process(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: w(:), cs0(:)
      real(dp) :: sigma
      integer :: i, k, nh, nw

      k = size(x, 2)
      call recursive_residuals(x, y, w, info)
      if (info /= 0) then
         allocate(process(0, 0))
         return
      end if
      nw = size(w)
      nh = floor(real(nw, dp) * h)
      if (nh < 1 .or. nh > nw .or. nw <= k) then
         allocate(process(0, 0))
         info = -1
         return
      end if
      sigma = sample_standard_deviation(w, nw - k)
      allocate(cs0(0:nw))
      cs0(0) = 0.0_dp
      do i = 1, nw
         cs0(i) = cs0(i - 1) + w(i)
      end do
      allocate(process(nw - nh + 1, 1))
      do i = 0, nw - nh
         process(i + 1, 1) = (cs0(i + nh) - cs0(i)) / &
            (sigma * sqrt(real(nw, dp)))
      end do
      info = 0
   end subroutine recursive_mosum

   subroutine recursive_estimates_process(x, y, process, info, rescale)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: process(:, :)
      integer, intent(out) :: info
      logical, intent(in), optional :: rescale
      real(dp), allocatable :: beta(:), e(:), beta_i(:), e_i(:)
      real(dp), allocatable :: qroot(:, :), q_i(:, :)
      real(dp) :: rss, rss_i, sigma, scale
      logical :: do_rescale
      integer :: i, k, n, rank, rank_i, ierr

      n = size(x, 1)
      k = size(x, 2)
      do_rescale = .true.
      if (present(rescale)) do_rescale = rescale
      call ols_fit(x, y, beta, e, rss, rank, info)
      if (info /= 0 .or. rank < k .or. n <= k) then
         allocate(process(0, 0))
         if (info == 0) info = -1
         return
      end if
      sigma = sqrt(rss / real(n - k, dp))
      call root_matrix(matmul(transpose(x), x), qroot, ierr)
      if (ierr /= 0) then
         allocate(process(0, 0))
         info = -2
         return
      end if
      qroot = qroot / (sigma * sqrt(real(n, dp)))
      allocate(process(n - k + 2, k))
      process = 0.0_dp
      do i = k, n - 1
         call ols_fit(x(1:i, :), y(1:i), beta_i, e_i, rss_i, &
            rank_i, ierr)
         if (ierr /= 0 .or. rank_i < k) then
            info = -3
            return
         end if
         if (do_rescale) then
            call root_matrix(matmul(transpose(x(1:i, :)), x(1:i, :)), &
               q_i, ierr)
            if (ierr /= 0) then
               info = -4
               return
            end if
            q_i = q_i / (sigma * sqrt(real(i, dp)))
            process(i - k + 2, :) = matmul(q_i, beta_i - beta)
         else
            process(i - k + 2, :) = matmul(qroot, beta_i - beta)
         end if
         scale = real(i, dp) / sqrt(real(n, dp))
         process(i - k + 2, :) = scale * process(i - k + 2, :)
      end do
      info = 0
   end subroutine recursive_estimates_process

   subroutine moving_estimates_process(x, y, h, process, info, rescale)
      real(dp), intent(in) :: x(:, :), y(:), h
      real(dp), allocatable, intent(out) :: process(:, :)
      integer, intent(out) :: info
      logical, intent(in), optional :: rescale
      real(dp), allocatable :: beta(:), e(:), beta_i(:), e_i(:)
      real(dp), allocatable :: qroot(:, :), q_i(:, :)
      real(dp) :: rss, rss_i, sigma
      logical :: do_rescale
      integer :: i, k, n, nh, rank, rank_i, ierr

      n = size(x, 1)
      k = size(x, 2)
      nh = floor(real(n, dp) * h)
      do_rescale = .true.
      if (present(rescale)) do_rescale = rescale
      if (nh <= k .or. nh > n) then
         allocate(process(0, 0))
         info = -1
         return
      end if
      call ols_fit(x, y, beta, e, rss, rank, info)
      if (info /= 0 .or. rank < k .or. n <= k) then
         allocate(process(0, 0))
         if (info == 0) info = -2
         return
      end if
      sigma = sqrt(rss / real(n - k, dp))
      call root_matrix(matmul(transpose(x), x), qroot, ierr)
      if (ierr /= 0) then
         allocate(process(0, 0))
         info = -3
         return
      end if
      qroot = qroot / (sigma * sqrt(real(n, dp)))
      allocate(process(n - nh + 1, k))
      do i = 0, n - nh
         call ols_fit(x(i + 1:i + nh, :), y(i + 1:i + nh), beta_i, &
            e_i, rss_i, rank_i, ierr)
         if (ierr /= 0 .or. rank_i < k) then
            info = -4
            return
         end if
         if (do_rescale) then
            call root_matrix(matmul(transpose(x(i + 1:i + nh, :)), &
               x(i + 1:i + nh, :)), q_i, ierr)
            if (ierr /= 0) then
               info = -5
               return
            end if
            q_i = q_i / (sigma * sqrt(real(nh, dp)))
            process(i + 1, :) = matmul(q_i, beta_i - beta)
         else
            process(i + 1, :) = matmul(qroot, beta_i - beta)
         end if
      end do
      process = real(nh, dp) * process / sqrt(real(n, dp))
      info = 0
   end subroutine moving_estimates_process

   subroutine score_cusum(x, y, process, info)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: process(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: beta(:), e(:), raw(:, :), qroot(:, :)
      real(dp), allocatable :: qinv(:, :), cs(:)
      real(dp) :: rss, sigma2
      integer :: i, j, k, n, rank, ierr

      n = size(x, 1)
      k = size(x, 2)
      call ols_fit(x, y, beta, e, rss, rank, info)
      if (info /= 0 .or. rank < k) then
         allocate(process(0, 0))
         if (info == 0) info = -1
         return
      end if
      sigma2 = rss / real(n, dp)
      allocate(raw(n, k + 1))
      do j = 1, k
         raw(:, j) = x(:, j) * e / sqrt(real(n, dp))
      end do
      raw(:, k + 1) = (e ** 2 - sigma2) / sqrt(real(n, dp))
      call root_matrix(matmul(transpose(raw), raw), qroot, ierr)
      if (ierr /= 0) then
         allocate(process(0, 0))
         info = -2
         return
      end if
      call inverse_matrix(qroot, qinv, ierr)
      if (ierr /= 0) then
         allocate(process(0, 0))
         info = -3
         return
      end if
      allocate(process(n + 1, k + 1), cs(k + 1))
      cs = 0.0_dp
      process(1, :) = 0.0_dp
      do i = 1, n
         cs = cs + raw(i, :)
         process(i + 1, :) = matmul(qinv, cs)
      end do
      info = 0
   end subroutine score_cusum

   subroutine score_mosum(x, y, h, process, info)
      real(dp), intent(in) :: x(:, :), y(:), h
      real(dp), allocatable, intent(out) :: process(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: beta(:), e(:), raw(:, :), qroot(:, :)
      real(dp), allocatable :: qinv(:, :), cs(:, :)
      real(dp) :: rss, sigma2
      integer :: i, j, k, n, nh, rank, ierr

      n = size(x, 1)
      k = size(x, 2)
      nh = floor(real(n, dp) * h)
      if (nh < 1 .or. nh > n) then
         allocate(process(0, 0))
         info = -1
         return
      end if
      call ols_fit(x, y, beta, e, rss, rank, info)
      if (info /= 0 .or. rank < k) then
         allocate(process(0, 0))
         if (info == 0) info = -2
         return
      end if
      sigma2 = rss / real(n, dp)
      allocate(raw(n, k + 1))
      do j = 1, k
         raw(:, j) = x(:, j) * e / sqrt(real(n, dp))
      end do
      raw(:, k + 1) = (e ** 2 - sigma2) / sqrt(real(n, dp))
      call root_matrix(matmul(transpose(raw), raw), qroot, ierr)
      if (ierr /= 0) then
         allocate(process(0, 0))
         info = -3
         return
      end if
      call inverse_matrix(qroot, qinv, ierr)
      if (ierr /= 0) then
         allocate(process(0, 0))
         info = -4
         return
      end if
      allocate(cs(0:n, k + 1))
      cs(0, :) = 0.0_dp
      do i = 1, n
         cs(i, :) = cs(i - 1, :) + raw(i, :)
      end do
      allocate(process(n - nh + 1, k + 1))
      do i = 0, n - nh
         process(i + 1, :) = matmul(qinv, cs(i + nh, :) - cs(i, :))
      end do
      info = 0
   end subroutine score_mosum

   pure real(dp) function process_max(process) result(value)
      real(dp), intent(in) :: process(:, :)
      if (size(process) == 0) then
         value = 0.0_dp
      else
         value = maxval(abs(process))
      end if
   end function process_max

   pure real(dp) function process_range(process) result(value)
      real(dp), intent(in) :: process(:, :)
      real(dp) :: one_range
      integer :: j

      value = 0.0_dp
      do j = 1, size(process, 2)
         one_range = maxval(process(:, j)) - minval(process(:, j))
         value = max(value, one_range)
      end do
   end function process_range

   pure real(dp) function process_max_l2(process) result(value)
      real(dp), intent(in) :: process(:, :)
      integer :: i

      value = 0.0_dp
      do i = 1, size(process, 1)
         value = max(value, sum(process(i, :) ** 2))
      end do
   end function process_max_l2

   pure real(dp) function process_mean_l2(process) result(value)
      real(dp), intent(in) :: process(:, :)
      real(dp) :: total
      integer :: i

      if (size(process, 1) == 0) then
         value = 0.0_dp
         return
      end if
      total = 0.0_dp
      do i = 1, size(process, 1)
         total = total + sum(process(i, :) ** 2)
      end do
      value = total / real(size(process, 1), dp)
   end function process_mean_l2
end module strucchange_efp
