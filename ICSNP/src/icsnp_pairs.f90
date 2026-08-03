! SPDX-License-Identifier: GPL-2.0-or-later
module icsnp_pairs
   use icsnp_kinds, only : dp
   use icsnp_status, only : icsnp_ok, icsnp_invalid_input
   implicit none
   private
   public :: pair_diff, pair_sum, pair_prod, row_norms
   public :: sum_sign_outers, sum_diff_sign_outers, spatial_ranks, signed_ranks

contains

   subroutine pair_diff(x, result, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: result(:,:)
      integer, intent(out), optional :: status
      integer :: n, p, i, j, r
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1) then
         allocate(result(0, max(0, p)))
         if (present(status)) status = icsnp_invalid_input
         return
      end if
      allocate(result(n * (n - 1) / 2, p))
      r = 0
      do i = 1, n - 1
         do j = i + 1, n
            r = r + 1
            result(r, :) = x(i, :) - x(j, :)
         end do
      end do
      if (present(status)) status = icsnp_ok
   end subroutine pair_diff

   subroutine pair_sum(x, result, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: result(:,:)
      integer, intent(out), optional :: status
      integer :: n, p, i, j, r
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1) then
         allocate(result(0, max(0, p)))
         if (present(status)) status = icsnp_invalid_input
         return
      end if
      allocate(result(n * (n - 1) / 2, p))
      r = 0
      do i = 1, n - 1
         do j = i + 1, n
            r = r + 1
            result(r, :) = x(i, :) + x(j, :)
         end do
      end do
      if (present(status)) status = icsnp_ok
   end subroutine pair_sum

   subroutine pair_prod(x, result, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: result(:,:)
      integer, intent(out), optional :: status
      integer :: n, p, i, j, r
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1) then
         allocate(result(0, max(0, p)))
         if (present(status)) status = icsnp_invalid_input
         return
      end if
      allocate(result(n * (n - 1) / 2, p))
      r = 0
      do i = 1, n - 1
         do j = i + 1, n
            r = r + 1
            result(r, :) = x(i, :) * x(j, :)
         end do
      end do
      if (present(status)) status = icsnp_ok
   end subroutine pair_prod

   pure function row_norms(x) result(norms)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: norms(size(x, 1))
      norms = sqrt(sum(x * x, dim=2))
   end function row_norms

   subroutine sum_sign_outers(x, result)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: result(size(x, 2), size(x, 2))
      real(dp) :: norm_value, u(size(x, 2))
      integer :: i
      result = 0.0_dp
      do i = 1, size(x, 1)
         norm_value = sqrt(dot_product(x(i, :), x(i, :)))
         if (norm_value > sqrt(tiny(1.0_dp))) then
            u = x(i, :) / norm_value
            result = result + spread(u, 2, size(u)) * spread(u, 1, size(u))
         end if
      end do
   end subroutine sum_sign_outers

   subroutine sum_diff_sign_outers(x, result, count)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: result(size(x, 2), size(x, 2))
      integer, intent(out), optional :: count
      real(dp) :: norm_value, u(size(x, 2))
      integer :: i, j, used
      result = 0.0_dp
      used = 0
      do i = 1, size(x, 1) - 1
         do j = i + 1, size(x, 1)
            u = x(i, :) - x(j, :)
            norm_value = sqrt(dot_product(u, u))
            if (norm_value > sqrt(tiny(1.0_dp))) then
               u = u / norm_value
               result = result + spread(u, 2, size(u)) * spread(u, 1, size(u))
               used = used + 1
            end if
         end do
      end do
      if (present(count)) count = used
   end subroutine sum_diff_sign_outers

   subroutine spatial_ranks(x, ranks)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: ranks(:,:)
      real(dp) :: difference(size(x, 2)), distance
      integer :: n, p, i, j
      n = size(x, 1)
      p = size(x, 2)
      allocate(ranks(n, p))
      ranks = 0.0_dp
      if (n < 1) return
      do i = 1, n
         do j = 1, n
            if (i == j) cycle
            difference = x(i, :) - x(j, :)
            distance = sqrt(dot_product(difference, difference))
            if (distance > sqrt(tiny(1.0_dp))) ranks(i, :) = ranks(i, :) + difference / distance
         end do
      end do
      ranks = ranks / real(n, dp)
   end subroutine spatial_ranks

   subroutine signed_ranks(x, ranks)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: ranks(:,:)
      real(dp) :: difference(size(x, 2)), sum_value(size(x, 2)), dm, dpv
      integer :: n, p, i, j
      n = size(x, 1)
      p = size(x, 2)
      allocate(ranks(n, p))
      ranks = 0.0_dp
      if (n < 1) return
      do i = 1, n
         do j = 1, n
            if (i == j) cycle
            difference = x(i, :) - x(j, :)
            sum_value = x(i, :) + x(j, :)
            dm = sqrt(dot_product(difference, difference))
            dpv = sqrt(dot_product(sum_value, sum_value))
            if (dm > sqrt(tiny(1.0_dp))) ranks(i, :) = ranks(i, :) + 0.5_dp * difference / dm
            if (dpv > sqrt(tiny(1.0_dp))) ranks(i, :) = ranks(i, :) + 0.5_dp * sum_value / dpv
         end do
      end do
      ranks = ranks / real(n, dp)
   end subroutine signed_ranks

end module icsnp_pairs
