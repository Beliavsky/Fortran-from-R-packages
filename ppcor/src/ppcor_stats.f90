! SPDX-License-Identifier: GPL-2.0-only
module ppcor_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ppcor_kinds, only : dp
   implicit none
   private
   integer, parameter, public :: ppcor_pearson = 1
   integer, parameter, public :: ppcor_kendall = 2
   integer, parameter, public :: ppcor_spearman = 3
   public :: association_matrix, method_name, method_from_name

contains

   subroutine association_matrix(x, method, a, info)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: method
      real(dp), allocatable, intent(out) :: a(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: ranks(:,:)
      integer :: j, n, p

      n = size(x,1)
      p = size(x,2)
      allocate(a(p,p))
      a = 0.0_dp
      info = 0
      if (n < 2 .or. p < 2) then
         info = 1
         return
      end if
      if (.not. all(ieee_is_finite(x))) then
         info = 2
         return
      end if

      select case (method)
      case (ppcor_pearson)
         call pearson_correlation(x, a, info)
      case (ppcor_spearman)
         allocate(ranks(n,p))
         do j = 1, p
            call average_ranks(x(:,j), ranks(:,j))
         end do
         call pearson_correlation(ranks, a, info)
      case (ppcor_kendall)
         call kendall_correlation(x, a)
      case default
         info = 3
      end select
   end subroutine association_matrix

   subroutine pearson_correlation(x, cor, info)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: cor(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: centered(:,:), sd(:)
      real(dp) :: denom
      integer :: i, j, n, p

      n = size(x,1)
      p = size(x,2)
      info = 0
      allocate(centered(n,p), sd(p))
      do j = 1, p
         centered(:,j) = x(:,j) - sum(x(:,j))/real(n,dp)
         sd(j) = sqrt(sum(centered(:,j)**2)/real(n-1,dp))
      end do
      if (any(sd <= 0.0_dp)) then
         info = 4
         cor = 0.0_dp
         return
      end if
      do j = 1, p
         do i = 1, j
            denom = real(n-1,dp)*sd(i)*sd(j)
            cor(i,j) = sum(centered(:,i)*centered(:,j))/denom
            cor(j,i) = cor(i,j)
         end do
      end do
      do i = 1, p
         cor(i,i) = 1.0_dp
      end do
   end subroutine pearson_correlation

   subroutine average_ranks(x, ranks)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: ranks(:)
      integer, allocatable :: idx(:)
      integer :: i, j, k, n, tmp
      real(dp) :: r

      n = size(x)
      allocate(idx(n))
      idx = [(i, i=1,n)]
      do i = 2, n
         tmp = idx(i)
         j = i-1
         do while (j >= 1)
            if (x(idx(j)) <= x(tmp)) exit
            idx(j+1) = idx(j)
            j = j-1
         end do
         idx(j+1) = tmp
      end do

      i = 1
      do while (i <= n)
         k = i
         do while (k < n)
            if (x(idx(k+1)) < x(idx(i)) .or. x(idx(k+1)) > x(idx(i))) exit
            k = k+1
         end do
         r = 0.5_dp*real(i+k,dp)
         do j = i, k
            ranks(idx(j)) = r
         end do
         i = k+1
      end do
   end subroutine average_ranks

   subroutine kendall_correlation(x, cor)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: cor(:,:)
      integer :: i, j, p

      p = size(x,2)
      cor = 0.0_dp
      do i = 1, p
         cor(i,i) = 1.0_dp
      end do
      do j = 2, p
         do i = 1, j-1
            cor(i,j) = kendall_tau_b(x(:,i), x(:,j))
            cor(j,i) = cor(i,j)
         end do
      end do
   end subroutine kendall_correlation

   real(dp) function kendall_tau_b(x, y) result(tau)
      real(dp), intent(in) :: x(:), y(:)
      integer :: i, j, n
      integer(kind=8) :: concordant, discordant, tie_x, tie_y, n0
      real(dp) :: dx, dy, denom

      n = size(x)
      concordant = 0_8
      discordant = 0_8
      tie_x = 0_8
      tie_y = 0_8
      do j = 2, n
         do i = 1, j-1
            dx = x(j)-x(i)
            dy = y(j)-y(i)
            if (.not. (dx > 0.0_dp .or. dx < 0.0_dp)) tie_x = tie_x + 1_8
            if (.not. (dy > 0.0_dp .or. dy < 0.0_dp)) tie_y = tie_y + 1_8
            if ((dx > 0.0_dp .and. dy > 0.0_dp) .or. &
                (dx < 0.0_dp .and. dy < 0.0_dp)) then
               concordant = concordant + 1_8
            else if ((dx > 0.0_dp .and. dy < 0.0_dp) .or. &
                     (dx < 0.0_dp .and. dy > 0.0_dp)) then
               discordant = discordant + 1_8
            end if
         end do
      end do
      n0 = int(n,8)*int(n-1,8)/2_8
      denom = sqrt(real(n0-tie_x,dp)*real(n0-tie_y,dp))
      if (denom > 0.0_dp) then
         tau = real(concordant-discordant,dp)/denom
      else
         tau = 0.0_dp
      end if
   end function kendall_tau_b

   pure function method_name(method) result(name)
      integer, intent(in) :: method
      character(len=8) :: name
      select case (method)
      case (ppcor_pearson)
         name = 'pearson '
      case (ppcor_kendall)
         name = 'kendall '
      case (ppcor_spearman)
         name = 'spearman'
      case default
         name = 'unknown '
      end select
   end function method_name

   pure integer function method_from_name(name) result(method)
      character(len=*), intent(in) :: name
      character(len=len(name)) :: lower
      integer :: i, c
      lower = name
      do i = 1, len(name)
         c = iachar(lower(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) lower(i:i) = achar(c+32)
      end do
      select case (trim(adjustl(lower)))
      case ('pearson')
         method = ppcor_pearson
      case ('kendall')
         method = ppcor_kendall
      case ('spearman')
         method = ppcor_spearman
      case default
         method = 0
      end select
   end function method_from_name

end module ppcor_stats
