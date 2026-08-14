! SPDX-License-Identifier: GPL-2.0-only
module mco_pareto
   use mco_kinds, only : dp
   implicit none
   private
   public :: dominates, constrained_dominates, nondominated_sort
   public :: crowding_distance, pareto_mask, pareto_filter
contains
   logical function dominates(a, b) result(ans)
      real(dp), intent(in) :: a(:), b(:)
      ans = all(a <= b) .and. any(a < b)
   end function dominates

   logical function constrained_dominates(a, av, b, bv) result(ans)
      real(dp), intent(in) :: a(:), b(:), av, bv
      real(dp), parameter :: tol = 100.0_dp*epsilon(1.0_dp)
      if (av <= tol .and. bv > tol) then
         ans = .true.
      else if (av > tol .and. bv <= tol) then
         ans = .false.
      else if (av > tol .and. bv > tol) then
         ans = av < bv
      else
         ans = dominates(a, b)
      end if
   end function constrained_dominates

   subroutine nondominated_sort(values, violation, rank, nfront)
      real(dp), intent(in) :: values(:,:), violation(:)
      integer, intent(out) :: rank(size(values,2))
      integer, intent(out), optional :: nfront
      integer :: n, i, j, r, left
      logical, allocatable :: assigned(:)
      n = size(values,2)
      if (size(violation) /= n) error stop "nondominated_sort: size mismatch"
      rank = 0
      allocate(assigned(n), source=.false.)
      r = 0
      do
         left = count(.not. assigned)
         if (left == 0) exit
         r = r + 1
         do i = 1, n
            if (assigned(i)) cycle
            do j = 1, n
               if (i == j .or. assigned(j)) cycle
               if (constrained_dominates(values(:,j), violation(j), values(:,i), violation(i))) exit
            end do
            if (j > n) rank(i) = -r
         end do
         do i = 1, n
            if (rank(i) == -r) then
               rank(i) = r
               assigned(i) = .true.
            end if
         end do
      end do
      if (present(nfront)) nfront = r
   end subroutine nondominated_sort

   subroutine crowding_distance(values, rank, distance)
      real(dp), intent(in) :: values(:,:)
      integer, intent(in) :: rank(:)
      real(dp), intent(out) :: distance(size(rank))
      integer :: n, m, r, i, j, k, nf
      integer, allocatable :: idx(:), work(:)
      real(dp) :: lo, hi
      n = size(values,2); m = size(values,1)
      if (size(rank) /= n) error stop "crowding_distance: size mismatch"
      distance = 0.0_dp
      do r = 1, maxval(rank)
         nf = count(rank == r)
         if (nf == 0) cycle
         allocate(idx(nf), work(nf))
         k = 0
         do i = 1, n
            if (rank(i) == r) then
               k = k + 1; idx(k) = i
            end if
         end do
         if (nf <= 2) then
            distance(idx) = huge(1.0_dp)
         else
            do j = 1, m
               work = idx
               call sort_indices_by_value(values(j,:), work)
               lo = values(j,work(1)); hi = values(j,work(nf))
               distance(work(1)) = huge(1.0_dp)
               distance(work(nf)) = huge(1.0_dp)
               if (hi > lo) then
                  do k = 2, nf-1
                     if (distance(work(k)) < huge(1.0_dp)/2) &
                        distance(work(k)) = distance(work(k)) + &
                        (values(j,work(k+1))-values(j,work(k-1)))/(hi-lo)
                  end do
               end if
            end do
         end if
         deallocate(idx, work)
      end do
   end subroutine crowding_distance

   subroutine sort_indices_by_value(x, idx)
      real(dp), intent(in) :: x(:)
      integer, intent(inout) :: idx(:)
      integer :: i, j, key
      do i = 2, size(idx)
         key = idx(i); j = i-1
         do while (j >= 1)
            if (x(idx(j)) <= x(key)) exit
            idx(j+1) = idx(j); j = j-1
         end do
         idx(j+1) = key
      end do
   end subroutine sort_indices_by_value

   function pareto_mask(front) result(mask)
      real(dp), intent(in) :: front(:,:)
      logical :: mask(size(front,2))
      integer :: i, j, n
      n = size(front,2); mask = .true.
      do i = 1, n
         if (.not. mask(i)) cycle
         do j = 1, n
            if (i == j) cycle
            if (dominates(front(:,j), front(:,i))) then
               mask(i) = .false.; exit
            end if
         end do
      end do
   end function pareto_mask

   function pareto_filter(front) result(filtered)
      real(dp), intent(in) :: front(:,:)
      real(dp), allocatable :: filtered(:,:)
      logical :: mask(size(front,2))
      integer :: i, k
      mask = pareto_mask(front)
      allocate(filtered(size(front,1), count(mask)))
      k = 0
      do i = 1, size(front,2)
         if (mask(i)) then
            k = k + 1; filtered(:,k) = front(:,i)
         end if
      end do
   end function pareto_filter
end module mco_pareto
