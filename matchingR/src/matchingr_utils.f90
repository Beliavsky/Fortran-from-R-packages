module matchingr_utils
   use matchingr_kinds, only : dp
   implicit none
   private
   public :: sort_index, sort_index_one_sided, rank_index
   public :: check_preferences, check_roommate_preferences
   public :: repeat_rows_real, repeat_cols_real

contains

   function sort_index(u) result(idx)
      real(dp), intent(in) :: u(:,:)
      integer, allocatable :: idx(:,:)
      integer :: i, j, k, n, m, key
      n = size(u,1); m = size(u,2)
      allocate(idx(n,m))
      do j = 1, m
         idx(:,j) = [(i, i=1,n)]
         do i = 2, n
            key = idx(i,j)
            k = i - 1
            do while (k >= 1)
               if (u(idx(k,j),j) > u(key,j)) exit
               if (.not.(u(idx(k,j),j) < u(key,j)) .and. &
                   .not.(u(idx(k,j),j) > u(key,j)) .and. idx(k,j) < key) exit
               idx(k+1,j) = idx(k,j)
               k = k - 1
            end do
            idx(k+1,j) = key
         end do
      end do
   end function sort_index

   function sort_index_one_sided(u) result(idx)
      real(dp), intent(in) :: u(:,:)
      integer, allocatable :: idx(:,:)
      integer :: j, k, n, m
      integer, allocatable :: full(:,:)
      n = size(u,1); m = size(u,2)
      full = sort_index(u)
      allocate(idx(n,m))
      idx = full
      ! matchingR's helper assumes u omits own utility and maps row labels
      ! around the diagonal: values >= the owner's 1-based id shift by one.
      do j = 1, m
         do k = 1, n
            if (idx(k,j) >= j) idx(k,j) = idx(k,j) + 1
         end do
      end do
   end function sort_index_one_sided

   function rank_index(sorted_idx) result(rank)
      integer, intent(in) :: sorted_idx(:,:)
      integer, allocatable :: rank(:,:)
      integer :: i, j, n, m, v, offset
      n = size(sorted_idx,1); m = size(sorted_idx,2)
      allocate(rank(n,m)); rank = 0
      offset = merge(1, 0, minval(sorted_idx) == 0)
      do j = 1, m
         do i = 1, n
            v = sorted_idx(i,j) + offset
            if (v >= 1 .and. v <= n) rank(v,j) = i - 1
         end do
      end do
   end function rank_index

   logical function check_preferences(pref, zero_based) result(ok)
      integer, intent(in) :: pref(:,:)
      logical, intent(out), optional :: zero_based
      integer :: j, n, i, v, base
      logical, allocatable :: seen(:)
      n = size(pref,1)
      ok = .true.
      if (minval(pref) == 0) then
         base = 0
      else
         base = 1
      end if
      if (present(zero_based)) zero_based = (base == 0)
      allocate(seen(n))
      do j = 1, size(pref,2)
         seen = .false.
         do i = 1, n
            v = pref(i,j) - base + 1
            if (v < 1 .or. v > n) then
               ok = .false.; return
            end if
            if (seen(v)) then
               ok = .false.; return
            end if
            seen(v) = .true.
         end do
      end do
   end function check_preferences

   logical function check_roommate_preferences(pref, zero_based) result(ok)
      integer, intent(in) :: pref(:,:)
      logical, intent(out), optional :: zero_based
      integer :: n, j, i, v, base
      logical, allocatable :: seen(:)
      n = size(pref,2)
      if (size(pref,1) /= n-1) then
         ok = .false.; return
      end if
      if (minval(pref) == 0) then
         base = 0
      else
         base = 1
      end if
      if (present(zero_based)) zero_based = (base == 0)
      allocate(seen(n))
      do j = 1, n
         seen = .false.; seen(j) = .true.
         do i = 1, n-1
            v = pref(i,j) - base + 1
            if (v < 1 .or. v > n .or. seen(v)) then
               ok = .false.; return
            end if
            seen(v) = .true.
         end do
      end do
      ok = .true.
   end function check_roommate_preferences

   function repeat_rows_real(x, reps) result(y)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: reps(:)
      real(dp), allocatable :: y(:,:)
      integer :: i, r, pos
      if (size(reps) /= size(x,1)) error stop "repeat_rows_real: reps size mismatch"
      if (any(reps < 0)) error stop "repeat_rows_real: negative repetition"
      allocate(y(sum(reps),size(x,2)))
      pos = 0
      do i = 1, size(x,1)
         do r = 1, reps(i)
            pos = pos + 1; y(pos,:) = x(i,:)
         end do
      end do
   end function repeat_rows_real

   function repeat_cols_real(x, reps) result(y)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: reps(:)
      real(dp), allocatable :: y(:,:)
      integer :: j, r, pos
      if (size(reps) /= size(x,2)) error stop "repeat_cols_real: reps size mismatch"
      if (any(reps < 0)) error stop "repeat_cols_real: negative repetition"
      allocate(y(size(x,1),sum(reps)))
      pos = 0
      do j = 1, size(x,2)
         do r = 1, reps(j)
            pos = pos + 1; y(:,pos) = x(:,j)
         end do
      end do
   end function repeat_cols_real

end module matchingr_utils
