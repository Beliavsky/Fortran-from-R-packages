program test_minkowski_pip
   use polyclip
   implicit none
   type(poly_set) :: pattern, path, c
   integer, allocatable :: ans(:)

   allocate(pattern%path(1), path%path(1))
   pattern%path(1) = make_path([-1._dp,1._dp,1._dp,-1._dp], [-1._dp,-1._dp,1._dp,1._dp])
   path%path(1) = make_path([0._dp,4._dp,4._dp,0._dp], [0._dp,0._dp,4._dp,4._dp])
   call polyminkowski(pattern, path, c, closed=.true., eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1 .and. c%path(1)%size() == 4, 'Minkowski path count')
   call check(abs(total_area(c) - 36._dp) < 1.0e-12_dp, 'Minkowski area')

   call pointinpolygon([2._dp,0._dp,-1._dp], [2._dp,2._dp,2._dp], path%path(1), ans, &
      eps=1._dp, x0=0._dp, y0=0._dp)
   call check(size(ans) == 3, 'point in polygon length')
   call check(all(ans == [1,-1,0]), 'point in polygon values')

   print *, 'test_minkowski_pip: ok'
contains
   subroutine check(ok, message)
      logical, intent(in) :: ok
      character(*), intent(in) :: message
      if (.not. ok) then
         write(*, '(a)') 'FAILED: '//message
         error stop 1
      end if
   end subroutine check
   pure real(dp) function total_area(s) result(aout)
      type(poly_set), intent(in) :: s
      integer :: i, j, k, n
      aout = 0._dp
      do k = 1, s%size()
         n = s%path(k)%size(); if (n < 3) cycle
         j = n
         do i = 1, n
            aout = aout + 0.5_dp * (s%path(k)%x(j) * s%path(k)%y(i) - &
               s%path(k)%x(i) * s%path(k)%y(j))
            j = i
         end do
      end do
   end function total_area
end program test_minkowski_pip
