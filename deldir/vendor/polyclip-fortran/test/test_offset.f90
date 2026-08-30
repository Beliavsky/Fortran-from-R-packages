program test_offset
   use polyclip
   implicit none
   type(poly_set) :: a, c, donut

   allocate(a%path(1), donut%path(2))
   a%path(1) = make_path([0._dp,4._dp,4._dp,0._dp], [0._dp,0._dp,4._dp,4._dp])
   call polyoffset(a, 1._dp, c, join_type=join_square, arc_tolerance=.01_dp, &
      eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1 .and. c%path(1)%size() == 8, 'square join geometry')
   call check(abs(total_area(c) - 34._dp) < 1.0e-12_dp, 'square join area')
   call polyoffset(a, 1._dp, c, join_type=join_miter, arc_tolerance=.01_dp, &
      eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1 .and. c%path(1)%size() == 4, 'miter join geometry')
   call check(abs(total_area(c) - 36._dp) < 1.0e-12_dp, 'miter join area')

   a%path(1) = make_path([0._dp,10._dp,10._dp,0._dp], [0._dp,0._dp,10._dp,10._dp])
   call polyoffset(a, 3._dp, c, join_type=join_round, arc_tolerance=.25_dp, &
      eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1 .and. c%path(1)%size() == 12, 'round join geometry')
   call check(abs(total_area(c) - 244._dp) < 1.0e-12_dp, 'round join area')

   a%path(1) = make_path([0._dp,20._dp,20._dp,0._dp], [0._dp,0._dp,20._dp,20._dp])
   call polyoffset(a, -3._dp, c, join_type=join_miter, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1, 'negative offset count')
   call check(abs(total_area(c) - 196._dp) < 1.0e-12_dp, 'negative offset area')

   donut%path(1) = a%path(1)
   donut%path(2) = make_path([5._dp,5._dp,15._dp,15._dp], [5._dp,15._dp,15._dp,5._dp])
   call polyoffset(donut, 2._dp, c, join_type=join_miter, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 2, 'donut positive path count')
   call check(abs(total_area(c) - 540._dp) < 1.0e-12_dp, 'donut positive area')
   call polyoffset(donut, -2._dp, c, join_type=join_miter, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 2, 'donut negative path count')
   call check(abs(total_area(c) - 60._dp) < 1.0e-12_dp, 'donut negative area')

   a%path(1) = make_path([0._dp,10._dp,10._dp], [0._dp,0._dp,10._dp])
   call polylineoffset(a, 2._dp, c, join_type=join_round, end_type=end_open_round, &
      arc_tolerance=.25_dp, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1 .and. c%path(1)%size() == 11, 'open round buffer geometry')
   call check(abs(total_area(c) - 91._dp) < 1.0e-12_dp, 'open round buffer area')

   print *, 'test_offset: ok'
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
end program test_offset
