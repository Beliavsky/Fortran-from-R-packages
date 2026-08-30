program test_boolean
   use polyclip
   implicit none
   type(poly_set) :: a, b, c
   real(dp) :: ar

   allocate(a%path(1), b%path(1))
   a%path(1) = make_path([0._dp,4._dp,4._dp,0._dp], [0._dp,0._dp,4._dp,4._dp])
   b%path(1) = make_path([2._dp,6._dp,6._dp,2._dp], [2._dp,2._dp,6._dp,6._dp])

   call polyclip_apply(a, b, c, op=clip_intersection, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1, 'intersection path count')
   call check(c%path(1)%size() == 4, 'intersection vertex count')
   call check(abs(total_area(c) - 4._dp) < 1.0e-12_dp, 'intersection area')

   call polyclip_apply(a, b, c, op=clip_union, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1, 'union path count')
   call check(c%path(1)%size() == 8, 'union vertex count')
   call check(abs(total_area(c) - 28._dp) < 1.0e-12_dp, 'union area')

   call polyclip_apply(a, b, c, op=clip_difference, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1, 'difference path count')
   call check(abs(total_area(c) - 12._dp) < 1.0e-12_dp, 'difference area')

   call polyclip_apply(a, b, c, op=clip_xor, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 2, 'xor path count')
   call check(abs(total_area(c) - 24._dp) < 1.0e-12_dp, 'xor area')

   a%path(1) = make_path([0._dp,6._dp,-10._dp,10._dp,-6._dp], &
      [10._dp,-8._dp,3._dp,3._dp,-8._dp])
   call polysimplify(a, c, fill_type=fill_evenodd, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 5, 'self-intersection evenodd path count')
   call check(abs(total_area(c) - 80._dp) < 1.0e-12_dp, 'self-intersection evenodd area')
   call polysimplify(a, c, fill_type=fill_nonzero, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1, 'self-intersection nonzero path count')
   call check(c%path(1)%size() == 10, 'self-intersection nonzero vertices')
   call check(abs(total_area(c) - 116._dp) < 1.0e-12_dp, 'self-intersection nonzero area')

   a%path(1) = make_path([-5._dp,5._dp,15._dp], [5._dp,5._dp,5._dp])
   b%path(1) = make_path([0._dp,10._dp,10._dp,0._dp], [0._dp,0._dp,10._dp,10._dp])
   call polyclip_apply(a, b, c, op=clip_intersection, closed=.false., &
      eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1, 'open intersection path count')
   call check(c%path(1)%size() == 3, 'open intersection preserves vertex')
   call check(all(abs(c%path(1)%x - [10._dp,5._dp,0._dp]) < 1.0e-12_dp), &
      'open intersection traversal')
   call polyclip_apply(a, b, c, op=clip_difference, closed=.false., &
      eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 2, 'open difference path count')

   a = b
   a%path(1)%x = a%path(1)%x(4:1:-1)
   a%path(1)%y = a%path(1)%y(4:1:-1)
   call polysimplify(a, c, fill_type=fill_positive, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 0, 'positive fill clockwise')
   call polysimplify(a, c, fill_type=fill_negative, eps=1._dp, x0=0._dp, y0=0._dp)
   call check(c%size() == 1, 'negative fill clockwise')
   ar = total_area(c)
   call check(abs(ar - 100._dp) < 1.0e-12_dp, 'negative fill area')

   print *, 'test_boolean: ok'
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
         n = s%path(k)%size()
         if (n < 3) cycle
         j = n
         do i = 1, n
            aout = aout + 0.5_dp * (s%path(k)%x(j) * s%path(k)%y(i) - &
               s%path(k)%x(i) * s%path(k)%y(j))
            j = i
         end do
      end do
   end function total_area
end program test_boolean
