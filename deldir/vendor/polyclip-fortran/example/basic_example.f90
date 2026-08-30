program basic_example
   use polyclip
   implicit none
   type(poly_set) :: a, b, intersection
   integer :: i

   allocate(a%path(1), b%path(1))
   a%path(1) = make_path([0._dp,4._dp,4._dp,0._dp], [0._dp,0._dp,4._dp,4._dp])
   b%path(1) = make_path([2._dp,6._dp,6._dp,2._dp], [2._dp,2._dp,6._dp,6._dp])

   call polyclip_apply(a, b, intersection, op=clip_intersection)
   print '(a,i0)', 'number of output polygons: ', intersection%size()
   if (intersection%size() > 0) then
      do i = 1, intersection%path(1)%size()
         print '(2f12.6)', intersection%path(1)%x(i), intersection%path(1)%y(i)
      end do
   end if
end program basic_example
