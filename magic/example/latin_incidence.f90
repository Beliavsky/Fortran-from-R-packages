program latin_incidence
   use magic, only : ik, integer_tensor, latin_square, incidence, &
                     unincidence, is_latin_square, is_incidence
   implicit none
   integer(ik), allocatable :: latin(:, :), recovered(:, :)
   type(integer_tensor) :: inc
   integer :: i

   latin = latin_square(5)
   inc = incidence(latin)
   recovered = unincidence(inc)

   do i = 1, size(latin, 1)
      write(*, '(*(i3))') latin(i, :)
   end do
   write(*, '(a,l1)') "Latin: ", is_latin_square(latin)
   write(*, '(a,l1)') "proper incidence tensor: ", is_incidence(inc, .false.)
   write(*, '(a,l1)') "round trip: ", all(latin == recovered)
end program latin_incidence
