program create_magic_square
   use magic, only : ik, magic_square_of_order, is_magic, magic_constant
   implicit none
   integer(ik), allocatable :: square(:, :)
   integer :: i

   square = magic_square_of_order(6)
   do i = 1, size(square, 1)
      write(*, '(*(i5))') square(i, :)
   end do
   write(*, '(a,l1)') "magic: ", is_magic(square)
   write(*, '(a,i0)') "magic constant: ", magic_constant(size(square, 1))
end program create_magic_square
