program demo_magic
   use magic
   implicit none
   integer(ik), allocatable :: square(:, :), hadamard(:, :), anti(:, :)
   type(integer_tensor) :: cube
   integer :: i

   square = magic_square_of_order(5)
   write(*, '(a)') "Normal magic square of order 5"
   do i = 1, 5
      write(*, '(*(i4))') square(i, :)
   end do
   write(*, '(a,l1)') "is_magic: ", is_magic(square)

   cube = magiccube_2np1(1)
   write(*, '(a,l1)') "order-3 cube is magic: ", is_magichypercube(cube)

   hadamard = sylvester_hadamard(3)
   write(*, '(a,i0,a,l1)') "Hadamard order ", size(hadamard, 1), ": ", is_hadamard(hadamard)

   anti = sam_square(5, 2)
   write(*, '(a,l1)') "SAM(5,2) is sparse antimagic: ", is_sam(anti)
end program demo_magic
