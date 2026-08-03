program panmagic_square
   use magic, only : ik, panmagic_6np1, is_panmagic, is_normal_square
   implicit none
   integer(ik), allocatable :: square(:, :)
   integer :: i

   square = panmagic_6np1(1)
   do i = 1, size(square, 1)
      write(*, '(*(i4))') square(i, :)
   end do
   write(*, '(a,l1)') "panmagic: ", is_panmagic(square)
   write(*, '(a,l1)') "normal: ", is_normal_square(square)
end program panmagic_square
