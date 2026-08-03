program magic_cube
   use magic, only : integer_tensor, magiccube_2np1, is_magichypercube
   implicit none
   type(integer_tensor) :: cube
   integer :: layer, row, column, n

   cube = magiccube_2np1(1)
   n = cube%shape(1)
   do layer = 1, n
      write(*, '(a,i0)') "layer ", layer
      do row = 1, n
         write(*, '(*(i4))') [(cube%get([row, column, layer]), column=1,n)]
      end do
   end do
   write(*, '(a,l1)') "magic hypercube: ", is_magichypercube(cube)
end program magic_cube
