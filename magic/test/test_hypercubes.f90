program test_hypercubes
   use magic
   implicit none
   type(integer_tensor) :: cube, shifted
   type(integer_tensor), allocatable :: subcubes(:)
   integer :: m, dimension

   do m = 1, 4
      cube = magiccube_2np1(m)
      call check(is_magichypercube(cube), "odd magic cube")
   end do

   do dimension = 2, 5
      cube = magichypercube_4n(1, dimension)
      call check(is_magichypercube(cube), "4n magic hypercube")
      shifted = tensor_shift(cube, [(m, m=1,dimension)])
      call check(is_semimagichypercube(shifted), "shifted semimagic hypercube")
   end do

   cube = magiccube_2np1(1)
   subcubes = diagonal_subhypercubes(cube)
   call check(size(subcubes) == 20, "rank-three diagonal subhypercube count")

   print '(a)', "test_hypercubes: PASS"
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') "FAIL: " // message
         error stop 1
      end if
   end subroutine check
end program test_hypercubes
