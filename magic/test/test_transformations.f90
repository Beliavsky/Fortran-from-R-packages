program test_transformations
   use magic
   implicit none
   integer(ik), allocatable :: a(:, :), b(:, :), values(:)
   integer :: k

   a = magic_square_of_order(7)
   do k = 0, 7
      b = transform_square(a, k)
      call check(is_magic(b), "dihedral transformation")
      call check(all(standardize_square(b) == standardize_square(a)), "canonical form")
   end do

   b = rotate_square(rotate_square(a))
   call check(all(b == rotate_square(a, 2)), "rotation composition")
   call check(all(rotate_square(a, 2) == reverse_square(a)), "half turn")

   allocate(values(10))
   values = [(int(k, ik), k=1,10)]
   call check(all(shift_vector(values, -2) == [3_ik,4_ik,5_ik,6_ik,7_ik,8_ik,9_ik,10_ik,1_ik,2_ik]), &
              "vector shift")

   a = circulant(values)
   call check(is_circulant_square(a, 1, 1), "circulant square")
   call check(is_persymmetric(matmul(transpose(a), a)), "circulant Gram persymmetry")

   print '(a)', "test_transformations: PASS"
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') "FAIL: " // message
         error stop 1
      end if
   end subroutine check
end program test_transformations
