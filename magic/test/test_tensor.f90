program test_tensor
   use magic
   implicit none
   type(integer_tensor) :: a, b, c, d
   integer, allocatable :: dims(:)

   a = sequence_tensor([2, 3, 4])
   b = tensor_reverse(a, [.true., .false., .true.])
   c = tensor_reverse(b, [.true., .false., .true.])
   call check(tensor_equal(a, c), "double reversal")

   b = tensor_shift(a, a%shape)
   call check(tensor_equal(a, b), "full-period shift")

   b = tensor_permute(a, [3, 1, 2])
   c = tensor_permute(b, [2, 3, 1])
   call check(tensor_equal(a, c), "inverse permutation")

   b = make_tensor([2, 2], 1_ik)
   c = make_tensor([1, 1], 2_ik)
   d = tensor_block_diag(b, c)
   call check(all(d%shape == [3, 3]), "block-diagonal shape")
   call check(d%get([1, 1]) == 1_ik .and. d%get([3, 3]) == 2_ik, "block-diagonal values")
   call check(d%get([1, 3]) == 0_ik, "block-diagonal padding")

   b = make_tensor([2, 2, 2], 1_ik)
   c = tensor_subsums(b, [2, 2, 2])
   call check(all(c%values == 8_ik), "wrapped subarray sums")

   b = sequence_tensor([2, 3])
   c = tensor_pad(b, [1, 2], PAD_EXTEND)
   d = tensor_drop(c, [-1, -2])
   call check(tensor_equal(b, d), "padding and dropping")

   d = tensor_take(c, [2, 3])
   call check(tensor_equal(b, d), "APL take")

   dims = first_nonsingleton_dimensions(make_tensor([1, 1, 2, 1, 3]), 2)
   call check(all(dims == [3, 5]), "first nonsingleton dimensions")

   print '(a)', "test_tensor: PASS"
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') "FAIL: " // message
         error stop 1
      end if
   end subroutine check
end program test_tensor
