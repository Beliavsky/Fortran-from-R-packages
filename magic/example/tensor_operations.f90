program tensor_operations
   use magic, only : integer_tensor, ik, sequence_tensor, tensor_shift, &
                     tensor_reverse, tensor_block_diag
   implicit none
   type(integer_tensor) :: a, shifted, reversed, block

   a = sequence_tensor([2, 3])
   shifted = tensor_shift(a, [1, -1])
   reversed = tensor_reverse(a, [.true., .false.])
   block = tensor_block_diag(a, shifted, -1_ik)

   write(*, '(a,*(i0,1x))') "a shape: ", a%shape
   write(*, '(a,*(i0,1x))') "shifted values: ", shifted%values
   write(*, '(a,*(i0,1x))') "reversed values: ", reversed%values
   write(*, '(a,*(i0,1x))') "block shape: ", block%shape
end program tensor_operations
