program test_sparse
   use sanic, only: dp, sparse_matrix, sparsify, sparse_from_triplets, sparse_to_dense, &
                    sparse_matvec, sparse_transpose_matvec, sparse_matmat, sparse_is_symmetric
   implicit none

   real(dp), parameter :: tol = 1.0e-12_dp
   real(dp) :: dense(3, 3), expected(3, 3), rhs(3, 2)
   real(dp), allocatable :: restored(:, :), y(:), yt(:), ym(:, :)
   type(sparse_matrix) :: a, b
   integer :: info

   dense = reshape([2.0_dp, 0.0_dp, 1.0e-14_dp, &
                    -1.0_dp, 3.0_dp, 0.0_dp, &
                    0.0_dp, 4.0_dp, 5.0_dp], [3, 3])
   a = sparsify(dense, tol=1.0e-13_dp)
   restored = sparse_to_dense(a)
   expected = dense
   expected(3, 1) = 0.0_dp
   call require(maxval(abs(restored - expected)) < tol, "sparsify tolerance")

   b = sparse_from_triplets(3, 3, [1, 1, 2, 3, 3], [1, 1, 2, 2, 3], &
                            [1.0_dp, 2.0_dp, 4.0_dp, -1.0_dp, 5.0_dp])
   restored = sparse_to_dense(b)
   expected = 0.0_dp
   expected(1, 1) = 3.0_dp
   expected(2, 2) = 4.0_dp
   expected(3, 2) = -1.0_dp
   expected(3, 3) = 5.0_dp
   call require(maxval(abs(restored - expected)) < tol, "triplet duplicate accumulation")

   call sparse_matvec(b, [1.0_dp, 2.0_dp, -1.0_dp], y, info)
   call require(info == 0 .and. maxval(abs(y - [3.0_dp, 8.0_dp, -7.0_dp])) < tol, "sparse matvec")

   call sparse_transpose_matvec(b, [1.0_dp, 2.0_dp, -1.0_dp], yt, info)
   call require(info == 0 .and. maxval(abs(yt - [3.0_dp, 9.0_dp, -5.0_dp])) < tol, "transpose matvec")

   rhs = reshape([1.0_dp, 0.0_dp, -1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], [3, 2])
   call sparse_matmat(b, rhs, ym, info)
   call require(info == 0 .and. maxval(abs(ym - matmul(expected, rhs))) < tol, "sparse matmat")

   a = sparsify(reshape([2.0_dp, 1.0_dp, 0.0_dp, &
                        1.0_dp, 3.0_dp, 4.0_dp, &
                        0.0_dp, 4.0_dp, 5.0_dp], [3, 3]))
   call require(sparse_is_symmetric(a), "sparse symmetry true")
   call require(.not. sparse_is_symmetric(b), "sparse symmetry false")
   a = sparse_from_triplets(2, 2, [1, 1, 2, 2, 1, 2], [1, 2, 1, 2, 2, 1], &
                            [2.0_dp, 1.0_dp, 4.0_dp, 5.0_dp, 2.0_dp, -1.0_dp])
   call require(sparse_is_symmetric(a), "sparse duplicate symmetry")

   print '(a)', 'test_sparse: PASS'

contains

   subroutine require(condition, label)
      logical, intent(in) :: condition !! True when the regression condition has passed.
      character(len=*), intent(in) :: label !! Short diagnostic label printed on failure.
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAILED:', trim(label)
         error stop 1
      end if
   end subroutine require

end program test_sparse
