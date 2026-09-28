program test_direct
   use sanic, only: dp, sparse_matrix, sparsify, solve2_dense, solve2_sparse, &
                    solve_chol_dense, solve_chol_sparse, solve_lu_dense, solve_lu_sparse, solve_qr_dense, solve_qr_sparse, &
                    dense_is_symmetric
   implicit none

   real(dp), parameter :: tol = 5.0e-10_dp
   real(dp) :: a(3, 3), aindef(3, 3), aspd(3, 3), arect(4, 3), rank_def(4, 3)
   real(dp) :: b(3, 1), b4(4, 1), b_rank(4, 1), truth(3, 1)
   real(dp) :: a4(4, 4), b_fill(4, 1), spd_fill(4, 4), truth4(4, 1)
   real(dp) :: wide_qr(2, 3), b_wide(2, 1), expected_wide(3, 1)
   real(dp), allocatable :: inv(:, :), x(:, :)
   type(sparse_matrix) :: sa
   integer :: info, i

   truth(:, 1) = [1.0_dp, -2.0_dp, 0.5_dp]
   a = reshape([3.0_dp, 1.0_dp, 2.0_dp, &
                2.0_dp, 4.0_dp, 0.0_dp, &
                -1.0_dp, 2.0_dp, 5.0_dp], [3, 3])
   b = matmul(a, truth)

   aspd = 0.0_dp
   aspd(1, 1) = 1.0_dp
   aspd(2, 2) = 1.0_dp
   aspd(3, 3) = 1.0_dp
   aspd(1, 2) = 1.0e-17_dp
   call require(dense_is_symmetric(aspd, tol=0.0_dp), "zero symmetry tolerance uses epsilon")

   call solve_lu_dense(a, x, info, b=b, pivot=1)
   call require(info == 0 .and. maxval(abs(x - truth)) < tol, "partial LU")
   call solve_lu_dense(a, x, info, b=b, pivot=2)
   call require(info == 0 .and. maxval(abs(x - truth)) < tol, "full LU")
   call solve_lu_dense(a, inv, info)
   call require(info == 0 .and. maxval(abs(matmul(a, inv) - identity3())) < 2.0e-9_dp, "LU inverse")
   call solve_lu_dense(a, x, info, b=b, pivot=3)
   call require(info == 1 .and. size(x) == 0, "invalid dense pivot rejected")

   aspd = matmul(transpose(a), a) + 0.5_dp*identity3()
   b = matmul(aspd, truth)
   call solve_chol_dense(aspd, x, info, b=b, pivot=0)
   call require(info == 0 .and. maxval(abs(x - truth)) < tol, "LLT")
   call solve_chol_dense(aspd, x, info, b=b, pivot=1)
   call require(info == 0 .and. maxval(abs(x - truth)) < tol, "LDLT")

   aindef = reshape([0.0_dp, 1.0_dp, 1.0_dp, &
                     1.0_dp, 2.0_dp, 0.0_dp, &
                     1.0_dp, 0.0_dp, 3.0_dp], [3, 3])
   b = matmul(aindef, truth)
   call solve_chol_dense(aindef, x, info, b=b, pivot=1)
   call require(info == 0 .and. maxval(abs(x - truth)) < 2.0e-9_dp, "pivoted indefinite LDLT")

   b = matmul(aspd, truth)
   call solve2_dense(aspd, x, info, b=b)
   call require(info == 0 .and. maxval(abs(x - truth)) < tol, "solve2 symmetric")

   arect = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                    0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
                    2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp], [4, 3])
   b4 = matmul(arect, truth)
   call solve_qr_dense(arect, x, info, b=b4, pivot=1)
   call require(info == 0 .and. maxval(abs(x - truth)) < 5.0e-9_dp, "QR least squares")
   call solve2_dense(arect, x, info, b=b4)
   call require(info == 0 .and. maxval(abs(x - truth)) < 5.0e-9_dp, "solve2 rectangular")
   sa = sparsify(arect)
   do i = 0, 2
      call solve_qr_sparse(sa, x, info, b=b4, ordering=i)
      call require(info == 0 .and. maxval(abs(x - truth)) < 5.0e-9_dp, "sparse QR ordering")
   end do

   rank_def(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   rank_def(:, 2) = 2.0_dp*rank_def(:, 1)
   rank_def(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
   b_rank(:, 1) = 3.0_dp*rank_def(:, 1) - rank_def(:, 3)
   call solve_qr_dense(rank_def, x, info, b=b_rank, pivot=0)
   call require(info == 2, "unpivoted rank-deficient QR status")
   do i = 1, 2
      call solve_qr_dense(rank_def, x, info, b=b_rank, pivot=i)
      call require(info == 0, "pivoted rank-deficient QR status")
      call require(maxval(abs(matmul(rank_def, x) - b_rank)) < 2.0e-9_dp, "rank-deficient QR residual")
   end do
   sa = sparsify(rank_def)
   do i = 0, 2
      call solve_qr_sparse(sa, x, info, b=b_rank, ordering=i)
      call require(info == 0, "sparse rank-deficient QR status")
      call require(maxval(abs(matmul(rank_def, x) - b_rank)) < 2.0e-9_dp, "sparse rank-deficient QR residual")
   end do

   wide_qr = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [2, 3])
   b_wide(:, 1) = [1.0_dp, 2.0_dp]
   expected_wide(:, 1) = [1.0_dp, 2.0_dp, 0.0_dp]
   call solve_qr_dense(wide_qr, x, info, b=b_wide, pivot=0)
   call require(info == 0 .and. maxval(abs(x - expected_wide)) < 2.0e-10_dp, "wide unpivoted QR basic solution")
   expected_wide(:, 1) = [0.0_dp, 1.0_dp, 1.0_dp]
   call solve_qr_dense(wide_qr, x, info, b=b_wide, pivot=1)
   call require(info == 0 .and. maxval(abs(x - expected_wide)) < 2.0e-10_dp, &
                "wide column-pivoted QR basic solution")
   call solve_qr_dense(wide_qr, x, info, b=b_wide, pivot=2)
   call require(info == 0 .and. maxval(abs(matmul(wide_qr, x) - b_wide)) < 2.0e-10_dp, &
                "wide full-pivot QR residual")
   sa = sparsify(wide_qr)
   do i = 0, 2
      call solve_qr_sparse(sa, x, info, b=b_wide, ordering=i)
      call require(info == 0 .and. maxval(abs(matmul(wide_qr, x) - b_wide)) < 2.0e-10_dp, &
                   "sparse wide QR residual")
   end do

   sa = sparsify(aspd)
   b = matmul(aspd, truth)
   call solve2_sparse(sa, x, info, b=b)
   call require(info == 0 .and. maxval(abs(x - truth)) < tol, "sparse direct")
   call solve_chol_sparse(sa, x, info, b=b, ordering=3)
   call require(info == 1 .and. size(x) == 0, "invalid sparse ordering rejected")

   truth4(:, 1) = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp]
   a4 = reshape([0.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, &
                 1.0_dp, 4.0_dp, 3.0_dp, 0.0_dp, &
                 0.0_dp, 1.0_dp, 4.0_dp, 2.0_dp, &
                 2.0_dp, 0.0_dp, 1.0_dp, 3.0_dp], [4, 4])
   b_fill = matmul(a4, truth4)
   sa = sparsify(a4)
   do i = 0, 2
      call solve_lu_sparse(sa, x, info, b=b_fill, ordering=i)
      call require(info == 0 .and. maxval(abs(x - truth4)) < 2.0e-9_dp, "sparse LU ordering")
   end do

   spd_fill = reshape([5.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
                       1.0_dp, 5.0_dp, 1.0_dp, 0.0_dp, &
                       0.0_dp, 1.0_dp, 5.0_dp, 1.0_dp, &
                       1.0_dp, 0.0_dp, 1.0_dp, 5.0_dp], [4, 4])
   b_fill = matmul(spd_fill, truth4)
   sa = sparsify(spd_fill)
   do i = 0, 2
      call solve_chol_sparse(sa, x, info, b=b_fill, pivot=0, ordering=i)
      call require(info == 0 .and. maxval(abs(x - truth4)) < 2.0e-9_dp, "sparse LLT ordering")
      call solve_chol_sparse(sa, x, info, b=b_fill, pivot=1, ordering=i)
      call require(info == 0 .and. maxval(abs(x - truth4)) < 2.0e-9_dp, "sparse LDLT ordering")
   end do

   print '(a)', 'test_direct: PASS'

contains

   pure function identity3() result(x)
      real(dp) :: x(3, 3)
      integer :: j
      x = 0.0_dp
      do j = 1, 3
         x(j, j) = 1.0_dp
      end do
   end function identity3

   subroutine require(condition, label)
      logical, intent(in) :: condition !! True when the regression condition has passed.
      character(len=*), intent(in) :: label !! Short diagnostic label printed on failure.
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAILED:', trim(label)
         error stop 1
      end if
   end subroutine require

end program test_direct
