program test_krylov
   use sanic, only: dp, sparse_matrix, sparsify, krylov_result, arnoldi_sparse, lanczos_sparse
   implicit none

   real(dp) :: a(3, 3), b(3)
   type(sparse_matrix) :: sa
   type(krylov_result) :: out

   a = 0.0_dp
   a(1, 1) = 1.0_dp
   a(2, 2) = 2.0_dp
   a(3, 3) = 4.0_dp
   b = [1.0_dp, 1.0_dp, 1.0_dp]
   sa = sparsify(a)
   call lanczos_sparse(sa, out, b=b, iter=3, tol=1.0e-14_dp, eigen=.true.)
   call require(out%info == 0 .and. out%iterations == 3, "Lanczos status")
   call require(maxval(abs(sort3(real(out%values, dp)) - [1.0_dp, 2.0_dp, 4.0_dp])) < 1.0e-8_dp, &
                "Lanczos eigenvalues")
   call lanczos_sparse(sa, out, b=b, iter=3, tol=-1.0_dp, eigen=.false.)
   call require(out%info == 1 .and. size(out%q) == 0, "negative Krylov tolerance rejected")

   a = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
                2.0_dp, 3.0_dp, 0.0_dp, &
                0.0_dp, 1.0_dp, 5.0_dp], [3, 3])
   sa = sparsify(a)
   call arnoldi_sparse(sa, out, b=b, symmetric=.false., iter=3, tol=1.0e-14_dp, eigen=.true.)
   call require(out%info == 0 .and. out%iterations == 3, "Arnoldi status")
   call require(maxval(abs(sort3(real(out%values, dp)) - [1.0_dp, 3.0_dp, 5.0_dp])) < 2.0e-6_dp, &
                "Arnoldi eigenvalues")

   a = 0.0_dp
   a(1, 1) = 1.0_dp
   a(2, 2) = 2.0_dp
   a(3, 3) = 3.0_dp
   b = [1.0_dp, 0.0_dp, 0.0_dp]
   sa = sparsify(a)
   call arnoldi_sparse(sa, out, b=b, symmetric=.false., iter=3, tol=1.0e-14_dp, eigen=.true.)
   call require(out%info == 0 .and. out%iterations == 1, "Arnoldi early-breakdown status")
   call require(all(shape(out%h) == [3, 3]) .and. all(shape(out%q) == [3, 3]), "Arnoldi requested shapes")
   call require(maxval(abs(out%h(2:3, :))) < 1.0e-14_dp .and. maxval(abs(out%q(:, 2:3))) < 1.0e-14_dp, &
                "Arnoldi zero-filled inactive projection")
   call require(maxval(abs(sort3(real(out%values, dp)) - [0.0_dp, 0.0_dp, 1.0_dp])) < 1.0e-10_dp, &
                "Arnoldi zero-padded eigenvalues")

   call lanczos_sparse(sa, out, b=b, iter=3, tol=1.0e-14_dp, eigen=.true.)
   call require(out%info == 0 .and. out%iterations == 1, "Lanczos early-breakdown status")
   call require(all(shape(out%q) == [3, 3]) .and. size(out%diagonal) == 3 .and. size(out%subdiagonal) == 2, &
                "Lanczos requested shapes")
   call require(maxval(abs(out%q(:, 2:3))) < 1.0e-14_dp .and. maxval(abs(out%diagonal(2:3))) < 1.0e-14_dp, &
                "Lanczos zero-filled inactive projection")
   call require(maxval(abs(sort3(real(out%values, dp)) - [0.0_dp, 0.0_dp, 1.0_dp])) < 1.0e-10_dp, &
                "Lanczos zero-padded eigenvalues")

   print '(a)', 'test_krylov: PASS'

contains

   pure function sort3(x) result(y)
      real(dp), intent(in) :: x(3) !! Three values to sort in ascending order.
      real(dp) :: y(3), t
      integer :: i, j
      y = x
      do i = 1, 2
         do j = i + 1, 3
            if (y(j) < y(i)) then
               t = y(i)
               y(i) = y(j)
               y(j) = t
            end if
         end do
      end do
   end function sort3

   subroutine require(condition, label)
      logical, intent(in) :: condition !! True when the regression condition has passed.
      character(len=*), intent(in) :: label !! Short diagnostic label printed on failure.
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAILED:', trim(label)
         error stop 1
      end if
   end subroutine require

end program test_krylov
