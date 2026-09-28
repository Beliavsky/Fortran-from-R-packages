program test_spectral
   use sanic, only: dp, eigen_result, svd_result, eigen2, svd2, sanic_svd_bdc, sanic_svd_jacobi
   implicit none

   real(dp) :: a2(2, 2), a3(3, 3), lower_defined(3, 3), rect(4, 3), wide(3, 4)
   real(dp) :: sigma_tall(3, 3), sigma_wide(3, 4)
   real(dp) :: left0(3, 3), right0(3, 3), diag0(3, 3), c, s
   real(dp), allocatable :: recon(:, :)
   type(eigen_result) :: eig
   type(svd_result) :: svd
   complex(dp) :: residual2(2), residual3(3)
   integer :: i

   a2 = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
   call eigen2(a2, eig, symmetric=.true.)
   call require(eig%info == 0, "symmetric eigen status")
   call require(maxval(abs(real(eig%values, dp) - [1.0_dp, 3.0_dp])) < 1.0e-10_dp, "symmetric eigenvalues")
   do i = 1, 2
      residual2 = matmul(cmplx(a2, 0.0_dp, kind=dp), eig%vectors(:, i)) - eig%values(i)*eig%vectors(:, i)
      call require(maxval(abs(residual2)) < 1.0e-9_dp, "symmetric eigenvectors")
   end do

   lower_defined = reshape([2.0_dp, 1.0_dp, 0.0_dp, &
                            99.0_dp, 2.0_dp, 1.0_dp, &
                            -77.0_dp, 88.0_dp, 2.0_dp], [3, 3])
   call eigen2(lower_defined, eig, symmetric=.true.)
   call require(eig%info == 0, "self-adjoint lower-triangle status")
   call require(maxval(abs(real(eig%values, dp) - [2.0_dp - sqrt(2.0_dp), 2.0_dp, &
                                                   2.0_dp + sqrt(2.0_dp)])) < 2.0e-10_dp, &
                "self-adjoint lower-triangle values")
   a3 = reshape([2.0_dp, 1.0_dp, 0.0_dp, &
                 1.0_dp, 2.0_dp, 1.0_dp, &
                 0.0_dp, 1.0_dp, 2.0_dp], [3, 3])
   do i = 1, 3
      residual3 = matmul(cmplx(a3, 0.0_dp, kind=dp), eig%vectors(:, i)) - eig%values(i)*eig%vectors(:, i)
      call require(maxval(abs(residual3)) < 2.0e-9_dp, "self-adjoint lower-triangle vectors")
   end do

   a2 = reshape([0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp], [2, 2])
   call eigen2(a2, eig, symmetric=.false., max_iter=2000)
   call require(eig%info == 0, "general eigen status")
   call require(maxval(abs(abs(aimag(eig%values)) - 1.0_dp)) < 2.0e-6_dp, "complex eigenvalues")
   do i = 1, 2
      residual2 = matmul(cmplx(a2, 0.0_dp, kind=dp), eig%vectors(:, i)) - eig%values(i)*eig%vectors(:, i)
      call require(maxval(abs(residual2)) < 1.0e-7_dp, "general eigenvectors")
   end do

   a3 = 0.0_dp
   a3(1, 1) = 1.0_dp
   a3(1, 2) = 1000.0_dp
   a3(2, 2) = 2.0_dp
   a3(2, 3) = 1000.0_dp
   a3(3, 3) = 3.0_dp
   call eigen2(a3, eig, symmetric=.false., max_iter=10000)
   call require(eig%info == 0, "nonnormal eigen status")
   do i = 1, 3
      call require(minval(abs(eig%values(i) - cmplx([1.0_dp, 2.0_dp, 3.0_dp], 0.0_dp, kind=dp))) < 1.0e-10_dp, &
                   "nonnormal eigenvalue")
      residual3 = matmul(cmplx(a3, 0.0_dp, kind=dp), eig%vectors(:, i)) - eig%values(i)*eig%vectors(:, i)
      call require(maxval(abs(residual3)) < 2.0e-10_dp, "nonnormal eigenvector")
   end do

   rect = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                   0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
                   2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp], [4, 3])
   call svd2(rect, svd, method=sanic_svd_jacobi, thin=.true.)
   call require(svd%info == 0, "thin SVD status")
   sigma_tall = 0.0_dp
   do i = 1, 3
      sigma_tall(i, i) = svd%d(i)
   end do
   recon = matmul(matmul(svd%u, sigma_tall), transpose(svd%v))
   call require(maxval(abs(recon - rect)) < 2.0e-10_dp, "thin SVD reconstruction")
   call require(all(shape(svd%u) == [4, 3]) .and. all(shape(svd%v) == [3, 3]), "thin SVD shapes")

   wide = transpose(rect)
   call svd2(wide, svd, method=sanic_svd_bdc, thin=.false.)
   call require(svd%info == 0, "full wide SVD status")
   sigma_wide = 0.0_dp
   do i = 1, 3
      sigma_wide(i, i) = svd%d(i)
   end do
   recon = matmul(matmul(svd%u, sigma_wide), transpose(svd%v))
   call require(maxval(abs(recon - wide)) < 2.0e-10_dp, "full wide SVD reconstruction")
   call require(all(shape(svd%u) == [3, 3]) .and. all(shape(svd%v) == [4, 4]), "full wide SVD shapes")

   call svd2(rect, svd, vectors=.false.)
   call require(svd%info == 0 .and. size(svd%d) == 3, "values-only SVD status")
   call require(size(svd%u) == 0 .and. size(svd%v) == 0, "values-only SVD vectors")

   left0 = 0.0_dp
   right0 = 0.0_dp
   diag0 = 0.0_dp
   c = cos(0.7_dp)
   s = sin(0.7_dp)
   left0(1, 1) = c
   left0(1, 2) = -s
   left0(2, 1) = s
   left0(2, 2) = c
   left0(3, 3) = 1.0_dp
   c = cos(0.4_dp)
   s = sin(0.4_dp)
   right0(1, 1) = c
   right0(1, 3) = -s
   right0(2, 2) = 1.0_dp
   right0(3, 1) = s
   right0(3, 3) = c
   diag0(1, 1) = 1.0_dp
   diag0(2, 2) = 1.0e-8_dp
   diag0(3, 3) = 1.0e-12_dp
   a3 = matmul(matmul(left0, diag0), transpose(right0))
   call svd2(a3, svd, method=sanic_svd_jacobi)
   call require(svd%info == 0, "ill-conditioned SVD status")
   call require(abs(svd%d(1) - 1.0_dp) < 1.0e-13_dp, "largest ill-conditioned singular value")
   call require(abs(svd%d(2) - 1.0e-8_dp) < 1.0e-18_dp, "middle ill-conditioned singular value")
   call require(abs(svd%d(3) - 1.0e-12_dp) < 1.0e-22_dp, "small ill-conditioned singular value")
   sigma_tall = 0.0_dp
   do i = 1, 3
      sigma_tall(i, i) = svd%d(i)
   end do
   recon = matmul(matmul(svd%u, sigma_tall), transpose(svd%v))
   call require(maxval(abs(recon - a3)) < 2.0e-15_dp, "ill-conditioned SVD reconstruction")

   call svd2(a3, svd, method=sanic_svd_bdc)
   call require(svd%info == 0, "ill-conditioned BDC-like SVD status")
   call require(abs(svd%d(1) - 1.0_dp) < 1.0e-13_dp, "BDC-like largest singular value")
   call require(abs(svd%d(2) - 1.0e-8_dp) < 1.0e-18_dp, "BDC-like middle singular value")
   call require(abs(svd%d(3) - 1.0e-12_dp) < 1.0e-22_dp, "BDC-like smallest singular value")

   print '(a)', 'test_spectral: PASS'

contains

   subroutine require(condition, label)
      logical, intent(in) :: condition !! True when the regression condition has passed.
      character(len=*), intent(in) :: label !! Short diagnostic label printed on failure.
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAILED:', trim(label)
         error stop 1
      end if
   end subroutine require

end program test_spectral
