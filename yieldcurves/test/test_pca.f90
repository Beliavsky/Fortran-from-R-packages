! SPDX-License-Identifier: MIT
program test_pca
  use yieldcurves
  implicit none
  real(dp), parameter :: x(6,3) = reshape([ &
    0.01_dp,0.02_dp,0.03_dp,0.04_dp,0.05_dp,0.06_dp, &
    0.02_dp,0.01_dp,0.04_dp,0.03_dp,0.06_dp,0.05_dp, &
    0.03_dp,0.04_dp,0.02_dp,0.05_dp,0.04_dp,0.07_dp], [6,3])
  real(dp), parameter :: variance_ref(3) = [0.757727354918054_dp, &
    0.230424844466192_dp, 0.011847800615754_dp]
  real(dp), parameter :: center_ref(3) = [0.035_dp, 0.035_dp, 0.0416666666666667_dp]
  type(pca_result_t) :: fit, scaled_fit
  real(dp) :: gram(3,3)
  integer :: i

  fit = yc_pca(x, 3)
  call require(fit%ok, 'PCA status')
  call close_vector(fit%center, center_ref, 1.0e-14_dp, 'PCA center')
  call close_vector(fit%variance_explained, variance_ref, 1.0e-8_dp, 'PCA variance')
  gram = matmul(transpose(fit%loadings),fit%loadings)
  do i=1,3
    call close_scalar(gram(i,i),1.0_dp,1.0e-12_dp,'loading norm')
  end do
  call close_scalar(sum(fit%variance_explained),1.0_dp,1.0e-12_dp,'variance sum')
  call close_scalar(fit%cumulative_variance(3),1.0_dp,1.0e-12_dp,'cumulative variance')
  scaled_fit = yc_pca(x, 2, .true.)
  call require(scaled_fit%ok, 'scaled PCA status')
  call require(all(scaled_fit%scale > 0.0_dp), 'scaled PCA standard deviations')
  print '(a)', 'test_pca: PASS'
contains
  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'failed: '//trim(label)
      error stop 1
    end if
  end subroutine require
  subroutine close_scalar(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual-expected) > tolerance) then
      write(*,'(a,3es24.15)') trim(label)//': ', actual, expected, abs(actual-expected)
      error stop 1
    end if
  end subroutine close_scalar
  subroutine close_vector(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: label
    if (size(actual) /= size(expected) .or. maxval(abs(actual-expected)) > tolerance) then
      write(*,'(a,es24.15)') trim(label)//' maximum error: ', maxval(abs(actual-expected))
      error stop 1
    end if
  end subroutine close_vector
end program test_pca
