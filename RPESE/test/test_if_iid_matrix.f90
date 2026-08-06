program test_if_iid_matrix
  use rpese, only : dp, pi, rpese_options, se_result, matrix_se_result, &
    rpese_success, se_if_iid, estimate_se, estimate_se_matrix
  implicit none
  integer, parameter :: n = 80
  real(dp) :: x(n), data(n, 2), expected
  type(rpese_options) :: options
  type(se_result) :: result
  type(matrix_se_result) :: matrix_result
  integer :: i

  do i = 1, n
    x(i) = 0.01_dp + 0.02_dp * sin(2.0_dp * pi * real(i, dp) / 13.0_dp) + &
      0.005_dp * cos(2.0_dp * pi * real(i, dp) / 7.0_dp)
  end do
  options = rpese_options()
  call estimate_se(x, 'mean', se_if_iid, result, options)
  call assert_true(result%status == rpese_success, 'mean IF iid status')
  expected = sqrt(sum((x - sum(x) / real(n, dp)) ** 2) / real(n * n, dp))
  call assert_close(result%standard_error, expected, 1.0e-12_dp, 'mean IF iid SE')
  call assert_true(abs(result%return_correlation) <= 1.0_dp, 'return correlation range')

  data(:, 1) = x
  data(:, 2) = 2.0_dp * x + 0.003_dp
  call estimate_se_matrix(data, 'sd', se_if_iid, matrix_result, options)
  call assert_true(matrix_result%status == rpese_success, 'matrix status')
  call assert_true(size(matrix_result%column) == 2, 'matrix column count')
  call assert_close(matrix_result%column(2)%estimate, 2.0_dp * matrix_result%column(1)%estimate, &
    1.0e-12_dp, 'matrix SD scaling')

  print '(a)', 'test_if_iid_matrix: PASS'
contains
  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      print '(a,2es24.14)', trim(label) // ' failed: ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', trim(label) // ' failed'
      error stop 1
    end if
  end subroutine assert_true
end program test_if_iid_matrix
