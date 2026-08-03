program test_simulation
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch
  implicit none
  type(mfgarch_model) :: model
  type(mfgarch_simulation) :: first, second, student
  integer :: status

  model%asymmetric = .true.
  model%k = 3
  model%alpha = 0.06_dp
  model%beta = 0.88_dp
  model%gamma = 0.06_dp
  model%m = 0.0_dp
  model%theta = 0.25_dp
  model%w1 = 1.0_dp
  model%w2 = 3.0_dp

  call simulate_mfgarch(120, model, 0.8_dp, 0.25_dp, 5, 48, first, status, &
    seed=123456_int64, correlation=0.2_dp)
  call assert_true(status == mfgarch_success, 'normal simulation status')
  call simulate_mfgarch(120, model, 0.8_dp, 0.25_dp, 5, 48, second, status, &
    seed=123456_int64, correlation=0.2_dp)
  call assert_true(maxval(abs(first%returns-second%returns)) <= tiny(1.0_dp), 'reproducibility')
  call assert_true(size(first%returns) == 120, 'daily dimensions')
  call assert_true(size(first%intraday_returns) == 120*48, 'intraday dimensions')
  call assert_true(all(first%tau > 0.0_dp) .and. all(first%g > 0.0_dp), 'positive components')
  call assert_true(all(first%realized_variance >= 0.0_dp), 'positive realized variance')

  call simulate_mfgarch(60, model, 0.8_dp, 0.25_dp, 5, 48, student, status, &
    seed=654321_int64, student_t_df=6.0_dp)
  call assert_true(status == mfgarch_success, 'student t simulation status')
  call assert_true(size(student%returns) == 60, 'student t dimensions')

  print '(a)', 'test_simulation: PASS'

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_simulation
