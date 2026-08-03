program test_alternative_simulations
  use, intrinsic :: iso_fortran_env, only : int64
  use mfgarch
  implicit none
  type(mfgarch_model) :: model
  type(mfgarch_simulation) :: rv_simulation, diffusion_simulation
  integer :: status

  model%asymmetric = .false.
  model%gamma = 0.0_dp
  model%k = 2
  model%alpha = 0.06_dp
  model%beta = 0.90_dp
  model%m = -0.2_dp
  model%theta = 0.001_dp
  model%w1 = 1.0_dp
  model%w2 = 3.0_dp

  call simulate_mfgarch_rv_dependent(40, model, 5, 48, .false., rv_simulation, status, &
    seed=777_int64)
  call assert_true(status == mfgarch_success, 'RV-dependent status')
  call assert_true(size(rv_simulation%returns) == 40, 'RV-dependent dimensions')
  call assert_true(all(rv_simulation%tau > 0.0_dp), 'RV-dependent tau positive')

  call simulate_mfgarch_diffusion(20, model, 0.7_dp, 0.15_dp, 5, 12, &
    diffusion_simulation, status, seed=888_int64)
  call assert_true(status == mfgarch_success, 'diffusion status')
  call assert_true(size(diffusion_simulation%returns) == 20, 'diffusion dimensions')
  call assert_true(all(diffusion_simulation%tau > 0.0_dp), 'diffusion tau positive')
  call assert_true(all(diffusion_simulation%realized_variance >= 0.0_dp), 'diffusion RV positive')

  print '(a)', 'test_alternative_simulations: PASS'

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_alternative_simulations
