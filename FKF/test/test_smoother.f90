! SPDX-License-Identifier: GPL-2.0-or-later
program test_smoother
  use fkf_module
  use fkf_linalg, only : spd_inverse_logdet, symmetrize
  use test_support
  implicit none

  type(fkf_model) :: model
  type(fkf_result) :: fit
  type(fks_result) :: smooth
  real(dp) :: y(2, 12), x(2), noise(2)
  real(dp) :: ahat(2, 12), vhat(2, 2, 12), jmat(2, 2), pinv(2, 2), logdet
  integer :: i, info

  model%a0 = [0.0_dp, 0.0_dp]
  model%p0 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
  model%dt = reshape([0.0_dp, 0.0_dp], [2, 1])
  model%ct = reshape([0.0_dp, 0.0_dp], [2, 1])
  allocate(model%tt(2, 2, 1), model%zt(2, 2, 1), model%hht(2, 2, 1), model%ggt(2, 2, 1))
  model%tt(:, :, 1) = reshape([1.0_dp, 0.0_dp, 0.4_dp, 0.8_dp], [2, 2])
  model%zt(:, :, 1) = reshape([1.0_dp, 0.2_dp, 0.0_dp, 1.0_dp], [2, 2])
  model%hht(:, :, 1) = reshape([0.05_dp, 0.01_dp, 0.01_dp, 0.08_dp], [2, 2])
  model%ggt(:, :, 1) = reshape([0.20_dp, 0.03_dp, 0.03_dp, 0.25_dp], [2, 2])

  x = [0.3_dp, -0.2_dp]
  do i = 1, 12
    noise = [0.07_dp * sin(real(i, dp)), 0.05_dp * cos(0.7_dp * real(i, dp))]
    y(:, i) = matmul(model%zt(:, :, 1), x) + noise
    x = matmul(model%tt(:, :, 1), x) + [0.02_dp * cos(real(i, dp)), -0.01_dp * sin(real(i, dp))]
  end do

  call kalman_filter(model, y, fit, .true.)
  call kalman_smooth(model, y, fit, smooth)
  call assert_true(smooth%status == fkf_success, 'smoother status')

  ahat(:, 12) = fit%att(:, 12)
  vhat(:, :, 12) = fit%ptt(:, :, 12)
  do i = 11, 1, -1
    call spd_inverse_logdet(fit%pt(:, :, i + 1), pinv, logdet, info)
    call assert_true(info == 0, 'RTS predicted covariance invertible')
    jmat = matmul(matmul(fit%ptt(:, :, i), transpose(model%tt(:, :, 1))), pinv)
    ahat(:, i) = fit%att(:, i) + matmul(jmat, ahat(:, i + 1) - fit%at(:, i + 1))
    vhat(:, :, i) = fit%ptt(:, :, i) + &
      matmul(matmul(jmat, vhat(:, :, i + 1) - fit%pt(:, :, i + 1)), transpose(jmat))
    call symmetrize(vhat(:, :, i))
  end do

  call assert_close(maxval(abs(smooth%ahatt - ahat)), 0.0_dp, 2.0e-11_dp, 'smoother state vs RTS')
  call assert_close(maxval(abs(smooth%vt - vhat)), 0.0_dp, 2.0e-11_dp, 'smoother covariance vs RTS')
  call assert_true(all([(smooth%vt(1, 1, i) <= fit%pt(1, 1, i) + 1.0e-12_dp, i=1,12)]), &
    'smoothing does not increase first-state variance')
  call finish_test('Durbin-Koopman smoother versus RTS smoother')
end program test_smoother
