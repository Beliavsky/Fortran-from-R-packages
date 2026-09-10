! SPDX-License-Identifier: GPL-2.0-or-later
program test_models
    use dlm, only : dp, dlm_model, dlm_success
    use dlm, only : ar_trans_pars, dlm_mod_arma, dlm_mod_poly, dlm_mod_reg, dlm_mod_seas, dlm_mod_trig
    implicit none

    type(dlm_model) :: model
    real(dp) :: raw(3), coeff(3), x(3, 2)
    integer :: info

    raw = 0.0_dp
    call ar_trans_pars(raw, coeff)
    if (maxval(abs(coeff)) > 1.0e-14_dp) error stop "ARtransPars zero map failed"

    call dlm_mod_poly(2, model, info, d_v=2.0_dp, d_w=[0.0_dp, 0.25_dp])
    if (info /= dlm_success) error stop "dlmModPoly construction failed"
    if (any(shape(model%gg) /= [2, 2])) error stop "dlmModPoly GG shape failed"
    if (abs(model%gg(1, 2) - 1.0_dp) > 1.0e-14_dp) error stop "dlmModPoly transition failed"
    if (abs(model%w(2, 2) - 0.25_dp) > 1.0e-14_dp) error stop "dlmModPoly W failed"

    call dlm_mod_seas(4, model, info)
    if (info /= dlm_success) error stop "dlmModSeas construction failed"
    if (any(shape(model%gg) /= [3, 3])) error stop "dlmModSeas GG shape failed"
    if (maxval(abs(model%gg(1, :) + 1.0_dp)) > 1.0e-14_dp) error stop "dlmModSeas first row failed"

    x = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [3, 2])
    call dlm_mod_reg(x, model, info)
    if (info /= dlm_success) error stop "dlmModReg construction failed"
    if (size(model%m0) /= 3) error stop "dlmModReg state size failed"
    if (any(model%jff /= reshape([0, 1, 2], [1, 3]))) error stop "dlmModReg JFF failed"

    call dlm_mod_arma(model, info, ar=[0.7_dp], ma=[0.2_dp], sigma2=1.5_dp)
    if (info /= dlm_success) error stop "dlmModARMA construction failed"
    if (any(shape(model%gg) /= [2, 2])) error stop "dlmModARMA state size failed"
    if (abs(model%gg(1, 1) - 0.7_dp) > 1.0e-14_dp) error stop "dlmModARMA AR failed"
    if (abs(model%gg(1, 2) - 1.0_dp) > 1.0e-14_dp) error stop "dlmModARMA shift failed"
    if (abs(model%w(2, 2) - 0.06_dp) > 1.0e-14_dp) error stop "dlmModARMA MA covariance failed"

    call dlm_mod_trig(model, info, period=12, q=2)
    if (info /= dlm_success) error stop "dlmModTrig construction failed"
    if (size(model%m0) /= 4) error stop "dlmModTrig state size failed"
    if (abs(model%ff(1, 1) - 1.0_dp) > 1.0e-14_dp) error stop "dlmModTrig FF failed"
    if (abs(model%ff(1, 3) - 1.0_dp) > 1.0e-14_dp) error stop "dlmModTrig harmonic FF failed"

    print *, "test_models: PASS"
end program test_models
