! SPDX-License-Identifier: GPL-2.0-or-later
program local_level
    use dlm, only : dp, dlm_model, dlm_filter_result, dlm_smooth_result
    use dlm, only : dlm_success, dlm_filter, dlm_mod_poly, dlm_smooth
    implicit none

    type(dlm_model) :: model
    type(dlm_filter_result) :: filtered
    type(dlm_smooth_result) :: smoothed
    real(dp) :: y(6, 1)
    integer :: info

    y(:, 1) = [0.9_dp, 1.2_dp, 1.0_dp, 1.4_dp, 1.3_dp, 1.6_dp]
    call dlm_mod_poly(1, model, info, d_v=0.25_dp, d_w=[0.05_dp], &
                      m0=[0.0_dp], c0=reshape([1.0_dp], [1, 1]))
    if (info /= dlm_success) error stop "model construction failed"
    call dlm_filter(y, model, filtered, info)
    if (info /= dlm_success) error stop "filter failed"
    call dlm_smooth(filtered, smoothed, info)
    if (info /= dlm_success) error stop "smoother failed"

    print '(a,f10.5)', "filtered final state = ", filtered%m(size(filtered%m, 1), 1)
    print '(a,f10.5)', "smoothed initial state = ", smoothed%s(1, 1)
end program local_level
