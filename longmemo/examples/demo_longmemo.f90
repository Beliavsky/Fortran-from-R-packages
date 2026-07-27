! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

program demo_longmemo
    use longmemo_kinds, only : dp
    use longmemo_stats, only : set_random_seed
    use longmemo
    implicit none

    real(dp), allocatable :: x(:)
    real(dp) :: eta(1)
    type(spectrum_result) :: theoretical
    type(whittle_result) :: fit
    type(fexp_result) :: fexp

    call set_random_seed(20260722)
    call sim_fgn_fft(1024, 0.75_dp, x)

    eta = [0.75_dp]
    call spec_fgn(eta, size(x), theoretical)
    call whittle_estimate(x, "fGn", fit, start_eta=[0.6_dp], covariance_m=2048)
    call fexp_estimate(x, order_poly=2, pvalmax=0.5_dp, result=fexp)

    print '(a,i0)', "observations:       ", size(x)
    print '(a,f10.6)', "true H:             ", eta(1)
    print '(a,f10.6)', "Whittle H:          ", fit%eta(1)
    print '(a,f10.6)', "Whittle standard error: ", fit%std_error(1)
    print '(a,f10.6)', "goodness-of-fit p:  ", fit%goodness_of_fit%p_value
    print '(a,f10.6)', "FEXP H:             ", fexp%hurst
    print '(a,i0)', "FEXP order:         ", fexp%order_poly
    print '(a,f10.6)', "first normalized fGn spectrum value: ", theoretical%spectrum(1)
end program demo_longmemo
