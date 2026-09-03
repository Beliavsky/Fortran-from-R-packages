program test_diagnostics
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use mice, only : dp, mice_ok, md_pairs_result, flux_result, md_pattern_result, md_pairs, flux, md_pattern, quickpred
    implicit none
    real(dp) :: data(4, 3), nan
    integer, allocatable :: pred(:, :)
    type(md_pairs_result) :: pairs
    type(flux_result) :: fl
    type(md_pattern_result) :: pat
    integer :: info

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    data(:, 1) = [1.0_dp, 2.0_dp, nan, 4.0_dp]
    data(:, 2) = [2.0_dp, nan, 6.0_dp, 8.0_dp]
    data(:, 3) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    call md_pairs(data, pairs, info)
    if (info /= mice_ok) error stop "md pairs status"
    if (pairs%rr(1, 2) /= 2 .or. pairs%mr(1, 2) /= 1 .or. pairs%rm(1, 2) /= 1) error stop "md pairs counts"
    call flux(data, fl, info)
    if (info /= mice_ok .or. size(fl%influx) /= 3) error stop "flux"
    call md_pattern(data, pat, info)
    if (info /= mice_ok .or. sum(pat%frequency) /= 4) error stop "md pattern"
    call quickpred(data, 0.1_dp, 0.0_dp, pred, info)
    if (info /= mice_ok .or. pred(1, 1) /= 0 .or. pred(2, 2) /= 0 .or. pred(3, 3) /= 0) error stop "quickpred diagonal"
    print *, "test_diagnostics: PASS"
end program test_diagnostics
