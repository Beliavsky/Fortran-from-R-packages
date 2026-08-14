program basic
    use desirability, only : dp, d_max, d_target, d_overall, d_overall_type, &
        hold, predict, predict_all
    implicit none

    type(d_overall_type) :: overall
    real(dp) :: outcomes(1, 2)
    real(dp), allocatable :: details(:, :), score(:)

    ! Example from the original package vignette.
    overall = d_overall([hold(d_max(80.0_dp, 97.0_dp)), &
        hold(d_target(55.0_dp, 57.5_dp, 60.0_dp))])

    outcomes(1, :) = [81.09_dp, 59.85_dp]
    details = predict_all(overall, outcomes)
    score = predict(overall, outcomes)

    print '(a,f10.6)', "conversion desirability: ", details(1, 1)
    print '(a,f10.6)', "activity desirability:   ", details(1, 2)
    print '(a,f10.6)', "overall desirability:    ", score(1)
end program basic
