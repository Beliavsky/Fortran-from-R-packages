program basic_survey
    use survey
    implicit none

    type(survey_design_t) :: d
    type(svystat_t) :: m
    real(dp) :: y(4,1), w(4)
    integer :: psu(4,1), strata(4,1)

    y(:,1) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
    w = [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp]
    psu(:,1) = [1, 2, 3, 4]
    strata(:,1) = 1

    call make_design(w, psu, d, strata=strata)
    m = svy_mean(y, d)

    print '(a,f12.6)', 'weighted mean = ', m%estimate(1)
    print '(a,f12.6)', 'SE            = ', sqrt(m%variance(1,1))
end program basic_survey
