program tree_fit_demo
    use clue, only: dp, fit_ultrametric_ip, non_ultrametricity, &
        fit_ultrametric_target_l2, is_ultrametric
    implicit none
    real(dp) :: d(4,4)
    real(dp), allocatable :: target(:,:)
    integer :: merge(3,2)

    d = reshape([ &
        0.0_dp,1.0_dp,4.0_dp,5.0_dp, &
        1.0_dp,0.0_dp,3.0_dp,4.0_dp, &
        4.0_dp,3.0_dp,0.0_dp,2.0_dp, &
        5.0_dp,4.0_dp,2.0_dp,0.0_dp ], [4,4])

    print '(a,es12.4)', 'Initial ultrametric deviation: ', non_ultrametricity(d)
    call fit_ultrametric_ip(d,maxiter=10000,tol=1.0e-10_dp)
    print '(a,l1)', 'IP fit is ultrametric: ', is_ultrametric(d,1.0e-8_dp)

    merge(1,:) = [-1,-2]
    merge(2,:) = [-3,-4]
    merge(3,:) = [1,2]
    target = fit_ultrametric_target_l2(d,merge)
    print '(a,l1)', 'Target fit is ultrametric: ', is_ultrametric(target,1.0e-8_dp)
end program tree_fit_demo
