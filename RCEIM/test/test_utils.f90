program test_utils
    use rceim, only : dp, enforce_domain, sort_population_by_score
    implicit none
    real(dp) :: p(3,2), lo(2), hi(2), s(3)

    lo = [-1.0_dp, 0.0_dp]
    hi = [ 1.0_dp, 2.0_dp]
    p(1,:) = [-2.0_dp, 1.0_dp]
    p(2,:) = [ 0.0_dp, 3.0_dp]
    p(3,:) = [ 2.0_dp,-1.0_dp]
    call enforce_domain(p, lo, hi)
    if (maxval(abs(p(1,:)-[-1.0_dp,1.0_dp])) > 1.0e-14_dp) error stop 1
    if (maxval(abs(p(2,:)-[ 0.0_dp,2.0_dp])) > 1.0e-14_dp) error stop 2
    if (maxval(abs(p(3,:)-[ 1.0_dp,0.0_dp])) > 1.0e-14_dp) error stop 3

    s = [3.0_dp,1.0_dp,2.0_dp]
    p(:,1) = [30.0_dp,10.0_dp,20.0_dp]
    p(:,2) = 0.0_dp
    call sort_population_by_score(p,s)
    if (maxval(abs(s-[1.0_dp,2.0_dp,3.0_dp])) > 1.0e-14_dp) error stop 4
    if (maxval(abs(p(:,1)-[10.0_dp,20.0_dp,30.0_dp])) > 1.0e-14_dp) error stop 5
end program test_utils
