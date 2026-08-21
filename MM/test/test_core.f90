program test_core
    use multiplicative_multinomial
    implicit none
    type(paras_type) :: par
    type(suffstats_type) :: ss, ess
    real(dp) :: th(3,3), prob(3), weights(3), ll, expected_ll, c
    integer :: y3(3), ymat(3,3), y2(2)

    call check_close(lmultinomial([2, 1]), log(3.0_dp), 1.0e-13_dp, "lmultinomial")
    call check_close(multinomial([2, 1]), 3.0_dp, 1.0e-13_dp, "multinomial")

    par = paras(2)
    y2 = [1, 1]
    th = 1.0_dp
    call check_close(mm_single([2, 0], par), 0.25_dp, 1.0e-13_dp, "ordinary multinomial")
    call check_close(normc(2, par), 1.0_dp, 1.0e-13_dp, "theta=1 normalizer")

    par = paras(2)
    call set_theta(par, reshape([1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], [2,2]))
    call check_close(dmm(y2, par), 2.0_dp / 3.0_dp, 1.0e-13_dp, "simple dMM")

    prob = [0.2_dp, 0.3_dp, 0.5_dp]
    th = 1.0_dp
    th(1,2) = 1.4_dp
    th(1,3) = 0.8_dp
    th(2,3) = 1.7_dp
    par = paras_from_p_theta(prob, th)
    call check_close(normc(3, par), 1.853920000000001_dp, 2.0e-13_dp, "reference normalizer")
    y3 = [1, 1, 1]
    call check_close(dmm(y3, par), 0.18486234573228613_dp, 2.0e-13_dp, "reference density")

    ymat(1,:) = [2, 1, 0]
    ymat(2,:) = [1, 1, 1]
    ymat(3,:) = [0, 2, 1]
    weights = [2.0_dp, 3.0_dp, 1.0_dp]
    ll = mm_loglik(ymat, par, weights)
    expected_ll = -13.16014610575452_dp
    call check_close(ll, expected_ll, 5.0e-12_dp, "reference loglik")

    ss = suffstats(ymat, weights)
    call check_vector(ss%row_sums, [7.0_dp, 7.0_dp, 4.0_dp], 1.0e-13_dp, "row sums")
    call check_close(ss%cross_prods(1,1), 11.0_dp, 1.0e-13_dp, "cross 11")
    call check_close(ss%cross_prods(1,2), 7.0_dp, 1.0e-13_dp, "cross 12")
    call check_close(ss%cross_prods(2,2), 9.0_dp, 1.0e-13_dp, "cross 22")
    c = 6.0_dp * log_gamma(4.0_dp)
    c = c - 2.0_dp * sum(log_gamma(real(ymat(1,:) + 1, dp)))
    c = c - 3.0_dp * sum(log_gamma(real(ymat(2,:) + 1, dp)))
    c = c - 1.0_dp * sum(log_gamma(real(ymat(3,:) + 1, dp)))
    call check_close(mm_support(par, ss) + c, ll, 5.0e-12_dp, "support/full likelihood relation")

    ess = expected_suffstats(par, 3)
    call check_vector(ess%row_sums, [0.42422543_dp, 1.15242837_dp, 1.42334621_dp], &
        5.0e-9_dp, "expected row sums")
    call check_close(ess%cross_prods(2,3), 1.30724087_dp, 5.0e-9_dp, "expected cross product")

    print '(a)', 'test_core: PASS'

contains

    subroutine check_close(x, ref, tol, label)
        real(dp), intent(in) :: x, ref, tol
        character(len=*), intent(in) :: label
        if (abs(x - ref) > tol * max(1.0_dp, abs(ref))) then
            write(*,'(a,2es24.14)') trim(label)//' failed: ', x, ref
            error stop 1
        end if
    end subroutine check_close

    subroutine check_vector(x, ref, tol, label)
        real(dp), intent(in) :: x(:), ref(:), tol
        character(len=*), intent(in) :: label
        if (size(x) /= size(ref) .or. maxval(abs(x-ref)) > tol * max(1.0_dp, maxval(abs(ref)))) then
            write(*,'(a)') trim(label)//' failed'
            error stop 1
        end if
    end subroutine check_vector

end program test_core
