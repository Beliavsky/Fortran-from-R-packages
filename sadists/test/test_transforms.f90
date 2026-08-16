program test_transforms
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use sadists
    implicit none
    integer :: failures
    real(dp) :: d1(3), d2(3), n1(3), n2(3)
    real(dp) :: xv(3), out(3)

    failures = 0

    call check(ddnbeta(0.4_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp), &
        6.314006043963866_dp, 2.0e-9_dp, 'ddnbeta', failures)
    call check(pdnbeta(0.4_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp), &
        0.5974452323308629_dp, 2.0e-9_dp, 'pdnbeta', failures)
    call check(qdnbeta(0.7_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp), &
        0.41711082301940877_dp, 2.0e-9_dp, 'qdnbeta', failures)

    call check(ddneta(0.2_dp, 40.0_dp, 1.2_dp, 2.0_dp), &
        2.6791125228912076_dp, 2.0e-10_dp, 'ddneta', failures)
    call check(pdneta(0.2_dp, 40.0_dp, 1.2_dp, 2.0_dp), &
        0.5450687330208478_dp, 2.0e-10_dp, 'pdneta', failures)
    call check(qdneta(0.7_dp, 40.0_dp, 1.2_dp, 2.0_dp), &
        0.2600581051107843_dp, 2.0e-10_dp, 'qdneta', failures)

    call check(dprodchisqpow(50.0_dp, [20.0_dp, 30.0_dp], &
        [2.0_dp, 3.0_dp], [1.0_dp, 0.5_dp]), &
        0.001394014881694272_dp, 3.0e-9_dp, 'dprodchisqpow', failures)
    call check(pprodchisqpow(50.0_dp, [20.0_dp, 30.0_dp], &
        [2.0_dp, 3.0_dp], [1.0_dp, 0.5_dp]), &
        0.011539357519794959_dp, 3.0e-9_dp, 'pprodchisqpow', failures)
    call check(qprodchisqpow(0.7_dp, [20.0_dp, 30.0_dp], &
        [2.0_dp, 3.0_dp], [1.0_dp, 0.5_dp]), &
        143.2559835673459_dp, 3.0e-9_dp, 'qprodchisqpow', failures)

    d1 = [10.0_dp, 20.0_dp, 30.0_dp]
    d2 = [100.0_dp, 80.0_dp, 60.0_dp]
    n1 = [1.0_dp, 0.5_dp, 2.0_dp]
    n2 = [0.2_dp, 1.0_dp, 0.0_dp]
    call check(dproddnf(0.02_dp, d1, d2, n1, n2), &
        4.740314150240827e-5_dp, 3.0e-10_dp, 'dproddnf', failures)
    call check(pproddnf(0.02_dp, d1, d2, n1, n2), &
        1.3043819867973512e-7_dp, 3.0e-10_dp, 'pproddnf', failures)
    call check(qproddnf(0.7_dp, d1, d2, n1, n2), &
        1.4861997730282381_dp, 3.0e-10_dp, 'qproddnf', failures)

    call check(qdnbeta(0.0_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp), &
        0.0_dp, 0.0_dp, 'qdnbeta p=0', failures)
    call check(qdnbeta(1.0_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp), &
        1.0_dp, 0.0_dp, 'qdnbeta p=1', failures)
    call check(qdneta(0.0_dp, 40.0_dp, 1.2_dp, 2.0_dp), &
        -1.0_dp, 0.0_dp, 'qdneta p=0', failures)
    call check(qdneta(1.0_dp, 40.0_dp, 1.2_dp, 2.0_dp), &
        1.0_dp, 0.0_dp, 'qdneta p=1', failures)

    call check(pdnf(1.1_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp, &
        lower_tail=.false.), 1.0_dp - 0.6433328491074706_dp, &
        2.0e-9_dp, 'pdnf upper', failures)
    call check(exp(pdnf(1.1_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp, &
        lower_tail=.false., log_p=.true.)), &
        1.0_dp - 0.6433328491074706_dp, 2.0e-9_dp, 'pdnf log upper', failures)

    xv = [0.8_dp, 1.0_dp, 1.2_dp]
    call ddnf_vec(xv, out, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp)
    if (.not. all(ieee_is_finite(out))) then
        print *, 'ddnf_vec produced non-finite values'
        failures = failures + 1
    end if

    if (failures == 0) then
        print *, 'test_transforms: PASS'
    else
        print *, 'test_transforms: FAIL', failures
        error stop 1
    end if

contains

    subroutine check(actual, expected, tol, name, failures)
        real(dp), intent(in) :: actual, expected, tol
        character(*), intent(in) :: name
        integer, intent(inout) :: failures
        if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            print '(a,2es24.16)', trim(name)//': ', actual, expected
            failures = failures + 1
        end if
    end subroutine check

end program test_transforms
