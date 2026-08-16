program test_reference
    use sadists
    implicit none
    integer :: failures
    real(dp) :: a(2), b(2)

    failures = 0

    call check(ddnf(1.1_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp), &
        1.3306960107003345_dp, 2.0e-9_dp, 'ddnf', failures)
    call check(pdnf(1.1_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp), &
        0.6433328491074706_dp, 2.0e-9_dp, 'pdnf', failures)
    call check(qdnf(0.7_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp), &
        1.1449471755302056_dp, 2.0e-9_dp, 'qdnf', failures)

    call check(ddnt(0.7_dp, 40.0_dp, 1.2_dp, 2.0_dp), &
        0.36134244687231126_dp, 2.0e-10_dp, 'ddnt', failures)
    call check(pdnt(0.7_dp, 40.0_dp, 1.2_dp, 2.0_dp), &
        0.3136305848750700_dp, 2.0e-10_dp, 'pdnt', failures)
    call check(qdnt(0.7_dp, 40.0_dp, 1.2_dp, 2.0_dp), &
        1.7033594122130098_dp, 2.0e-10_dp, 'qdnt', failures)

    a = [10.0_dp, 30.0_dp]
    b = [0.5_dp, -0.2_dp]
    call check(dupsilon(0.8_dp, a, b), 0.3484975192972518_dp, &
        2.0e-12_dp, 'dupsilon', failures)
    call check(pupsilon(0.8_dp, a, b), 0.6940833519024893_dp, &
        2.0e-12_dp, 'pupsilon', failures)
    call check(qupsilon(0.7_dp, a, b), 0.8170515059749764_dp, &
        2.0e-12_dp, 'qupsilon', failures)

    call check(dlambdap(1.0_dp, 25.0_dp, 1.3_dp), &
        0.37713705728165425_dp, 2.0e-12_dp, 'dlambdap', failures)
    call check(plambdap(1.0_dp, 25.0_dp, 1.3_dp), &
        0.3888744067371140_dp, 2.0e-12_dp, 'plambdap', failures)
    call check(qlambdap(0.7_dp, 25.0_dp, 1.3_dp), &
        1.8200648148751357_dp, 2.0e-12_dp, 'qlambdap', failures)

    call check(dkprime(1.0_dp, 50.0_dp, 60.0_dp, 1.1_dp, 0.8_dp), &
        0.4858419836785802_dp, 2.0e-10_dp, 'dkprime', failures)
    call check(pkprime(1.0_dp, 50.0_dp, 60.0_dp, 1.1_dp, 0.8_dp), &
        0.45167383607661166_dp, 2.0e-10_dp, 'pkprime', failures)
    call check(qkprime(0.7_dp, 50.0_dp, 60.0_dp, 1.1_dp, 0.8_dp), &
        1.5306735019490592_dp, 2.0e-10_dp, 'qkprime', failures)

    a = [1.0_dp, -0.4_dp]
    b = [30.0_dp, 20.0_dp]
    call check(dsumchisqpow(25.0_dp, a, b, [2.0_dp, 1.0_dp], &
        [1.0_dp, 0.5_dp]), 0.04591335810491421_dp, 3.0e-11_dp, &
        'dsumchisqpow', failures)
    call check(psumchisqpow(25.0_dp, a, b, [2.0_dp, 1.0_dp], &
        [1.0_dp, 0.5_dp]), 0.2811274064982566_dp, 3.0e-11_dp, &
        'psumchisqpow', failures)
    call check(qsumchisqpow(0.7_dp, a, b, [2.0_dp, 1.0_dp], &
        [1.0_dp, 0.5_dp]), 33.95470477814333_dp, 3.0e-9_dp, &
        'qsumchisqpow', failures)

    a = [1.0_dp, -0.5_dp]
    b = [20.0_dp, 30.0_dp]
    call check(dsumlogchisq(1.2_dp, a, b, [2.0_dp, 3.0_dp]), &
        1.056717787506512_dp, 3.0e-9_dp, 'dsumlogchisq', failures)
    call check(psumlogchisq(1.2_dp, a, b, [2.0_dp, 3.0_dp]), &
        0.36250226446532785_dp, 3.0e-9_dp, 'psumlogchisq', failures)
    call check(qsumlogchisq(0.7_dp, a, b, [2.0_dp, 3.0_dp]), &
        1.5001181983507512_dp, 3.0e-9_dp, 'qsumlogchisq', failures)

    a = [1.2_dp, -0.7_dp]
    b = [0.5_dp, 0.3_dp]
    call check(dprodnormal(-0.7_dp, a, b), 0.8144698077341044_dp, &
        2.0e-12_dp, 'dprodnormal', failures)
    call check(pprodnormal(-0.7_dp, a, b), 0.5556987369016666_dp, &
        2.0e-12_dp, 'pprodnormal', failures)
    call check(qprodnormal(0.7_dp, a, b), -0.5229375053057881_dp, &
        2.0e-12_dp, 'qprodnormal', failures)

    if (failures == 0) then
        print *, 'test_reference: PASS'
    else
        print *, 'test_reference: FAIL', failures
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

end program test_reference
