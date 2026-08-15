program test_delaporte
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use delaporte, only : dp, ddelap, pdelap, qdelap, ddelap_vec, pdelap_vec, &
        qdelap_vec, rdelap, seed_delaporte, momdelap
    implicit none

    integer :: failures, status
    real(dp) :: params(3), x(30000), mean_x
    real(dp) :: out3(3)

    failures = 0

    call check_close(ddelap(4.0_dp, 0.5_dp, 4.0_dp, 0.2_dp), &
        0.0547024400602606_dp, 2.0e-14_dp, 'ddelap alpha < 0.8')
    call check_close(ddelap(1.0_dp, 1.0_dp, 1.0e-10_dp, 2.0_dp), &
        0.270670566459692_dp, 2.0e-14_dp, 'ddelap small beta')
    call check_close(ddelap(0.0_dp, 1.0_dp, 4.0_dp, 2.0_dp), &
        0.0270670566473225_dp, 2.0e-14_dp, 'ddelap x=0')
    call check_true(ddelap(1.1_dp, 1.0_dp, 2.0_dp, 3.0_dp) == 0.0_dp, &
        'ddelap noninteger is zero')
    call check_true(ddelap(huge(1.0_dp), 1.0_dp, 2.0_dp, 3.0_dp) == 0.0_dp, &
        'ddelap +large is zero')
    call check_true(ieee_is_nan(ddelap(1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp)), &
        'ddelap invalid parameter NaN')

    call check_close(pdelap(6.0_dp, 2.9647_dp, 0.005_dp / 2.9647_dp, &
        0.0057_dp, lower_tail=.false.), 0.0_dp, 2.0e-14_dp, &
        'pdelap upper tail near one')
    call check_close(qdelap(0.4_dp, 1.0_dp, 4.0_dp, 2.0_dp), 4.0_dp, &
        0.0_dp, 'qdelap basic')
    call check_close(qdelap(0.49_dp, 4.0_dp, 6.0_dp, 3.0_dp), 25.0_dp, &
        0.0_dp, 'qdelap lower tail')
    call check_close(qdelap(0.4_dp, 1.0_dp, 4.0_dp, 2.0_dp, &
        lower_tail=.false.), 6.0_dp, 0.0_dp, 'qdelap upper tail')
    call check_close(pdelap(huge(1.0_dp), 1.0_dp, 2.0_dp, 3.0_dp), &
        1.0_dp, 0.0_dp, 'pdelap huge q')
    call check_close(qdelap(-0.255_dp, 20.0_dp, 15.0_dp, 50.0_dp, &
        log_p=.true.), 400.0_dp, 0.0_dp, 'qdelap log p')
    call check_true(.not. ieee_is_finite(qdelap(1.0_dp, 3.0_dp, 1.0_dp, &
        2.0_dp)), 'qdelap p=1 is infinity')

    call ddelap_vec([0.0_dp, 1.0_dp, 2.0_dp], [1.0_dp, 2.0_dp, 3.0_dp], &
        [4.0_dp, 1.0_dp, 2.0_dp], [2.0_dp, 5.0_dp, 7.0_dp], out3)
    call check_close(out3(1), 0.0270670566473225_dp, 2.0e-14_dp, &
        'ddelap_vec 1')
    call check_close(out3(2), 0.0101069204986282_dp, 2.0e-14_dp, &
        'ddelap_vec 2')
    call check_close(out3(3), 0.0013903385524195401_dp, 2.0e-14_dp, &
        'ddelap_vec 3')

    call pdelap_vec([0.0_dp, 1.0_dp, 2.0_dp], [1.0_dp, 2.0_dp, 3.0_dp], &
        [4.0_dp, 1.0_dp, 2.0_dp], [2.0_dp, 5.0_dp, 7.0_dp], out3)
    call check_close(out3(1), 0.0270670566473225_dp, 2.0e-14_dp, &
        'pdelap_vec 1')
    call check_close(out3(2), 0.0117914072483995674_dp, 2.0e-14_dp, &
        'pdelap_vec 2')
    call check_close(out3(3), 0.0017280726137360276_dp, 2.0e-14_dp, &
        'pdelap_vec 3')

    call qdelap_vec([0.4_dp, 0.07_dp, 0.4_dp], [1.0_dp, 2.0_dp], &
        [4.0_dp, 1.0_dp], [2.0_dp, 5.0_dp], out3)
    call check_close(out3(1), 4.0_dp, 0.0_dp, 'qdelap_vec 1')
    call check_close(out3(2), 3.0_dp, 0.0_dp, 'qdelap_vec 2')

    call momdelap([5.0_dp, 7.0_dp, 9.0_dp, 9.0_dp, 10.0_dp, 11.0_dp, &
        11.0_dp, 13.0_dp, 17.0_dp, 24.0_dp], params, 2, status)
    call check_true(status == 0, 'momdelap status')
    call check_close(params(1), 0.88342721893491116_dp, 1.0e-12_dp, &
        'momdelap alpha')
    call check_close(params(2), 4.51388888888888928_dp, 1.0e-12_dp, &
        'momdelap beta')
    call check_close(params(3), 7.61230769230769155_dp, 1.0e-12_dp, &
        'momdelap lambda')

    call momdelap([5.0_dp, 7.0_dp, 9.0_dp, 9.0_dp, 10.0_dp, 11.0_dp, &
        11.0_dp, 13.0_dp, 17.0_dp, 24.0_dp], params, 1, status)
    call check_close(params(1), 1.4520314593282133_dp, 1.0e-12_dp, &
        'momdelap type1 alpha')
    call momdelap([5.0_dp, 7.0_dp, 9.0_dp, 9.0_dp, 10.0_dp, 11.0_dp, &
        11.0_dp, 13.0_dp, 17.0_dp, 24.0_dp], params, 3, status)
    call check_close(params(1), 2.3979594095667354_dp, 1.0e-12_dp, &
        'momdelap type3 alpha')

    call seed_delaporte(4175)
    call rdelap(size(x), 10.0_dp, 2.0_dp, 10.0_dp, x, exact=.false.)
    mean_x = sum(x) / real(size(x), dp)
    call check_true(abs(mean_x - 30.0_dp) < 0.35_dp, &
        'rdelap mixture sample mean')

    if (failures == 0) then
        print '(a)', 'test_delaporte: PASS'
    else
        print '(a,i0)', 'test_delaporte: FAIL ', failures
        error stop 1
    end if

contains

    subroutine check_close(actual, expected, tol, label)
        real(dp), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            failures = failures + 1
            print '(a,2(1x,es24.16))', 'FAIL '//trim(label)//':', actual, expected
        end if
    end subroutine check_close

    subroutine check_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL '//trim(label)
        end if
    end subroutine check_true

end program test_delaporte
