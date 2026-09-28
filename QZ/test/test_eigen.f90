program test_eigen
    use qz_api, only: dp, qz_dgeev, qz_zgeev, qz_dggev, qz_zggev, &
        qz_eigen_real_result, qz_eigen_complex_result, qz_geigen_real_result, qz_geigen_complex_result
    implicit none

    real(dp) :: a(2, 2), b(2, 2)
    complex(dp) :: ac(2, 2), bc(2, 2)
    type(qz_eigen_real_result) :: er
    type(qz_eigen_complex_result) :: ec
    type(qz_geigen_real_result) :: ger
    type(qz_geigen_complex_result) :: gec
    integer :: j

    a = reshape([0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp], [2, 2])
    call qz_dgeev(a, er)
    call check(er%info == 0, 'qz_dgeev info')
    call check(maxval(abs(matmul(cmplx(a, 0.0_dp, kind=dp), er%v) - &
        er%v * spread(er%w, 1, 2))) < 1.0e-11_dp, 'qz_dgeev residual')
    call check(abs(sum(er%w)) < 1.0e-12_dp, 'qz_dgeev eigenvalue sum')
    call check(abs(abs(aimag(er%w(1))) - 1.0_dp) < 1.0e-12_dp, 'qz_dgeev complex pair')

    ac = cmplx(0.0_dp, 0.0_dp, kind=dp)
    ac(1, 1) = cmplx(1.0_dp, 2.0_dp, kind=dp)
    ac(2, 2) = cmplx(3.0_dp, -1.0_dp, kind=dp)
    call qz_zgeev(ac, ec)
    call check(ec%info == 0, 'qz_zgeev info')
    do j = 1, 2
        call check(maxval(abs(matmul(ac, ec%v(:, j)) - ec%w(j) * ec%v(:, j))) < 1.0e-11_dp, &
            'qz_zgeev residual')
    end do

    a = 0.0_dp
    b = 0.0_dp
    a(1, 1) = 2.0_dp
    a(2, 2) = 3.0_dp
    b(1, 1) = 1.0_dp
    b(2, 2) = 2.0_dp
    call qz_dggev(a, b, ger)
    call check(ger%info == 0, 'qz_dggev info')
    do j = 1, 2
        call check(maxval(abs(matmul(cmplx(a, 0.0_dp, kind=dp), ger%v(:, j)) - &
            ger%lambda(j) * matmul(cmplx(b, 0.0_dp, kind=dp), ger%v(:, j)))) < 1.0e-11_dp, &
            'qz_dggev residual')
    end do
    call check(any(abs(ger%lambda - cmplx(2.0_dp, 0.0_dp, kind=dp)) < 1.0e-12_dp), &
        'qz_dggev lambda 2')
    call check(any(abs(ger%lambda - cmplx(1.5_dp, 0.0_dp, kind=dp)) < 1.0e-12_dp), &
        'qz_dggev lambda 1.5')

    ac = cmplx(a, 0.0_dp, kind=dp)
    bc = cmplx(b, 0.0_dp, kind=dp)
    ac(1, 2) = cmplx(0.0_dp, 0.25_dp, kind=dp)
    call qz_zggev(ac, bc, gec)
    call check(gec%info == 0, 'qz_zggev info')
    do j = 1, 2
        call check(maxval(abs(matmul(ac, gec%v(:, j)) - gec%lambda(j) * matmul(bc, gec%v(:, j)))) < 1.0e-10_dp, &
            'qz_zggev residual')
    end do

    print '(a)', 'test_eigen: PASS'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition !! Condition that must hold for the deterministic regression test.
        character(len=*), intent(in) :: message !! Failure label printed by ERROR STOP context.
        if (.not. condition) then
            print '(a)', trim(message)
            error stop 1
        end if
    end subroutine check
end program test_eigen
