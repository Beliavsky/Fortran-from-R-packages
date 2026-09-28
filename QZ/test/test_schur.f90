program test_schur
    use qz_api, only: dp, qz_dgees, qz_zgees, qz_dgges, qz_zgges, &
        qz_schur_real_result, qz_schur_complex_result, qz_gschur_real_result, qz_gschur_complex_result
    implicit none

    real(dp) :: a(3, 3), ar(3, 3), br(3, 3)
    complex(dp) :: ac(3, 3), bc(3, 3)
    type(qz_schur_real_result) :: sr
    type(qz_schur_complex_result) :: sc
    type(qz_gschur_real_result) :: gr
    type(qz_gschur_complex_result) :: gc

    a = reshape([1.0_dp, 2.0_dp, 0.0_dp, -1.0_dp, 3.0_dp, 1.0_dp, &
        0.5_dp, -2.0_dp, 4.0_dp], [3, 3])
    call qz_dgees(a, sr)
    call check(sr%info == 0, 'qz_dgees info')
    call check(maxval(abs(a - matmul(sr%q, matmul(sr%t, transpose(sr%q))))) < 1.0e-10_dp, &
        'qz_dgees reconstruction')
    call check(maxval(abs(matmul(transpose(sr%q), sr%q) - eye3())) < 1.0e-11_dp, &
        'qz_dgees orthogonality')

    ac = cmplx(a, 0.0_dp, kind=dp)
    ac(1, 3) = cmplx(0.5_dp, 0.75_dp, kind=dp)
    call qz_zgees(ac, sc)
    call check(sc%info == 0, 'qz_zgees info')
    call check(maxval(abs(ac - matmul(sc%q, matmul(sc%t, conjg(transpose(sc%q)))))) < 1.0e-10_dp, &
        'qz_zgees reconstruction')

    ar = a
    br = 0.0_dp
    br(1, 1) = 2.0_dp
    br(2, 2) = 1.5_dp
    br(3, 3) = 1.0_dp
    br(1, 2) = 0.1_dp
    call qz_dgges(ar, br, gr)
    call check(gr%info == 0, 'qz_dgges info')
    call check(maxval(abs(ar - matmul(gr%q, matmul(gr%s, transpose(gr%z))))) < 1.0e-10_dp, &
        'qz_dgges A reconstruction')
    call check(maxval(abs(br - matmul(gr%q, matmul(gr%t, transpose(gr%z))))) < 1.0e-10_dp, &
        'qz_dgges B reconstruction')

    ac = cmplx(ar, 0.0_dp, kind=dp)
    bc = cmplx(br, 0.0_dp, kind=dp)
    ac(2, 3) = ac(2, 3) + cmplx(0.0_dp, 0.4_dp, kind=dp)
    bc(1, 3) = cmplx(0.0_dp, -0.15_dp, kind=dp)
    call qz_zgges(ac, bc, gc)
    call check(gc%info == 0, 'qz_zgges info')
    call check(maxval(abs(ac - matmul(gc%q, matmul(gc%s, conjg(transpose(gc%z)))))) < 1.0e-10_dp, &
        'qz_zgges A reconstruction')
    call check(maxval(abs(bc - matmul(gc%q, matmul(gc%t, conjg(transpose(gc%z)))))) < 1.0e-10_dp, &
        'qz_zgges B reconstruction')

    print '(a)', 'test_schur: PASS'
contains
    pure function eye3() result(x)
        real(dp) :: x(3, 3)
        integer :: i
        x = 0.0_dp
        do i = 1, 3
            x(i, i) = 1.0_dp
        end do
    end function eye3

    subroutine check(condition, message)
        logical, intent(in) :: condition !! Condition that must hold for the deterministic regression test.
        character(len=*), intent(in) :: message !! Failure label printed by ERROR STOP context.
        if (.not. condition) then
            print '(a)', trim(message)
            error stop 1
        end if
    end subroutine check
end program test_schur
