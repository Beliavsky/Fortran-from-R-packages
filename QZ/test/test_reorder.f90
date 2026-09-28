program test_reorder
    use qz_api, only: dp, qz_dtrsen, qz_ztrsen, qz_dtgsen, qz_ztgsen, ordqz, &
        qz_reorder_real_result, qz_reorder_complex_result, &
        qz_greorder_real_result, qz_greorder_complex_result
    implicit none

    real(dp) :: t(3, 3), q(3, 3), b(3, 3)
    complex(dp) :: tc(3, 3), qc(3, 3), bc(3, 3)
    logical :: select(3)
    type(qz_reorder_real_result) :: rr, orr
    type(qz_reorder_complex_result) :: rc, orc
    type(qz_greorder_real_result) :: gr, ogr
    type(qz_greorder_complex_result) :: gc, ogc

    t = diag3(2.0_dp, -1.0_dp, 0.5_dp)
    q = diag3(1.0_dp, 1.0_dp, 1.0_dp)
    select = [.false., .true., .false.]
    call qz_dtrsen(t, q, select, rr)
    call check(rr%info == 0, 'qz_dtrsen info')
    call check(rr%m == 1, 'qz_dtrsen selected dimension')
    call check(abs(real(rr%w(1), dp) + 1.0_dp) < 1.0e-12_dp, 'qz_dtrsen ordering')

    tc = cmplx(t, 0.0_dp, kind=dp)
    qc = cmplx(q, 0.0_dp, kind=dp)
    tc(1, 1) = cmplx(2.0_dp, 0.5_dp, kind=dp)
    call qz_ztrsen(tc, qc, select, rc)
    call check(rc%info == 0, 'qz_ztrsen info')
    call check(rc%m == 1, 'qz_ztrsen selected dimension')
    call check(abs(real(rc%w(1), dp) + 1.0_dp) < 1.0e-12_dp, 'qz_ztrsen ordering')

    b = diag3(1.0_dp, 1.0_dp, 1.0_dp)
    call qz_dtgsen(t, b, q, q, select, gr)
    call check(gr%info == 0, 'qz_dtgsen info')
    call check(gr%m == 1, 'qz_dtgsen selected dimension')
    call check(abs(real(gr%lambda(1), dp) + 1.0_dp) < 1.0e-12_dp, 'qz_dtgsen ordering')

    bc = cmplx(b, 0.0_dp, kind=dp)
    call qz_ztgsen(tc, bc, qc, qc, select, gc)
    call check(gc%info == 0, 'qz_ztgsen info')
    call check(gc%m == 1, 'qz_ztgsen selected dimension')
    call check(abs(real(gc%lambda(1), dp) + 1.0_dp) < 1.0e-12_dp, 'qz_ztgsen ordering')

    call ordqz(t, orr, keyword='lhp')
    call check(orr%info == 0 .and. orr%m == 1, 'ordqz real standard status')
    call check(real(orr%w(1), dp) < 0.0_dp, 'ordqz lhp standard')

    call ordqz(tc, orc, keyword='udi')
    call check(orc%info == 0 .and. orc%m == 1, 'ordqz complex standard status')
    call check(abs(orc%w(1)) < 1.0_dp, 'ordqz udi standard')

    call ordqz(t, b, ogr, keyword='rhp')
    call check(ogr%info == 0 .and. ogr%m == 2, 'ordqz real generalized status')
    call check(all(real(ogr%lambda(1:2), dp) >= 0.0_dp), 'ordqz rhp generalized')

    call ordqz(tc, bc, ogc, keyword='lhp')
    call check(ogc%info == 0 .and. ogc%m == 1, 'ordqz complex generalized status')
    call check(real(ogc%lambda(1), dp) < 0.0_dp, 'ordqz lhp complex generalized')

    t = diag3(-1.0_dp, 2.0_dp, 0.5_dp)
    b = diag3(1.0_dp, 0.0_dp, 1.0_dp)
    call ordqz(t, b, ogr, keyword='rhp.fo')
    call check(ogr%info == 0 .and. ogr%m == 1, 'ordqz finite-only generalized status')
    call check(abs(real(ogr%lambda(1), dp) - 0.5_dp) < 1.0e-12_dp, 'ordqz excludes infinite generalized root')

    print '(a)', 'test_reorder: PASS'
contains
    pure function diag3(a, b, c) result(x)
        real(dp), intent(in) :: a !! First diagonal element.
        real(dp), intent(in) :: b !! Second diagonal element.
        real(dp), intent(in) :: c !! Third diagonal element.
        real(dp) :: x(3, 3)
        x = 0.0_dp
        x(1, 1) = a
        x(2, 2) = b
        x(3, 3) = c
    end function diag3

    subroutine check(condition, message)
        logical, intent(in) :: condition !! Condition that must hold for the deterministic regression test.
        character(len=*), intent(in) :: message !! Failure label printed by ERROR STOP context.
        if (.not. condition) then
            print '(a)', trim(message)
            error stop 1
        end if
    end subroutine check
end program test_reorder
