program test_differsums
    use multiplicative_multinomial
    implicit none
    type(paras_type) :: par
    type(mm_fit_type) :: fit
    real(dp) :: th(2,2), prob(2), ll, manual
    integer :: y(4,2), i

    prob = [0.35_dp, 0.65_dp]
    th = 1.0_dp
    th(1,2) = 1.25_dp
    par = paras_from_p_theta(prob, th)
    y(1,:) = [1, 0]
    y(2,:) = [1, 1]
    y(3,:) = [0, 3]
    y(4,:) = [2, 1]
    ll = mm_loglik(y, par)
    manual = 0.0_dp
    do i = 1, size(y,1)
        manual = manual + log(dmm(y(i,:), par))
    end do
    call check_close(ll, manual, 5.0e-13_dp, "different-sum likelihood")

    call optimizer_differsums(y, fit, start=par, max_iter=200, tol=1.0e-7_dp)
    if (fit%loglik + 1.0e-8_dp < ll) error stop "different-sum optimizer decreased likelihood"

    print '(a)', 'test_differsums: PASS'

contains

    subroutine check_close(x, ref, tol, label)
        real(dp), intent(in) :: x, ref, tol
        character(len=*), intent(in) :: label
        if (abs(x-ref) > tol * max(1.0_dp, abs(ref))) then
            write(*,'(a,2es24.14)') trim(label)//' failed: ', x, ref
            error stop 1
        end if
    end subroutine check_close

end program test_differsums
