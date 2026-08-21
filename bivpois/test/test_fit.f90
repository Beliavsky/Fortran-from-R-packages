program test_fit
    use bivpois, only : dp, bp_mle2_result, bp_mle_result, bp_profile_result, &
                        bp_mle2, bp_mle, lambda3_profile, profile_loglik
    implicit none
    integer, parameter :: n = 80
    integer, parameter :: x1(n) = [ &
        5,11,4,6,3,8,7,4,3,5,0,13,3,4,6,3,3,7,7,4, &
        3,6,4,3,3,6,3,5,3,5,5,5,2,7,7,1,2,3,4,6, &
        5,7,3,7,2,3,4,5,4,7,8,3,3,3,9,5,5,10,8,2, &
        5,4,3,3,7,4,2,7,6,3,4,5,5,4,6,4,7,2,10,6 ]
    integer, parameter :: x2(n) = [ &
        6,10,7,7,3,7,6,5,4,8,4,12,4,10,5,6,3,6,9,7, &
        3,3,7,6,9,5,10,7,10,7,7,6,4,11,6,6,3,10,11,7, &
        5,5,9,3,6,9,5,8,5,9,10,7,8,2,9,6,8,13,12,6, &
        2,5,5,5,10,4,6,8,5,8,6,12,6,10,8,7,6,3,9,7 ]
    type(bp_mle2_result) :: fit2
    type(bp_mle_result) :: fit
    type(bp_profile_result) :: prof
    real(dp) :: ll0

    fit2 = bp_mle2(x1, x2, tol=1.0e-10_dp)
    call assert_close(fit2%lambda(3), 2.48070493735253_dp, 3.0e-6_dp, "lambda3")
    call assert_close(fit2%lambda(1), 2.38179506264747_dp, 3.0e-6_dp, "lambda1")
    call assert_close(fit2%lambda(2), 4.31929506264747_dp, 3.0e-6_dp, "lambda2")
    call assert_close(fit2%loglik, -355.3619376814744_dp, 3.0e-9_dp, "loglik")

    ll0 = profile_loglik(x1, x2, 0.0_dp)
    call assert_close(ll0, -363.6849481214147_dp, 3.0e-9_dp, "independent loglik")

    fit = bp_mle(x1, x2, tol=1.0e-10_dp)
    call assert_close(fit%lambda(3), fit2%lambda(3), 1.0e-10_dp, "full/fast lambda3")
    if (.not. fit%converged) error stop "full fit convergence"
    if (fit%rho <= 0.0_dp .or. fit%rho >= 1.0_dp) error stop "rho bounds"
    if (fit%var_observed <= 0.0_dp) error stop "observed variance"
    if (fit%var_asymptotic <= 0.0_dp) error stop "asymptotic variance"
    call assert_close(fit%var_observed, 0.2396492644_dp, 5.0e-4_dp, "observed variance reference")
    call assert_close(fit%var_asymptotic, 0.2535201175_dp, 1.0e-7_dp, "asymptotic variance reference")
    if (any(fit%pvalue < 0.0_dp) .or. any(fit%pvalue > 1.0_dp)) error stop "pvalue bounds"

    prof = lambda3_profile(x1, x2)
    call assert_close(prof%mle, 2.48_dp, 1.0e-14_dp, "profile grid MLE")
    call assert_close(prof%ci(1), 1.39_dp, 1.0e-14_dp, "profile lower CI")
    call assert_close(prof%ci(2), 3.30_dp, 1.0e-14_dp, "profile upper CI")

    print '(a)', 'test_fit: PASS'

contains

    subroutine assert_close(a, b, tol, label)
        real(dp), intent(in) :: a, b, tol
        character(*), intent(in) :: label
        if (abs(a - b) > tol) then
            print '(a,2es24.15)', trim(label)//": ", a, b
            error stop 1
        end if
    end subroutine assert_close

end program test_fit
