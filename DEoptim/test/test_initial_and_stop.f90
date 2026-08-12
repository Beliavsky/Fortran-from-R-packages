program test_initial_and_stop
    use deoptim, only : dp, i8, de_control, de_result, deoptim_solve, de_success, de_unsupported
    implicit none
    type(de_control) :: ctrl
    type(de_result) :: res
    real(dp) :: lo(2), hi(2), init(10,2)
    integer :: i

    lo = -2.0_dp
    hi = 2.0_dp
    do i = 1, 10
        init(i,1) = -1.0_dp + 0.2_dp * real(i-1,dp)
        init(i,2) =  1.0_dp - 0.2_dp * real(i-1,dp)
    end do
    init(1,:) = 0.0_dp

    ctrl = de_control()
    ctrl%itermax = 100
    ctrl%vtr = 0.0_dp
    ctrl%seed = 3_i8
    call deoptim_solve(sphere, lo, hi, res, ctrl, initialpop=init)
    if (res%status /= de_success) error stop "initial population solve failed"
    if (res%iter /= 0) error stop "VTR should stop after initial population"
    if (res%nfeval /= 10_i8) error stop "initial evaluation count mismatch"
    if (abs(res%bestval) > 0.0_dp) error stop "initial optimum not found"

    ctrl%bs = .true.
    call deoptim_solve(sphere, lo, hi, res, ctrl, initialpop=init)
    if (res%status /= de_unsupported) error stop "bs compatibility behavior mismatch"

    print *, "test_initial_and_stop: PASS"
contains
    function sphere(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = sum(x*x)
    end function sphere
end program test_initial_and_stop
