program test_statuses
    use lpsolve
    implicit none
    type(lp_result) :: r
    real(dp) :: c(1), a(1,1), b(1)
    integer :: sense(1)
    real(dp), allocatable :: a0(:,:), b0(:)
    integer, allocatable :: s0(:)

    c = 1.0_dp
    a(1,1) = 1.0_dp
    b = -1.0_dp
    sense = LP_LE
    call solve_lp(LP_MAX, c, a, sense, b, r)
    if (r%status /= LP_INFEASIBLE) error stop 'infeasible status'

    allocate(a0(0,1), b0(0), s0(0))
    call solve_lp(LP_MAX, c, a0, s0, b0, r)
    if (r%status /= LP_UNBOUNDED) error stop 'unbounded status'
end program test_statuses
