program test_sparse_q8
    use lpsolve
    implicit none
    type(q8_triplets) :: q8
    type(sparse_constraints) :: sp
    type(lp_result) :: r
    type(lp_control) :: ctl
    real(dp) :: c(64), rhs(42)
    integer :: sense(42), bins(64), i

    call make_q8(q8)
    if (size(q8%value) /= 252) error stop 'q8 nonzero count'
    if (maxval(q8%constraint) /= 42) error stop 'q8 constraint count'

    sp%nrow = 42
    sp%ncol = 64
    allocate(sp%row(size(q8%value)), sp%col(size(q8%value)), sp%val(size(q8%value)))
    sp%row = q8%constraint
    sp%col = q8%variable
    sp%val = q8%value
    c = 1.0_dp
    rhs = 1.0_dp
    sense(1:16) = LP_EQ
    sense(17:42) = LP_LE
    bins = [(i,i=1,64)]
    ctl%max_nodes = 50000
    call solve_lp_sparse(LP_MAX, c, sp, sense, rhs, r, ctl, binary_variables=bins)
    if (r%status /= LP_OPTIMAL) error stop '8 queens status'
    if (abs(r%objective - 8.0_dp) > 1.0e-8_dp) error stop '8 queens objective'
    if (abs(sum(r%solution) - 8.0_dp) > 1.0e-8_dp) error stop '8 queens count'
end program test_sparse_q8
