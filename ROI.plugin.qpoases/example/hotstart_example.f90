program hotstart_example
    use qpoases
    implicit none

    type(qpoases_model) :: model
    real(dp) :: h(2,2), g(2), lb(2), ub(2)
    real(dp), allocatable :: x(:)
    integer :: status

    h = 0.0_dp
    h(1,1) = 1.0_dp
    h(2,2) = 1.0_dp
    lb = 0.0_dp
    ub = 10.0_dp

    g = [-1.0_dp,-2.0_dp]
    call init_qproblemb(model,h,g,lb,ub,200,status,hessian_type=hst_identity)
    call get_primal_solution(model,x)
    print '(a,2f12.6)', "initial x = ", x

    g = [-3.0_dp,-4.0_dp]
    call hotstart_qproblemb(model,g,lb,ub,200,status)
    call get_primal_solution(model,x)
    print '(a,2f12.6)', "hotstart x = ", x
end program hotstart_example
