module quadratic_example_objective
    use rceim, only : dp
    implicit none
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = (x(1)+1.0_dp)**2 + (x(2)+2.0_dp)**2
    end function objective
end module quadratic_example_objective

program quadratic_example
    use rceim, only : dp, rceim_options, rceim_result, ceim_optimize
    use quadratic_example_objective, only : objective
    implicit none
    type(rceim_options) :: opt
    type(rceim_result) :: res
    real(dp) :: lower(2), upper(2)

    lower = [-10.0_dp,-10.0_dp]
    upper = [ 10.0_dp, 10.0_dp]
    opt%n_total = 500
    opt%n_elite = 125
    opt%n_super = 1
    opt%epsilon = 0.01_dp
    opt%max_iter = 60
    opt%seed = 2026
    call ceim_optimize(objective, lower, upper, res, opt)
    write(*,'(a,2f12.6)') 'x = ', res%x
    write(*,'(a,es14.6)') 'f = ', res%value
    write(*,'(a,a)') 'stop = ', res%criterion
end program quadratic_example
