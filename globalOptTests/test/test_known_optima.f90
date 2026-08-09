program test_known_optima
    use global_opt_tests, only : dp, ackleys, beckerlago, goldprice, griewank, rastrigin, rosenbrock, salomon, schaffer1
    implicit none
    real(dp) :: x10(10), x2(2), x5(5)

    x10 = 0.0_dp
    if (abs(ackleys(x10)) > 1.0e-12_dp) error stop 'Ackley optimum failed'
    if (abs(griewank(x10)) > 1.0e-12_dp) error stop 'Griewank optimum failed'
    if (abs(rastrigin(x10)) > 1.0e-12_dp) error stop 'Rastrigin optimum failed'

    x10 = 1.0_dp
    if (abs(rosenbrock(x10)) > 1.0e-12_dp) error stop 'Rosenbrock optimum failed'

    x2 = [5.0_dp, -5.0_dp]
    if (abs(beckerlago(x2)) > 1.0e-12_dp) error stop 'Becker-Lago optimum failed'

    x2 = [0.0_dp, -1.0_dp]
    if (abs(goldprice(x2)-3.0_dp) > 1.0e-12_dp) error stop 'Goldstein-Price optimum failed'

    x2 = 0.0_dp
    if (abs(schaffer1(x2)) > 1.0e-12_dp) error stop 'Schaffer1 optimum failed'

    x5 = 0.0_dp
    if (abs(salomon(x5)) > 1.0e-12_dp) error stop 'Salomon optimum failed'

    print '(a)', 'test_known_optima: PASS'
end program test_known_optima
