program constrained_least_squares
    use limsolve
    implicit none
    type(solve_result) :: fit
    real(dp) :: a(2,2),b(2),e(1,2),f(1),g(2,2),h(2)
    a=0.0_dp; a(1,1)=1.0_dp; a(2,2)=1.0_dp
    b=[0.0_dp,0.0_dp]
    e(1,:)=[1.0_dp,1.0_dp]; f=[1.0_dp]
    g=0.0_dp; g(1,1)=1.0_dp; g(2,2)=1.0_dp; h=0.0_dp
    call lsei(a,b,e,f,g,h,fit)
    print '(a,2f14.8)', 'x = ',fit%x
    print '(a,es14.6)', '||A*x-b||^2 = ',fit%solution_norm
end program constrained_least_squares
