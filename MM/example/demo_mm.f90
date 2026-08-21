program demo_mm
    use multiplicative_multinomial
    implicit none
    type(paras_type) :: par, fitted
    type(glm_fit_type) :: gl
    real(dp) :: th(2,2), prob(2)
    integer :: y(2), obs(3,2), n(3)

    par = paras(2)
    th = 1.0_dp
    th(1,2) = 2.0_dp
    call set_theta(par, th)
    y = [1, 1]
    print '(a,f12.8)', 'P([1,1]) with theta12=2: ', dmm(y, par)

    obs(1,:) = [2, 0]
    obs(2,:) = [1, 1]
    obs(3,:) = [0, 2]
    n = [16, 48, 36]
    call lindsey_fit(obs, fitted, gl, real(n,dp))
    prob = p(fitted)
    th = theta(fitted)
    print '(a,2f12.8)', 'Lindsey p:     ', prob
    print '(a,f12.8)', 'Lindsey theta: ', th(1,2)
end program demo_mm
