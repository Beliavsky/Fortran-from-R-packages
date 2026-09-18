program otrimle_example
    use otrimle_mod
    implicit none

    real(dp) :: x(10, 2)
    integer :: initial(10)
    type(otrimle_fit) :: fit

    x(:, 1) = [-2.0_dp, -1.8_dp, -2.2_dp, -1.9_dp, -2.1_dp, 2.0_dp, 2.2_dp, 1.8_dp, 2.1_dp, 8.0_dp]
    x(:, 2) = [-2.1_dp, -1.9_dp, -2.0_dp, -2.2_dp, -1.8_dp, 2.1_dp, 1.9_dp, 2.0_dp, 2.2_dp, -7.0_dp]
    initial = [1, 1, 1, 1, 1, 2, 2, 2, 2, 0]

    call rimle(x, 2, fit, initial=initial, logicd=-8.0_dp, npr_max=0.30_dp, erc=20.0_dp)
    write (*, '(a,i0)') 'status code: ', fit%code
    write (*, '(a,f10.6)') 'improper log likelihood: ', fit%iloglik
    write (*, '(a,3(f10.6,1x))') 'component proportions: ', fit%exproportion
    write (*, '(a,2(f10.4,1x))') 'cluster 1 mean: ', fit%mean(:, 1)
    write (*, '(a,2(f10.4,1x))') 'cluster 2 mean: ', fit%mean(:, 2)
end program otrimle_example
