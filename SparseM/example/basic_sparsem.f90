program basic_sparsem
    use sparsem, only: dp, csr_matrix, slm_result
    use sparsem, only: csr_from_dense, csr_matvec, slm_fit_csr
    implicit none

    real(dp) :: x_dense(4, 2)
    real(dp) :: y(4)
    real(dp), allocatable :: fitted(:)
    type(csr_matrix) :: x
    type(slm_result) :: fit
    integer :: info

    x_dense(:, 1) = 1.0_dp
    x_dense(:, 2) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    y = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp]

    x = csr_from_dense(x_dense)
    fit = slm_fit_csr(x, y, info)
    if (info /= 0) error stop 'least-squares fit failed'

    fitted = csr_matvec(x, fit%coefficients, info)
    print '(a,2(1x,f10.6))', 'coefficients:', fit%coefficients
    print '(a,4(1x,f10.6))', 'fitted:', fitted
end program basic_sparsem
