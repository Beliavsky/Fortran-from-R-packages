program test_highlevel
    use qz_api, only: dp, qz_geigen, geigen, qz_decompose, fda_geigen, h_transpose, &
        qz_eigen_real_result, qz_reorder_real_result, qz_fda_geigen_result
    implicit none

    real(dp) :: a(2, 2), b(2, 2), c(2, 2), hr(2, 2)
    complex(dp) :: x(2, 2), hx(2, 2)
    type(qz_eigen_real_result) :: e1, e2
    type(qz_reorder_real_result) :: schur
    type(qz_fda_geigen_result) :: fda

    a = 0.0_dp
    a(1, 1) = 2.0_dp
    a(2, 2) = 1.0_dp
    call qz_geigen(a, e1)
    call geigen(a, e2)
    call check(e1%info == 0 .and. e2%info == 0, 'high-level geigen status')
    call check(maxval(abs(e1%w - e2%w)) < 1.0e-12_dp, 'geigen alias')

    call qz_decompose(a, schur)
    call check(schur%info == 0, 'qz_decompose status')
    call check(maxval(abs(a - matmul(schur%q, matmul(schur%t, transpose(schur%q))))) < 1.0e-12_dp, &
        'qz_decompose reconstruction')

    b = 0.0_dp
    c = 0.0_dp
    b(1, 1) = 1.0_dp
    b(2, 2) = 1.0_dp
    c = b
    call fda_geigen(a, b, c, fda)
    call check(fda%info == 0, 'fda_geigen status')
    call check(size(fda%values) == 2, 'fda_geigen value count')
    call check(abs(fda%values(1) - 2.0_dp) < 1.0e-11_dp, 'fda_geigen first value')
    call check(abs(fda%values(2) - 1.0_dp) < 1.0e-11_dp, 'fda_geigen second value')
    call check(abs(dot_product(fda%lmat(:, 1), matmul(b, fda%lmat(:, 1))) - 1.0_dp) < 1.0e-11_dp, &
        'fda_geigen L normalization')
    call check(abs(dot_product(fda%mmat(:, 1), matmul(c, fda%mmat(:, 1))) - 1.0_dp) < 1.0e-11_dp, &
        'fda_geigen M normalization')

    x(1, 1) = cmplx(1.0_dp, 2.0_dp, kind=dp)
    x(2, 1) = cmplx(3.0_dp, -1.0_dp, kind=dp)
    x(1, 2) = cmplx(-2.0_dp, 4.0_dp, kind=dp)
    x(2, 2) = cmplx(0.5_dp, -0.25_dp, kind=dp)
    hx = h_transpose(x)
    call check(maxval(abs(hx - conjg(transpose(x)))) < 1.0e-15_dp, 'Hermitian transpose complex')
    hr = h_transpose(a)
    call check(maxval(abs(hr - transpose(a))) < 1.0e-15_dp, 'Hermitian transpose real')

    print '(a)', 'test_highlevel: PASS'
contains
    subroutine check(condition, message)
        logical, intent(in) :: condition !! Condition that must hold for the deterministic regression test.
        character(len=*), intent(in) :: message !! Failure label printed by ERROR STOP context.
        if (.not. condition) then
            print '(a)', trim(message)
            error stop 1
        end if
    end subroutine check
end program test_highlevel
