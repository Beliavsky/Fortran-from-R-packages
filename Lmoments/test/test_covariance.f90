program test_covariance
    use lmoments, only: dp, lmom_cov, shifted_legendre, covnormpoly4
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp), parameter :: x(10) = [3.2_dp, -1.1_dp, 0.4_dp, 5.7_dp, 2.2_dp, &
        1.8_dp, 9.1_dp, -3.0_dp, 4.4_dp, 0.0_dp]
    real(dp), parameter :: ref(4,4) = reshape([ &
        1.242455555555558_dp, 0.224372222222220_dp, 0.119053174603176_dp, 0.049749206349215_dp, &
        0.224372222222220_dp, 0.259015520282186_dp, 0.072983068783077_dp, 0.072070634920573_dp, &
        0.119053174603176_dp, 0.072983068783077_dp, 0.109552380952347_dp, 0.036942063492170_dp, &
        0.049749206349215_dp, 0.072070634920574_dp, 0.036942063492172_dp, 0.072970294784531_dp ], [4,4])
    real(dp) :: cov(4,4), c(6,6), pcov(4,4)
    integer :: info

    call lmom_cov(x, cov, info)
    call assert_true(info == 0, 'cov info')
    call assert_mat(cov, ref, 3.0e-12_dp, 'covariance reference')

    call shifted_legendre(6, c, info)
    call assert_true(info == 0, 'legendre info')
    call assert_close(c(6,1), -1.0_dp, 1.0e-14_dp, 'P5 c0')
    call assert_close(c(6,2), 30.0_dp, 1.0e-14_dp, 'P5 c1')
    call assert_close(c(6,6), 252.0_dp, 1.0e-14_dp, 'P5 c5')

    pcov = covnormpoly4(x)
    call assert_true(all(ieee_is_finite(pcov)), 'parameter covariance finite')
    call assert_true(maxval(abs(pcov - transpose(pcov))) < 1.0e-9_dp, 'parameter covariance symmetric')

    print '(a)', 'test_covariance: PASS'
contains
    subroutine assert_mat(a, b, tol, label)
        real(dp), intent(in) :: a(:,:), b(:,:), tol
        character(*), intent(in) :: label
        if (maxval(abs(a - b)) > tol) then
            print *, label, maxval(abs(a-b))
            error stop 1
        end if
    end subroutine
    subroutine assert_close(a, b, tol, label)
        real(dp), intent(in) :: a, b, tol
        character(*), intent(in) :: label
        if (abs(a - b) > tol) then
            print *, label, a, b
            error stop 1
        end if
    end subroutine
    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, label
            error stop 1
        end if
    end subroutine
end program test_covariance
