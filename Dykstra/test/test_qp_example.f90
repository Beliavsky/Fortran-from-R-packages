program test_qp_example
    use dykstra, only : dp, dykstra_result, dykstra_solve
    implicit none

    real(dp) :: dmat(3,3), dvec(3), amat(3,4), bvec(4)
    type(dykstra_result) :: result

    dmat = 0.0_dp
    dmat(1,1) = 1.0_dp
    dmat(2,2) = 1.0_dp
    dmat(3,3) = 1.0_dp
    dvec = [1.0_dp, 1.5_dp, 1.0_dp]

    amat = 0.0_dp
    amat(:,1) = 1.0_dp
    amat(1,2) = 1.0_dp
    amat(2,3) = 1.0_dp
    amat(3,4) = 1.0_dp
    bvec = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

    call dykstra_solve(dmat, dvec, amat, result, bvec=bvec, meq=1)

    call assert_true(result%converged, "QP example did not converge")
    call assert_close(result%solution, [1.0_dp/6.0_dp, 2.0_dp/3.0_dp, 1.0_dp/6.0_dp], &
        5.0e-13_dp, "QP example solution")
    call assert_scalar(result%value, -13.0_dp/12.0_dp, 5.0e-13_dp, "QP example objective")

contains

    subroutine assert_close(x, y, tol, label)
        real(dp), intent(in) :: x(:), y(:), tol
        character(len=*), intent(in) :: label
        if (size(x) /= size(y)) error stop label // ": size mismatch"
        if (maxval(abs(x - y)) > tol) error stop label // ": numerical mismatch"
    end subroutine assert_close

    subroutine assert_scalar(x, y, tol, label)
        real(dp), intent(in) :: x, y, tol
        character(len=*), intent(in) :: label
        if (abs(x - y) > tol) error stop label // ": numerical mismatch"
    end subroutine assert_scalar

    subroutine assert_true(value, label)
        logical, intent(in) :: value
        character(len=*), intent(in) :: label
        if (.not. value) error stop label
    end subroutine assert_true

end program test_qp_example
