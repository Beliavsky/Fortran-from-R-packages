program test_mapping_storage
    use deoptim, only : dp, de_control, de_result, deoptim_solve, de_success
    implicit none
    type(de_control) :: ctrl
    type(de_result) :: res
    real(dp) :: lo(3), hi(3)
    integer :: i, j

    lo = -4.0_dp
    hi = 4.0_dp
    ctrl = de_control()
    ctrl%np = 40
    ctrl%itermax = 80
    ctrl%strategy = 2
    ctrl%seed = 9217
    ctrl%storepopfrom = 1
    ctrl%storepopfreq = 5
    ctrl%reltol = 0.0_dp
    ctrl%steptol = ctrl%itermax

    call deoptim_solve(sphere, lo, hi, res, ctrl, map=half_grid)
    if (res%status /= de_success) error stop "mapped solve failed"
    if (res%nstore /= 16) error stop "stored population count mismatch"
    if (size(res%storepop,1) /= ctrl%np .or. size(res%storepop,2) /= 3) &
        error stop "stored population shape mismatch"

    do i = 1, size(res%pop,1)
        do j = 1, size(res%pop,2)
            if (abs(2.0_dp * res%pop(i,j) - real(nint(2.0_dp * res%pop(i,j)),dp)) > 1.0e-12_dp) &
                error stop "fnMap-style mapping was not retained"
        end do
    end do
    if (res%bestval > 1.0e-12_dp) error stop "mapped sphere should find zero"
    print *, "test_mapping_storage: PASS"
contains
    function sphere(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = sum(x*x)
    end function sphere

    subroutine half_grid(x)
        real(dp), intent(inout) :: x(:)
        x = 0.5_dp * real(nint(2.0_dp*x), dp)
    end subroutine half_grid
end program test_mapping_storage
