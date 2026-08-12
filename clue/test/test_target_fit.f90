program test_target_fit
    use clue_kinds, only: dp
    use clue_target_fit, only: fit_ultrametric_target_l2, fit_ultrametric_target_l1
    use clue_trees, only: is_ultrametric
    implicit none
    real(dp) :: d(4,4)
    real(dp), allocatable :: u2(:,:), u1(:,:)
    integer :: merge(3,2)

    d = reshape([ &
        0.0_dp, 1.0_dp, 5.0_dp, 5.0_dp, &
        1.0_dp, 0.0_dp, 6.0_dp, 6.0_dp, &
        5.0_dp, 6.0_dp, 0.0_dp, 2.0_dp, &
        5.0_dp, 6.0_dp, 2.0_dp, 0.0_dp ], [4,4])
    merge(1,:) = [-1,-2]
    merge(2,:) = [-3,-4]
    merge(3,:) = [ 1, 2]

    u2 = fit_ultrametric_target_l2(d, merge)
    call check(abs(u2(1,2)-1.0_dp) < 1e-12_dp, 'L2 first merge')
    call check(abs(u2(3,4)-2.0_dp) < 1e-12_dp, 'L2 second merge')
    call check(abs(u2(1,3)-5.5_dp) < 1e-12_dp, 'L2 root merge')
    call check(is_ultrametric(u2, 1e-12_dp), 'L2 target is ultrametric')

    u1 = fit_ultrametric_target_l1(d, merge)
    call check(abs(u1(1,2)-1.0_dp) < 1e-12_dp, 'L1 first merge')
    call check(abs(u1(3,4)-2.0_dp) < 1e-12_dp, 'L1 second merge')
    call check(abs(u1(1,3)-5.0_dp) < 1e-12_dp, 'L1 root merge')
    call check(is_ultrametric(u1, 1e-12_dp), 'L1 target is ultrametric')

    print '(a)', 'test_target_fit: PASS'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print '(a)', 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine
end program test_target_fit
