program test_optimizer
    use caramel, only: dp, caramel_options, caramel_result, caramel_optimize, seed_random, pareto
    implicit none
    type(caramel_options) :: opt
    type(caramel_result) :: res
    real(dp) :: bounds(1,2), prec(2)
    logical :: minmax(2)
    integer, allocatable :: front(:)

    call seed_random(424242)
    bounds(1,:) = [-5.0_dp, 10.0_dp]
    prec = [1.0e-3_dp,1.0e-3_dp]
    minmax = [.false.,.false.]
    opt%popsize = 24
    opt%archsize = 30
    opt%maxrun = 120
    opt%repart_gene = [4,4,4,4]
    opt%gpp = 3
    opt%sensitivity = .true.

    call caramel_optimize(2, 1, minmax, bounds, prec, schaffer, res, opt)
    call check(res%success, 'optimizer success')
    call check(size(res%parameters,1) > 0, 'nonempty archive')
    call check(all(res%parameters(:,1) >= bounds(1,1) .and. res%parameters(:,1) <= bounds(1,2)), 'bounds')
    call check(size(res%objectives,2) == 2, 'objective columns')
    call check(size(res%save_crit,1) > 0 .and. size(res%save_crit,2) == 3, 'history')
    call check(size(res%derivatives,1) == size(res%parameters,1), 'sensitivity rows')

    allocate(front(size(res%objectives,1)))
    call pareto(-res%objectives, front)
    call check(all(front == 1), 'archive pareto')
    print '(a,i0,a,i0,a)', 'test_optimizer: PASS (front=', size(res%parameters,1), ', calls=', res%nrun, ')'
contains
    subroutine schaffer(x, values)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:)
        if (x(1) <= 1.0_dp) then
            values(1) = -x(1)
        else if (x(1) <= 3.0_dp) then
            values(1) = x(1) - 2.0_dp
        else if (x(1) <= 4.0_dp) then
            values(1) = 4.0_dp - x(1)
        else
            values(1) = x(1) - 4.0_dp
        end if
        values(2) = (x(1)-5.0_dp)**2
    end subroutine schaffer

    subroutine check(condition, name)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: name
        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(name)
            if (allocated(res%message)) write(*,'(a)') trim(res%message)
            error stop 1
        end if
    end subroutine check
end program test_optimizer
