program schaffer_example
    use caramel, only: dp, caramel_options, caramel_result, caramel_optimize, seed_random
    implicit none
    type(caramel_options) :: options
    type(caramel_result) :: result
    real(dp) :: bounds(1,2), prec(2)
    logical :: minimize(2)

    call seed_random(20260812)
    bounds(1,:) = [-5.0_dp, 10.0_dp]
    prec = 1.0e-3_dp
    minimize = .false.

    options%popsize = 40
    options%archsize = 50
    options%maxrun = 300
    options%repart_gene = [5,5,5,5]

    call caramel_optimize(2, 1, minimize, bounds, prec, schaffer, result, options)
    if (.not. result%success) then
        write(*,'(a)') 'caRamel failed: '//result%message
        error stop 1
    end if

    write(*,'(a,i0)') 'Pareto points: ', size(result%parameters,1)
    write(*,'(a,i0)') 'Objective calls: ', result%nrun
    write(*,'(a)') 'First few points: x, f1, f2'
    call print_points(min(8,size(result%parameters,1)))
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

    subroutine print_points(n)
        integer, intent(in) :: n
        integer :: i
        do i = 1, n
            write(*,'(3(es14.6,1x))') result%parameters(i,1), result%objectives(i,1), result%objectives(i,2)
        end do
    end subroutine print_points
end program schaffer_example
