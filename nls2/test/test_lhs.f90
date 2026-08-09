program test_lhs
    use nls2
    implicit none
    integer, parameter :: n = 25, p = 3
    integer :: i, j
    real(dp) :: lower(p), upper(p), points(n,p), s(n)

    lower = [-2.0_dp, 5.0_dp, 10.0_dp]
    upper = [ 3.0_dp, 9.0_dp, 20.0_dp]
    call seed_rng(9876)
    call latin_hypercube(lower, upper, points)
    call assert_true(all(points >= spread(lower,1,n)), 'lhs lower bounds')
    call assert_true(all(points <= spread(upper,1,n)), 'lhs upper bounds')
    do j = 1, p
        s = (points(:,j) - lower(j)) / (upper(j)-lower(j))
        call sort_real(s)
        do i = 1, n
            call assert_true(s(i) >= real(i-1,dp)/real(n,dp), 'lhs stratum low')
            call assert_true(s(i) <= real(i,dp)/real(n,dp), 'lhs stratum high')
        end do
    end do
contains
    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j
        real(dp) :: t
        do i = 2, size(x)
            t = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= t) exit
                x(j+1) = x(j)
                j = j - 1
            end do
            x(j+1) = t
        end do
    end subroutine sort_real
    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) error stop 'test_lhs: '//label
    end subroutine assert_true
end program test_lhs
