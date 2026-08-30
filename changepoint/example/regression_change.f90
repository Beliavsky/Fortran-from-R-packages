! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
program regression_change
use r_kinds, only : dp
use changepoint, only : changepoint_result, cp_regression_pelt
implicit none
integer, parameter :: n = 12
real(dp) :: y(n), x(n, 2)
type(changepoint_result) :: result
integer :: i

do i = 1, n
    x(i, 1) = 1.0_dp
    x(i, 2) = real(i, dp)
    if (i <= 6) then
        y(i) = 1.0_dp + 0.5_dp * x(i, 2)
    else
        y(i) = -2.0_dp + 1.5_dp * x(i, 2)
    end if
end do
call cp_regression_pelt(y, x, 1.0_dp, 3, -1.0_dp, result)
print '(a,i0)', 'number of changepoints: ', result%ncpts
if (result%ncpts > 0) print '(a,*(1x,i0))', 'locations:', result%cpts
end program regression_change
