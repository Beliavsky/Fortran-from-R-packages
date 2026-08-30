! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
program mean_changes
use r_kinds, only : dp
use changepoint, only : changepoint_result, cp_pelt, cp_cost_mean_normal
implicit none
real(dp) :: x(12)
type(changepoint_result) :: result
integer :: i

x = [0.0_dp, 0.1_dp, -0.1_dp, 0.05_dp, 5.0_dp, 5.1_dp, 4.9_dp, 5.05_dp, &
    -2.0_dp, -2.1_dp, -1.9_dp, -2.05_dp]
call cp_pelt(x, cp_cost_mean_normal, 1.0_dp, 2, result)
print '(a,i0)', 'number of changepoints: ', result%ncpts
if (result%ncpts > 0) then
    write(*, '(a)', advance='no') 'locations:'
    do i = 1, result%ncpts
        write(*, '(1x,i0)', advance='no') result%cpts(i)
    end do
    print *
end if
end program mean_changes
