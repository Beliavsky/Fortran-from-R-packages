program demo_leaps

use leaps, only: dp, regsubsets_result, regsubsets_fit, get_model, &
    model_coefficients

implicit none

integer, parameter :: n = 10, p = 4
real(dp) :: x(n,p), y(n), t, rss
real(dp), allocatable :: beta(:)
integer, allocatable :: ids(:)
type(regsubsets_result) :: fit
integer :: i, k, ier

do i = 1, n
    t = real(i, dp)
    x(i,1) = t
    x(i,2) = t*t
    x(i,3) = sin(t)
    x(i,4) = real(mod(i,2), dp)
end do

y = 2.0_dp + 0.5_dp*x(:,1) - 1.25_dp*x(:,3) + &
    [0.02_dp, -0.03_dp, 0.01_dp, 0.04_dp, -0.02_dp, &
     0.00_dp, 0.03_dp, -0.01_dp, 0.02_dp, -0.04_dp]

call regsubsets_fit(x, y, fit, nvmax=3, nbest=2, &
                    method='exhaustive', ier=ier)
if (ier /= 0) error stop 'regsubsets_fit failed'

write(*,'(a)') 'Best subset of each size:'
do k = 1, fit%nvmax
    call get_model(fit, k, 1, ids, rss, ier)
    if (ier /= 0) cycle
    write(*,'(a,i0,a,*(i0,1x))') '  predictors=', k, ' ids=', ids
    write(*,'(a,es16.8)') '    rss = ', rss
end do

call model_coefficients(fit, 2, 1, beta, ids, ier)
if (ier == 0) then
    write(*,'(/,a)') 'Coefficients for the best two-predictor model:'
    do i = 1, size(beta)
        if (ids(i) == 0) then
            write(*,'(a,es16.8)') '  intercept: ', beta(i)
        else
            write(*,'(a,i0,a,es16.8)') '  x', ids(i), ': ', beta(i)
        end if
    end do
end if

end program demo_leaps
