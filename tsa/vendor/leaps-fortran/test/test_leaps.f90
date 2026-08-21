program test_leaps

use leaps, only: dp, regsubsets_result, regsubsets_fit, get_model, &
    model_coefficients

implicit none

integer, parameter :: n = 12, p = 5
real(dp) :: x(n,p), y(n), t, z, noise(n), rss
real(dp), allocatable :: beta(:), vc(:,:)
integer, allocatable :: ids(:)
logical :: fi(p), fo(p)
type(regsubsets_result) :: fit
integer :: i, ier

noise = [0.10_dp, -0.08_dp, 0.04_dp, 0.03_dp, -0.06_dp, 0.09_dp, &
         -0.02_dp, 0.05_dp, -0.07_dp, 0.01_dp, 0.06_dp, -0.05_dp]

do i = 1, n
    t = real(i, dp)
    z = (t - 6.5_dp) / 3.0_dp
    x(i,1) = z
    x(i,2) = z**2
    x(i,3) = sin(t)
    x(i,4) = cos(0.7_dp * t)
    x(i,5) = real(mod(i,3) - 1, dp)
end do

y = 1.2_dp + 2.0_dp*x(:,1) - 1.4_dp*x(:,3) + 0.45_dp*x(:,5) + noise

call regsubsets_fit(x, y, fit, nvmax=4, nbest=3, method='exhaustive', ier=ier)
call check(ier == 0, 'exhaustive fit status')

call get_model(fit, 1, 1, ids, rss, ier)
call check(ier == 0, 'get 1-predictor model')
call check(all(ids == [0,1]), 'best 1-predictor model')
call check(close(rss, 12.150138283273455_dp, 2.0e-11_dp), 'best 1-predictor rss')

call get_model(fit, 2, 1, ids, rss, ier)
call check(all(ids == [0,1,3]), 'best 2-predictor model')
call check(close(rss, 1.5906048548576883_dp, 2.0e-11_dp), 'best 2-predictor rss')

call get_model(fit, 2, 2, ids, rss, ier)
call check(all(ids == [0,1,4]), 'second 2-predictor model')
call check(close(rss, 7.16062401240969_dp, 2.0e-11_dp), 'second 2-predictor rss')

call get_model(fit, 3, 1, ids, rss, ier)
call check(all(ids == [0,1,3,5]), 'best 3-predictor model')
call check(close(rss, 0.04172906573251695_dp, 2.0e-11_dp), 'best 3-predictor rss')

call model_coefficients(fit, 3, 1, beta, ids, ier, vc)
call check(ier == 0, 'coefficient recovery status')
call check(all(ids == [0,1,3,5]), 'coefficient ids')
call check(close(beta(1), 1.20824952_dp, 2.0e-8_dp), 'intercept coefficient')
call check(close(beta(2), 1.98813816_dp, 2.0e-8_dp), 'x1 coefficient')
call check(close(beta(3), -1.40802201_dp, 2.0e-8_dp), 'x3 coefficient')
call check(close(beta(4), 0.44322608_dp, 2.0e-8_dp), 'x5 coefficient')
call check(maxval(abs(vc - transpose(vc))) < 1.0e-13_dp, 'vcov symmetry')

call check_method('forward')
call check_method('backward')
call check_method('seqrep')

fi = .false.
fo = .false.
fi(3) = .true.
fo(5) = .true.
call regsubsets_fit(x, y, fit, nvmax=3, nbest=1, method='exhaustive', &
                    force_in=fi, force_out=fo, ier=ier)
call check(ier == 0, 'forced-variable fit status')
call get_model(fit, 2, 1, ids, rss, ier)
call check(all(ids == [0,1,3]), 'forced-in/out model')



call regsubsets_fit(x, y, fit, nvmax=2, nbest=1, method='exhaustive', &
                    intercept=.false., ier=ier)
call check(ier == 0, 'no-intercept fit status')
call get_model(fit, 2, 1, ids, rss, ier)
call check(ier == 0, 'no-intercept model retrieval')
call check(all(ids > 0), 'no-intercept model identifiers')

! Rank-deficient predictor: x2 is exactly twice x1.  The translation should
! detect and exclude the redundant column rather than failing the search.
x(:,2) = 2.0_dp * x(:,1)
fi = .false.
fo = .false.
call regsubsets_fit(x, y, fit, nvmax=3, nbest=1, method='exhaustive', ier=ier)
call check(ier == 0, 'rank-deficient fit status')
call check(count(fit%dependent) == 1, 'rank-deficient column detection')
call get_model(fit, 2, 1, ids, rss, ier)
call check(ier == 0, 'rank-deficient best model retrieval')
call check(.not. (any(ids == 1) .and. any(ids == 2)), &
           'dependent pair not simultaneously selected')

print '(a)', 'test_leaps: PASS'

contains

subroutine check_method(name)
    character(len=*), intent(in) :: name
    type(regsubsets_result) :: local_fit
    integer, allocatable :: local_ids(:)
    integer :: local_ier

    call regsubsets_fit(x, y, local_fit, nvmax=4, nbest=1, method=name, ier=local_ier)
    call check(local_ier == 0, trim(name)//' fit status')
    call get_model(local_fit, 1, 1, local_ids, ier=local_ier)
    call check(all(local_ids == [0,1]), trim(name)//' size 1')
    call get_model(local_fit, 2, 1, local_ids, ier=local_ier)
    call check(all(local_ids == [0,1,3]), trim(name)//' size 2')
    call get_model(local_fit, 3, 1, local_ids, ier=local_ier)
    call check(all(local_ids == [0,1,3,5]), trim(name)//' size 3')
end subroutine check_method


logical function close(a, b, tol)
    real(dp), intent(in) :: a, b, tol
    close = abs(a-b) <= tol * max(1.0_dp, abs(b))
end function close


subroutine check(ok, label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    if (.not. ok) then
        write(*,'(a)') 'FAIL: '//trim(label)
        error stop 1
    end if
end subroutine check

end program test_leaps
