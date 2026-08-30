! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_fit
use r_kinds, only : dp
use changepoint_regression, only : cp_regression_segment_fit
implicit none
private
public :: cp_segment_means
public :: cp_segment_variances_mle
public :: cp_segment_scales
public :: cp_segment_trend_fits
public :: cp_segment_regression_fits

contains

subroutine cp_segment_means(data, cpts, means, status)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cpts(:)
    real(dp), allocatable, intent(out) :: means(:)
    integer, intent(out) :: status
    integer, allocatable :: b(:)
    integer :: i, nseg

    call make_boundaries(size(data), cpts, b, status)
    if (status /= 0) then
        allocate(means(0))
        return
    end if
    nseg = size(b) - 1
    allocate(means(nseg))
    do i = 1, nseg
        means(i) = sum(data(b(i - 1) + 1:b(i))) / real(b(i) - b(i - 1), dp)
    end do
end subroutine cp_segment_means

subroutine cp_segment_variances_mle(data, cpts, variances, status)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cpts(:)
    real(dp), allocatable, intent(out) :: variances(:)
    integer, intent(out) :: status
    integer, allocatable :: b(:)
    integer :: i, nseg, nn
    real(dp) :: mu

    call make_boundaries(size(data), cpts, b, status)
    if (status /= 0) then
        allocate(variances(0))
        return
    end if
    nseg = size(b) - 1
    allocate(variances(nseg))
    do i = 1, nseg
        nn = b(i) - b(i - 1)
        mu = sum(data(b(i - 1) + 1:b(i))) / real(nn, dp)
        variances(i) = sum((data(b(i - 1) + 1:b(i)) - mu)**2) / real(nn, dp)
    end do
end subroutine cp_segment_variances_mle

subroutine cp_segment_scales(data, cpts, shape, scales, status)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cpts(:)
    real(dp), intent(in) :: shape
    real(dp), allocatable, intent(out) :: scales(:)
    integer, intent(out) :: status
    integer, allocatable :: b(:)
    integer :: i, nseg, nn

    if (shape <= 0.0_dp) then
        status = 1
        allocate(scales(0))
        return
    end if
    call make_boundaries(size(data), cpts, b, status)
    if (status /= 0) then
        allocate(scales(0))
        return
    end if
    nseg = size(b) - 1
    allocate(scales(nseg))
    do i = 1, nseg
        nn = b(i) - b(i - 1)
        scales(i) = sum(data(b(i - 1) + 1:b(i))) / (real(nn, dp) * shape)
    end do
end subroutine cp_segment_scales

subroutine cp_segment_trend_fits(data, cpts, theta_start, theta_end, status)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cpts(:)
    real(dp), allocatable, intent(out) :: theta_start(:), theta_end(:)
    integer, intent(out) :: status
    integer, allocatable :: b(:)
    integer :: i, j, nseg, nn
    real(dp) :: sum_y, sum_ty, denom

    call make_boundaries(size(data), cpts, b, status)
    if (status /= 0) then
        allocate(theta_start(0), theta_end(0))
        return
    end if
    nseg = size(b) - 1
    allocate(theta_start(nseg), theta_end(nseg))
    do i = 1, nseg
        nn = b(i) - b(i - 1)
        if (nn < 2) then
            status = 2
            return
        end if
        sum_y = 0.0_dp
        sum_ty = 0.0_dp
        do j = 1, nn
            sum_y = sum_y + data(b(i - 1) + j)
            sum_ty = sum_ty + real(j, dp) * data(b(i - 1) + j)
        end do
        denom = 2.0_dp * real(nn, dp) * real(2 * nn + 1, dp) - &
            3.0_dp * real(nn * (nn + 1), dp)
        theta_start(i) = (2.0_dp * sum_y * real(2 * nn + 1, dp) - 6.0_dp * sum_ty) / denom
        theta_end(i) = 6.0_dp * sum_ty / real((nn + 1) * (2 * nn + 1), dp) + &
            theta_start(i) * (1.0_dp - 3.0_dp * real(nn, dp) / real(2 * nn + 1, dp))
    end do
    status = 0
end subroutine cp_segment_trend_fits

subroutine cp_segment_regression_fits(y, x, cpts, beta, sigma2, status, rcond)
    real(dp), intent(in) :: y(:), x(:, :)
    integer, intent(in) :: cpts(:)
    real(dp), allocatable, intent(out) :: beta(:, :), sigma2(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: rcond
    integer, allocatable :: b(:)
    real(dp), allocatable :: btmp(:)
    real(dp) :: rss
    integer :: i, nseg, nn, p, st

    if (size(x, 1) /= size(y)) then
        status = 1
        allocate(beta(0, 0), sigma2(0))
        return
    end if
    call make_boundaries(size(y), cpts, b, status)
    if (status /= 0) then
        allocate(beta(0, 0), sigma2(0))
        return
    end if
    p = size(x, 2)
    nseg = size(b) - 1
    allocate(beta(nseg, p), sigma2(nseg))
    do i = 1, nseg
        nn = b(i) - b(i - 1)
        if (present(rcond)) then
            call cp_regression_segment_fit(y(b(i - 1) + 1:b(i)), x(b(i - 1) + 1:b(i), :), btmp, rss, st, rcond)
        else
            call cp_regression_segment_fit(y(b(i - 1) + 1:b(i)), x(b(i - 1) + 1:b(i), :), btmp, rss, st)
        end if
        if (st /= 0) then
            status = st
            return
        end if
        beta(i, :) = btmp
        if (nn > p) then
            sigma2(i) = rss / real(nn - p, dp)
        else
            sigma2(i) = 0.0_dp
        end if
    end do
    status = 0
end subroutine cp_segment_regression_fits

subroutine make_boundaries(n, cpts, b, status)
    integer, intent(in) :: n, cpts(:)
    integer, allocatable, intent(out) :: b(:)
    integer, intent(out) :: status
    integer :: i, previous_cpt

    status = 0
    allocate(b(0:size(cpts) + 1))
    b(0) = 0
    previous_cpt = 0
    do i = 1, size(cpts)
        if (cpts(i) <= 0 .or. cpts(i) >= n) then
            status = 1
            return
        end if
        if (cpts(i) <= previous_cpt) then
            status = 1
            return
        end if
        b(i) = cpts(i)
        previous_cpt = cpts(i)
    end do
    b(size(cpts) + 1) = n
end subroutine make_boundaries

end module changepoint_fit
