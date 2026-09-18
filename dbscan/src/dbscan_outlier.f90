module dbscan_outlier
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_is_nan
    use dbscan_kinds, only : dp
    use dbscan_types, only : knn_result, frnn_result
    use dbscan_nn, only : knn, frnn
    implicit none
    private

    real(dp), parameter :: pi_dp = acos(-1.0_dp)

    public :: lof, pointdensity

contains

    function lof(x, min_pts, status) result(score)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        integer, intent(in) :: min_pts !! Neighborhood size including the point itself; at least two.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        real(dp), allocatable :: score(:)
        type(knn_result) :: nn
        real(dp), allocatable :: kdist(:)
        real(dp), allocatable :: lrd(:)
        integer, allocatable :: count_n(:)
        integer, allocatable :: offsets(:)
        integer, allocatable :: ids(:)
        real(dp), allocatable :: dists(:)
        integer :: i
        integer :: j
        integer :: n
        integer :: p
        integer :: pos
        integer :: stat
        real(dp) :: reach_sum
        real(dp) :: kth
        real(dp) :: tol

        n = size(x, 1)
        if (min_pts < 2 .or. min_pts > n) then
            if (present(status)) status = 1
            allocate(score(0))
            return
        end if
        call knn(x, min_pts - 1, nn, stat)
        if (stat /= 0) then
            if (present(status)) status = stat
            allocate(score(0))
            return
        end if
        allocate(kdist(n), count_n(n), offsets(n + 1))
        kdist = nn%dist(:, min_pts - 1)
        do i = 1, n
            kth = kdist(i)
            tol = epsilon(1.0_dp) * max(1.0_dp, kth)
            count_n(i) = 0
            do j = 1, n
                if (j == i) cycle
                if (sqrt(sum((x(i, :) - x(j, :)) ** 2)) <= kth + tol) count_n(i) = count_n(i) + 1
            end do
        end do
        offsets(1) = 1
        do i = 1, n
            offsets(i + 1) = offsets(i) + count_n(i)
        end do
        allocate(ids(sum(count_n)), dists(sum(count_n)))
        do i = 1, n
            pos = offsets(i)
            kth = kdist(i)
            tol = epsilon(1.0_dp) * max(1.0_dp, kth)
            do j = 1, n
                if (j == i) cycle
                if (sqrt(sum((x(i, :) - x(j, :)) ** 2)) <= kth + tol) then
                    ids(pos) = j
                    dists(pos) = sqrt(sum((x(i, :) - x(j, :)) ** 2))
                    pos = pos + 1
                end if
            end do
        end do

        allocate(lrd(n), score(n))
        do i = 1, n
            reach_sum = 0.0_dp
            do p = offsets(i), offsets(i + 1) - 1
                j = ids(p)
                reach_sum = reach_sum + max(kdist(j), dists(p))
            end do
            if (abs(reach_sum) <= tiny(1.0_dp)) then
                lrd(i) = ieee_value(1.0_dp, ieee_positive_inf)
            else
                lrd(i) = real(count_n(i), dp) / reach_sum
            end if
        end do
        do i = 1, n
            score(i) = 0.0_dp
            do p = offsets(i), offsets(i + 1) - 1
                score(i) = score(i) + lrd(ids(p))
            end do
            score(i) = score(i) / real(count_n(i), dp) / lrd(i)
            if (ieee_is_nan(score(i))) score(i) = 1.0_dp
        end do
        if (present(status)) status = 0
    end function lof

    function pointdensity(x, eps, density_type, status) result(value)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
        real(dp), intent(in) :: eps !! Neighborhood radius or Gaussian kernel standard deviation.
        character(len=*), intent(in), optional :: density_type !! frequency, density, or gaussian.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input.
        real(dp), allocatable :: value(:)
        type(frnn_result) :: nn
        character(len=:), allocatable :: kind
        integer :: i
        integer :: p
        integer :: stat
        real(dp) :: sigma2

        kind = 'frequency'
        if (present(density_type)) kind = trim(adjustl(density_type))
        if (eps <= 0.0_dp) then
            if (present(status)) status = 1
            allocate(value(0))
            return
        end if
        if (kind == 'gaussian') then
            call frnn(x, 3.0_dp * eps, nn, stat)
        else
            call frnn(x, eps, nn, stat)
        end if
        if (stat /= 0) then
            if (present(status)) status = stat
            allocate(value(0))
            return
        end if
        allocate(value(nn%n))
        select case (kind)
        case ('frequency')
            do i = 1, nn%n
                value(i) = real(1 + nn%offset(i + 1) - nn%offset(i), dp)
            end do
        case ('density')
            do i = 1, nn%n
                value(i) = real(1 + nn%offset(i + 1) - nn%offset(i), dp)
            end do
            value = value / (2.0_dp * eps * real(nn%n, dp))
        case ('gaussian')
            sigma2 = eps ** 2
            do i = 1, nn%n
                value(i) = 0.0_dp
                do p = nn%offset(i), nn%offset(i + 1) - 1
                    value(i) = value(i) + exp(-(nn%dist(p) ** 2) / (2.0_dp * sigma2))
                end do
                value(i) = value(i) / (real(nn%n, dp) * eps * 2.0_dp * pi_dp)
            end do
        case default
            deallocate(value)
            allocate(value(0))
            if (present(status)) status = 2
            return
        end select
        if (present(status)) status = 0
    end function pointdensity

end module dbscan_outlier
