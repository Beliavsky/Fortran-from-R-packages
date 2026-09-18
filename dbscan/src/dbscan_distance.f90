module dbscan_distance
    use dbscan_kinds, only : dp
    use dbscan_utils, only : contains_nan
    implicit none
    private

    public :: euclidean_distance, pairwise_distance_matrix, distance_matrix

contains

    pure real(dp) function euclidean_distance(x, y) result(d)
        real(dp), intent(in) :: x(:) !! Coordinates of the first point.
        real(dp), intent(in) :: y(:) !! Coordinates of the second point; same length as x.
        d = sqrt(sum((x - y) ** 2))
    end function euclidean_distance

    subroutine pairwise_distance_matrix(x, d, metric, status)
        real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows and features in columns.
        real(dp), allocatable, intent(out) :: d(:, :) !! Symmetric pairwise distance matrix.
        character(len=*), intent(in), optional :: metric !! Metric: euclidean, sqeuclidean, or manhattan.
        integer, intent(out), optional :: status !! Zero on success; nonzero for invalid input or metric.
        character(len=:), allocatable :: met
        integer :: i
        integer :: j
        integer :: n
        real(dp) :: v

        if (present(status)) status = 0
        if (size(x, 1) < 1 .or. size(x, 2) < 1) then
            allocate(d(0, 0))
            if (present(status)) status = 1
            return
        end if
        if (contains_nan(x)) then
            allocate(d(0, 0))
            if (present(status)) status = 2
            return
        end if

        met = 'euclidean'
        if (present(metric)) met = trim(adjustl(metric))
        n = size(x, 1)
        allocate(d(n, n))
        d = 0.0_dp
        do i = 1, n
            do j = i + 1, n
                select case (met)
                case ('euclidean')
                    v = sqrt(sum((x(i, :) - x(j, :)) ** 2))
                case ('sqeuclidean')
                    v = sum((x(i, :) - x(j, :)) ** 2)
                case ('manhattan')
                    v = sum(abs(x(i, :) - x(j, :)))
                case default
                    d = 0.0_dp
                    if (present(status)) status = 3
                    return
                end select
                d(i, j) = v
                d(j, i) = v
            end do
        end do
    end subroutine pairwise_distance_matrix

    pure real(dp) function distance_matrix(d, i, j) result(value)
        real(dp), intent(in) :: d(:, :) !! Full symmetric distance matrix.
        integer, intent(in) :: i !! First one-based observation index.
        integer, intent(in) :: j !! Second one-based observation index.
        value = d(i, j)
    end function distance_matrix

end module dbscan_distance
