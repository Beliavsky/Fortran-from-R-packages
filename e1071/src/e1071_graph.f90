module e1071_graph
    use e1071_kinds, only: dp
    implicit none
    private

    type, public :: shortest_paths_result
        real(dp), allocatable :: distance(:, :)
        integer, allocatable :: middle(:, :)
    end type shortest_paths_result

    public :: all_shortest_paths, extract_path

contains

    subroutine all_shortest_paths(cost, result)
        real(dp), intent(in) :: cost(:, :) !! Square edge-cost matrix; absent edges should use a sufficiently large or infinite
        !! cost.
        type(shortest_paths_result), intent(out) :: result !! Floyd-Warshall distances and one-based intermediate-vertex matrix.
        integer :: n
        integer :: i
        integer :: j
        integer :: k

        if (size(cost, 1) /= size(cost, 2)) error stop "all_shortest_paths: cost matrix must be square"
        n = size(cost, 1)
        allocate(result%distance(n, n), result%middle(n, n))
        result%distance = cost
        result%middle = 0
        do i = 1, n
            result%distance(i, i) = 0.0_dp
        end do
        do k = 1, n
            do i = 1, n
                do j = 1, n
                    if (result%distance(i, k) + result%distance(k, j) < result%distance(i, j)) then
                        result%distance(i, j) = result%distance(i, k) + result%distance(k, j)
                        result%middle(i, j) = k
                    end if
                end do
            end do
        end do
    end subroutine all_shortest_paths

    function extract_path(result, from, to) result(path)
        type(shortest_paths_result), intent(in) :: result !! Floyd-Warshall result whose middle matrix defines shortest-path
        !! decomposition.
        integer, intent(in) :: from !! One-based source vertex index.
        integer, intent(in) :: to !! One-based destination vertex index.
        integer, allocatable :: path(:)
        integer, allocatable :: work(:)
        integer :: used
        integer :: n

        n = size(result%middle, 1)
        if (from < 1 .or. from > n .or. to < 1 .or. to > n) error stop "extract_path: vertex index out of range"
        allocate(work(max(2, 2 * n + 2)))
        used = 0
        call append_path(result%middle, from, to, work, used, .true.)
        allocate(path(used))
        path = work(:used)
    end function extract_path

    recursive subroutine append_path(middle, from, to, work, used, include_from)
        integer, intent(in) :: middle(:, :) !! One-based intermediate-vertex matrix produced by all_shortest_paths; zero means
        !! direct path.
        integer, intent(in) :: from !! Source vertex for this recursive path segment.
        integer, intent(in) :: to !! Destination vertex for this recursive path segment.
        integer, intent(inout) :: work(:) !! Preallocated path workspace receiving vertex indices in traversal order.
        integer, intent(inout) :: used !! Number of workspace entries already filled; incremented by this routine.
        logical, intent(in) :: include_from !! If true, append the segment source before its interior and destination vertices.
        integer :: k

        if (include_from) then
            used = used + 1
            if (used > size(work)) error stop "append_path: internal path workspace exhausted"
            work(used) = from
        end if
        k = middle(from, to)
        if (k == 0) then
            used = used + 1
            if (used > size(work)) error stop "append_path: internal path workspace exhausted"
            work(used) = to
        else
            call append_path(middle, from, k, work, used, .false.)
            call append_path(middle, k, to, work, used, .false.)
        end if
    end subroutine append_path

end module e1071_graph
