module proxy_string_measures
    use proxy_kinds, only: dp
    implicit none
    private

    public :: levenshtein_distance, levenshtein_auto, levenshtein_cross

contains

    pure function levenshtein_distance(x, y) result(distance)
        character(len=*), intent(in) :: x !! First string; distance counts single-character insertions, deletions, and
        !! substitutions.
        character(len=*), intent(in) :: y !! Second string compared with `x`; trailing Fortran padding is ignored with `len_trim`.
        real(dp) :: distance
        integer, allocatable :: previous(:)
        integer, allocatable :: current(:)
        integer :: cost
        integer :: i
        integer :: j
        integer :: nx
        integer :: ny

        nx = len_trim(x)
        ny = len_trim(y)
        allocate(previous(0:ny), current(0:ny))
        do j = 0, ny
            previous(j) = j
        end do
        do i = 1, nx
            current(0) = i
            do j = 1, ny
                cost = 1
                if (x(i:i) == y(j:j)) cost = 0
                current(j) = min(previous(j) + 1, current(j - 1) + 1, previous(j - 1) + cost)
            end do
            previous = current
        end do
        distance = real(previous(ny), dp)
    end function levenshtein_distance

    subroutine levenshtein_auto(x, distance)
        character(len=*), intent(in) :: x(:) !! Array of strings for which all pairwise Levenshtein distances are requested.
        real(dp), allocatable, intent(out) :: distance(:, :) !! Allocated symmetric square matrix with zero diagonal and edit
        !! distances off diagonal.
        integer :: i
        integer :: j

        allocate(distance(size(x), size(x)))
        do j = 1, size(x)
            do i = 1, size(x)
                distance(i, j) = levenshtein_distance(x(i), x(j))
            end do
        end do
    end subroutine levenshtein_auto

    subroutine levenshtein_cross(x, y, distance)
        character(len=*), intent(in) :: x(:) !! First string collection defining rows of the cross-distance matrix.
        character(len=*), intent(in) :: y(:) !! Second string collection defining columns of the cross-distance matrix.
        real(dp), allocatable, intent(out) :: distance(:, :) !! Allocated `size(x)` by `size(y)` matrix of Levenshtein edit
        !! distances.
        integer :: i
        integer :: j

        allocate(distance(size(x), size(y)))
        do j = 1, size(y)
            do i = 1, size(x)
                distance(i, j) = levenshtein_distance(x(i), y(j))
            end do
        end do
    end subroutine levenshtein_cross

end module proxy_string_measures
