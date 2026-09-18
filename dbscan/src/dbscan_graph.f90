module dbscan_graph
    use dbscan_types, only : frnn_result, knn_result
    implicit none
    private

    public :: connected_components_frnn, connected_components_knn
    public :: relabel_components

contains

    subroutine connected_components_frnn(nn, labels, mutual)
        type(frnn_result), intent(in) :: nn !! Adjacency-list graph in fixed-radius CSR representation.
        integer, allocatable, intent(out) :: labels(:) !! One-based connected-component labels.
        logical, intent(in), optional :: mutual !! Require reciprocal adjacency when true.
        integer, allocatable :: parent(:)
        integer :: i
        integer :: p
        integer :: j
        integer :: n
        logical :: is_mutual

        n = nn%n
        allocate(parent(n), labels(n))
        do i = 1, n
            parent(i) = i
        end do
        is_mutual = .false.
        if (present(mutual)) is_mutual = mutual

        do i = 1, n
            do p = nn%offset(i), nn%offset(i + 1) - 1
                j = nn%id(p)
                if (j < 1 .or. j > n) cycle
                if (is_mutual) then
                    if (.not. edge_exists(nn, j, i)) cycle
                end if
                call union_sets(parent, i, j)
            end do
        end do
        do i = 1, n
            labels(i) = find_root(parent, i)
        end do
        call relabel_components(labels)
    end subroutine connected_components_frnn

    subroutine connected_components_knn(nn, labels, mutual)
        type(knn_result), intent(in) :: nn !! k-nearest-neighbor graph.
        integer, allocatable, intent(out) :: labels(:) !! One-based connected-component labels.
        logical, intent(in), optional :: mutual !! Require reciprocal nearest-neighbor links when true.
        integer, allocatable :: parent(:)
        integer :: i
        integer :: j
        integer :: p
        integer :: n
        integer :: k
        logical :: is_mutual

        n = size(nn%id, 1)
        k = size(nn%id, 2)
        allocate(parent(n), labels(n))
        do i = 1, n
            parent(i) = i
        end do
        is_mutual = .false.
        if (present(mutual)) is_mutual = mutual

        do i = 1, n
            do p = 1, k
                j = nn%id(i, p)
                if (j < 1 .or. j > n) cycle
                if (is_mutual) then
                    if (.not. any(nn%id(j, :) == i)) cycle
                end if
                call union_sets(parent, i, j)
            end do
        end do
        do i = 1, n
            labels(i) = find_root(parent, i)
        end do
        call relabel_components(labels)
    end subroutine connected_components_knn

    pure subroutine relabel_components(labels)
        integer, intent(inout) :: labels(:) !! Arbitrary positive component labels to compact in encounter order.
        integer, allocatable :: roots(:)
        integer :: i
        integer :: j
        integer :: nr
        logical :: found

        allocate(roots(size(labels)))
        nr = 0
        do i = 1, size(labels)
            found = .false.
            do j = 1, nr
                if (roots(j) == labels(i)) then
                    labels(i) = j
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                nr = nr + 1
                roots(nr) = labels(i)
                labels(i) = nr
            end if
        end do
    end subroutine relabel_components

    pure logical function edge_exists(nn, i, j) result(found)
        type(frnn_result), intent(in) :: nn !! CSR graph whose row i is searched.
        integer, intent(in) :: i !! One-based source row.
        integer, intent(in) :: j !! One-based target vertex.
        integer :: p

        found = .false.
        do p = nn%offset(i), nn%offset(i + 1) - 1
            if (nn%id(p) == j) then
                found = .true.
                return
            end if
        end do
    end function edge_exists

    recursive integer function find_root(parent, x) result(root)
        integer, intent(inout) :: parent(:) !! Disjoint-set parent array with one-based indices.
        integer, intent(in) :: x !! Vertex whose representative is requested.

        if (parent(x) == x) then
            root = x
        else
            parent(x) = find_root(parent, parent(x))
            root = parent(x)
        end if
    end function find_root

    subroutine union_sets(parent, a, b)
        integer, intent(inout) :: parent(:) !! Disjoint-set parent array.
        integer, intent(in) :: a !! First vertex to merge.
        integer, intent(in) :: b !! Second vertex to merge.
        integer :: ra
        integer :: rb

        ra = find_root(parent, a)
        rb = find_root(parent, b)
        if (ra == rb) return
        if (ra < rb) then
            parent(rb) = ra
        else
            parent(ra) = rb
        end if
    end subroutine union_sets

end module dbscan_graph
