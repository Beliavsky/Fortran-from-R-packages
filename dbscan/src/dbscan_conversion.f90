module dbscan_conversion
    use dbscan_kinds, only : dp
    use dbscan_types, only : optics_result, dendrogram_result
    implicit none
    private

    public :: optics_to_dendrogram, dendrogram_to_reachability

contains

    subroutine optics_to_dendrogram(object, tree)
        type(optics_result), intent(in) :: object !! OPTICS ordering to convert to a reachability dendrogram.
        type(dendrogram_result), intent(out) :: tree !! Chain-based dendrogram preserving OPTICS leaf order.
        integer :: i
        integer :: n
        real(dp) :: h

        n = size(object%order)
        if (n < 2) then
            allocate(tree%merge(0, 2), tree%height(0), tree%order(n))
            if (n == 1) tree%order = object%order
            return
        end if
        allocate(tree%merge(n - 1, 2), tree%height(n - 1), tree%order(n))
        tree%order = object%order
        tree%merge(1, 1) = -object%order(1)
        tree%merge(1, 2) = -object%order(2)
        h = object%reachdist(object%order(2))
        if (h >= huge(1.0_dp) / 2.0_dp) h = object%coredist(object%order(1))
        tree%height(1) = h
        do i = 2, n - 1
            tree%merge(i, 1) = i - 1
            tree%merge(i, 2) = -object%order(i + 1)
            h = object%reachdist(object%order(i + 1))
            if (h >= huge(1.0_dp) / 2.0_dp) h = tree%height(i - 1)
            tree%height(i) = max(tree%height(i - 1), h)
        end do
    end subroutine optics_to_dendrogram

    subroutine dendrogram_to_reachability(tree, object)
        type(dendrogram_result), intent(in) :: tree !! Dendrogram with a leaf ordering and merge heights.
        type(optics_result), intent(out) :: object !! Reachability-style representation reconstructed from the tree.
        integer :: i
        integer :: n

        n = size(tree%order)
        allocate(object%order(n), object%reachdist(n), object%coredist(n), object%predecessor(n))
        object%order = tree%order
        object%reachdist = huge(1.0_dp)
        object%coredist = huge(1.0_dp)
        object%predecessor = 0
        if (n >= 1) object%reachdist(object%order(1)) = huge(1.0_dp)
        do i = 2, n
            object%reachdist(object%order(i)) = tree%height(i - 1)
            object%predecessor(object%order(i)) = object%order(i - 1)
        end do
    end subroutine dendrogram_to_reachability

end module dbscan_conversion
