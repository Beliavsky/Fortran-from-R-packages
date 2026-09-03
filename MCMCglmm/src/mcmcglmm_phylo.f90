! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_phylo
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use ape_types, only : phylo_tree
    use ape_tree_algorithms, only : node_depth_edgelength, mrca
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_matrix, only : sample_mvn_covariance
    implicit none
    private

    public :: phylogenetic_precision
    public :: breeding_values_phylo

contains

    pure subroutine phylogenetic_precision(tree, nodes, scale_tree, inverse_relationship, relationship, info)
        type(phylo_tree), intent(in) :: tree !! Rooted ape-style phylogenetic tree with branch lengths.
        integer, intent(in) :: nodes(:) !! One-based ape node numbers to include in the Brownian relationship matrix.
        logical, intent(in) :: scale_tree !! If true, divide depths by the common root-to-tip height of an ultrametric tree.
        real(dp), allocatable, intent(out) :: inverse_relationship(:, :) !! Dense inverse Brownian relationship matrix.
        real(dp), allocatable, intent(out) :: relationship(:, :) !! Brownian relationship matrix from root-to-MRCA depth.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid tree/nodes, scaling, or inversion failure.
        real(dp), allocatable :: depth(:)
        real(dp) :: height
        real(dp) :: tolerance
        integer :: ancestor
        integer :: i
        integer :: j
        integer :: total_nodes

        info = 0
        if (.not. tree%valid() .or. .not. tree%has_lengths()) then
            allocate(inverse_relationship(0, 0), relationship(0, 0))
            info = 1
            return
        end if
        total_nodes = tree%total_nodes()
        if (size(nodes) < 1 .or. any(nodes < 1) .or. any(nodes > total_nodes)) then
            allocate(inverse_relationship(0, 0), relationship(0, 0))
            info = 2
            return
        end if
        call node_depth_edgelength(tree, depth, info)
        if (info /= 0) then
            allocate(inverse_relationship(0, 0), relationship(0, 0))
            return
        end if
        height = 1.0_dp
        if (scale_tree) then
            height = depth(1)
            tolerance = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(height))
            if (height <= 0.0_dp .or. any(abs(depth(1:tree%n_tip) - height) > tolerance)) then
                allocate(inverse_relationship(0, 0), relationship(0, 0))
                info = 3
                return
            end if
        end if
        allocate(relationship(size(nodes), size(nodes)))
        do j = 1, size(nodes)
            relationship(j, j) = depth(nodes(j)) / height
            do i = 1, j - 1
                ancestor = mrca(tree, nodes(i), nodes(j))
                if (ancestor == 0) then
                    allocate(inverse_relationship(0, 0))
                    relationship = 0.0_dp
                    info = 4
                    return
                end if
                relationship(i, j) = depth(ancestor) / height
                relationship(j, i) = relationship(i, j)
            end do
        end do
        call inverse_matrix(relationship, inverse_relationship, info)
    end subroutine phylogenetic_precision

    pure subroutine breeding_values_phylo(state, tree, covariance, scale_tree, values, info)
        type(rng_state), intent(inout) :: state !! Generator state consumed by Brownian branch innovations.
        type(phylo_tree), intent(in) :: tree !! Rooted ape-style phylogenetic tree with one branch length per edge.
        real(dp), intent(in) :: covariance(:, :) !! Positive-definite trait covariance matrix G for unit branch length.
        logical, intent(in) :: scale_tree !! If true, scale all branches by the common ultrametric root-to-tip height.
        real(dp), allocatable, intent(out) :: values(:, :) !! All-node by trait Brownian breeding values; root is zero.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid inputs, scaling, or unresolved traversal.
        real(dp), allocatable :: depth(:)
        real(dp), allocatable :: innovation(:)
        real(dp), allocatable :: scaled_covariance(:, :)
        real(dp), allocatable :: zero_mean(:)
        real(dp) :: branch_length
        real(dp) :: height
        real(dp) :: tolerance
        logical, allocatable :: assigned(:)
        logical :: progress
        integer :: child
        integer :: edge_index
        integer :: parent
        integer :: pass
        integer :: root_node
        integer :: total_nodes
        integer :: traits

        info = 0
        traits = size(covariance, 1)
        if (.not. tree%valid() .or. .not. tree%has_lengths() .or. size(covariance, 2) /= traits .or. traits < 1) then
            allocate(values(0, 0))
            info = 1
            return
        end if
        if (any(tree%edge_length < 0.0_dp)) then
            allocate(values(0, 0))
            info = 2
            return
        end if
        call node_depth_edgelength(tree, depth, info)
        if (info /= 0) then
            allocate(values(0, 0))
            return
        end if
        height = 1.0_dp
        if (scale_tree) then
            height = depth(1)
            tolerance = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(height))
            if (height <= 0.0_dp .or. any(abs(depth(1:tree%n_tip) - height) > tolerance)) then
                allocate(values(0, 0))
                info = 3
                return
            end if
        end if
        total_nodes = tree%total_nodes()
        root_node = tree%root()
        if (root_node == 0) then
            allocate(values(0, 0))
            info = 4
            return
        end if
        allocate(values(total_nodes, traits), assigned(total_nodes), zero_mean(traits), scaled_covariance(traits, traits))
        values = 0.0_dp
        assigned = .false.
        zero_mean = 0.0_dp
        assigned(root_node) = .true.
        do pass = 1, total_nodes
            progress = .false.
            do edge_index = 1, tree%nedge()
                parent = tree%edge(edge_index, 1)
                child = tree%edge(edge_index, 2)
                if (assigned(child) .or. .not. assigned(parent)) cycle
                branch_length = tree%edge_length(edge_index) / height
                scaled_covariance = branch_length * covariance
                if (branch_length > 0.0_dp) then
                    call sample_mvn_covariance(state, zero_mean, scaled_covariance, innovation, info)
                    if (info /= 0) return
                    values(child, :) = values(parent, :) + innovation
                else
                    values(child, :) = values(parent, :)
                end if
                assigned(child) = .true.
                progress = .true.
            end do
            if (all(assigned)) exit
            if (.not. progress) exit
        end do
        if (.not. all(assigned)) info = 5
    end subroutine breeding_values_phylo

end module mcmcglmm_phylo
