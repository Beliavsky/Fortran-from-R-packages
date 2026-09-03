! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program phylo_example
    use mcmcglmm, only : dp, phylogenetic_precision
    use ape_types, only : phylo_tree
    implicit none

    integer :: info
    real(dp), allocatable :: inverse_relationship(:, :)
    real(dp), allocatable :: relationship(:, :)
    type(phylo_tree) :: tree

    tree%n_tip = 4
    tree%n_node = 3
    allocate(tree%edge(6, 2), tree%edge_length(6))
    tree%edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
    tree%edge_length = 0.5_dp

    call phylogenetic_precision(tree, [1, 2, 3, 4], .true., inverse_relationship, relationship, info)
    if (info /= 0) error stop 'phylogenetic precision failed'

    print '(a,f8.4)', 'relationship between sister tips 1 and 2: ', relationship(1, 2)
    print '(a,f8.4)', 'relationship between tips 1 and 3: ', relationship(1, 3)
end program phylo_example
