! Elimination-tree and fundamental-supernode analysis for SparseM.
!
! The algorithms in this module are modern free-form translations of the
! symbolic initialization stage distributed with SparseM's Ng-Peyton sparse
! Cholesky implementation. Source/algorithm credits preserved from upstream
! include Joseph W. H. Liu, Esmond G. Ng, Barry W. Peyton, and the Oak Ridge
! National Laboratory Mathematical Sciences Section. See NOTICE.md and LICENSE.
!
module sparsem_ng_symbolic
    use sparsem_types, only: csr_matrix
    use sparsem_ng_order, only: csr_symmetric_adjacency
    implicit none
    private

    public :: ng_symbolic_analysis, ng_supernodal_structure

contains

    subroutine ng_symbolic_analysis(a, perm, inv_perm, parent, col_count, nsuper, supernode, xsuper, &
            nnzl, nsub, info)
        type(csr_matrix), intent(in) :: a !! Square sparse matrix whose structural graph is analyzed.
        integer, intent(inout) :: perm(:) !! New-to-old ordering; postordered and child-reordered on return.
        integer, intent(inout) :: inv_perm(:) !! Old-to-new inverse ordering updated consistently with perm.
        integer, allocatable, intent(out) :: parent(:) !! Elimination-tree parent vector in final symbolic ordering.
        integer, allocatable, intent(out) :: col_count(:) !! Number of structural factor nonzeros in each column.
        integer, intent(out) :: nsuper !! Number of fundamental supernodes in the final symbolic ordering.
        integer, allocatable, intent(out) :: supernode(:) !! Supernode membership index for every ordered column.
        integer, allocatable, intent(out) :: xsuper(:) !! Supernode partition pointers; first nsuper+1 entries are active.
        integer, intent(out) :: nnzl !! Total structural nonzeros in the lower Cholesky factor, including the diagonal.
        integer, intent(out) :: nsub !! Number of compressed supernodal row subscripts estimated by FSUP1.
        integer, intent(out) :: info !! Zero on success; negative for nonsquare, permutation, or adjacency errors.
        integer, allocatable :: adjacency(:)
        integer, allocatable :: brother(:)
        integer, allocatable :: fdesc(:)
        integer, allocatable :: first_son(:)
        integer, allocatable :: invpos(:)
        integer, allocatable :: level(:)
        integer, allocatable :: nchild(:)
        integer, allocatable :: previous_leaf(:)
        integer, allocatable :: previous_neighbor(:)
        integer, allocatable :: row_count(:)
        integer, allocatable :: set_parent(:)
        integer, allocatable :: stack(:)
        integer, allocatable :: weight(:)
        integer, allocatable :: xadj(:)
        integer :: n
        integer :: stat

        info = 0
        nsuper = 0
        nnzl = 0
        nsub = 0
        n = a%nrow
        if (a%nrow /= a%ncol) then
            call allocate_empty_outputs(parent, col_count, supernode, xsuper)
            info = -1
            return
        end if
        if (size(perm) /= n .or. size(inv_perm) /= n) then
            call allocate_empty_outputs(parent, col_count, supernode, xsuper)
            info = -2
            return
        end if
        if (n == 0) then
            allocate(parent(0), col_count(0), supernode(0), xsuper(1))
            xsuper(1) = 1
            return
        end if

        call csr_symmetric_adjacency(a, xadj, adjacency, stat)
        if (stat /= 0) then
            call allocate_empty_outputs(parent, col_count, supernode, xsuper)
            info = stat
            return
        end if

        allocate(parent(n), col_count(n), supernode(n), xsuper(n + 1))
        allocate(first_son(n), brother(n), invpos(n), row_count(n), set_parent(n))
        allocate(previous_leaf(n), previous_neighbor(n), stack(n))
        allocate(level(0:n), weight(0:n), fdesc(0:n), nchild(0:n))

        call etordr_core(n, xadj, adjacency, perm, inv_perm, parent, first_son, brother, invpos)
        call fcnthn_core(n, xadj, adjacency, perm, inv_perm, parent, row_count, col_count, nnzl, &
            set_parent, previous_leaf, level, weight, fdesc, nchild, previous_neighbor)
        call chordr_core(n, perm, inv_perm, col_count, parent, first_son, brother, invpos, stack)
        call fsup1_core(n, parent, col_count, nsub, nsuper, supernode)
        call fsup2_core(n, nsuper, supernode, xsuper)
    end subroutine ng_symbolic_analysis

    subroutine allocate_empty_outputs(parent, col_count, supernode, xsuper)
        integer, allocatable, intent(out) :: parent(:) !! Elimination-tree output to initialize as empty.
        integer, allocatable, intent(out) :: col_count(:) !! Column-count output to initialize as empty.
        integer, allocatable, intent(out) :: supernode(:) !! Supernode-membership output to initialize as empty.
        integer, allocatable, intent(out) :: xsuper(:) !! Supernode-pointer output initialized to a single one.

        allocate(parent(0), col_count(0), supernode(0), xsuper(1))
        xsuper(1) = 1
    end subroutine allocate_empty_outputs

    subroutine etree_core(neqns, xadj, adjacency, perm, inv_perm, parent, ancestor)
        integer, intent(in) :: neqns !! Number of equations/vertices in the ordered graph.
        integer, intent(in) :: xadj(:) !! One-based adjacency pointers in original vertex numbering.
        integer, intent(in) :: adjacency(:) !! Undirected adjacency indices in original vertex numbering.
        integer, intent(in) :: perm(:) !! New-to-old ordering used to traverse graph vertices.
        integer, intent(in) :: inv_perm(:) !! Old-to-new inverse ordering used to compare neighbors.
        integer, intent(out) :: parent(:) !! Elimination-tree parent in ordered numbering; zero marks a root.
        integer, intent(out) :: ancestor(:) !! Path-compression workspace for elimination-tree construction.
        integer :: i
        integer :: j
        integer :: jstop
        integer :: jstart
        integer :: neighbor
        integer :: next_node
        integer :: node

        if (neqns <= 0) return
        do i = 1, neqns
            parent(i) = 0
            ancestor(i) = 0
            node = perm(i)
            jstart = xadj(node)
            jstop = xadj(node + 1) - 1
            if (jstart <= jstop) then
                do j = jstart, jstop
                    neighbor = inv_perm(adjacency(j))
                    if (neighbor < i) then
100                     continue
                        if (ancestor(neighbor) == i) cycle
                        if (ancestor(neighbor) > 0) then
                            next_node = ancestor(neighbor)
                            ancestor(neighbor) = i
                            neighbor = next_node
                            go to 100
                        end if
                        parent(neighbor) = i
                        ancestor(neighbor) = i
                    end if
                end do
            end if
        end do
    end subroutine etree_core

    subroutine betree_core(neqns, parent, first_son, brother)
        integer, intent(in) :: neqns !! Number of nodes in the elimination tree or forest.
        integer, intent(in) :: parent(:) !! Parent vector whose positive entries exceed their children in postorder.
        integer, intent(out) :: first_son(:) !! First-child links of the binary tree representation.
        integer, intent(out) :: brother(:) !! Next-sibling links of the binary tree representation.
        integer :: last_root
        integer :: node
        integer :: node_parent

        if (neqns <= 0) return
        first_son = 0
        brother = 0
        last_root = neqns
        if (neqns <= 1) return
        do node = neqns - 1, 1, -1
            node_parent = parent(node)
            if (node_parent <= 0 .or. node_parent == node) then
                brother(last_root) = node
                last_root = node
            else
                brother(node) = first_son(node_parent)
                first_son(node_parent) = node
            end if
        end do
        brother(last_root) = 0
    end subroutine betree_core

    subroutine etpost_core(root, first_son, brother, invpos, parent, stack)
        integer, intent(in) :: root !! Root index of the binary elimination-tree representation.
        integer, intent(in) :: first_son(:) !! First-child links used for depth-first postorder traversal.
        integer, intent(inout) :: brother(:) !! Sibling links, reused temporarily to hold reordered parents.
        integer, intent(out) :: invpos(:) !! Old-tree-position to new-postorder-position inverse permutation.
        integer, intent(inout) :: parent(:) !! Parent vector rewritten into the new postorder numbering.
        integer, intent(inout) :: stack(:) !! Integer traversal stack of length at least root.
        integer :: itop
        integer :: node_parent
        integer :: node
        integer :: num
        integer :: new_node

        num = 0
        itop = 0
        node = root
100     continue
        itop = itop + 1
        stack(itop) = node
        node = first_son(node)
        if (node > 0) go to 100
200     continue
        if (itop <= 0) go to 300
        node = stack(itop)
        itop = itop - 1
        num = num + 1
        invpos(node) = num
        node = brother(node)
        if (node <= 0) go to 200
        go to 100
300     continue
        do node = 1, num
            new_node = invpos(node)
            node_parent = parent(node)
            if (node_parent > 0) node_parent = invpos(node_parent)
            brother(new_node) = node_parent
        end do
        do new_node = 1, num
            parent(new_node) = brother(new_node)
        end do
    end subroutine etpost_core

    subroutine invinv_core(neqns, inv_perm, invpos, perm)
        integer, intent(in) :: neqns !! Number of equations/vertices in the composed permutations.
        integer, intent(inout) :: inv_perm(:) !! Existing inverse permutation composed with invpos on return.
        integer, intent(in) :: invpos(:) !! Inverse permutation of the new equivalent postordering.
        integer, intent(out) :: perm(:) !! New-to-old permutation inverse to the updated inv_perm.
        integer :: i
        integer :: intermediate
        integer :: node

        do i = 1, neqns
            intermediate = inv_perm(i)
            inv_perm(i) = invpos(intermediate)
        end do
        do i = 1, neqns
            node = inv_perm(i)
            perm(node) = i
        end do
    end subroutine invinv_core

    subroutine etordr_core(neqns, xadj, adjacency, perm, inv_perm, parent, first_son, brother, invpos)
        integer, intent(in) :: neqns !! Number of equations/vertices whose elimination tree is postordered.
        integer, intent(in) :: xadj(:) !! One-based adjacency pointers in original vertex numbering.
        integer, intent(in) :: adjacency(:) !! Undirected adjacency indices in original vertex numbering.
        integer, intent(inout) :: perm(:) !! New-to-old ordering updated to an elimination-tree postorder.
        integer, intent(inout) :: inv_perm(:) !! Old-to-new inverse ordering updated consistently with perm.
        integer, intent(out) :: parent(:) !! Parent vector of the postordered elimination tree.
        integer, intent(inout) :: first_son(:) !! First-child workspace for tree postordering.
        integer, intent(inout) :: brother(:) !! Sibling workspace for tree postordering.
        integer, intent(inout) :: invpos(:) !! Postordering inverse-permutation and traversal-stack workspace.

        call etree_core(neqns, xadj, adjacency, perm, inv_perm, parent, invpos)
        call betree_core(neqns, parent, first_son, brother)
        call etpost_core(neqns, first_son, brother, invpos, parent, perm)
        call invinv_core(neqns, inv_perm, invpos, perm)
    end subroutine etordr_core

    subroutine fcnthn_core(neqns, xadj, adjacency, perm, inv_perm, etree_parent, row_count, col_count, &
            nnzl, set_parent, previous_leaf, level, weight, first_descendant, nchild, previous_neighbor)
        integer, intent(in) :: neqns !! Number of equations/vertices in the postordered graph.
        integer, intent(in) :: xadj(:) !! One-based adjacency pointers in original vertex numbering.
        integer, intent(in) :: adjacency(:) !! Undirected adjacency indices in original vertex numbering.
        integer, intent(in) :: perm(:) !! New-to-old postordering used to inspect original graph rows.
        integer, intent(in) :: inv_perm(:) !! Old-to-new inverse postordering for graph neighbors.
        integer, intent(in) :: etree_parent(:) !! Elimination-tree parent vector in postorder numbering.
        integer, intent(out) :: row_count(:) !! Structural lower-factor nonzero counts by row, including diagonal.
        integer, intent(out) :: col_count(:) !! Structural lower-factor nonzero counts by column, including diagonal.
        integer, intent(out) :: nnzl !! Total number of structural nonzeros in the lower Cholesky factor.
        integer, intent(out) :: set_parent(:) !! Disjoint-set parent workspace for least-common-ancestor searches.
        integer, intent(out) :: previous_leaf(:) !! Previous low-neighbor leaf recorded for each high neighbor.
        integer, intent(out) :: level(0:) !! Elimination-tree levels, with element zero used for the artificial root.
        integer, intent(out) :: weight(0:) !! Column-count accumulation weights, including artificial-root slot zero.
        integer, intent(out) :: first_descendant(0:) !! Lowest-numbered descendant for each tree node and root zero.
        integer, intent(out) :: nchild(0:) !! Number of elimination-tree children for each node and root zero.
        integer, intent(out) :: previous_neighbor(:) !! Previous low neighbor processed for each high neighbor.
        integer :: high_neighbor
        integer :: first_desc
        integer :: j
        integer :: jstop
        integer :: jstart
        integer :: k
        integer :: last1
        integer :: last2
        integer :: lca
        integer :: leaf_flag
        integer :: low_neighbor
        integer :: old_neighbor
        integer :: parent
        integer :: previous
        integer :: temp
        integer :: xsuper_root

        xsuper_root = 1
        level(0) = 0
        do k = neqns, 1, -1
            row_count(k) = 1
            col_count(k) = 0
            set_parent(k) = k
            previous_leaf(k) = 0
            level(k) = level(etree_parent(k)) + 1
            weight(k) = 1
            first_descendant(k) = k
            nchild(k) = 0
            previous_neighbor(k) = 0
        end do
        nchild(0) = 0
        first_descendant(0) = 0

        do k = 1, neqns
            parent = etree_parent(k)
            weight(parent) = 0
            nchild(parent) = nchild(parent) + 1
            first_desc = first_descendant(k)
            if (first_desc < first_descendant(parent)) first_descendant(parent) = first_desc
        end do

        do low_neighbor = 1, neqns
            leaf_flag = 0
            first_desc = first_descendant(low_neighbor)
            old_neighbor = perm(low_neighbor)
            jstart = xadj(old_neighbor)
            jstop = xadj(old_neighbor + 1) - 1
            do j = jstart, jstop
                high_neighbor = inv_perm(adjacency(j))
                if (high_neighbor > low_neighbor) then
                    if (first_desc > previous_neighbor(high_neighbor)) then
                        weight(low_neighbor) = weight(low_neighbor) + 1
                        previous = previous_leaf(high_neighbor)
                        if (previous == 0) then
                            row_count(high_neighbor) = row_count(high_neighbor) + &
                                level(low_neighbor) - level(high_neighbor)
                        else
                            last1 = previous
                            last2 = set_parent(last1)
                            lca = set_parent(last2)
300                         continue
                            if (lca /= last2) then
                                set_parent(last1) = lca
                                last1 = lca
                                last2 = set_parent(last1)
                                lca = set_parent(last2)
                                go to 300
                            end if
                            row_count(high_neighbor) = row_count(high_neighbor) + &
                                level(low_neighbor) - level(lca)
                            weight(lca) = weight(lca) - 1
                        end if
                        previous_leaf(high_neighbor) = low_neighbor
                        leaf_flag = 1
                    end if
                    previous_neighbor(high_neighbor) = low_neighbor
                end if
            end do
            parent = etree_parent(low_neighbor)
            weight(parent) = weight(parent) - 1
            if (leaf_flag == 1 .or. nchild(low_neighbor) >= 2) xsuper_root = low_neighbor
            set_parent(xsuper_root) = parent
        end do

        nnzl = 0
        do k = 1, neqns
            temp = col_count(k) + weight(k)
            col_count(k) = temp
            nnzl = nnzl + temp
            parent = etree_parent(k)
            if (parent /= 0) col_count(parent) = col_count(parent) + temp
        end do
    end subroutine fcnthn_core

    subroutine btree2_core(neqns, parent, col_count, first_son, brother, last_son)
        integer, intent(in) :: neqns !! Number of nodes in the elimination tree or forest.
        integer, intent(in) :: parent(:) !! Parent vector in the current postorder numbering.
        integer, intent(in) :: col_count(:) !! Factor column counts used to select the last child of each parent.
        integer, intent(out) :: first_son(:) !! First-child links after child ordering by column count.
        integer, intent(out) :: brother(:) !! Sibling links after child ordering by column count.
        integer, intent(out) :: last_son(:) !! Last-child workspace for each parent.
        integer :: last_child
        integer :: last_root
        integer :: node
        integer :: node_parent

        if (neqns <= 0) return
        first_son = 0
        brother = 0
        last_son = 0
        last_root = neqns
        if (neqns <= 1) return
        do node = neqns - 1, 1, -1
            node_parent = parent(node)
            if (node_parent <= 0 .or. node_parent == node) then
                brother(last_root) = node
                last_root = node
            else
                last_child = last_son(node_parent)
                if (last_child /= 0) then
                    if (col_count(node) >= col_count(last_child)) then
                        brother(node) = first_son(node_parent)
                        first_son(node_parent) = node
                    else
                        brother(last_child) = node
                        last_son(node_parent) = node
                    end if
                else
                    first_son(node_parent) = node
                    last_son(node_parent) = node
                end if
            end if
        end do
        brother(last_root) = 0
    end subroutine btree2_core

    subroutine epost2_core(root, first_son, brother, invpos, parent, col_count, stack)
        integer, intent(in) :: root !! Root index of the binary elimination-tree representation.
        integer, intent(in) :: first_son(:) !! First-child links for the child-reordered tree.
        integer, intent(inout) :: brother(:) !! Sibling links, reused temporarily for reordered parents.
        integer, intent(out) :: invpos(:) !! Old-tree-position to new-postorder-position inverse permutation.
        integer, intent(inout) :: parent(:) !! Parent vector rewritten into the child-reordered postorder.
        integer, intent(inout) :: col_count(:) !! Factor column counts permuted into the new postorder.
        integer, intent(inout) :: stack(:) !! Traversal and temporary column-count workspace.
        integer :: itop
        integer :: node_parent
        integer :: node
        integer :: num
        integer :: new_node

        num = 0
        itop = 0
        node = root
100     continue
        itop = itop + 1
        stack(itop) = node
        node = first_son(node)
        if (node > 0) go to 100
200     continue
        if (itop <= 0) go to 300
        node = stack(itop)
        itop = itop - 1
        num = num + 1
        invpos(node) = num
        node = brother(node)
        if (node <= 0) go to 200
        go to 100

300     continue
        do node = 1, num
            new_node = invpos(node)
            node_parent = parent(node)
            if (node_parent > 0) node_parent = invpos(node_parent)
            brother(new_node) = node_parent
        end do
        do new_node = 1, num
            parent(new_node) = brother(new_node)
        end do
        do node = 1, num
            new_node = invpos(node)
            stack(new_node) = col_count(node)
        end do
        do node = 1, num
            col_count(node) = stack(node)
        end do
    end subroutine epost2_core

    subroutine chordr_core(neqns, perm, inv_perm, col_count, parent, first_son, brother, invpos, stack)
        integer, intent(in) :: neqns !! Number of equations/vertices in the elimination tree.
        integer, intent(inout) :: perm(:) !! New-to-old ordering updated after child-aware postordering.
        integer, intent(inout) :: inv_perm(:) !! Old-to-new inverse ordering updated consistently with perm.
        integer, intent(inout) :: col_count(:) !! Factor column counts permuted with the child-aware postorder.
        integer, intent(inout) :: parent(:) !! Elimination-tree parent vector updated to the new postorder.
        integer, intent(inout) :: first_son(:) !! First-child workspace for the reordered tree.
        integer, intent(inout) :: brother(:) !! Sibling workspace for the reordered tree.
        integer, intent(inout) :: invpos(:) !! Inverse postordering workspace.
        integer, intent(inout) :: stack(:) !! Traversal and temporary column-count workspace.

        call btree2_core(neqns, parent, col_count, first_son, brother, invpos)
        call epost2_core(neqns, first_son, brother, invpos, parent, col_count, stack)
        call invinv_core(neqns, inv_perm, invpos, perm)
    end subroutine chordr_core

    subroutine fsup1_core(neqns, etree_parent, col_count, nsub, nsuper, supernode)
        integer, intent(in) :: neqns !! Number of columns in the postordered factor pattern.
        integer, intent(in) :: etree_parent(:) !! Elimination-tree parent vector in final postorder numbering.
        integer, intent(in) :: col_count(:) !! Structural nonzero count for each lower-factor column.
        integer, intent(out) :: nsub !! Number of row subscripts in fundamental-supernode compressed storage.
        integer, intent(out) :: nsuper !! Number of fundamental supernodes.
        integer, intent(out) :: supernode(:) !! Supernode membership index for each ordered column.
        integer :: column

        if (neqns <= 0) then
            nsuper = 0
            nsub = 0
            return
        end if
        nsuper = 1
        supernode(1) = 1
        nsub = col_count(1)
        do column = 2, neqns
            if (etree_parent(column - 1) == column) then
                if (col_count(column - 1) == col_count(column) + 1) then
                    supernode(column) = nsuper
                    cycle
                end if
            end if
            nsuper = nsuper + 1
            supernode(column) = nsuper
            nsub = nsub + col_count(column)
        end do
    end subroutine fsup1_core

    subroutine fsup2_core(neqns, nsuper, supernode, xsuper)
        integer, intent(in) :: neqns !! Number of columns represented by the supernode membership vector.
        integer, intent(in) :: nsuper !! Number of fundamental supernodes.
        integer, intent(in) :: supernode(:) !! Supernode membership index for each ordered column.
        integer, intent(out) :: xsuper(:) !! Supernode starting-column pointers, with xsuper(nsuper+1)=neqns+1.
        integer :: column
        integer :: last_super
        integer :: super

        if (neqns <= 0) then
            xsuper(1) = 1
            return
        end if
        last_super = nsuper + 1
        do column = neqns, 1, -1
            super = supernode(column)
            if (super /= last_super) xsuper(last_super) = column + 1
            last_super = super
        end do
        xsuper(1) = 1
    end subroutine fsup2_core


    subroutine ng_supernodal_structure(a, perm, inv_perm, col_count, nsuper, xsuper, supernode, nsub, &
            xlindx, lindx, xlnz, info)
        type(csr_matrix), intent(in) :: a !! Square sparse matrix whose reordered graph defines the factor structure.
        integer, intent(in) :: perm(:) !! New-to-old postordering produced by ng_symbolic_analysis.
        integer, intent(in) :: inv_perm(:) !! Old-to-new inverse postordering corresponding to perm.
        integer, intent(in) :: col_count(:) !! Structural lower-factor nonzero count for every reordered column.
        integer, intent(in) :: nsuper !! Number of fundamental supernodes in the symbolic factorization.
        integer, intent(in) :: xsuper(:) !! Supernode partition pointers, with active extent nsuper+1.
        integer, intent(in) :: supernode(:) !! Supernode membership number for every reordered column.
        integer, intent(in) :: nsub !! Required length of the compressed supernodal row-subscript vector.
        integer, allocatable, intent(out) :: xlindx(:) !! Pointers into lindx for each supernode, length nsuper+1.
        integer, allocatable, intent(out) :: lindx(:) !! Compressed row indices of the supernodal factor structure.
        integer, allocatable, intent(out) :: xlnz(:) !! Column pointers for numerical lower-factor values, length n+1.
        integer, intent(out) :: info !! Zero on success; negative for inconsistent dimensions or symbolic structure.
        integer, allocatable :: adjacency(:)
        integer, allocatable :: marker(:)
        integer, allocatable :: merge_link(:)
        integer, allocatable :: reach_link(:)
        integer, allocatable :: xadj(:)
        integer :: n
        integer :: stat

        info = 0
        n = a%nrow
        if (a%nrow /= a%ncol) then
            allocate(xlindx(1), lindx(0), xlnz(1))
            xlindx(1) = 1
            xlnz(1) = 1
            info = -1
            return
        end if
        if (size(perm) /= n .or. size(inv_perm) /= n .or. size(col_count) /= n .or. &
                size(supernode) /= n .or. size(xsuper) < nsuper + 1 .or. nsub < 0) then
            allocate(xlindx(1), lindx(0), xlnz(1))
            xlindx(1) = 1
            xlnz(1) = 1
            info = -2
            return
        end if
        if (n == 0) then
            allocate(xlindx(1), lindx(0), xlnz(1))
            xlindx(1) = 1
            xlnz(1) = 1
            return
        end if

        call csr_symmetric_adjacency(a, xadj, adjacency, stat)
        if (stat /= 0) then
            allocate(xlindx(1), lindx(0), xlnz(1))
            xlindx(1) = 1
            xlnz(1) = 1
            info = stat
            return
        end if

        allocate(xlindx(nsuper + 1), lindx(nsub), xlnz(n + 1))
        allocate(merge_link(max(1, nsuper)), reach_link(0:n + 1), marker(n))
        call symfc2_core(n, xadj, adjacency, perm, inv_perm, col_count, nsuper, xsuper, supernode, &
            nsub, xlindx, lindx, xlnz, merge_link, reach_link, marker, info)
    end subroutine ng_supernodal_structure

    subroutine symfc2_core(neqns, xadj, adjacency, perm, inv_perm, col_count, nsuper, xsuper, supernode, &
            nofsub, xlindx, lindx, xlnz, merge_link, reach_link, marker, flag)
        integer, intent(in) :: neqns !! Number of equations in the postordered sparse system.
        integer, intent(in) :: xadj(:) !! One-based pointers for the symmetric graph adjacency structure.
        integer, intent(in) :: adjacency(:) !! Symmetric graph adjacency indices in original numbering.
        integer, intent(in) :: perm(:) !! New-to-old postordering used by symbolic factorization.
        integer, intent(in) :: inv_perm(:) !! Old-to-new inverse postordering used for graph indices.
        integer, intent(in) :: col_count(:) !! Structural lower-factor nonzero count for each ordered column.
        integer, intent(in) :: nsuper !! Number of fundamental supernodes.
        integer, intent(in) :: xsuper(:) !! First-column pointers for the fundamental supernodes.
        integer, intent(in) :: supernode(:) !! Supernode membership of each ordered column.
        integer, intent(in) :: nofsub !! Exact available length of the compressed row-subscript vector.
        integer, intent(out) :: xlindx(:) !! Output pointers into lindx for each supernode.
        integer, intent(out) :: lindx(:) !! Output compressed row indices for each supernode.
        integer, intent(out) :: xlnz(:) !! Output column pointers for numerical factor storage.
        integer, intent(inout) :: merge_link(:) !! Workspace linking children in the supernodal elimination tree.
        integer, intent(inout) :: reach_link(0:) !! Workspace linked list holding each supernode reach set.
        integer, intent(inout) :: marker(:) !! Workspace marking indices already inserted in the current reach set.
        integer, intent(out) :: flag !! Zero on success; -2 if structural counts are inconsistent.
        integer :: first_column
        integer :: head
        integer :: i
        integer :: j_begin
        integer :: j_end
        integer :: j_pointer
        integer :: j_super
        integer :: j_width
        integer :: k_nnz
        integer :: k_begin
        integer :: k_end
        integer :: k_pointer
        integer :: k_super
        integer :: length
        integer :: last_column
        integer :: new_index
        integer :: next_index
        integer :: node
        integer :: nz_begin
        integer :: nz_end
        integer :: parent_column
        integer :: parent_super
        integer :: point
        integer :: tail
        integer :: width

        flag = 0
        if (neqns <= 0) return
        if (nofsub < 0) then
            flag = -2
            return
        end if

        nz_end = 0
        head = 0
        tail = neqns + 1
        point = 1
        do i = 1, neqns
            marker(i) = 0
            xlnz(i) = point
            point = point + col_count(i)
        end do
        xlnz(neqns + 1) = point

        point = 1
        do k_super = 1, nsuper
            merge_link(k_super) = 0
            first_column = xsuper(k_super)
            xlindx(k_super) = point
            point = point + col_count(first_column)
        end do
        xlindx(nsuper + 1) = point
        if (point - 1 /= nofsub) then
            flag = -2
            return
        end if

        do k_super = 1, nsuper
            first_column = xsuper(k_super)
            last_column = xsuper(k_super + 1) - 1
            width = last_column - first_column + 1
            length = col_count(first_column)
            k_nnz = 0
            reach_link(head) = tail
            j_super = merge_link(k_super)

            if (j_super > 0) then
                j_width = xsuper(j_super + 1) - xsuper(j_super)
                j_begin = xlindx(j_super) + j_width
                j_end = xlindx(j_super + 1) - 1
                do j_pointer = j_end, j_begin, -1
                    new_index = lindx(j_pointer)
                    k_nnz = k_nnz + 1
                    marker(new_index) = k_super
                    reach_link(new_index) = reach_link(head)
                    reach_link(head) = new_index
                end do

                j_super = merge_link(j_super)
                do while (j_super /= 0 .and. k_nnz < length)
                    j_width = xsuper(j_super + 1) - xsuper(j_super)
                    j_begin = xlindx(j_super) + j_width
                    j_end = xlindx(j_super + 1) - 1
                    next_index = head
                    do j_pointer = j_begin, j_end
                        new_index = lindx(j_pointer)
                        do
                            i = next_index
                            next_index = reach_link(i)
                            if (new_index <= next_index) exit
                        end do
                        if (new_index < next_index) then
                            k_nnz = k_nnz + 1
                            reach_link(i) = new_index
                            reach_link(new_index) = next_index
                            marker(new_index) = k_super
                            next_index = new_index
                        end if
                    end do
                    j_super = merge_link(j_super)
                end do
            end if

            if (k_nnz < length) then
                node = perm(first_column)
                k_begin = xadj(node)
                k_end = xadj(node + 1) - 1
                do k_pointer = k_begin, k_end
                    new_index = inv_perm(adjacency(k_pointer))
                    if (new_index > first_column .and. marker(new_index) /= k_super) then
                        next_index = head
                        do
                            i = next_index
                            next_index = reach_link(i)
                            if (new_index <= next_index) exit
                        end do
                        k_nnz = k_nnz + 1
                        reach_link(i) = new_index
                        reach_link(new_index) = next_index
                        marker(new_index) = k_super
                    end if
                end do
            end if

            if (reach_link(head) /= first_column) then
                reach_link(first_column) = reach_link(head)
                reach_link(head) = first_column
                k_nnz = k_nnz + 1
            end if

            nz_begin = nz_end + 1
            nz_end = nz_end + k_nnz
            if (nz_end + 1 /= xlindx(k_super + 1)) then
                flag = -2
                return
            end if
            i = head
            do k_pointer = nz_begin, nz_end
                i = reach_link(i)
                lindx(k_pointer) = i
            end do

            if (length > width) then
                parent_column = lindx(xlindx(k_super) + width)
                parent_super = supernode(parent_column)
                merge_link(k_super) = merge_link(parent_super)
                merge_link(parent_super) = k_super
            end if
        end do
    end subroutine symfc2_core

end module sparsem_ng_symbolic
