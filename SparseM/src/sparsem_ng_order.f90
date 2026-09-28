! Multiple-minimum-degree ordering for the SparseM translation.
!
! The MMD core is a modern free-form Fortran translation of the ordering
! routines distributed with SparseM's Ng-Peyton sparse Cholesky code.
! Algorithm/source credits preserved from upstream include Joseph W. H. Liu,
! Esmond G. Ng, Barry W. Peyton, the University of Waterloo SPARSPAK-A
! release (January 1984), and the Oak Ridge National Laboratory Mathematical
! Sciences Section. See NOTICE.md and LICENSE for package-level provenance.
!
module sparsem_ng_order
    use sparsem_kinds, only: dp
    use sparsem_types, only: csr_matrix
    implicit none
    private

    public :: minimum_degree_order
    public :: csr_symmetric_permute
    public :: csr_symmetric_adjacency

contains

    subroutine minimum_degree_order(a, perm, inv_perm, nofsub, info)
        type(csr_matrix), intent(in) :: a !! Square sparse matrix whose undirected structural graph is ordered.
        integer, allocatable, intent(out) :: perm(:) !! New-to-old permutation; perm(k) is the original index at position k.
        integer, allocatable, intent(out) :: inv_perm(:) !! Old-to-new inverse permutation satisfying inv_perm(perm(k)) = k.
        integer, intent(out), optional :: nofsub !! Ng-Peyton estimate of factor subscripts accumulated by multiple minimum degree.
        integer, intent(out), optional :: info !! Zero on success; -1 for nonsquare input; -2 for invalid structural column indices.
        integer, allocatable :: xadj(:)
        integer, allocatable :: adjacency(:)
        integer, allocatable :: work(:)
        integer :: stat
        integer :: nsub

        if (present(info)) info = 0
        if (present(nofsub)) nofsub = 0
        if (a%nrow /= a%ncol) then
            allocate(perm(0), inv_perm(0))
            if (present(info)) info = -1
            return
        end if
        if (a%nrow == 0) then
            allocate(perm(0), inv_perm(0))
            return
        end if

        call csr_symmetric_adjacency(a, xadj, adjacency, stat)
        if (stat /= 0) then
            allocate(perm(0), inv_perm(0))
            if (present(info)) info = stat
            return
        end if

        allocate(perm(a%nrow), inv_perm(a%nrow), work(4 * a%nrow))
        call ordmmd_core(a%nrow, xadj, adjacency, inv_perm, perm, work, nsub, stat)
        if (present(nofsub)) nofsub = nsub
        if (present(info)) info = stat
    end subroutine minimum_degree_order

    function csr_symmetric_permute(a, perm, inv_perm, info) result(b)
        type(csr_matrix), intent(in) :: a !! Square CSR matrix to permute symmetrically as A(perm,perm).
        integer, intent(in) :: perm(:) !! New-to-old row/column permutation of length n.
        integer, intent(in) :: inv_perm(:) !! Old-to-new inverse permutation of length n.
        integer, intent(out), optional :: info !! Zero on success; -1 for shape/permutation mismatch; -2 for invalid indices.
        type(csr_matrix) :: b
        integer :: i
        integer :: old_col
        integer :: old_row
        integer :: pos
        integer :: start_pos
        integer :: stop_pos

        if (present(info)) info = 0
        if (a%nrow /= a%ncol .or. size(perm) /= a%nrow .or. size(inv_perm) /= a%nrow) then
            call initialize_empty_csr(b)
            if (present(info)) info = -1
            return
        end if
        if (.not. valid_permutation(perm, inv_perm)) then
            call initialize_empty_csr(b)
            if (present(info)) info = -2
            return
        end if

        b%nrow = a%nrow
        b%ncol = a%ncol
        allocate(b%values(size(a%values)))
        allocate(b%col_ind(size(a%col_ind)))
        allocate(b%row_ptr(a%nrow + 1))
        b%row_ptr(1) = 1
        pos = 1
        do i = 1, a%nrow
            old_row = perm(i)
            start_pos = pos
            do old_col = a%row_ptr(old_row), a%row_ptr(old_row + 1) - 1
                if (a%col_ind(old_col) < 1 .or. a%col_ind(old_col) > a%ncol) then
                    call initialize_empty_csr(b)
                    if (present(info)) info = -2
                    return
                end if
                b%values(pos) = a%values(old_col)
                b%col_ind(pos) = inv_perm(a%col_ind(old_col))
                pos = pos + 1
            end do
            stop_pos = pos - 1
            call sort_csr_row(b%col_ind, b%values, start_pos, stop_pos)
            b%row_ptr(i + 1) = pos
        end do
    end function csr_symmetric_permute

    subroutine csr_symmetric_adjacency(a, xadj, adjacency, info)
        type(csr_matrix), intent(in) :: a !! Square CSR matrix providing the graph structure; diagonal entries are ignored.
        integer, allocatable, intent(out) :: xadj(:) !! One-based row pointers for the undirected adjacency structure.
        integer, allocatable, intent(out) :: adjacency(:) !! Unique one-based undirected neighbors, grouped by row.
        integer, intent(out) :: info !! Zero on success or -2 when a stored column index is out of range.
        integer, allocatable :: candidate_ptr(:)
        integer, allocatable :: candidates(:)
        integer, allocatable :: counts(:)
        integer, allocatable :: marker(:)
        integer, allocatable :: next_pos(:)
        integer :: candidate_nnz
        integer :: col
        integer :: i
        integer :: k
        integer :: pos
        integer :: total

        info = 0
        allocate(counts(a%nrow), source=0)
        candidate_nnz = 0
        do i = 1, a%nrow
            do pos = a%row_ptr(i), a%row_ptr(i + 1) - 1
                col = a%col_ind(pos)
                if (col < 1 .or. col > a%ncol) then
                    allocate(xadj(1), adjacency(0))
                    xadj(1) = 1
                    info = -2
                    return
                end if
                if (col /= i) then
                    counts(i) = counts(i) + 1
                    counts(col) = counts(col) + 1
                    candidate_nnz = candidate_nnz + 2
                end if
            end do
        end do

        allocate(candidate_ptr(a%nrow + 1))
        candidate_ptr(1) = 1
        do i = 1, a%nrow
            candidate_ptr(i + 1) = candidate_ptr(i) + counts(i)
        end do
        allocate(candidates(candidate_nnz))
        allocate(next_pos(a%nrow))
        next_pos = candidate_ptr(1:a%nrow)
        do i = 1, a%nrow
            do pos = a%row_ptr(i), a%row_ptr(i + 1) - 1
                col = a%col_ind(pos)
                if (col == i) cycle
                candidates(next_pos(i)) = col
                next_pos(i) = next_pos(i) + 1
                candidates(next_pos(col)) = i
                next_pos(col) = next_pos(col) + 1
            end do
        end do

        counts = 0
        allocate(marker(a%nrow), source=0)
        do i = 1, a%nrow
            do k = candidate_ptr(i), candidate_ptr(i + 1) - 1
                col = candidates(k)
                if (marker(col) /= i) then
                    marker(col) = i
                    counts(i) = counts(i) + 1
                end if
            end do
        end do

        allocate(xadj(a%nrow + 1))
        xadj(1) = 1
        do i = 1, a%nrow
            xadj(i + 1) = xadj(i) + counts(i)
        end do
        total = xadj(a%nrow + 1) - 1
        allocate(adjacency(total))
        marker = 0
        next_pos = xadj(1:a%nrow)
        do i = 1, a%nrow
            do k = candidate_ptr(i), candidate_ptr(i + 1) - 1
                col = candidates(k)
                if (marker(col) /= i) then
                    marker(col) = i
                    adjacency(next_pos(i)) = col
                    next_pos(i) = next_pos(i) + 1
                end if
            end do
            call sort_integer_segment(adjacency, xadj(i), xadj(i + 1) - 1)
        end do
    end subroutine csr_symmetric_adjacency

    pure logical function valid_permutation(perm, inv_perm) result(valid)
        integer, intent(in) :: perm(:) !! Candidate new-to-old permutation.
        integer, intent(in) :: inv_perm(:) !! Candidate old-to-new inverse permutation.
        integer :: i

        valid = size(perm) == size(inv_perm)
        if (.not. valid) return
        do i = 1, size(perm)
            if (perm(i) < 1 .or. perm(i) > size(perm)) then
                valid = .false.
                return
            end if
            if (inv_perm(perm(i)) /= i) then
                valid = .false.
                return
            end if
        end do
    end function valid_permutation

    pure subroutine sort_integer_segment(values, first, last)
        integer, intent(inout) :: values(:) !! Integer array containing the segment to sort in ascending order.
        integer, intent(in) :: first !! First one-based position of the segment.
        integer, intent(in) :: last !! Last one-based position of the segment; may be less than first for an empty segment.
        integer :: i
        integer :: j
        integer :: key

        do i = first + 1, last
            key = values(i)
            j = i - 1
            do while (j >= first)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_integer_segment

    pure subroutine sort_csr_row(columns, values, first, last)
        integer, intent(inout) :: columns(:) !! CSR column indices whose selected row segment is sorted.
        real(dp), intent(inout) :: values(:) !! CSR numerical values moved with their column indices.
        integer, intent(in) :: first !! First one-based position of the CSR row segment.
        integer, intent(in) :: last !! Last one-based position of the CSR row segment; may be less than first when empty.
        integer :: i
        integer :: j
        integer :: key_column
        real(dp) :: key_value

        do i = first + 1, last
            key_column = columns(i)
            key_value = values(i)
            j = i - 1
            do while (j >= first)
                if (columns(j) <= key_column) exit
                columns(j + 1) = columns(j)
                values(j + 1) = values(j)
                j = j - 1
            end do
            columns(j + 1) = key_column
            values(j + 1) = key_value
        end do
    end subroutine sort_csr_row

    subroutine initialize_empty_csr(a)
        type(csr_matrix), intent(out) :: a !! CSR object initialized to an allocated zero-by-zero matrix.

        a%nrow = 0
        a%ncol = 0
        allocate(a%values(0), a%col_ind(0), a%row_ptr(1))
        a%row_ptr(1) = 1
    end subroutine initialize_empty_csr

    subroutine ordmmd_core(neqns, xadj, adjacency, inv_perm, perm, work, nofsub, info)
        integer, intent(in) :: neqns !! Number of graph vertices/equations to order.
        integer, intent(in) :: xadj(:) !! One-based adjacency pointers of length neqns+1.
        integer, intent(inout) :: adjacency(:) !! Adjacency workspace; quotient-graph updates overwrite entries.
        integer, intent(out) :: inv_perm(:) !! Old-to-new inverse ordering produced by multiple minimum degree.
        integer, intent(out) :: perm(:) !! New-to-old ordering produced by multiple minimum degree.
        integer, intent(inout) :: work(:) !! Integer workspace of length at least 4*neqns.
        integer, intent(out) :: nofsub !! Estimated number of factor subscripts accumulated during ordering.
        integer, intent(out) :: info !! Zero on success or -1 when work is too small.
        integer :: delta
        integer :: max_integer

        info = 0
        nofsub = 0
        if (size(work) < 4 * neqns) then
            info = -1
            return
        end if
        if (neqns <= 0) return

        delta = 0
        max_integer = 32767
        call genmmd_core(neqns, xadj, adjacency, inv_perm, perm, delta, &
            work(1:neqns), work(neqns + 1:2 * neqns), &
            work(2 * neqns + 1:3 * neqns), work(3 * neqns + 1:4 * neqns), &
            max_integer, nofsub)
    end subroutine ordmmd_core

    subroutine genmmd_core(neqns, xadj, adjacency, inv_perm, perm, delta, dhead, qsize, llist, marker, &
            max_integer, nofsub)
        integer, intent(in) :: neqns !! Number of equations/vertices in the quotient graph.
        integer, intent(in) :: xadj(:) !! One-based adjacency pointers for the graph.
        integer, intent(inout) :: adjacency(:) !! Quotient-graph adjacency storage overwritten during elimination.
        integer, intent(inout) :: inv_perm(:) !! Inverse-order workspace, finalized as old-to-new numbering.
        integer, intent(inout) :: perm(:) !! Permutation/back-link workspace, finalized as new-to-old numbering.
        integer, intent(in) :: delta !! Multiple-elimination degree tolerance; SparseM uses zero.
        integer, intent(inout) :: dhead(:) !! Heads of the linked degree lists.
        integer, intent(inout) :: qsize(:) !! Current quotient/supernode sizes.
        integer, intent(inout) :: llist(:) !! Linked-list workspace for eliminated quotient elements.
        integer, intent(inout) :: marker(:) !! Marker tags used to avoid repeated graph scans.
        integer, intent(in) :: max_integer !! Large positive sentinel reserved for permanent marks.
        integer, intent(out) :: nofsub !! Accumulated estimate of factor subscripts.
        integer :: ehead
        integer :: i
        integer :: mdeg
        integer :: mdlmt
        integer :: mdnode
        integer :: nextmd
        integer :: num
        integer :: tag

        if (neqns <= 0) return
        nofsub = 0
        call mmdint_core(neqns, xadj, dhead, inv_perm, perm, qsize, llist, marker)
        num = 1
        nextmd = dhead(1)

100     continue
        if (nextmd <= 0) go to 200
        mdnode = nextmd
        nextmd = inv_perm(mdnode)
        marker(mdnode) = max_integer
        inv_perm(mdnode) = -num
        num = num + 1
        go to 100

200     continue
        if (num > neqns) go to 1000
        tag = 1
        dhead(1) = 0
        mdeg = 2

300     continue
        if (dhead(mdeg) > 0) go to 400
        mdeg = mdeg + 1
        go to 300

400     continue
        mdlmt = mdeg + delta
        ehead = 0

500     continue
        mdnode = dhead(mdeg)
        if (mdnode > 0) go to 600
        mdeg = mdeg + 1
        if (mdeg > mdlmt) go to 900
        go to 500

600     continue
        nextmd = inv_perm(mdnode)
        dhead(mdeg) = nextmd
        if (nextmd > 0) perm(nextmd) = -mdeg
        inv_perm(mdnode) = -num
        nofsub = nofsub + mdeg + qsize(mdnode) - 2
        if (num + qsize(mdnode) > neqns) go to 1000

        tag = tag + 1
        if (tag < max_integer) go to 800
        tag = 1
        do i = 1, neqns
            if (marker(i) < max_integer) marker(i) = 0
        end do

800     continue
        call mmdelm_core(mdnode, xadj, adjacency, dhead, inv_perm, perm, qsize, llist, marker, &
            max_integer, tag)
        num = num + qsize(mdnode)
        llist(mdnode) = ehead
        ehead = mdnode
        if (delta >= 0) go to 500

900     continue
        if (num > neqns) go to 1000
        call mmdupd_core(ehead, neqns, xadj, adjacency, delta, mdeg, dhead, inv_perm, perm, qsize, &
            llist, marker, max_integer, tag)
        go to 300

1000    continue
        call mmdnum_core(neqns, perm, inv_perm, qsize)
    end subroutine genmmd_core

    subroutine mmdint_core(neqns, xadj, dhead, dforw, dbakw, qsize, llist, marker)
        integer, intent(in) :: neqns !! Number of equations/vertices in the graph.
        integer, intent(in) :: xadj(:) !! One-based adjacency pointers used to compute initial degrees.
        integer, intent(out) :: dhead(:) !! Heads of degree lists, initialized for every possible degree.
        integer, intent(out) :: dforw(:) !! Forward links in degree lists.
        integer, intent(out) :: dbakw(:) !! Backward links or negative degree-list heads.
        integer, intent(out) :: qsize(:) !! Quotient-node sizes, initialized to one.
        integer, intent(out) :: llist(:) !! Quotient-element linked-list workspace, initialized to zero.
        integer, intent(out) :: marker(:) !! Marker workspace, initialized to zero.
        integer :: fnode
        integer :: ndeg
        integer :: node

        do node = 1, neqns
            dhead(node) = 0
            qsize(node) = 1
            marker(node) = 0
            llist(node) = 0
        end do
        do node = 1, neqns
            ndeg = xadj(node + 1) - xadj(node) + 1
            fnode = dhead(ndeg)
            dforw(node) = fnode
            dhead(ndeg) = node
            if (fnode > 0) dbakw(fnode) = node
            dbakw(node) = -ndeg
        end do
    end subroutine mmdint_core

    subroutine mmdnum_core(neqns, perm, inv_perm, qsize)
        integer, intent(in) :: neqns !! Number of equations/vertices being numbered.
        integer, intent(inout) :: perm(:) !! Merge-tree workspace converted to the final new-to-old permutation.
        integer, intent(inout) :: inv_perm(:) !! Elimination labels converted to the final old-to-new permutation.
        integer, intent(in) :: qsize(:) !! Quotient/supernode sizes at elimination; zero marks merged nodes.
        integer :: father
        integer :: next_father
        integer :: node
        integer :: node_qsize
        integer :: num
        integer :: root

        do node = 1, neqns
            node_qsize = qsize(node)
            if (node_qsize <= 0) perm(node) = inv_perm(node)
            if (node_qsize > 0) perm(node) = -inv_perm(node)
        end do

        do node = 1, neqns
            if (perm(node) > 0) cycle
            father = node
200         continue
            if (perm(father) > 0) go to 300
            father = -perm(father)
            go to 200
300         continue
            root = father
            num = perm(root) + 1
            inv_perm(node) = -num
            perm(root) = num
            father = node
400         continue
            next_father = -perm(father)
            if (next_father <= 0) cycle
            perm(father) = -root
            father = next_father
            go to 400
        end do

        do node = 1, neqns
            num = -inv_perm(node)
            inv_perm(node) = num
            perm(num) = node
        end do
    end subroutine mmdnum_core

    subroutine mmdelm_core(mdnode, xadj, adjacency, dhead, dforw, dbakw, qsize, llist, marker, &
            max_integer, tag)
        integer, intent(in) :: mdnode !! Current minimum-degree quotient node to eliminate.
        integer, intent(in) :: xadj(:) !! One-based adjacency pointers for quotient-graph storage.
        integer, intent(inout) :: adjacency(:) !! Quotient-graph adjacency storage updated in place.
        integer, intent(inout) :: dhead(:) !! Heads of degree linked lists updated as reachable nodes are removed.
        integer, intent(inout) :: dforw(:) !! Forward degree-list links and quotient-node state codes.
        integer, intent(inout) :: dbakw(:) !! Backward degree-list links and quotient-node state codes.
        integer, intent(inout) :: qsize(:) !! Quotient/supernode sizes updated when indistinguishable nodes merge.
        integer, intent(inout) :: llist(:) !! Linked list of eliminated quotient elements.
        integer, intent(inout) :: marker(:) !! Marker tags used for reachable-set construction and purging.
        integer, intent(in) :: max_integer !! Large positive sentinel used to permanently mark merged nodes.
        integer, intent(in) :: tag !! Current marker generation for this elimination step.
        integer :: element
        integer :: i
        integer :: istop
        integer :: istart
        integer :: j
        integer :: jstop
        integer :: jstart
        integer :: link
        integer :: neighbor
        integer :: node
        integer :: npv
        integer :: nq_neighbors
        integer :: next_node
        integer :: previous_node
        integer :: reachable_limit
        integer :: reachable_location
        integer :: reachable_node
        integer :: quotient_neighbor

        marker(mdnode) = tag
        istart = xadj(mdnode)
        istop = xadj(mdnode + 1) - 1
        element = 0
        reachable_location = istart
        reachable_limit = istop

        do i = istart, istop
            neighbor = adjacency(i)
            if (neighbor == 0) go to 300
            if (marker(neighbor) >= tag) cycle
            marker(neighbor) = tag
            if (dforw(neighbor) < 0) go to 100
            adjacency(reachable_location) = neighbor
            reachable_location = reachable_location + 1
            cycle
100         continue
            llist(neighbor) = element
            element = neighbor
        end do

300     continue
        if (element <= 0) go to 1000
        adjacency(reachable_limit) = -element
        link = element
400     continue
        jstart = xadj(link)
        jstop = xadj(link + 1) - 1
        do j = jstart, jstop
            node = adjacency(j)
            link = -node
            if (node < 0) then
                go to 400
            else if (node == 0) then
                go to 900
            end if
            if (marker(node) >= tag .or. dforw(node) < 0) cycle
            marker(node) = tag
600         continue
            if (reachable_location < reachable_limit) go to 700
            link = -adjacency(reachable_limit)
            reachable_location = xadj(link)
            reachable_limit = xadj(link + 1) - 1
            go to 600
700         continue
            adjacency(reachable_location) = node
            reachable_location = reachable_location + 1
        end do
900     continue
        element = llist(element)
        go to 300

1000    continue
        if (reachable_location <= reachable_limit) adjacency(reachable_location) = 0
        link = mdnode
1100    continue
        istart = xadj(link)
        istop = xadj(link + 1) - 1
        do i = istart, istop
            reachable_node = adjacency(i)
            link = -reachable_node
            if (reachable_node < 0) then
                go to 1100
            else if (reachable_node == 0) then
                go to 1800
            end if

            previous_node = dbakw(reachable_node)
            if (previous_node /= 0 .and. previous_node /= -max_integer) then
                next_node = dforw(reachable_node)
                if (next_node > 0) dbakw(next_node) = previous_node
                if (previous_node > 0) dforw(previous_node) = next_node
                npv = -previous_node
                if (previous_node < 0) dhead(npv) = next_node
            end if

            jstart = xadj(reachable_node)
            jstop = xadj(reachable_node + 1) - 1
            quotient_neighbor = jstart
            do j = jstart, jstop
                neighbor = adjacency(j)
                if (neighbor == 0) go to 1500
                if (marker(neighbor) >= tag) cycle
                adjacency(quotient_neighbor) = neighbor
                quotient_neighbor = quotient_neighbor + 1
            end do
1500        continue
            nq_neighbors = quotient_neighbor - jstart
            if (nq_neighbors > 0) go to 1600

            qsize(mdnode) = qsize(mdnode) + qsize(reachable_node)
            qsize(reachable_node) = 0
            marker(reachable_node) = max_integer
            dforw(reachable_node) = -mdnode
            dbakw(reachable_node) = -max_integer
            cycle

1600        continue
            dforw(reachable_node) = nq_neighbors + 1
            dbakw(reachable_node) = 0
            adjacency(quotient_neighbor) = mdnode
            quotient_neighbor = quotient_neighbor + 1
            if (quotient_neighbor <= jstop) adjacency(quotient_neighbor) = 0
        end do
1800    continue
    end subroutine mmdelm_core

    subroutine mmdupd_core(ehead, neqns, xadj, adjacency, delta, mdeg, dhead, dforw, dbakw, qsize, &
            llist, marker, max_integer, tag)
        integer, intent(in) :: ehead !! Head of the linked list of quotient elements formed in the elimination step.
        integer, intent(in) :: neqns !! Number of equations/vertices in the quotient graph.
        integer, intent(in) :: xadj(:) !! One-based adjacency pointers for quotient-graph storage.
        integer, intent(in) :: adjacency(:) !! Quotient-graph adjacency storage created by the elimination routine.
        integer, intent(in) :: delta !! Multiple-elimination degree tolerance; SparseM uses zero.
        integer, intent(inout) :: mdeg !! Current minimum degree, reduced when an updated node has a smaller degree.
        integer, intent(inout) :: dhead(:) !! Heads of the degree linked lists.
        integer, intent(inout) :: dforw(:) !! Forward degree-list links and quotient-node state codes.
        integer, intent(inout) :: dbakw(:) !! Backward degree-list links and quotient-node state codes.
        integer, intent(inout) :: qsize(:) !! Quotient/supernode sizes, including indistinguishable-node merges.
        integer, intent(inout) :: llist(:) !! Linked lists for quotient elements and nodes awaiting degree updates.
        integer, intent(inout) :: marker(:) !! Marker tags used while counting external degrees.
        integer, intent(in) :: max_integer !! Large positive sentinel reserved for permanently marked merged nodes.
        integer, intent(inout) :: tag !! Marker generation advanced during degree updates.
        integer :: degree
        integer :: degree0
        integer :: element
        integer :: enode
        integer :: first_node
        integer :: i
        integer :: iq2
        integer :: istop
        integer :: istart
        integer :: j
        integer :: jstop
        integer :: jstart
        integer :: link
        integer :: mdeg0
        integer :: mtag
        integer :: neighbor
        integer :: node
        integer :: q2head
        integer :: qxhead

        mdeg0 = mdeg + delta
        element = ehead

100     continue
        if (element <= 0) return
        mtag = tag + mdeg0
        if (mtag < max_integer) go to 300
        tag = 1
        do i = 1, neqns
            if (marker(i) < max_integer) marker(i) = 0
        end do
        mtag = tag + mdeg0

300     continue
        q2head = 0
        qxhead = 0
        degree0 = 0
        link = element

400     continue
        istart = xadj(link)
        istop = xadj(link + 1) - 1
        do i = istart, istop
            enode = adjacency(i)
            link = -enode
            if (enode < 0) then
                go to 400
            else if (enode == 0) then
                go to 800
            end if

            if (qsize(enode) == 0) cycle
            degree0 = degree0 + qsize(enode)
            marker(enode) = mtag
            if (dbakw(enode) /= 0) cycle
            if (dforw(enode) == 2) go to 600
            llist(enode) = qxhead
            qxhead = enode
            cycle
600         continue
            llist(enode) = q2head
            q2head = enode
        end do

800     continue
        enode = q2head
        iq2 = 1

900     continue
        if (enode <= 0) go to 1500
        if (dbakw(enode) /= 0) go to 2200
        tag = tag + 1
        degree = degree0
        istart = xadj(enode)
        neighbor = adjacency(istart)
        if (neighbor == element) neighbor = adjacency(istart + 1)
        link = neighbor
        if (dforw(neighbor) < 0) go to 1000
        degree = degree + qsize(neighbor)
        go to 2100

1000    continue
        istart = xadj(link)
        istop = xadj(link + 1) - 1
        do i = istart, istop
            node = adjacency(i)
            link = -node
            if (node == enode) cycle
            if (node < 0) then
                go to 1000
            else if (node == 0) then
                go to 2100
            end if

            if (qsize(node) == 0) cycle
            if (marker(node) >= tag) go to 1200
            marker(node) = tag
            degree = degree + qsize(node)
            cycle

1200        continue
            if (dbakw(node) /= 0) cycle
            if (dforw(node) /= 2) go to 1300
            qsize(enode) = qsize(enode) + qsize(node)
            qsize(node) = 0
            marker(node) = max_integer
            dforw(node) = -enode
            dbakw(node) = -max_integer
            cycle

1300        continue
            if (dbakw(node) == 0) dbakw(node) = -max_integer
        end do
        go to 2100

1500    continue
        enode = qxhead
        iq2 = 0

1600    continue
        if (enode <= 0) go to 2300
        if (dbakw(enode) /= 0) go to 2200
        tag = tag + 1
        degree = degree0
        istart = xadj(enode)
        istop = xadj(enode + 1) - 1
        do i = istart, istop
            neighbor = adjacency(i)
            if (neighbor == 0) go to 2100
            if (marker(neighbor) >= tag) cycle
            marker(neighbor) = tag
            link = neighbor
            if (dforw(neighbor) < 0) go to 1700
            degree = degree + qsize(neighbor)
            cycle

1700        continue
            jstart = xadj(link)
            jstop = xadj(link + 1) - 1
            do j = jstart, jstop
                node = adjacency(j)
                link = -node
                if (node < 0) then
                    go to 1700
                else if (node == 0) then
                    go to 2000
                end if
                if (marker(node) >= tag) cycle
                marker(node) = tag
                degree = degree + qsize(node)
            end do
2000        continue
        end do

2100    continue
        degree = degree - qsize(enode) + 1
        first_node = dhead(degree)
        dforw(enode) = first_node
        dbakw(enode) = -degree
        if (first_node > 0) dbakw(first_node) = enode
        dhead(degree) = enode
        if (degree < mdeg) mdeg = degree

2200    continue
        enode = llist(enode)
        if (iq2 == 1) go to 900
        go to 1600

2300    continue
        tag = mtag
        element = llist(element)
        go to 100
    end subroutine mmdupd_core

end module sparsem_ng_order
