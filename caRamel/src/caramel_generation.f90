module caramel_generation
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use caramel_kinds, only: dp
    use caramel_delaunay, only: delaunay_nd
    use caramel_linalg, only: cholesky_lower
    use caramel_pareto, only: dominate
    use caramel_random, only: random_normal, random_uniform, random_index
    use caramel_utils, only: rselect, val2rank, matvcov, vol_splx, dimprove, column_mean, column_sd
    implicit none
    private

    type, public :: index_block
        integer, allocatable :: index(:)
    end type index_block

    public :: cinterp, cextrap, crecombination, cusecovar, new_xval

contains

    subroutine cinterp(param, crit, simplices, volume, nnew, xnew, pcrit)
        real(dp), intent(in) :: param(:,:), crit(:,:)
        integer, intent(in) :: simplices(:,:)
        real(dp), intent(in) :: volume(:)
        integer, intent(in) :: nnew
        real(dp), allocatable, intent(out) :: xnew(:,:), pcrit(:,:)
        integer, allocatable :: choice(:)
        real(dp), allocatable :: bary(:), weights(:)
        integer :: i, j, s, nv

        if (size(volume) /= size(simplices,1)) error stop "cinterp: inconsistent volume dimension"
        allocate(xnew(max(0,nnew),size(param,2)), pcrit(max(0,nnew),size(crit,2)))
        if (nnew <= 0) return
        if (size(simplices,1) == 0) then
            xnew = ieee_value(0.0_dp, ieee_quiet_nan)
            pcrit = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        weights = max(volume, 0.0_dp)
        if (sum(weights) <= 0.0_dp) weights = 1.0_dp
        call rselect(nnew, weights, choice)
        nv = size(simplices,2)
        allocate(bary(nv))
        do i = 1, nnew
            s = choice(i)
            do j = 1, nv
                bary(j) = random_uniform()
            end do
            if (sum(bary) <= 0.0_dp) bary = 1.0_dp
            bary = bary / sum(bary)
            xnew(i,:) = matmul(bary, param(simplices(s,:),:))
            pcrit(i,:) = matmul(bary, crit(simplices(s,:),:))
        end do
    end subroutine cinterp

    subroutine cextrap(param, crit, directions, lengths, nnew, xnew, pcrit)
        real(dp), intent(in) :: param(:,:), crit(:,:)
        integer, intent(in) :: directions(:,:)
        real(dp), intent(in) :: lengths(:)
        integer, intent(in) :: nnew
        real(dp), allocatable, intent(out) :: xnew(:,:), pcrit(:,:)
        integer, allocatable :: choice(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: boost, lambda, mean_length, u
        integer :: i, e, a, b

        if (size(directions,2) /= 2 .or. size(lengths) /= size(directions,1)) then
            error stop "cextrap: inconsistent edge dimensions"
        end if
        allocate(xnew(max(0,nnew),size(param,2)), pcrit(max(0,nnew),size(crit,2)))
        if (nnew <= 0) return
        if (size(directions,1) == 0) then
            xnew = ieee_value(0.0_dp, ieee_quiet_nan)
            pcrit = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        weights = max(lengths, 0.0_dp)
        if (sum(weights) <= 0.0_dp) weights = 1.0_dp
        call rselect(nnew, weights, choice)
        mean_length = sum(lengths) / real(size(lengths),dp)
        do i = 1, nnew
            e = choice(i)
            a = directions(e,1)
            b = directions(e,2)
            if (lengths(e) > tiny(1.0_dp)) then
                boost = mean_length / lengths(e)
            else
                boost = 1.0_dp
            end if
            u = min(random_uniform(), 1.0_dp - epsilon(1.0_dp))
            lambda = -log(max(tiny(1.0_dp), 1.0_dp - u))
            xnew(i,:) = param(b,:) + boost * lambda * (param(b,:) - param(a,:))
            pcrit(i,:) = crit(b,:) + boost * lambda * (crit(b,:) - crit(a,:))
        end do
    end subroutine cextrap

    subroutine crecombination(param, blocks, nnew, xnew)
        real(dp), intent(in) :: param(:,:)
        type(index_block), intent(in), optional :: blocks(:)
        integer, intent(in) :: nnew
        real(dp), allocatable, intent(out) :: xnew(:,:)
        type(index_block), allocatable :: all_blocks(:)
        logical, allocatable :: covered(:)
        integer :: npar, nblock, extra, i, j, b, parent, pos

        npar = size(param,2)
        allocate(xnew(max(0,nnew),npar))
        if (nnew <= 0) return
        if (size(param,1) == 0) then
            xnew = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        allocate(covered(npar))
        covered = .false.
        nblock = 0
        if (present(blocks)) nblock = size(blocks)
        do b = 1, nblock
            do j = 1, size(blocks(b)%index)
                if (blocks(b)%index(j) < 1 .or. blocks(b)%index(j) > npar) then
                    error stop "crecombination: block index out of range"
                end if
                covered(blocks(b)%index(j)) = .true.
            end do
        end do
        extra = count(.not. covered)
        allocate(all_blocks(nblock+extra))
        do b = 1, nblock
            all_blocks(b)%index = blocks(b)%index
        end do
        pos = nblock
        do j = 1, npar
            if (.not. covered(j)) then
                pos = pos + 1
                allocate(all_blocks(pos)%index(1))
                all_blocks(pos)%index(1) = j
            end if
        end do

        do i = 1, nnew
            xnew(i,:) = 0.0_dp
            do b = 1, size(all_blocks)
                parent = random_index(size(param,1))
                xnew(i,all_blocks(b)%index) = param(parent,all_blocks(b)%index)
            end do
        end do
    end subroutine crecombination

    subroutine cusecovar(xref, amplif, nnew, xnew)
        real(dp), intent(in) :: xref(:,:)
        real(dp), intent(in) :: amplif
        integer, intent(in) :: nnew
        real(dp), allocatable, intent(out) :: xnew(:,:)
        real(dp), allocatable :: rr(:,:), cov(:,:), l(:,:), sd(:), g(:), z(:), noise(:)
        integer, allocatable :: pvar(:)
        logical :: ok
        integer :: p, nv, i, j, k, attempt
        real(dp) :: jitter

        p = size(xref,2)
        allocate(xnew(max(0,nnew),p))
        if (nnew <= 0) return
        if (size(xref,1) == 0) then
            xnew = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        g = column_mean(xref)
        sd = column_sd(xref)
        pvar = pack([(j,j=1,p)], sd > 0.0_dp)
        xnew = spread(g,1,nnew)
        nv = size(pvar)
        if (nv == 0) return

        allocate(rr(p,p))
        call matvcov(xref, g, rr)
        allocate(cov(nv,nv), l(nv,nv), z(nv), noise(nv))
        do i = 1, nv
            do j = 1, nv
                cov(i,j) = amplif * sd(pvar(i)) * rr(pvar(i),pvar(j)) * amplif * sd(pvar(j))
            end do
        end do
        jitter = 0.0_dp
        do attempt = 1, 8
            call cholesky_lower(cov, l, ok)
            if (ok) exit
            if (jitter <= 0.0_dp) then
                jitter = 1.0e-12_dp * max(1.0_dp, maxval(abs(cov)))
            else
                jitter = 10.0_dp * jitter
            end if
            do k = 1, nv
                cov(k,k) = cov(k,k) + jitter
            end do
        end do
        if (.not. ok) then
            l = 0.0_dp
            do k = 1, nv
                l(k,k) = amplif * sd(pvar(k))
            end do
        end if

        do i = 1, nnew
            do j = 1, nv
                z(j) = random_normal()
            end do
            noise = matmul(l,z)
            do j = 1, nv
                xnew(i,pvar(j)) = g(pvar(j)) + noise(j)
            end do
        end do
    end subroutine cusecovar

    subroutine new_xval(param, crit, isperf, sp, bounds, repart_gene, fireworks, xnew, project_crit, blocks)
        real(dp), intent(in) :: param(:,:), crit(:,:), sp(:), bounds(:,:)
        logical, intent(in) :: isperf(:), fireworks
        integer, intent(in) :: repart_gene(4)
        real(dp), allocatable, intent(out) :: xnew(:,:), project_crit(:,:)
        type(index_block), intent(in), optional :: blocks(:)
        real(dp), allocatable :: obj(:,:), obj_unique(:,:), ranks(:), volume(:), edge_length(:)
        real(dp), allocatable :: xi(:,:), ci(:,:), xe(:,:), ce(:,:), xc(:,:), xr(:,:)
        real(dp), allocatable :: rr(:), obj_arch(:,:), crit_arch(:,:), xref(:,:)
        integer, allocatable :: front_rank(:), archive_idx(:), unique_map(:), simplices(:,:), kept_simp(:,:), orig_simp(:,:)
        integer, allocatable :: local_edges(:,:), edges(:,:), tmp_edges(:,:), reference_idx(:), ipp(:)
        logical :: delaunay_ok
        integer :: nvec, nobj, npar, n_inter, n_extra, n_cov, n_recomb
        integer :: i, j, s, row, total_cap, nfire, maximin, idx
        real(dp), parameter :: a = 3.0_dp / 8.0_dp
        real(dp) :: nanv, best

        nvec = size(param,1)
        npar = size(param,2)
        nobj = size(crit,2)
        if (size(crit,1) /= nvec .or. size(isperf) /= nobj .or. size(sp) /= npar) then
            error stop "new_xval: inconsistent dimensions"
        end if
        if (size(bounds,1) /= npar .or. size(bounds,2) /= 2) error stop "new_xval: invalid bounds"
        nanv = ieee_value(0.0_dp, ieee_quiet_nan)

        allocate(obj(nvec,nobj), ranks(nvec), front_rank(nvec))
        do j = 1, nobj
            call val2rank(crit(:,j), 3, ranks)
            obj(:,j) = (ranks - a) / (real(nvec,dp) + 1.0_dp - 2.0_dp*a)
            if (.not. isperf(j)) obj(:,j) = 1.0_dp - obj(:,j)
        end do
        call dominate(obj, front_rank)
        archive_idx = pack([(i,i=1,nvec)], front_rank == 1)
        obj_arch = obj(archive_idx,:)
        crit_arch = crit(archive_idx,:)

        n_inter = max(0,repart_gene(1))
        n_extra = max(0,repart_gene(2))
        n_cov = max(0,repart_gene(3))
        n_recomb = max(0,repart_gene(4))

        call unique_rows_first(obj, obj_unique, unique_map)
        if ((n_inter > 0 .or. n_extra > 0 .or. n_cov > 0) .and. size(obj_unique,1) >= nobj + 1) then
            call delaunay_nd(obj_unique, simplices, delaunay_ok)
        else
            allocate(simplices(0,nobj+1))
            delaunay_ok = .false.
        end if

        if (delaunay_ok) then
            call filter_frontal_simplices(simplices, unique_map, front_rank, kept_simp)
        else
            allocate(kept_simp(0,nobj+1))
        end if

        allocate(orig_simp(size(kept_simp,1),size(kept_simp,2)))
        do s = 1, size(kept_simp,1)
            do j = 1, size(kept_simp,2)
                orig_simp(s,j) = unique_map(kept_simp(s,j))
            end do
        end do

        allocate(volume(size(kept_simp,1)))
        allocate(edges(0,2), edge_length(0))
        do s = 1, size(kept_simp,1)
            volume(s) = vol_splx(obj(unique_map(kept_simp(s,:)),:))
            call dimprove(obj(unique_map(kept_simp(s,:)),:), &
                          front_rank(unique_map(kept_simp(s,:))), local_edges, rr)
            if (size(local_edges,1) > 0) then
                allocate(tmp_edges(size(local_edges,1),2))
                do i = 1, size(local_edges,1)
                    tmp_edges(i,1) = unique_map(kept_simp(s,local_edges(i,1)))
                    tmp_edges(i,2) = unique_map(kept_simp(s,local_edges(i,2)))
                end do
                call append_edges_unique(edges, edge_length, tmp_edges, rr)
                deallocate(tmp_edges)
            end if
        end do

        if (n_inter > 0 .and. size(kept_simp,1) > 0) then
            call cinterp(param, crit, orig_simp, volume, n_inter, xi, ci)
        else
            allocate(xi(0,npar), ci(0,nobj))
        end if
        if (n_extra > 0 .and. size(edges,1) > 0) then
            call cextrap(param, crit, edges, edge_length, n_extra, xe, ce)
        else
            allocate(xe(0,npar), ce(0,nobj))
        end if

        if (size(kept_simp,1) > 0) then
            reference_idx = unique_int(reshape(orig_simp, [size(orig_simp)]))
        else if (size(archive_idx) > 0) then
            reference_idx = archive_idx
        else
            reference_idx = [(i,i=1,nvec)]
        end if
        xref = param(reference_idx,:)
        if (n_cov > 0) then
            call cusecovar(xref, sqrt(2.0_dp), n_cov, xc)
        else
            allocate(xc(0,npar))
        end if

        if (n_recomb > 0 .and. size(archive_idx) > 1) then
            if (present(blocks)) then
                call crecombination(param(archive_idx,:), blocks, n_recomb, xr)
            else
                call crecombination(param(archive_idx,:), nnew=n_recomb, xnew=xr)
            end if
        else
            allocate(xr(0,npar))
        end if

        nfire = 0
        if (fireworks .and. size(archive_idx) > 0) nfire = (nobj + 1) * count(sp > 0.0_dp)
        total_cap = size(xi,1) + size(xe,1) + size(xc,1) + size(xr,1) + nfire
        allocate(xnew(total_cap,npar), project_crit(total_cap,nobj))
        row = 0
        call append_pair(xi, ci, xnew, project_crit, row)
        call append_pair(xe, ce, xnew, project_crit, row)
        if (size(xc,1) > 0) then
            xnew(row+1:row+size(xc,1),:) = xc
            project_crit(row+1:row+size(xc,1),:) = nanv
            row = row + size(xc,1)
        end if
        if (size(xr,1) > 0) then
            xnew(row+1:row+size(xr,1),:) = xr
            project_crit(row+1:row+size(xr,1),:) = nanv
            row = row + size(xr,1)
        end if

        if (fireworks .and. size(archive_idx) > 0) then
            allocate(ipp(nobj+1))
            do j = 1, nobj
                ipp(j) = maxloc(obj_arch(:,j),dim=1)
            end do
            maximin = 1
            best = minval(obj_arch(1,:))
            do i = 2, size(obj_arch,1)
                if (minval(obj_arch(i,:)) > best) then
                    best = minval(obj_arch(i,:))
                    maximin = i
                end if
            end do
            ipp(nobj+1) = maximin
            do i = 1, size(ipp)
                idx = archive_idx(ipp(i))
                do j = 1, npar
                    if (sp(j) <= 0.0_dp) cycle
                    row = row + 1
                    xnew(row,:) = param(idx,:)
                    xnew(row,j) = xnew(row,j) + random_normal() * sp(j)
                    project_crit(row,:) = crit(idx,:)
                end do
            end do
        end if

        if (row < total_cap) then
            xnew = xnew(1:row,:)
            project_crit = project_crit(1:row,:)
        end if
        do i = 1, size(xnew,1)
            do j = 1, npar
                xnew(i,j) = min(bounds(j,2), max(bounds(j,1), xnew(i,j)))
            end do
        end do
    end subroutine new_xval

    subroutine unique_rows_first(x, ux, map)
        real(dp), intent(in) :: x(:,:)
        real(dp), allocatable, intent(out) :: ux(:,:)
        integer, allocatable, intent(out) :: map(:)
        real(dp), allocatable :: tmp(:,:)
        integer, allocatable :: imap(:)
        integer :: i, j, nuniq, d
        logical :: duplicate

        d = size(x,2)
        allocate(tmp(size(x,1),d), imap(size(x,1)))
        nuniq = 0
        do i = 1, size(x,1)
            duplicate = .false.
            do j = 1, nuniq
                if (maxval(abs(x(i,:) - tmp(j,:))) <= 0.0_dp) then
                    duplicate = .true.
                    exit
                end if
            end do
            if (.not. duplicate) then
                nuniq = nuniq + 1
                tmp(nuniq,:) = x(i,:)
                imap(nuniq) = i
            end if
        end do
        allocate(ux(nuniq,d), map(nuniq))
        if (nuniq > 0) then
            ux = tmp(1:nuniq,:)
            map = imap(1:nuniq)
        end if
    end subroutine unique_rows_first

    subroutine filter_frontal_simplices(simplices, map, rank, kept)
        integer, intent(in) :: simplices(:,:), map(:), rank(:)
        integer, allocatable, intent(out) :: kept(:,:)
        integer, allocatable :: tmp(:,:)
        integer :: i, n

        allocate(tmp(size(simplices,1),size(simplices,2)))
        n = 0
        do i = 1, size(simplices,1)
            if (any(rank(map(simplices(i,:))) == 1)) then
                n = n + 1
                tmp(n,:) = simplices(i,:)
            end if
        end do
        allocate(kept(n,size(simplices,2)))
        if (n > 0) kept = tmp(1:n,:)
    end subroutine filter_frontal_simplices

    subroutine append_edges_unique(edges, lengths, new_edges, new_lengths)
        integer, allocatable, intent(inout) :: edges(:,:)
        real(dp), allocatable, intent(inout) :: lengths(:)
        integer, intent(in) :: new_edges(:,:)
        real(dp), intent(in) :: new_lengths(:)
        integer, allocatable :: e2(:,:)
        real(dp), allocatable :: l2(:)
        integer :: i, j, n, add
        logical :: exists

        n = size(edges,1)
        add = 0
        allocate(e2(n+size(new_edges,1),2), l2(n+size(new_edges,1)))
        if (n > 0) then
            e2(1:n,:) = edges
            l2(1:n) = lengths
        end if
        do i = 1, size(new_edges,1)
            exists = .false.
            do j = 1, n + add
                if (all(e2(j,:) == new_edges(i,:))) then
                    exists = .true.
                    exit
                end if
            end do
            if (.not. exists) then
                add = add + 1
                e2(n+add,:) = new_edges(i,:)
                l2(n+add) = new_lengths(i)
            end if
        end do
        deallocate(edges, lengths)
        allocate(edges(n+add,2), lengths(n+add))
        if (n+add > 0) then
            edges = e2(1:n+add,:)
            lengths = l2(1:n+add)
        end if
    end subroutine append_edges_unique

    function unique_int(x) result(u)
        integer, intent(in) :: x(:)
        integer, allocatable :: u(:)
        integer, allocatable :: tmp(:)
        integer :: i, j, n
        logical :: found

        allocate(tmp(size(x)))
        n = 0
        do i = 1, size(x)
            found = .false.
            do j = 1, n
                if (tmp(j) == x(i)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                n = n + 1
                tmp(n) = x(i)
            end if
        end do
        allocate(u(n))
        if (n > 0) u = tmp(1:n)
    end function unique_int

    subroutine append_pair(x, c, xall, call, row)
        real(dp), intent(in) :: x(:,:), c(:,:)
        real(dp), intent(inout) :: xall(:,:), call(:,:)
        integer, intent(inout) :: row
        integer :: n

        n = size(x,1)
        if (n == 0) return
        xall(row+1:row+n,:) = x
        call(row+1:row+n,:) = c
        row = row + n
    end subroutine append_pair

end module caramel_generation
