module mstate_crprep
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use mstate_kinds, only : dp
    use mstate_types, only : crprep_type
    implicit none
    private
    public :: crprep

contains

    subroutine crprep(tstop, status, trans_codes, cens_code, result, tstart, id, strata, keep, shorten, origin, &
                      prec_factor, info)
        real(dp), intent(in) :: tstop(:)
        integer, intent(in) :: status(:), trans_codes(:), cens_code
        type(crprep_type), intent(out) :: result
        real(dp), intent(in), optional :: tstart(:), keep(:, :)
        integer, intent(in), optional :: id(:), strata(:)
        logical, intent(in), optional :: shorten
        real(dp), intent(in), optional :: origin, prec_factor
        integer, intent(out), optional :: info

        real(dp), allocatable :: ts(:), te(:), kv(:, :)
        integer, allocatable :: ids(:), str(:), nw(:)
        integer :: n, nf, nkeep, f, i, j, nr, pos, nall, st
        integer :: failcode, n_event_times, row0, row1
        real(dp), allocatable :: event_times(:)
        real(dp) :: org, pf, prec, g0, h0
        logical :: do_shorten, calc_trunc
        type(crprep_type) :: full

        if (present(info)) info = 0
        n = size(tstop)
        if (size(status) /= n .or. size(trans_codes) == 0) then
            if (present(info)) info = 1
            call empty_result(result)
            return
        end if
        if (present(tstart)) then
            if (size(tstart) /= n) then
                if (present(info)) info = 2
                call empty_result(result)
                return
            end if
        end if
        if (present(id)) then
            if (size(id) /= n) then
                if (present(info)) info = 3
                call empty_result(result)
                return
            end if
        end if
        if (present(strata)) then
            if (size(strata) /= n) then
                if (present(info)) info = 4
                call empty_result(result)
                return
            end if
        end if
        nkeep = 0
        if (present(keep)) then
            if (size(keep, 1) /= n) then
                if (present(info)) info = 5
                call empty_result(result)
                return
            end if
            nkeep = size(keep, 2)
        end if

        org = 0.0_dp
        if (present(origin)) org = origin
        pf = 1000.0_dp
        if (present(prec_factor)) pf = prec_factor
        prec = epsilon(1.0_dp) * pf
        do_shorten = .true.
        if (present(shorten)) do_shorten = shorten

        allocate(ts(n), te(n), ids(n), str(n))
        ts = -org
        if (present(tstart)) ts = tstart - org
        te = tstop - org
        ids = [(i, i=1,n)]
        if (present(id)) ids = id
        str = 1
        if (present(strata)) str = strata
        if (nkeep > 0) then
            allocate(kv(n, nkeep)); kv = keep
        else
            allocate(kv(n, 0))
        end if

        if (any(.not. ieee_is_finite(ts)) .or. any(.not. ieee_is_finite(te))) then
            if (present(info)) info = 6
            call empty_result(result)
            return
        end if
        if (any(te < ts)) then
            if (present(info)) info = 7
            call empty_result(result)
            return
        end if
        calc_trunc = any(abs(ts + org) > 0.0_dp)

        ! First count rows exactly as create.wData.omega does, for every failcode.
        nf = size(trans_codes)
        allocate(nw(n))
        nall = 0
        do f = 1, nf
            failcode = trans_codes(f)
            do i = 1, n
                call fail_event_times_stratum(te, status, str, str(i), failcode, event_times)
                nw(i) = 1
                if (status(i) /= cens_code .and. status(i) /= failcode) then
                    nw(i) = 1 + count(event_times > te(i))
                end if
                nall = nall + nw(i)
                if (allocated(event_times)) deallocate(event_times)
            end do
        end do

        full%n = nall
        full%nkeep = nkeep
        full%has_truncation = calc_trunc
        allocate(full%id(nall), full%status(nall), full%strata(nall), full%count(nall), full%failcode(nall))
        allocate(full%tstart(nall), full%tstop(nall), full%weight_cens(nall), full%weight_trunc(nall))
        allocate(full%keep(nall, nkeep))
        full%count = 0
        full%weight_cens = 1.0_dp
        full%weight_trunc = 1.0_dp
        if (nkeep > 0) full%keep = 0.0_dp

        pos = 0
        do f = 1, nf
            failcode = trans_codes(f)
            do i = 1, n
                call fail_event_times_stratum(te, status, str, str(i), failcode, event_times)
                n_event_times = size(event_times)
                nr = 1
                if (status(i) /= cens_code .and. status(i) /= failcode) nr = 1 + count(event_times > te(i))
                row0 = pos + 1
                do j = 1, nr
                    pos = pos + 1
                    full%id(pos) = ids(i)
                    full%status(pos) = status(i)
                    full%strata(pos) = str(i)
                    full%failcode(pos) = failcode
                    if (nkeep > 0) full%keep(pos,:) = kv(i,:)
                    if (j == 1) then
                        full%tstart(pos) = ts(i)
                        full%tstop(pos) = te(i)
                    else
                        if (j == 2) then
                            full%tstart(pos) = te(i)
                        else
                            full%tstart(pos) = event_times(n_event_times - nr + j - 1)
                        end if
                        full%tstop(pos) = event_times(n_event_times - nr + j)
                    end if
                end do
                row1 = pos

                ! Censoring product-limit weight, normalized to the first interval for the subject.
                st = str(i)
                g0 = censor_survival_at(ts, te, status, str, st, cens_code, full%tstop(row0)-prec, prec)
                if (nr == 1 .and. ieee_is_finite(g0)) then
                    full%weight_cens(row0) = 1.0_dp
                else
                    do j = row0, row1
                        if (g0 > 0.0_dp) then
                            full%weight_cens(j) = censor_survival_at(ts, te, status, str, st, cens_code, &
                                                                    full%tstop(j)-prec, prec) / g0
                        else
                            full%weight_cens(j) = 0.0_dp
                        end if
                    end do
                end if

                if (calc_trunc) then
                    h0 = trunc_survival_at(ts, te, str, st, -full%tstop(row0), prec)
                    if (nr == 1 .and. ieee_is_finite(h0)) then
                        full%weight_trunc(row0) = 1.0_dp
                    else
                        do j = row0, row1
                            if (h0 > 0.0_dp) then
                                full%weight_trunc(j) = trunc_survival_at(ts, te, str, st, -full%tstop(j), prec) / h0
                            else
                                full%weight_trunc(j) = 0.0_dp
                            end if
                        end do
                    end if
                end if
                if (allocated(event_times)) deallocate(event_times)
            end do
        end do

        ! Count rows within each subject/failcode block before optional shortening.
        call assign_counts(full)
        if (do_shorten) then
            call shorten_crprep(full, result)
        else
            result = full
        end if
    end subroutine crprep

    subroutine fail_event_times_stratum(tstop, status, strata, stratum, failcode, event_times)
        real(dp), intent(in) :: tstop(:)
        integer, intent(in) :: status(:), strata(:), stratum, failcode
        real(dp), allocatable, intent(out) :: event_times(:)
        real(dp), allocatable :: tmp(:)
        integer :: i, n
        n = count((status == failcode) .and. (strata == stratum))
        allocate(tmp(n)); n = 0
        do i = 1, size(tstop)
            if (status(i) == failcode .and. strata(i) == stratum) then
                n = n + 1; tmp(n) = tstop(i)
            end if
        end do
        call sort_unique(tmp, event_times)
    end subroutine fail_event_times_stratum

    real(dp) function censor_survival_at(tstart, tstop, status, strata, stratum, cens_code, query, prec) result(surv)
        real(dp), intent(in) :: tstart(:), tstop(:), query, prec
        integer, intent(in) :: status(:), strata(:), stratum, cens_code
        real(dp), allocatable :: ev(:), u(:)
        real(dp) :: tt, risk, deaths
        integer :: i, k, n
        n = count((strata == stratum) .and. (status == cens_code))
        if (n == 0) then
            surv = 1.0_dp; return
        end if
        allocate(ev(n)); n = 0
        do i = 1, size(tstop)
            if (strata(i) == stratum .and. status(i) == cens_code) then
                n = n + 1; ev(n) = tstop(i) + prec
            end if
        end do
        call sort_unique(ev, u)
        surv = 1.0_dp
        do k = 1, size(u)
            tt = u(k)
            if (tt > query) exit
            risk = 0.0_dp; deaths = 0.0_dp
            do i = 1, size(tstop)
                if (strata(i) /= stratum) cycle
                if (tstart(i) < tt .and. (tstop(i) + merge(prec,0.0_dp,status(i)==cens_code)) >= tt) &
                    risk = risk + 1.0_dp
                if (status(i) == cens_code .and. same_time(tstop(i)+prec, tt)) deaths = deaths + 1.0_dp
            end do
            if (risk > 0.0_dp .and. deaths > 0.0_dp) surv = surv * max(0.0_dp, 1.0_dp - deaths/risk)
        end do
    end function censor_survival_at

    real(dp) function trunc_survival_at(tstart, tstop, strata, stratum, query, prec) result(surv)
        real(dp), intent(in) :: tstart(:), tstop(:), query, prec
        integer, intent(in) :: strata(:), stratum
        real(dp), allocatable :: ev(:), u(:)
        real(dp) :: aa, bb, tt, risk, deaths
        integer :: i, k, n
        n = count(strata == stratum)
        allocate(ev(n)); n = 0
        do i = 1, size(tstop)
            if (strata(i) == stratum) then
                n = n + 1
                ev(n) = -(tstart(i) + 2.0_dp*prec)
            end if
        end do
        call sort_unique(ev, u)
        surv = 1.0_dp
        do k = 1, size(u)
            tt = u(k)
            if (tt > query) exit
            risk = 0.0_dp; deaths = 0.0_dp
            do i = 1, size(tstop)
                if (strata(i) /= stratum) cycle
                aa = -tstop(i)
                bb = -(tstart(i) + 2.0_dp*prec)
                if (aa < tt .and. bb >= tt) risk = risk + 1.0_dp
                if (same_time(bb,tt)) deaths = deaths + 1.0_dp
            end do
            if (risk > 0.0_dp .and. deaths > 0.0_dp) surv = surv * max(0.0_dp, 1.0_dp - deaths/risk)
        end do
    end function trunc_survival_at

    subroutine assign_counts(x)
        type(crprep_type), intent(inout) :: x
        integer :: i, c
        if (x%n == 0) return
        c = 0
        do i = 1, x%n
            if (i == 1) then
                c = 0
            else if (x%id(i) /= x%id(i-1) .or. x%failcode(i) /= x%failcode(i-1)) then
                c = 0
            end if
            c = c + 1
            x%count(i) = c
        end do
    end subroutine assign_counts

    subroutine shorten_crprep(src, dst)
        type(crprep_type), intent(in) :: src
        type(crprep_type), intent(out) :: dst
        logical, allocatable :: keep_row(:)
        integer :: i, j, nout, first, last, block_start
        if (src%n == 0) then
            call empty_result(dst); return
        end if
        allocate(keep_row(src%n)); keep_row = .false.
        i = 1
        do while (i <= src%n)
            first = i
            do while (i < src%n)
                if (src%id(i+1) /= src%id(first)) exit
                if (src%failcode(i+1) /= src%failcode(first)) exit
                i = i + 1
            end do
            last = i
            keep_row(first) = .true.
            if (first < last) then
                do j = first+1, last-1
                    if (.not. same_time(src%weight_cens(j),src%weight_cens(j+1))) keep_row(j)=.true.
                    if (src%has_truncation) then
                        if (.not. same_time(src%weight_trunc(j),src%weight_trunc(j+1))) keep_row(j)=.true.
                    end if
                end do
                keep_row(last) = .true.
            end if
            i = i + 1
        end do
        nout = count(keep_row)
        call allocate_crprep(dst,nout,src%nkeep,src%has_truncation)
        j = 0; i = 1
        do while (i <= src%n)
            first = i
            do while (i < src%n)
                if (src%id(i+1) /= src%id(first)) exit
                if (src%failcode(i+1) /= src%failcode(first)) exit
                i = i + 1
            end do
            last = i
            block_start = first
            do while (block_start <= last)
                if (keep_row(block_start)) then
                    j = j + 1
                    call copy_row(src,block_start,dst,j)
                    if (block_start == first) then
                        dst%tstart(j) = src%tstart(first)
                    else if (j > 1) then
                        dst%tstart(j) = dst%tstop(j-1)
                    end if
                end if
                block_start = block_start + 1
            end do
            i = i + 1
        end do
        call assign_counts(dst)
    end subroutine shorten_crprep

    subroutine allocate_crprep(x,n,nkeep,has_truncation)
        type(crprep_type), intent(out) :: x
        integer, intent(in) :: n,nkeep
        logical, intent(in) :: has_truncation
        x%n=n; x%nkeep=nkeep; x%has_truncation=has_truncation
        allocate(x%id(n),x%status(n),x%strata(n),x%count(n),x%failcode(n))
        allocate(x%tstart(n),x%tstop(n),x%weight_cens(n),x%weight_trunc(n),x%keep(n,nkeep))
    end subroutine allocate_crprep

    subroutine copy_row(src,isrc,dst,idst)
        type(crprep_type), intent(in) :: src
        integer, intent(in) :: isrc,idst
        type(crprep_type), intent(inout) :: dst
        dst%id(idst)=src%id(isrc); dst%status(idst)=src%status(isrc); dst%strata(idst)=src%strata(isrc)
        dst%count(idst)=src%count(isrc); dst%failcode(idst)=src%failcode(isrc)
        dst%tstart(idst)=src%tstart(isrc); dst%tstop(idst)=src%tstop(isrc)
        dst%weight_cens(idst)=src%weight_cens(isrc); dst%weight_trunc(idst)=src%weight_trunc(isrc)
        if (src%nkeep>0) dst%keep(idst,:)=src%keep(isrc,:)
    end subroutine copy_row

    subroutine empty_result(x)
        type(crprep_type), intent(out) :: x
        call allocate_crprep(x,0,0,.false.)
    end subroutine empty_result

    subroutine sort_unique(x,u)
        real(dp),intent(in)::x(:)
        real(dp),allocatable,intent(out)::u(:)
        real(dp),allocatable::y(:)
        real(dp)::key
        integer::i,j,n
        if(size(x)==0)then;allocate(u(0));return;end if
        allocate(y(size(x)));y=x
        do i=2,size(y)
            key=y(i);j=i-1
            do while(j>=1)
                if(y(j)<=key)exit
                y(j+1)=y(j);j=j-1
            end do
            y(j+1)=key
        end do
        n=1
        do i=2,size(y)
            if(.not.same_time(y(i),y(n)))then;n=n+1;y(n)=y(i);end if
        end do
        allocate(u(n));u=y(1:n)
    end subroutine sort_unique

    pure logical function same_time(a,b) result(eq)
        real(dp),intent(in)::a,b
        real(dp)::sc
        sc=max(1.0_dp,abs(a),abs(b))
        eq=abs(a-b)<=16.0_dp*epsilon(1.0_dp)*sc
    end function same_time

end module mstate_crprep
