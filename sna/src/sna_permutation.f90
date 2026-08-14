! Label permutation and optimization routines translated from R/sna permutation.R.
! Upstream copyright (C) Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_permutation
    use sna_kinds, only : dp, i8
    implicit none
    private

    abstract interface
        function graph_objective(a, b) result(v)
            import dp
            real(dp), intent(in) :: a(:,:), b(:,:)
            real(dp) :: v
        end function graph_objective
    end interface

    public :: graph_objective, numperm, rperm, rmperm, permute_matrix
    public :: lab_optimize_mc, lab_optimize_hillclimb, lab_optimize_anneal
    public :: lab_optimize_exhaustive, lab_optimize_gumbel

contains

    function numperm(n, permnum) result(o)
        integer, intent(in) :: n
        integer(i8), intent(in) :: permnum
        integer, allocatable :: o(:)
        integer(i8) :: pnum, relpos
        integer :: i, p, left
        allocate(o(n))
        o = 0
        pnum = permnum
        do i = 1, n
            left = n-i+1
            relpos = modulo(pnum, int(left,i8))
            p = 1
            do
                if (o(p) == 0) then
                    if (relpos == 0_i8) then
                        o(p) = i
                        exit
                    else
                        relpos = relpos - 1_i8
                    end if
                end if
                p = p + 1
            end do
            pnum = pnum / int(left,i8)
        end do
    end function numperm

    function rperm(exchange_list) result(o)
        integer, intent(in) :: exchange_list(:)
        integer, allocatable :: o(:)
        integer, allocatable :: ind(:)
        integer :: n, i, j, k, ng
        integer :: cls
        real(dp) :: u
        logical :: seen
        n = size(exchange_list)
        allocate(o(n))
        o = [(i,i=1,n)]
        do i = 1, n
            cls = exchange_list(i)
            seen = .false.
            do j = 1, i-1
                if (exchange_list(j) == cls) then
                    seen = .true.
                    exit
                end if
            end do
            if (seen) cycle
            ng = count(exchange_list == cls)
            allocate(ind(ng))
            k = 0
            do j = 1, n
                if (exchange_list(j) == cls) then
                    k = k+1
                    ind(k)=j
                end if
            end do
            do j = ng, 2, -1
                call random_number(u)
                k = 1 + int(u*real(j,dp))
                if (k > j) k=j
                call swap_int(o(ind(j)), o(ind(k)))
            end do
            deallocate(ind)
        end do
    end function rperm

    function rmperm(m, exchange_list) result(p)
        real(dp), intent(in) :: m(:,:)
        integer, intent(in), optional :: exchange_list(:)
        real(dp), allocatable :: p(:,:)
        integer, allocatable :: o(:), el(:)
        integer :: n, i
        n=size(m,1)
        allocate(el(n))
        if (present(exchange_list)) then
            el=exchange_list
        else
            el=0
        end if
        o=rperm(el)
        allocate(p(n,n))
        do i=1,n
            p(i,:)=m(o(i),o)
        end do
    end function rmperm

    function permute_matrix(m,o) result(p)
        real(dp), intent(in) :: m(:,:)
        integer, intent(in) :: o(:)
        real(dp), allocatable :: p(:,:)
        integer :: i,n
        n=size(o)
        allocate(p(n,n))
        do i=1,n
            p(i,:)=m(o(i),o)
        end do
    end function permute_matrix

    function lab_optimize_mc(d1,d2,fun,exchange_list,seek_max,draws) result(best)
        real(dp), intent(in) :: d1(:,:), d2(:,:)
        procedure(graph_objective) :: fun
        integer, intent(in), optional :: exchange_list(:)
        logical, intent(in), optional :: seek_max
        integer, intent(in), optional :: draws
        real(dp) :: best, v
        integer :: n, k, nd
        integer, allocatable :: el(:), o(:)
        real(dp), allocatable :: p(:,:)
        logical :: mx
        n=size(d1,1)
        nd=1000
        if(present(draws)) nd=draws
        mx=.false.
        if(present(seek_max)) mx=seek_max
        allocate(el(n))
        if(present(exchange_list)) then
        el=exchange_list
        else
        el=0
        end if
        best=fun(d1,d2)
        do k=1,nd
            o=rperm(el)
            p=permute_matrix(d2,o)
            v=fun(d1,p)
            if(mx) then
            best=max(best,v)
            else
            best=min(best,v)
            end if
        end do
    end function lab_optimize_mc

    function lab_optimize_hillclimb(d1,d2,fun,exchange_list,seek_max) result(best)
        real(dp), intent(in) :: d1(:,:), d2(:,:)
        procedure(graph_objective) :: fun
        integer, intent(in), optional :: exchange_list(:)
        logical, intent(in), optional :: seek_max
        real(dp) :: best, v, bestv
        integer :: n,i,j,bi,bj
        integer, allocatable :: el(:), o(:), cand(:)
        real(dp), allocatable :: p(:,:)
        logical :: mx, improved
        n=size(d1,1)
        mx=.false.
        if(present(seek_max)) mx=seek_max
        allocate(el(n),o(n),cand(n))
        if(present(exchange_list)) then
        el=exchange_list
        else
        el=0
        end if
        o=[(i,i=1,n)]
        best=fun(d1,d2)
        do
            improved=.false.
            bestv=best
            bi=0
            bj=0
            do i=1,n-1
                do j=i+1,n
                    if(el(i)/=el(j)) cycle
                    cand=o
                    call swap_int(cand(i),cand(j))
                    p=permute_matrix(d2,cand)
                    v=fun(d1,p)
                    if((mx .and. v>bestv) .or. (.not.mx .and. v<bestv)) then
                        bestv=v
                        bi=i
                        bj=j
                    end if
                end do
            end do
            if(bi>0) then
                call swap_int(o(bi),o(bj))
                best=bestv
                improved=.true.
            end if
            if(.not.improved) exit
        end do
    end function lab_optimize_hillclimb

    function lab_optimize_anneal(d1,d2,fun,exchange_list,seek_max,prob_init,prob_decay,freeze_time,full_neighborhood) &
        & result(global_best)
        real(dp), intent(in) :: d1(:,:), d2(:,:)
        procedure(graph_objective) :: fun
        integer, intent(in), optional :: exchange_list(:)
        logical, intent(in), optional :: seek_max, full_neighborhood
        real(dp), intent(in), optional :: prob_init, prob_decay
        integer, intent(in), optional :: freeze_time
        real(dp) :: global_best,best,v,prob,pdec,u,bestcand
        integer :: n,i,j,bi,bj,ftime,iter,nlegal,kk
        integer, allocatable :: el(:),o(:),cand(:),li(:),lj(:)
        real(dp), allocatable :: p(:,:)
        logical :: mx,full,improved,accept
        n=size(d1,1)
        mx=.false.
        if(present(seek_max)) mx=seek_max
        full=.true.
        if(present(full_neighborhood)) full=full_neighborhood
        prob=1.0_dp
        if(present(prob_init)) prob=prob_init
        pdec=0.99_dp
        if(present(prob_decay)) pdec=prob_decay
        ftime=1000
        if(present(freeze_time)) ftime=freeze_time
        allocate(el(n),o(n),cand(n))
        if(present(exchange_list)) then
        el=exchange_list
        else
        el=0
        end if
        o=[(i,i=1,n)]
        best=fun(d1,d2)
        global_best=best
        nlegal=0
        do i=1,n-1
        do j=i+1,n
        if(el(i)==el(j)) nlegal=nlegal+1
        end do
        end do
        if(nlegal==0) return
        allocate(li(nlegal),lj(nlegal))
        kk=0
        do i=1,n-1
        do j=i+1,n
        if(el(i)==el(j)) then
        kk=kk+1
        li(kk)=i
        lj(kk)=j
        end if
        end do
        end do
        iter=0
        do
            iter=iter+1
            improved=.false.
            if(full) then
                bestcand=best
                bi=0
                bj=0
                do kk=1,nlegal
                    cand=o
                    call swap_int(cand(li(kk)),cand(lj(kk)))
                    p=permute_matrix(d2,cand)
                    v=fun(d1,p)
                    if((mx .and. v>bestcand) .or. (.not.mx .and. v<bestcand)) then
                        bestcand=v
                        bi=li(kk)
                        bj=lj(kk)
                    end if
                end do
                if(bi>0) then
                    call swap_int(o(bi),o(bj))
                    best=bestcand
                    improved=.true.
                else if(ftime>0) then
                    call random_number(u)
                    if(u<prob) then
                        call random_number(u)
                        kk=1+int(u*real(nlegal,dp))
                        if(kk>nlegal)kk=nlegal
                        call swap_int(o(li(kk)),o(lj(kk)))
                        p=permute_matrix(d2,o)
                        best=fun(d1,p)
                    end if
                end if
            else
                call random_number(u)
                kk=1+int(u*real(nlegal,dp))
                if(kk>nlegal)kk=nlegal
                cand=o
                call swap_int(cand(li(kk)),cand(lj(kk)))
                p=permute_matrix(d2,cand)
                v=fun(d1,p)
                accept=(mx.and.v>best).or.((.not.mx).and.v<best)
                if(.not.accept .and. ftime>0) then
                call random_number(u)
                accept=u<prob
                end if
                if(accept) then
                o=cand
                improved=((mx.and.v>best).or.((.not.mx).and.v<best))
                best=v
                end if
            end if
            if((mx.and.best>global_best).or.((.not.mx).and.best<global_best)) global_best=best
            ftime=ftime-1
            prob=prob*pdec
            if(ftime<=0 .and. .not.improved) exit
            if(iter>1000000) exit
        end do
    end function lab_optimize_anneal

    function lab_optimize_exhaustive(d1,d2,fun,exchange_list,seek_max,max_permutations) result(best)
        real(dp), intent(in) :: d1(:,:), d2(:,:)
        procedure(graph_objective) :: fun
        integer, intent(in), optional :: exchange_list(:)
        logical, intent(in), optional :: seek_max
        integer(i8), intent(in), optional :: max_permutations
        real(dp) :: best,v
        integer :: n,i
        integer(i8) :: k,nperm,cap
        integer, allocatable :: el(:),o(:)
        real(dp), allocatable :: p(:,:)
        logical :: mx,legal
        n=size(d1,1)
        mx=.false.
        if(present(seek_max)) mx=seek_max
        allocate(el(n))
        if(present(exchange_list)) then
        el=exchange_list
        else
        el=0
        end if
        nperm=factorial_i8(n)
        cap=nperm
        if(present(max_permutations)) cap=min(cap,max_permutations)
        best=fun(d1,d2)
        do k=0_i8,cap-1_i8
            o=numperm(n,k)
            legal=.true.
            do i=1,n
                if(el(i)/=el(o(i))) then
                legal=.false.
                exit
                end if
            end do
            if(.not.legal) cycle
            p=permute_matrix(d2,o)
            v=fun(d1,p)
            if(mx) then
            best=max(best,v)
            else
            best=min(best,v)
            end if
        end do
    end function lab_optimize_exhaustive

    pure function factorial_i8(n) result(v)
        integer,intent(in)::n
        integer(i8)::v
        integer::i
        v=1_i8
        do i=2,n
            if(v>huge(v)/int(i,i8)) then
            v=huge(v)
            return
            end if
            v=v*int(i,i8)
        end do
    end function factorial_i8

    pure subroutine swap_int(a,b)
        integer,intent(inout)::a,b
        integer::t
        t=a
        a=b
        b=t
    end subroutine swap_int
    function lab_optimize_gumbel(d1,d2,fun,exchange_list,draws,tol,estimator) result(v)
        ! Monte Carlo/Gumbel extreme-value approximation from sna::lab.optimize.gumbel.
        real(dp), intent(in) :: d1(:,:),d2(:,:)
        procedure(graph_objective) :: fun
        integer, intent(in), optional :: exchange_list(:),draws
        real(dp), intent(in), optional :: tol
        character(len=*), intent(in), optional :: estimator
        real(dp) :: v,b,b_old,bdiff,mfg,a,atol
        real(dp), allocatable :: fg(:),pdat(:,:)
        integer, allocatable :: ex(:),ord(:)
        integer :: n,nd,i
        character(len=12) :: est
        n=size(d1,1)
        if(size(d1,2)/=n.or.any(shape(d2)/=shape(d1))) error stop 'lab_optimize_gumbel: graph orders differ'
        nd=500
        if(present(draws))nd=draws
        atol=1.0e-5_dp
        if(present(tol))atol=tol
        est='median'
        if(present(estimator))est=trim(estimator)
        allocate(ex(n))
        if(present(exchange_list))then
        if(size(exchange_list)/=n)error stop 'lab_optimize_gumbel: exchange-list size mismatch'
        ex=exchange_list
        else
        ex=0
        end if
        allocate(fg(nd))
        do i=1,nd
            ord=rperm(ex)
            pdat=permute_matrix(d2,ord)
            fg(i)=fun(d1,pdat)
        end do
        mfg=sum(fg)/real(nd,dp)
        b=1.0_dp
        bdiff=huge(1.0_dp)
        do i=1,10000
            b_old=b
            ! Shift before exponentiation; the ratio is invariant to the common shift.
            block
                real(dp) :: zmax,den
                real(dp), allocatable :: w(:)
                zmax=maxval(-fg/max(abs(b),sqrt(tiny(1.0_dp))))
                w=exp(-fg/max(abs(b),sqrt(tiny(1.0_dp)))-zmax)
                den=sum(w)
                if(den<=tiny(1.0_dp))exit
                b=mfg-sum(fg*w)/den
            end block
            if(abs(b)<=sqrt(tiny(1.0_dp)))b=sign(sqrt(tiny(1.0_dp)),merge(b,1.0_dp,b/=0.0_dp))
            bdiff=abs(b_old-b)
            if(bdiff<=atol)exit
        end do
        block
            real(dp) :: zmax
            zmax=maxval(-fg/b)
            a=-b*(zmax+log(sum(exp(-fg/b-zmax))/real(nd,dp)))
        end block
        select case(trim(est))
        case('mean')
        v=a+0.5772156649015328606_dp*b
        case('mode')
        v=a
        case default
        v=a-b*log(log(2.0_dp))
        end select
    end function lab_optimize_gumbel

end module sna_permutation
