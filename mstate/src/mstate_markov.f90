module mstate_markov
    use mstate_kinds, only : dp
    use mstate_types, only : msdata_type, etmdata_type, transition_map, markov_test_result
    use mstate_transitions, only : trans2q
    use survival_types, only : coxph_result
    use survival_cox, only : coxph_fit_counting
    use survival_linalg, only : invert_matrix
    implicit none
    private
    public :: markov_test

contains

    subroutine markov_test(ms, tr, transition, grid, b, result, covariates, offset, dist, seed, info)
        type(msdata_type), intent(in) :: ms
        type(transition_map), intent(in) :: tr
        integer, intent(in) :: transition, b
        real(dp), intent(in) :: grid(:)
        type(markov_test_result), intent(out) :: result
        real(dp), intent(in), optional :: covariates(:, :), offset(:)
        character(len=*), intent(in), optional :: dist
        integer, intent(in), optional :: seed
        integer, intent(out), optional :: info
        type(etmdata_type) :: etm
        type(coxph_result) :: fit
        integer, allocatable :: q(:, :), qual(:), relpat(:), ridx(:), status(:), nobs(:)
        real(dp), allocatable :: cov_etm(:, :), off_etm(:), zmat(:, :), roff(:)
        real(dp), allocatable :: t0(:), t1(:), lp(:), explp(:), s0(:), incr(:), resid(:, :)
        real(dp), allocatable :: z3(:, :, :), zbar(:, :, :), zbar0(:, :), obs(:, :)
        real(dp), allocatable :: hmat(:, :, :), multiplier(:, :, :), estcov(:, :, :), estvar(:, :)
        real(dp), allocatable :: wb(:, :, :), nwb(:, :, :), nch(:, :), gvec(:), imul(:)
        real(dp), allocatable :: subvec(:), inv(:, :), vals(:)
        integer :: froms, tos, nq, ng, nc, nr, i, j, k, v, x, y, qidx, indx1, indx2, wbidx
        integer :: nrel, nsubdim
        real(dp) :: den, term1, term2, accum, trace1, rtmp, lower, upper
        character(len=8) :: distribution
        logical :: ok

        if (present(info)) info = 0
        if (transition < 1 .or. transition > tr%ntrans .or. b < 1 .or. size(grid) < 1) then
            if (present(info)) info = 1
            return
        end if
        if (present(covariates)) then
            if (size(covariates,1) /= ms%n) then
                if (present(info)) info = 2
                return
            end if
        end if
        if (present(offset)) then
            if (size(offset) /= ms%n) then
                if (present(info)) info = 3
                return
            end if
        end if
        if (present(seed)) call set_seed_scalar(seed)
        distribution = 'poisson'
        if (present(dist)) distribution = adjustl(dist)
        if (index(distribution,'poisson') /= 1 .and. index(distribution,'normal') /= 1) then
            if (present(info)) info = 4
            return
        end if

        froms = tr%from(transition)
        tos = tr%to(transition)
        call trans2q(tr,q)
        allocate(qual(tr%nstate)); nq = 1; qual(1) = froms
        do i = 1, tr%nstate
            if (i == froms) cycle
            if (q(i,froms) > 0) then
                nq = nq + 1; qual(nq) = i
            end if
        end do
        call sort_int(qual(1:nq))
        qual = qual

        call msdata_to_etm_cov(ms, etm, covariates, cov_etm, offset, off_etm)
        call relevant_ids(etm, froms, relpat)
        nrel = size(relpat)
        nr = count(etm%from == froms)
        if (nr == 0) then
            if (present(info)) info = 5
            return
        end if
        allocate(ridx(nr), t0(nr), t1(nr), status(nr), roff(nr))
        nc = 0
        if (present(covariates)) nc = size(covariates,2)
        allocate(zmat(nr,nc))
        k = 0
        do i = 1, etm%n
            if (etm%from(i) /= froms) cycle
            k = k + 1
            ridx(k) = i; t0(k) = etm%entry(i); t1(k) = etm%exit(i)
            status(k) = merge(1,0,etm%to(i) == tos)
            if (nc > 0) zmat(k,:) = cov_etm(i,:)
            roff(k) = off_etm(i)
        end do

        allocate(lp(nr), explp(nr))
        if (nc > 0) then
            call coxph_fit_counting(t0, t1, status, zmat, fit, 'efron', offset=roff)
            lp = matmul(zmat,fit%coef) + roff
            allocate(result%beta(nc), result%beta_vcov(nc,nc))
            result%beta = fit%coef; result%beta_vcov = fit%var
        else
            lp = roff
            allocate(result%beta(0), result%beta_vcov(0,0))
        end if
        explp = exp(min(lp,700.0_dp))
        ng = size(grid)
        allocate(nobs(ng)); nobs = 0
        do j = 1, ng
            nobs(j) = count(status == 1 .and. t1 > grid(j))
        end do

        allocate(z3(nr,ng,nq)); z3 = 0.0_dp
        do i = 1, nr
            do j = 1, ng
                do qidx = 1, nq
                    z3(i,j,qidx) = state_at_time(etm, etm%id(ridx(i)), grid(j), qual(qidx))
                end do
            end do
        end do

        allocate(s0(nr), incr(nr), resid(nr,ng))
        do i = 1, nr
            s0(i) = sum(explp, mask=(t0 < t1(i) .and. t1 >= t1(i)))
            if (s0(i) > 0.0_dp) then
                incr(i) = real(status(i),dp)/s0(i)
            else
                incr(i) = 0.0_dp
            end if
        end do
        do i = 1, nr
            do j = 1, ng
                resid(i,j) = real(status(i),dp)*merge(1.0_dp,0.0_dp,t1(i)>grid(j)) - explp(i)* &
                    (cumhaz_at(incr,t1,max(grid(j),t1(i))) - cumhaz_at(incr,t1,max(grid(j),t0(i))))
            end do
        end do

        allocate(obs(ng,nq)); obs = 0.0_dp
        do qidx = 1, nq
            do j = 1, ng
                obs(j,qidx) = sum(resid(:,j)*z3(:,j,qidx), mask=t1>grid(j))
            end do
        end do

        allocate(zbar(nr,ng,nq)); zbar = 0.0_dp
        allocate(zbar0(nr,nc)); zbar0 = 0.0_dp
        do i = 1, nr
            den = sum(explp, mask=(t0 < t1(i) .and. t1 >= t1(i)))
            if (den <= 0.0_dp) cycle
            if (nc > 0) then
                do j = 1, nc
                    zbar0(i,j) = sum(zmat(:,j)*explp, mask=(t0<t1(i) .and. t1>=t1(i)))/den
                end do
            end if
            do y = 1, ng
                do qidx = 1, nq
                    zbar(i,y,qidx) = sum(z3(:,y,qidx)*explp, mask=(t0<t1(i) .and. t1>=t1(i)))/den
                end do
            end do
        end do

        allocate(hmat(ng,nc,nq), multiplier(ng,nc,nq)); hmat = 0.0_dp; multiplier = 0.0_dp
        if (nc > 0) then
            do j = 1, nc
                do qidx = 1, nq
                    do y = 1, ng
                        accum = 0.0_dp
                        do x = 1, nr
                            if (t1(x) <= grid(y)) cycle
                            do i = 1, nr
                                if (.not. (t1(i) > t0(x) .and. t1(i) <= t1(x))) cycle
                                accum = accum + explp(x)*(zmat(x,j)-zbar0(i,j))* &
                                    (z3(x,y,qidx)-zbar(i,y,qidx))*incr(i)
                            end do
                        end do
                        hmat(y,j,qidx) = accum
                    end do
                end do
            end do
            do qidx = 1, nq
                multiplier(:,:,qidx) = matmul(hmat(:,:,qidx),fit%var)
            end do
        end if

        allocate(estcov(ng,nq,nq)); estcov = 0.0_dp
        do indx1 = 1, nq
            do indx2 = indx1, nq
                do y = 1, ng
                    accum = 0.0_dp
                    do v = 1, nr
                        do i = 1, nr
                            if (.not. (t0(v) < t1(i) .and. t1(v) >= t1(i))) cycle
                            if (nc > 0) then
                                term1 = (z3(v,y,indx1)-zbar(i,y,indx1))*merge(1.0_dp,0.0_dp,t1(i)>grid(y))
                                term2 = (z3(v,y,indx2)-zbar(i,y,indx1))*merge(1.0_dp,0.0_dp,t1(i)>grid(y))
                                do j = 1, nc
                                    term1 = term1 - multiplier(y,j,indx1)*(zmat(v,j)-zbar0(i,j))
                                    term2 = term2 - multiplier(y,j,indx2)*(zmat(v,j)-zbar0(i,j))
                                end do
                            else
                                term1 = z3(v,y,indx1)-zbar(i,y,indx1)
                                term2 = z3(v,y,indx2)-zbar(i,y,indx2)
                            end if
                            accum = accum + term1*term2*explp(v)*merge(1.0_dp,0.0_dp,t1(i)>grid(y))*incr(i)
                        end do
                    end do
                    estcov(y,indx1,indx2) = accum
                    estcov(y,indx2,indx1) = accum
                end do
            end do
        end do

        allocate(estvar(ng,nq)); estvar = 0.0_dp
        allocate(result%zbar(ng,nq), result%obs_chisq_trace(ng), result%est_cov(ng,nq,nq))
        result%est_cov = estcov
        do qidx = 1, nq
            estvar(:,qidx) = estcov(:,qidx,qidx)
            do j = 1, ng
                if (estvar(j,qidx) > 0.0_dp) result%zbar(j,qidx) = obs(j,qidx)/sqrt(estvar(j,qidx))
            end do
        end do
        result%obs_chisq_trace = 0.0_dp
        nsubdim = nq-1
        if (nsubdim > 0) then
            allocate(subvec(nsubdim), inv(nsubdim,nsubdim))
            do y = 1, ng
                call invert_matrix(estcov(y,2:nq,2:nq),inv,ok)
                if (ok) then
                    subvec = obs(y,2:nq)
                    result%obs_chisq_trace(y) = dot_product(subvec,matmul(inv,subvec))
                end if
            end do
        end if

        allocate(wb(b,ng,nq),nwb(b,ng,nq),nch(b,ng),gvec(nr),imul(nc))
        wb = 0.0_dp; nwb = 0.0_dp; nch = 0.0_dp
        do wbidx = 1, b
            call wild_weights(distribution,gvec)
            do qidx = 1, nq
                do y = 1, ng
                    rtmp = sum(real(status,dp)*(z3(:,y,qidx)-zbar(:,y,qidx))* &
                               merge(1.0_dp,0.0_dp,t1>grid(y))*gvec)
                    trace1 = 0.0_dp
                    if (nc > 0) then
                        do j = 1, nc
                            imul(j) = sum(real(status,dp)*(zmat(:,j)-zbar0(:,j))*gvec)
                        end do
                        trace1 = dot_product(hmat(y,:,qidx),matmul(fit%var,imul))
                    end if
                    wb(wbidx,y,qidx) = rtmp-trace1
                    if (estvar(y,qidx) > 0.0_dp) nwb(wbidx,y,qidx) = wb(wbidx,y,qidx)/sqrt(estvar(y,qidx))
                end do
            end do
            if (nsubdim > 0) then
                do y = 1, ng
                    call invert_matrix(estcov(y,2:nq,2:nq),inv,ok)
                    if (ok) then
                        subvec = wb(wbidx,y,2:nq)
                        nch(wbidx,y) = dot_product(subvec,matmul(inv,subvec))
                    end if
                end do
            end if
        end do

        allocate(result%n_wb_trace(b,ng,nq), result%nch_wb_trace(b,ng))
        result%n_wb_trace = nwb; result%nch_wb_trace = nch
        allocate(result%orig_stat(nq), result%p_stat_wb(nq), result%est_quant(2,ng,nq))
        allocate(vals(b))
        do qidx = 1, nq
            result%orig_stat(qidx) = sum(abs(result%zbar(:,qidx)))/real(ng,dp)
            do wbidx = 1, b
                vals(wbidx) = sum(abs(nwb(wbidx,:,qidx)))/real(ng,dp)
            end do
            result%p_stat_wb(qidx) = real(count(vals>result%orig_stat(qidx)),dp)/real(b,dp)
            do y = 1, ng
                vals = nwb(:,y,qidx)
                call quantile_pair(vals,0.025_dp,0.975_dp,lower,upper)
                result%est_quant(1,y,qidx)=lower; result%est_quant(2,y,qidx)=upper
            end do
        end do
        result%orig_ch_stat = sum(result%obs_chisq_trace)/real(ng,dp)
        do wbidx = 1, b
            vals(wbidx) = sum(nch(wbidx,:))/real(ng,dp)
        end do
        result%p_ch_stat_wb = real(count(vals>result%orig_ch_stat),dp)/real(b,dp)

        result%transition = transition; result%from_state = froms; result%to_state = tos
        result%b = b; result%nsub = nrel; result%dist = distribution
        allocate(result%qualset(nq),result%nobs_grid(ng),result%grid(ng))
        result%qualset = qual(1:nq); result%nobs_grid = nobs; result%grid = grid
    end subroutine markov_test

    subroutine msdata_to_etm_cov(ms, etm, cov_ms, cov_etm, off_ms, off_etm)
        type(msdata_type), intent(in) :: ms
        type(etmdata_type), intent(out) :: etm
        real(dp), intent(in), optional :: cov_ms(:, :), off_ms(:)
        real(dp), allocatable, intent(out) :: cov_etm(:, :), off_etm(:)
        integer :: i, j, n, row, nc, tostate
        logical :: same
        nc=0; if(present(cov_ms))nc=size(cov_ms,2)
        n=0;i=1
        do while(i<=ms%n)
            n=n+1;j=i+1
            do while(j<=ms%n)
                same=ms%id(j)==ms%id(i).and.ms%from(j)==ms%from(i).and. &
                     ms%tstart(j)==ms%tstart(i).and.ms%tstop(j)==ms%tstop(i)
                if(.not.same)exit
                j=j+1
            end do
            i=j
        end do
        etm%n=n
        allocate(etm%id(n),etm%from(n),etm%to(n),etm%entry(n),etm%exit(n),cov_etm(n,nc),off_etm(n))
        cov_etm=0.0_dp;off_etm=0.0_dp;row=0;i=1
        do while(i<=ms%n)
            row=row+1;j=i;tostate=99
            do while(j<=ms%n)
                same=ms%id(j)==ms%id(i).and.ms%from(j)==ms%from(i).and. &
                     ms%tstart(j)==ms%tstart(i).and.ms%tstop(j)==ms%tstop(i)
                if(.not.same)exit
                if(ms%status(j)==1)tostate=ms%to(j)
                j=j+1
            end do
            etm%id(row)=ms%id(i);etm%from(row)=ms%from(i);etm%to(row)=tostate
            etm%entry(row)=ms%tstart(i);etm%exit(row)=ms%tstop(i)
            if(present(cov_ms).and.nc>0)cov_etm(row,:)=cov_ms(i,:)
            if(present(off_ms))off_etm(row)=off_ms(i)
            i=j
        end do
    end subroutine msdata_to_etm_cov

    subroutine relevant_ids(etm, froms, ids)
        type(etmdata_type), intent(in) :: etm
        integer, intent(in) :: froms
        integer, allocatable, intent(out) :: ids(:)
        integer, allocatable :: tmp(:)
        integer :: i,n
        allocate(tmp(etm%n));n=0
        do i=1,etm%n
            if(etm%from(i)/=froms)cycle
            if(n==0.or..not.any(tmp(1:n)==etm%id(i)))then
                n=n+1;tmp(n)=etm%id(i)
            end if
        end do
        allocate(ids(n));if(n>0)ids=tmp(1:n)
        call sort_int(ids)
    end subroutine relevant_ids

    real(dp) function state_at_time(etm,id,t,state) result(v)
        type(etmdata_type),intent(in)::etm
        integer,intent(in)::id,state
        real(dp),intent(in)::t
        integer::i
        v=0.0_dp
        do i=1,etm%n
            if(etm%id(i)==id.and.etm%entry(i)<t.and.etm%exit(i)>=t)then
                if(etm%from(i)==state)v=1.0_dp
                return
            end if
        end do
    end function state_at_time

    real(dp) function cumhaz_at(incr,times,t) result(h)
        real(dp),intent(in)::incr(:),times(:),t
        h=sum(incr,mask=times<=t)
    end function cumhaz_at

    subroutine wild_weights(dist,g)
        character(len=*),intent(in)::dist
        real(dp),intent(out)::g(:)
        integer::i,k
        real(dp)::u,p,l,u1,u2,pi
        pi=acos(-1.0_dp)
        if(index(dist,'poisson')==1)then
            l=exp(-1.0_dp)
            do i=1,size(g)
                k=0;p=1.0_dp
                do while(p>l)
                    k=k+1;call random_number(u);p=p*u
                end do
                g(i)=real(k-1,dp)-1.0_dp
            end do
        else
            i=1
            do while(i<=size(g))
                call random_number(u1);call random_number(u2);u1=max(u1,tiny(1.0_dp))
                g(i)=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
                if(i+1<=size(g))g(i+1)=sqrt(-2.0_dp*log(u1))*sin(2.0_dp*pi*u2)
                i=i+2
            end do
        end if
    end subroutine wild_weights

    subroutine quantile_pair(x,p1,p2,q1,q2)
        real(dp),intent(in)::x(:),p1,p2
        real(dp),intent(out)::q1,q2
        real(dp),allocatable::y(:)
        integer::i,j
        real(dp)::key
        allocate(y(size(x)));y=x
        do i=2,size(y)
            key=y(i);j=i-1
            do while(j>=1)
                if(y(j)<=key)exit
                y(j+1)=y(j);j=j-1
            end do
            y(j+1)=key
        end do
        q1=type7_quantile(y,p1);q2=type7_quantile(y,p2)
    end subroutine quantile_pair

    real(dp) function type7_quantile(y,p) result(q)
        real(dp),intent(in)::y(:),p
        real(dp)::h,frac
        integer::lo,hi,n
        n=size(y);if(n==1)then;q=y(1);return;end if
        h=1.0_dp+(real(n-1,dp))*p;lo=floor(h);hi=ceiling(h);frac=h-real(lo,dp)
        q=(1.0_dp-frac)*y(lo)+frac*y(hi)
    end function type7_quantile

    subroutine sort_int(x)
        integer,intent(inout)::x(:)
        integer::i,j,key
        do i=2,size(x)
            key=x(i);j=i-1
            do while(j>=1)
                if(x(j)<=key)exit
                x(j+1)=x(j);j=j-1
            end do
            x(j+1)=key
        end do
    end subroutine sort_int

    subroutine set_seed_scalar(seed)
        integer,intent(in)::seed
        integer::n,i
        integer,allocatable::put(:)
        call random_seed(size=n);allocate(put(n))
        do i=1,n
            put(i)=mod(abs(seed)+104729*i,2147483646)+1
        end do
        call random_seed(put=put)
    end subroutine set_seed_scalar

end module mstate_markov
