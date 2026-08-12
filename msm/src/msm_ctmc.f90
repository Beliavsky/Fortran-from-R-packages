! SPDX-License-Identifier: GPL-2.0-or-later
! Continuous-time Markov computations translated from msm pijt.c and R summaries.
module msm_ctmc
    use msm_kinds, only : dp
    use msm_linalg, only : expm, expm_frechet, solve_linear, eye
    implicit none
    private
    integer, parameter, public :: obs_panel = 1
    integer, parameter, public :: obs_exact = 2
    integer, parameter, public :: obs_death = 3
    public :: make_generator, valid_generator, transition_matrix, transition_derivatives
    public :: observation_kernel, death_transition_density, ctmc_minus2loglik, ctmc_gradient
    public :: ctmc_aggregate_minus2loglik, ctmc_censored_minus2loglik
    public :: piecewise_transition_matrix, qmatrix_covariates
    public :: sojourn_times, next_state_probabilities, absorbing_states, transient_states
    public :: expected_total_time, expected_first_passage, eventual_passage_probability
    public :: passage_probability_matrix, expected_visits, state_prevalence
contains
    function make_generator(offdiag) result(q)
        real(dp), intent(in) :: offdiag(:,:)
        real(dp), allocatable :: q(:,:)
        integer :: n, i
        n = size(offdiag,1)
        if (size(offdiag,2) /= n) error stop "make_generator: matrix must be square"
        q = offdiag
        do i = 1, n
            q(i,i) = 0.0_dp
            if (any(q(i,:) < 0.0_dp)) error stop "make_generator: negative off-diagonal intensity"
            q(i,i) = -sum(q(i,:))
        end do
    end function make_generator

    pure function valid_generator(q, tol) result(ok)
        real(dp), intent(in) :: q(:,:)
        real(dp), intent(in), optional :: tol
        logical :: ok
        real(dp) :: eps
        integer :: i, n
        eps = 1.0e-10_dp
        if (present(tol)) eps = tol
        n = size(q,1)
        if (size(q,2) /= n) then
            ok = .false.
            return
        end if
        ok = .true.
        do i = 1, n
            if (q(i,i) > eps) ok = .false.
            if (abs(sum(q(i,:))) > eps*(1.0_dp+sum(abs(q(i,:))))) ok = .false.
        end do
        do i = 1, n
            if (any(q(i,1:i-1) < -eps)) ok = .false.
            if (i < n) then
                if (any(q(i,i+1:n) < -eps)) ok = .false.
            end if
        end do
    end function valid_generator

    function transition_matrix(q, t, clean) result(p)
        real(dp), intent(in) :: q(:,:), t
        logical, intent(in), optional :: clean
        real(dp), allocatable :: p(:,:)
        logical :: do_clean
        integer :: i, j, n
        if (t < 0.0_dp) error stop "transition_matrix: t must be nonnegative"
        n = size(q,1)
        if (size(q,2) /= n) error stop "transition_matrix: q must be square"
        p = expm(q*t)
        do_clean = .true.
        if (present(clean)) do_clean = clean
        if (do_clean) then
            do j = 1, n
                do i = 1, n
                    if (p(i,j) < epsilon(1.0_dp) .and. p(i,j) > -1.0e-12_dp) p(i,j) = 0.0_dp
                    if (p(i,j) > 1.0_dp-epsilon(1.0_dp) .and. p(i,j) < 1.0_dp+1.0e-12_dp) p(i,j) = 1.0_dp
                end do
            end do
        end if
    end function transition_matrix

    function observation_kernel(q, t, obstype) result(kern)
        real(dp), intent(in) :: q(:,:), t
        integer, intent(in) :: obstype
        real(dp), allocatable :: kern(:,:)
        real(dp) :: pii
        integer :: i, j, n
        n = size(q,1)
        allocate(kern(n,n))
        if (obstype == obs_exact) then
            do i = 1, n
                pii = exp(q(i,i)*t)
                do j = 1, n
                    if (i == j) then
                        kern(i,j) = pii
                    else
                        kern(i,j) = pii*q(i,j)
                    end if
                end do
            end do
        else
            kern = transition_matrix(q,t)
        end if
    end function observation_kernel

    function death_transition_density(from_state, to_state, p, q) result(v)
        integer, intent(in) :: from_state, to_state
        real(dp), intent(in) :: p(:,:), q(:,:)
        real(dp) :: v
        integer :: j, n
        n = size(q,1)
        if (from_state == to_state) then
            v = 1.0_dp
            return
        end if
        v = 0.0_dp
        do j = 1, n
            if (j /= to_state) v = v + p(from_state,j)*q(j,to_state)
        end do
    end function death_transition_density

    subroutine transition_derivatives(q, dq, t, dpmat, obstype)
        ! dq(:,:,k) = derivative of Q with respect to parameter k.
        ! For panel data use a block-exponential Frechet derivative, robust to repeated eigenvalues.
        real(dp), intent(in) :: q(:,:), dq(:,:,:), t
        real(dp), allocatable, intent(out) :: dpmat(:,:,:)
        integer, intent(in), optional :: obstype
        real(dp), allocatable :: f(:,:), l(:,:)
        integer :: n, np, k, i, j, ot
        n = size(q,1); np = size(dq,3)
        if (size(q,2) /= n .or. size(dq,1) /= n .or. size(dq,2) /= n) &
            error stop "transition_derivatives: shape mismatch"
        allocate(dpmat(n,n,np))
        ot = obs_panel
        if (present(obstype)) ot = obstype
        if (ot == obs_exact) then
            do k = 1, np
                do i = 1, n
                    do j = 1, n
                        if (i == j) then
                            dpmat(i,j,k) = dq(i,i,k)*t*exp(q(i,i)*t)
                        else
                            dpmat(i,j,k) = exp(q(i,i)*t)*(dq(i,j,k) + dq(i,i,k)*q(i,j)*t)
                        end if
                    end do
                end do
            end do
        else
            do k = 1, np
                call expm_frechet(q*t, dq(:,:,k)*t, f, l)
                dpmat(:,:,k) = l
            end do
        end if
    end subroutine transition_derivatives

    function ctmc_aggregate_minus2loglik(from_state,to_state,timelag,nocc,q,obstype) result(nll2)
        integer, intent(in) :: from_state(:),to_state(:),nocc(:)
        real(dp), intent(in) :: timelag(:),q(:,:,:)
        integer, intent(in), optional :: obstype(:)
        real(dp) :: nll2,pr
        real(dp), allocatable :: p(:,:)
        integer :: i,k,ot,nrec
        nrec=size(from_state)
        if(size(to_state)/=nrec.or.size(timelag)/=nrec.or.size(nocc)/=nrec) error stop "ctmc_aggregate_minus2loglik: size mismatch"
        if(size(q,3)/=1.and.size(q,3)/=nrec) error stop "ctmc_aggregate_minus2loglik: q third dimension"
        nll2=0.0_dp
        do i=1,nrec
            k=merge(1,i,size(q,3)==1); ot=obs_panel; if(present(obstype)) ot=obstype(i)
            p=observation_kernel(q(:,:,k),timelag(i),merge(obs_exact,obs_panel,ot==obs_exact))
            if(ot==obs_death) then
                pr=death_transition_density(from_state(i),to_state(i),p,q(:,:,k))
            else
                pr=p(from_state(i),to_state(i))
            end if
            if(pr<=0.0_dp) then; nll2=huge(1.0_dp); return; end if
            nll2=nll2-2.0_dp*real(nocc(i),dp)*log(pr)
        end do
    end function ctmc_aggregate_minus2loglik

    function ctmc_censored_minus2loglik(times,q,allowed_state,obstype,death_state) result(nll2)
        ! allowed_state(i,k) says true state i is compatible with observation k.
        real(dp), intent(in) :: times(:),q(:,:,:)
        logical, intent(in) :: allowed_state(:,:)
        integer, intent(in), optional :: obstype(:),death_state(:)
        real(dp) :: nll2,sc
        real(dp), allocatable :: a(:),an(:),p(:,:)
        integer :: n,nobs,k,i,j,qk,ot,ds
        n=size(q,1); nobs=size(times)
        if(size(allowed_state,1)/=n.or.size(allowed_state,2)/=nobs) error stop "ctmc_censored_minus2loglik: shape mismatch"
        allocate(a(n),an(n)); a=0.0_dp
        where(allowed_state(:,1)) a=1.0_dp
        sc=sum(a); if(sc<=0.0_dp) then; nll2=huge(1.0_dp); return; end if
        a=a/sc; nll2=0.0_dp
        do k=2,nobs
            qk=min(k-1,size(q,3)); ot=obs_panel; if(present(obstype)) ot=obstype(k)
            p=observation_kernel(q(:,:,qk),times(k)-times(k-1),merge(obs_exact,obs_panel,ot==obs_exact))
            an=0.0_dp
            if(ot==obs_death) then
                ds=0; if(present(death_state)) ds=death_state(k)
                if(ds<=0) error stop "ctmc_censored_minus2loglik: death_state required"
                do j=1,n
                    if(.not.allowed_state(j,k)) cycle
                    do i=1,n
                        an(j)=an(j)+a(i)*p(i,j)*q(j,ds,qk)
                    end do
                end do
            else
                an=matmul(transpose(p),a)
                where(.not.allowed_state(:,k)) an=0.0_dp
            end if
            sc=sum(an); if(sc<=0.0_dp) then; nll2=huge(1.0_dp); return; end if
            nll2=nll2-2.0_dp*log(sc); a=an/sc
        end do
    end function ctmc_censored_minus2loglik

    function ctmc_gradient(states,times,q,dq,obstype) result(g)
        integer, intent(in) :: states(:)
        real(dp), intent(in) :: times(:),q(:,:),dq(:,:,:)
        integer, intent(in), optional :: obstype(:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: p(:,:),dpm(:,:,:)
        real(dp) :: pr,dpr,dt
        integer :: np,k,r,s,j,ot
        np=size(dq,3); allocate(g(np)); g=0.0_dp
        do k=2,size(states)
            r=states(k-1); s=states(k); dt=times(k)-times(k-1); ot=obs_panel; if(present(obstype)) ot=obstype(k)
            p=observation_kernel(q,dt,merge(obs_exact,obs_panel,ot==obs_exact))
            call transition_derivatives(q,dq,dt,dpm,merge(obs_exact,obs_panel,ot==obs_exact))
            if(ot==obs_death) then
                pr=death_transition_density(r,s,p,q)
                do j=1,np
                    dpr=death_derivative(r,s,p,q,dpm(:,:,j),dq(:,:,j))
                    g(j)=g(j)-2.0_dp*dpr/pr
                end do
            else
                pr=p(r,s)
                do j=1,np; g(j)=g(j)-2.0_dp*dpm(r,s,j)/pr; end do
            end if
        end do
    end function ctmc_gradient

    pure function death_derivative(r,s,p,q,dpq,dq) result(v)
        integer,intent(in)::r,s
        real(dp),intent(in)::p(:,:),q(:,:),dpq(:,:),dq(:,:)
        real(dp)::v
        integer::j
        if(r==s) then; v=0.0_dp; return; end if
        v=0.0_dp
        do j=1,size(q,1)
            if(j/=s) v=v+dpq(r,j)*q(j,s)+p(r,j)*dq(j,s)
        end do
    end function death_derivative

    function ctmc_minus2loglik(states, times, q, obstype) result(nll2)
        integer, intent(in) :: states(:)
        real(dp), intent(in) :: times(:), q(:,:)
        integer, intent(in), optional :: obstype(:)
        real(dp) :: nll2, pr, dt
        real(dp), allocatable :: p(:,:)
        integer :: i, ot
        if (size(states) /= size(times)) error stop "ctmc_minus2loglik: states/times length mismatch"
        nll2 = 0.0_dp
        do i = 2, size(states)
            dt = times(i)-times(i-1)
            if (dt < 0.0_dp) error stop "ctmc_minus2loglik: times must be ordered"
            ot = obs_panel
            if (present(obstype)) ot = obstype(i)
            p = observation_kernel(q,dt,merge(obs_exact,obs_panel,ot==obs_exact))
            if (ot == obs_death) then
                pr = death_transition_density(states(i-1),states(i),p,q)
            else
                pr = p(states(i-1),states(i))
            end if
            if (pr <= 0.0_dp) then
                nll2 = huge(1.0_dp)
                return
            end if
            nll2 = nll2 - 2.0_dp*log(pr)
        end do
    end function ctmc_minus2loglik

    function piecewise_transition_matrix(q, change_times, from_time, to_time) result(p)
        ! q(:,:,k) applies from change_times(k) until change_times(k+1).
        real(dp), intent(in) :: q(:,:,:), change_times(:), from_time, to_time
        real(dp), allocatable :: p(:,:), pk(:,:)
        real(dp) :: left, right
        integer :: n, k
        n = size(q,1)
        if (size(q,2) /= n .or. size(q,3) /= size(change_times)) error stop "piecewise_transition_matrix: shape mismatch"
        if (to_time < from_time) error stop "piecewise_transition_matrix: to_time < from_time"
        p = eye(n)
        do k = 1, size(change_times)
            left = max(from_time,change_times(k))
            if (k < size(change_times)) then
                right = min(to_time,change_times(k+1))
            else
                right = to_time
            end if
            if (right > left) then
                pk = transition_matrix(q(:,:,k),right-left)
                p = matmul(p,pk)
            end if
            if (right >= to_time) exit
        end do
    end function piecewise_transition_matrix

    function qmatrix_covariates(base_q, beta, x, transition_index) result(q)
        ! Log-linear covariate model used by msm: q_ij(x)=q_ij(0)*exp(beta_k' x).
        real(dp), intent(in) :: base_q(:,:), beta(:,:), x(:)
        integer, intent(in) :: transition_index(:,:)
        real(dp), allocatable :: q(:,:)
        integer :: n, i, j, k
        n = size(base_q,1)
        if (size(base_q,2) /= n .or. any(shape(transition_index) /= [n,n])) &
            error stop "qmatrix_covariates: shape mismatch"
        if (size(beta,2) /= size(x)) error stop "qmatrix_covariates: beta/x mismatch"
        q = 0.0_dp
        do i = 1, n
            do j = 1, n
                if (i /= j .and. transition_index(i,j) > 0) then
                    k = transition_index(i,j)
                    if (k > size(beta,1)) error stop "qmatrix_covariates: invalid transition index"
                    q(i,j) = base_q(i,j)*exp(dot_product(beta(k,:),x))
                end if
            end do
            q(i,i) = -sum(q(i,:))
        end do
    end function qmatrix_covariates

    function sojourn_times(q) result(tau)
        real(dp), intent(in) :: q(:,:)
        real(dp), allocatable :: tau(:)
        integer :: i, n
        n = size(q,1); allocate(tau(n))
        do i = 1, n
            if (q(i,i) < 0.0_dp) then
                tau(i) = -1.0_dp/q(i,i)
            else
                tau(i) = huge(1.0_dp)
            end if
        end do
    end function sojourn_times

    function next_state_probabilities(q) result(r)
        real(dp), intent(in) :: q(:,:)
        real(dp), allocatable :: r(:,:)
        integer :: i, n
        n = size(q,1); allocate(r(n,n)); r = 0.0_dp
        do i = 1, n
            if (q(i,i) < 0.0_dp) r(i,:) = q(i,:)/(-q(i,i))
            r(i,i) = 0.0_dp
        end do
    end function next_state_probabilities

    function absorbing_states(q, tol) result(mask)
        real(dp), intent(in) :: q(:,:)
        real(dp), intent(in), optional :: tol
        logical, allocatable :: mask(:)
        real(dp) :: eps
        integer :: i, n
        eps = 1.0e-12_dp; if (present(tol)) eps = tol
        n = size(q,1); allocate(mask(n))
        do i = 1, n
            mask(i) = sum(abs(q(i,:))) <= eps
        end do
    end function absorbing_states

    function transient_states(q, tol) result(mask)
        real(dp), intent(in) :: q(:,:)
        real(dp), intent(in), optional :: tol
        logical, allocatable :: mask(:), ab(:)
        allocate(ab(size(q,1)))
        if (present(tol)) then
            ab = absorbing_states(q,tol)
        else
            ab = absorbing_states(q)
        end if
        allocate(mask(size(ab))); mask = .not. ab
    end function transient_states

    function expected_total_time(q, horizon) result(et)
        ! et(i,j) = expected time in state j over [0,horizon], starting in i.
        ! If horizon absent, returns fundamental matrix for transient states and zero elsewhere.
        real(dp), intent(in) :: q(:,:)
        real(dp), intent(in), optional :: horizon
        real(dp), allocatable :: et(:,:), z(:,:), ez(:,:), qt(:,:), rhs(:,:), nt(:,:)
        logical, allocatable :: tr(:)
        integer, allocatable :: ind(:)
        integer :: n, ntc, i, a, b
        n = size(q,1)
        allocate(et(n,n)); et = 0.0_dp
        if (present(horizon)) then
            if (horizon < 0.0_dp) error stop "expected_total_time: negative horizon"
            allocate(z(2*n,2*n)); z = 0.0_dp
            z(1:n,1:n) = q
            z(1:n,n+1:2*n) = eye(n)
            ez = expm(z*horizon)
            et = ez(1:n,n+1:2*n)
            return
        end if
        tr = transient_states(q); ntc = count(tr)
        if (ntc == 0) return
        allocate(ind(ntc)); a = 0
        do i = 1, n
            if (tr(i)) then; a=a+1; ind(a)=i; end if
        end do
        allocate(qt(ntc,ntc),rhs(ntc,ntc)); rhs = eye(ntc)
        do a=1,ntc; do b=1,ntc; qt(a,b)=-q(ind(a),ind(b)); end do; end do
        nt = solve_linear(qt,rhs)
        do a=1,ntc; do b=1,ntc; et(ind(a),ind(b))=nt(a,b); end do; end do
    end function expected_total_time

    function expected_first_passage(q, target) result(tau)
        ! Unconditional mean first-passage time.  States whose eventual
        ! passage probability is less than one have infinite mean passage
        ! time and are reported as HUGE().
        real(dp), intent(in) :: q(:,:)
        logical, intent(in) :: target(:)
        real(dp), allocatable :: tau(:), qt(:,:), rhs(:,:), sol(:,:), hit(:)
        integer, allocatable :: ind(:)
        logical, allocatable :: sure(:)
        integer :: n, m, i, a, b
        n=size(q,1)
        if(size(q,2)/=n .or. size(target)/=n) &
            error stop "expected_first_passage: dimension mismatch"
        allocate(tau(n)); tau=0.0_dp
        hit=eventual_passage_probability(q,target)
        allocate(sure(n)); sure=target .or. hit >= 1.0_dp-1.0e-8_dp
        where(.not.sure) tau=huge(1.0_dp)
        m=count(sure .and. .not.target)
        if(m==0) return
        allocate(ind(m)); a=0
        do i=1,n
            if(sure(i) .and. .not.target(i)) then
                a=a+1; ind(a)=i
            end if
        end do
        allocate(qt(m,m),rhs(m,1)); rhs(:,1)=1.0_dp
        do a=1,m
            do b=1,m
                qt(a,b)=-q(ind(a),ind(b))
            end do
        end do
        sol=solve_linear(qt,rhs)
        do a=1,m
            tau(ind(a))=sol(a,1)
        end do
    end function expected_first_passage

    function eventual_passage_probability(q, target) result(h)
        ! Hitting probabilities for a target set.  Restrict the linear
        ! system to states from which a directed path to the target exists;
        ! this avoids a singular system when the chain contains unrelated
        ! absorbing or recurrent classes.
        real(dp), intent(in) :: q(:,:)
        logical, intent(in) :: target(:)
        real(dp), allocatable :: h(:), qt(:,:), rhs(:,:), sol(:,:)
        integer, allocatable :: ind(:)
        logical, allocatable :: reach(:)
        logical :: changed
        integer :: n,m,i,j,a,b
        real(dp), parameter :: edge_tol=100.0_dp*epsilon(1.0_dp)
        n=size(q,1)
        if(size(q,2)/=n .or. size(target)/=n) &
            error stop "eventual_passage_probability: dimension mismatch"
        allocate(h(n)); h=0.0_dp; where(target) h=1.0_dp
        allocate(reach(n)); reach=target
        do
            changed=.false.
            do i=1,n
                if(reach(i)) cycle
                do j=1,n
                    if(i/=j .and. q(i,j)>edge_tol .and. reach(j)) then
                        reach(i)=.true.; changed=.true.; exit
                    end if
                end do
            end do
            if(.not.changed) exit
        end do
        m=count(reach .and. .not.target)
        if(m==0) return
        allocate(ind(m)); a=0
        do i=1,n
            if(reach(i) .and. .not.target(i)) then
                a=a+1; ind(a)=i
            end if
        end do
        allocate(qt(m,m),rhs(m,1)); rhs=0.0_dp
        do a=1,m
            do b=1,m
                qt(a,b)=-q(ind(a),ind(b))
            end do
            do j=1,n
                if(target(j)) rhs(a,1)=rhs(a,1)+q(ind(a),j)
            end do
        end do
        sol=solve_linear(qt,rhs)
        do a=1,m
            h(ind(a))=min(1.0_dp,max(0.0_dp,sol(a,1)))
        end do
    end function eventual_passage_probability

    function passage_probability_matrix(q,t) result(pp)
        real(dp), intent(in) :: q(:,:),t
        real(dp), allocatable :: pp(:,:),qr(:,:),p(:,:)
        integer :: n,j
        n=size(q,1); allocate(pp(n,n))
        do j=1,n
            qr=q; qr(j,:)=0.0_dp
            p=transition_matrix(qr,t)
            pp(:,j)=p(:,j)
        end do
    end function passage_probability_matrix

    function expected_visits(q,horizon) result(ev)
        real(dp), intent(in) :: q(:,:)
        real(dp), intent(in), optional :: horizon
        real(dp), allocatable :: ev(:,:),et(:,:),r(:,:)
        integer :: i,n
        n=size(q,1); allocate(et(n,n))
        if(present(horizon)) then
            et=expected_total_time(q,horizon)
        else
            et=expected_total_time(q)
        end if
        r=q
        do i=1,n; r(i,i)=0.0_dp; end do
        ev=matmul(et,r)
    end function expected_visits

    function state_prevalence(initp,q,times) result(prev)
        real(dp), intent(in) :: initp(:),q(:,:),times(:)
        real(dp), allocatable :: prev(:,:),p(:,:)
        integer :: k,n
        n=size(q,1); if(size(initp)/=n) error stop "state_prevalence: initp size"
        allocate(prev(n,size(times)))
        do k=1,size(times)
            p=transition_matrix(q,times(k)); prev(:,k)=matmul(transpose(p),initp)
        end do
    end function state_prevalence
end module msm_ctmc
