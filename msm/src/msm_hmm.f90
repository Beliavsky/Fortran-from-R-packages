! SPDX-License-Identifier: GPL-2.0-or-later
! Forward/backward likelihood and Viterbi algorithms translated from msm lik.c.
module msm_hmm
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use msm_kinds, only : dp
    use msm_ctmc, only : obs_panel, obs_exact, obs_death, observation_kernel, transition_derivatives
    use msm_emissions, only : emission_model, emission_probability
    implicit none
    private
    type, public :: hmm_result
        real(dp) :: minus2loglik = 0.0_dp
        real(dp), allocatable :: filtered(:,:)
        real(dp), allocatable :: smoothed(:,:)
        integer, allocatable :: viterbi(:)
    end type hmm_result
    public :: hmm_forward_backward, hmm_minus2loglik, hmm_q_gradient
contains
    function observation_probabilities(obscol, models, true_state, allowed) result(pout)
        real(dp), intent(in) :: obscol(:)
        type(emission_model), intent(in) :: models(:,:)
        integer, intent(in) :: true_state
        logical, intent(in), optional :: allowed(:)
        real(dp), allocatable :: pout(:)
        integer :: nout,nst,i,k
        nout=size(models,1); nst=size(models,2)
        if(size(obscol)/=nout) error stop "observation_probabilities: outcome count mismatch"
        allocate(pout(nst)); pout=1.0_dp
        do i=1,nst
            if(true_state>0 .and. i/=true_state) then
                pout(i)=0.0_dp
                cycle
            end if
            if(present(allowed)) then
                if(.not.allowed(i)) then
                    pout(i)=0.0_dp
                    cycle
                end if
            end if
            do k=1,nout
                if(.not.ieee_is_nan(obscol(k))) pout(i)=pout(i)*emission_probability(models(k,i),obscol(k))
            end do
        end do
    end function observation_probabilities

    subroutine hmm_forward_backward(q, times, obs, models, initp, result, obstype, true_state, allowed_state, death_state)
        ! q(:,:,k) applies on interval times(k) -> times(k+1). A singleton third dimension is reused.
        real(dp), intent(in) :: q(:,:,:), times(:), obs(:,:), initp(:)
        type(emission_model), intent(in) :: models(:,:)
        type(hmm_result), intent(out) :: result
        integer, intent(in), optional :: obstype(:), true_state(:), death_state(:)
        logical, intent(in), optional :: allowed_state(:,:)
        real(dp), allocatable :: alpha(:,:), beta(:,:), scale(:), p(:,:), pout(:), tmp(:), logdelta(:,:), best(:)
        integer, allocatable :: psi(:,:)
        real(dp) :: s, val, loglik, bestval
        integer :: nst,nobs,nout,i,j,k,ot,ts,ds,qk,arg
        nst=size(q,1); nobs=size(times); nout=size(obs,1)
        if(size(q,2)/=nst .or. size(q,3)<1) error stop "hmm_forward_backward: q shape"
        if(size(obs,2)/=nobs .or. size(models,1)/=nout .or. size(models,2)/=nst) error stop "hmm_forward_backward: shape"
        if(size(initp)/=nst) error stop "hmm_forward_backward: initp size"
        allocate(alpha(nst,nobs),beta(nst,nobs),scale(nobs),psi(nst,nobs),logdelta(nst,nobs),pout(nst))
        ts=0; if(present(true_state)) ts=true_state(1)
        if(present(allowed_state)) then
            pout=observation_probabilities(obs(:,1),models,ts,allowed_state(:,1))
        else
            pout=observation_probabilities(obs(:,1),models,ts)
        end if
        alpha(:,1)=initp*pout
        scale(1)=sum(alpha(:,1))
        if(scale(1)<=0.0_dp) then; result%minus2loglik=huge(1.0_dp); return; end if
        alpha(:,1)=alpha(:,1)/scale(1); loglik=log(scale(1))
        do i=1,nst
            if(initp(i)>0.0_dp .and. pout(i)>0.0_dp) then
                logdelta(i,1)=log(initp(i))+log(pout(i))
            else
                logdelta(i,1)=-huge(1.0_dp)
            end if
            psi(i,1)=0
        end do
        do k=2,nobs
            qk=min(k-1,size(q,3)); ot=obs_panel; if(present(obstype)) ot=obstype(k)
            p=observation_kernel(q(:,:,qk),times(k)-times(k-1),merge(obs_exact,obs_panel,ot==obs_exact))
            ts=0; if(present(true_state)) ts=true_state(k)
            if(present(allowed_state)) then
                pout=observation_probabilities(obs(:,k),models,ts,allowed_state(:,k))
            else
                pout=observation_probabilities(obs(:,k),models,ts)
            end if
            allocate(tmp(nst)); tmp=0.0_dp
            if(ot==obs_death) then
                ds=0; if(present(death_state)) ds=death_state(k)
                if(ds<=0 .and. ts>0) ds=ts
                if(ds<=0) error stop "hmm_forward_backward: death observation needs death_state or true_state"
                do j=1,nst
                    do i=1,nst
                        tmp(j)=tmp(j)+alpha(i,k-1)*p(i,j)*q(j,ds,qk)
                    end do
                end do
            else
                tmp=matmul(transpose(p),alpha(:,k-1))*pout
            end if
            scale(k)=sum(tmp)
            if(scale(k)<=0.0_dp) then; result%minus2loglik=huge(1.0_dp); return; end if
            alpha(:,k)=tmp/scale(k); loglik=loglik+log(scale(k)); deallocate(tmp)
            do j=1,nst
                bestval=-huge(1.0_dp); arg=1
                do i=1,nst
                    if(p(i,j)>0.0_dp .and. logdelta(i,k-1)>-huge(1.0_dp)/2) then
                        val=logdelta(i,k-1)+log(p(i,j))
                        if(ot==obs_death) then
                            ds=0; if(present(death_state)) ds=death_state(k); if(ds<=0 .and. ts>0) ds=ts
                            if(q(j,ds,qk)>0.0_dp) val=val+log(q(j,ds,qk)); if(q(j,ds,qk)<=0.0_dp) val=-huge(1.0_dp)
                        else
                            if(pout(j)>0.0_dp) val=val+log(pout(j)); if(pout(j)<=0.0_dp) val=-huge(1.0_dp)
                        end if
                        if(val>bestval) then; bestval=val; arg=i; end if
                    end if
                end do
                logdelta(j,k)=bestval; psi(j,k)=arg
            end do
        end do
        result%minus2loglik=-2.0_dp*loglik
        result%filtered=alpha
        allocate(result%smoothed(nst,nobs)); beta(:,nobs)=1.0_dp; result%smoothed(:,nobs)=alpha(:,nobs)
        do k=nobs-1,1,-1
            qk=min(k,size(q,3)); ot=obs_panel; if(present(obstype)) ot=obstype(k+1)
            p=observation_kernel(q(:,:,qk),times(k+1)-times(k),merge(obs_exact,obs_panel,ot==obs_exact))
            ts=0; if(present(true_state)) ts=true_state(k+1)
            if(present(allowed_state)) then
                pout=observation_probabilities(obs(:,k+1),models,ts,allowed_state(:,k+1))
            else
                pout=observation_probabilities(obs(:,k+1),models,ts)
            end if
            if(ot==obs_death) then
                ds=0; if(present(death_state)) ds=death_state(k+1); if(ds<=0 .and. ts>0) ds=ts
                do i=1,nst
                    beta(i,k)=0.0_dp
                    do j=1,nst
                        beta(i,k)=beta(i,k)+p(i,j)*q(j,ds,qk)*beta(j,k+1)
                    end do
                end do
            else
                beta(:,k)=matmul(p,pout*beta(:,k+1))
            end if
            beta(:,k)=beta(:,k)/scale(k+1)
            result%smoothed(:,k)=alpha(:,k)*beta(:,k)
            s=sum(result%smoothed(:,k)); if(s>0.0_dp) result%smoothed(:,k)=result%smoothed(:,k)/s
        end do
        allocate(result%viterbi(nobs)); best=logdelta(:,nobs); arg=maxloc(best,dim=1); result%viterbi(nobs)=arg
        do k=nobs,2,-1
            result%viterbi(k-1)=psi(result%viterbi(k),k)
        end do
    end subroutine hmm_forward_backward

    function hmm_minus2loglik(q,times,obs,models,initp) result(v)
        real(dp), intent(in) :: q(:,:),times(:),obs(:,:),initp(:)
        type(emission_model), intent(in) :: models(:,:)
        real(dp) :: v
        real(dp), allocatable :: q3(:,:,:)
        type(hmm_result) :: r
        allocate(q3(size(q,1),size(q,2),1)); q3(:,:,1)=q
        call hmm_forward_backward(q3,times,obs,models,initp,r)
        v=r%minus2loglik
    end function hmm_minus2loglik

    function hmm_q_gradient(q,dq,times,obs,models,initp) result(g)
        ! Analytic score for a constant Q and perturbations dq(:,:,p).
        ! This implements the forward sensitivity recursion underlying msm's hidden-model derivative.
        real(dp), intent(in) :: q(:,:),dq(:,:,:),times(:),obs(:,:),initp(:)
        type(emission_model), intent(in) :: models(:,:)
        real(dp), allocatable :: g(:)
        real(dp), allocatable :: a(:),an(:),da(:,:),dan(:,:),p(:,:),dp3(:,:,:),pout(:),num(:)
        real(dp) :: sc, dsc, loglik
        integer :: nst,nobs,np,k,j
        nst=size(q,1); nobs=size(times); np=size(dq,3)
        allocate(g(np),a(nst),an(nst),da(nst,np),dan(nst,np),num(np)); g=0.0_dp; da=0.0_dp
        pout=observation_probabilities(obs(:,1),models,0); a=initp*pout; sc=sum(a)
        if(sc<=0.0_dp) then; g=huge(1.0_dp); return; end if
        a=a/sc; loglik=log(sc)
        do k=2,nobs
            p=observation_kernel(q,times(k)-times(k-1),obs_panel)
            call transition_derivatives(q,dq,times(k)-times(k-1),dp3,obs_panel)
            pout=observation_probabilities(obs(:,k),models,0)
            an=matmul(transpose(p),a)*pout
            do j=1,np
                dan(:,j)=(matmul(transpose(dp3(:,:,j)),a)+matmul(transpose(p),da(:,j)))*pout
            end do
            sc=sum(an); if(sc<=0.0_dp) then; g=huge(1.0_dp); return; end if
            do j=1,np
                dsc=sum(dan(:,j)); g(j)=g(j)+dsc/sc
                dan(:,j)=(dan(:,j)*sc-an*dsc)/(sc*sc)
            end do
            a=an/sc; da=dan; loglik=loglik+log(sc)
        end do
        g=-2.0_dp*g
    end function hmm_q_gradient
end module msm_hmm
