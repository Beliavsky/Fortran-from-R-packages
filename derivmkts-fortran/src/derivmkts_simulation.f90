! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_simulation
    use derivmkts_kinds, only: dp
    use derivmkts_math, only: cholesky_lower
    use derivmkts_rng, only: rng_state, seed_rng, normal_rng, poisson_rng
    use derivmkts_types, only: simulation_result
    implicit none
    private
    public :: simprice

    interface simprice
        module procedure simprice_scalar
        module procedure simprice_multi
    end interface simprice

contains

    function simprice_scalar(s0,v,r,tt,d,trials,periods,jump,lambda,alphaj,vj,seed,scalar_v_is_stddev) result(out)
        real(dp), intent(in) :: s0,v,r,tt,d
        integer, intent(in) :: trials,periods
        logical, intent(in), optional :: jump,scalar_v_is_stddev
        real(dp), intent(in), optional :: lambda,alphaj,vj
        integer, intent(in), optional :: seed
        type(simulation_result) :: out
        real(dp) :: cov(1,1),ra(1),da(1),la(1),aa(1),ja(1)
        logical :: isstd
        isstd=.true.
        if(present(scalar_v_is_stddev))isstd=scalar_v_is_stddev
        cov(1,1)=merge(v*v,v,isstd)
        ra(1)=r
        da(1)=d
        la(1)=0.0_dp
        aa(1)=0.0_dp
        ja(1)=0.0_dp
        if(present(lambda))la(1)=lambda
        if(present(alphaj))aa(1)=alphaj
        if(present(vj))ja(1)=vj
        out=simprice_multi(s0,cov,ra,tt,da,trials,periods,jump,la,aa,ja,seed)
    end function simprice_scalar

    function simprice_multi(s0,covariance,r,tt,d,trials,periods,jump,lambda,alphaj,vj,seed) result(out)
        real(dp), intent(in) :: s0,covariance(:,:),r(:),tt,d(:)
        integer, intent(in) :: trials,periods
        logical, intent(in), optional :: jump
        real(dp), intent(in), optional :: lambda(:),alphaj(:),vj(:)
        integer, intent(in), optional :: seed
        type(simulation_result) :: out
        type(rng_state) :: rng
        real(dp), allocatable :: l(:,:),z(:),corr(:),la(:),aa(:),ja(:)
        real(dp) :: h,drift,kappa,jumpfactor
        integer :: nasset,i,j,a,nj,seedval
        logical :: dojump,ok
        nasset=size(covariance,1)
        if(size(covariance,2)/=nasset .or. size(r)/=nasset .or. size(d)/=nasset) return
        if(trials<1 .or. periods<1 .or. tt<0.0_dp .or. s0<=0.0_dp) return
        allocate(l(nasset,nasset),z(nasset),corr(nasset),la(nasset),aa(nasset),ja(nasset))
        call cholesky_lower(covariance,l,ok)
        if(.not.ok) return
        la=0.0_dp
        aa=0.0_dp
        ja=0.0_dp
        if(present(lambda)) then
            if(size(lambda)/=nasset)return
            la=lambda
        end if
        if(present(alphaj)) then
            if(size(alphaj)/=nasset)return
            aa=alphaj
        end if
        if(present(vj)) then
            if(size(vj)/=nasset)return
            ja=vj
        end if
        dojump=.false.
        if(present(jump))dojump=jump
        seedval=1234567
        if(present(seed))seedval=seed
        call seed_rng(rng,seedval)
        allocate(out%price(trials,periods+1,nasset),out%jumps(trials,periods+1,nasset))
        out%price(:,1,:)=s0
        out%jumps=0
        h=tt/real(periods,dp)
        do i=1,trials
            do j=1,periods
                do a=1,nasset
                    z(a)=normal_rng(rng)
                end do
                corr=matmul(l,z)
                do a=1,nasset
                    kappa=exp(aa(a))-1.0_dp
                    drift=(r(a)-d(a)-0.5_dp*covariance(a,a))*h
                    if(dojump)drift=drift-kappa*la(a)*h
                    jumpfactor=0.0_dp
                    nj=0
                    if(dojump) then
                        nj=poisson_rng(rng,la(a)*h)
                        if(nj>0) jumpfactor=(aa(a)-0.5_dp*ja(a)*ja(a))*real(nj,dp)+ &
                            ja(a)*sqrt(real(nj,dp))*normal_rng(rng)
                    end if
                    out%price(i,j+1,a)=out%price(i,j,a)*exp(drift+sqrt(h)*corr(a)+jumpfactor)
                    out%jumps(i,j+1,a)=nj
                end do
            end do
        end do
        out%valid=.true.
    end function simprice_multi

end module derivmkts_simulation
