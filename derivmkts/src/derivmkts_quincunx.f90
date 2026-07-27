! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_quincunx
    use derivmkts_kinds, only: dp
    use derivmkts_rng, only: rng_state, seed_rng, binomial_rng
    use derivmkts_types, only: quincunx_result
    implicit none
    private
    public :: quincunx
contains
    function quincunx(n,numballs,probright,seed) result(out)
        integer, intent(in), optional :: n,numballs,seed
        real(dp), intent(in), optional :: probright
        type(quincunx_result) :: out
        type(rng_state) :: rng
        integer :: levels,balls,seedval,i,k
        real(dp) :: p
        levels=3
        if(present(n))levels=n
        balls=20
        if(present(numballs))balls=numballs
        p=0.5_dp
        if(present(probright))p=probright
        seedval=24680
        if(present(seed))seedval=seed
        if(levels<0 .or. balls<0 .or. p<0.0_dp .or. p>1.0_dp)return
        call seed_rng(rng,seedval)
        allocate(out%counts(0:levels),out%expected(0:levels))
        out%counts=0
        do i=1,balls
            k=binomial_rng(rng,levels,p)
            out%counts(k)=out%counts(k)+1
        end do
        do k=0,levels
            out%expected(k)=real(balls,dp)*binomial_probability(levels,k,p)
        end do
    end function quincunx

    pure real(dp) function binomial_probability(n,k,p) result(prob)
        integer,intent(in)::n,k
        real(dp),intent(in)::p
        integer::i,kk
        real(dp)::c
        kk=min(k,n-k)
        c=1.0_dp
        do i=1,kk
            c=c*real(n-kk+i,dp)/real(i,dp)
        end do
        prob=c*p**k*(1.0_dp-p)**(n-k)
    end function binomial_probability
end module derivmkts_quincunx
