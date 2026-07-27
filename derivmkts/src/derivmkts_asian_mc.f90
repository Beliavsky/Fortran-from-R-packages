! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_asian_mc
    use derivmkts_kinds, only: dp
    use derivmkts_rng, only: rng_state, seed_rng, normal_rng
    use derivmkts_math, only: mean_value, sample_sd, covariance, variance
    use derivmkts_types, only: asian_mc_result, option_pair
    use derivmkts_asian_analytic, only: geomavgprice, geomavgstrike
    use derivmkts_black_scholes, only: bscall, bsput
    implicit none
    private
    public :: arithasianmc, geomasianmc, arithavgpricecv
contains

    function arithasianmc(s,k,v,r,tt,d,m,numsim,seed) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d
        integer, intent(in) :: m,numsim
        integer, intent(in), optional :: seed
        type(asian_mc_result) :: out
        real(dp), allocatable :: paths(:,:),avg(:),terminal(:),payoff(:)
        real(dp) :: disc
        call make_paths(s,v,r,tt,d,m,numsim,seed,paths)
        if(.not.allocated(paths))return
        allocate(avg(numsim),terminal(numsim),payoff(numsim))
        avg=sum(paths,dim=2)/real(m,dp)
        terminal=paths(:,m)
        disc=exp(-r*tt)
        payoff=max(avg-k,0.0_dp)
        out%avg_price_call=disc*mean_value(payoff)
        out%sd_avg_price_call=disc*sample_sd(payoff)
        payoff=max(k-avg,0.0_dp)
        out%avg_price_put=disc*mean_value(payoff)
        out%sd_avg_price_put=disc*sample_sd(payoff)
        payoff=max(terminal-avg,0.0_dp)
        out%avg_strike_call=disc*mean_value(payoff)
        out%sd_avg_strike_call=disc*sample_sd(payoff)
        payoff=max(avg-terminal,0.0_dp)
        out%avg_strike_put=disc*mean_value(payoff)
        out%sd_avg_strike_put=disc*sample_sd(payoff)
        payoff=max(terminal-k,0.0_dp)
        out%vanilla_call=disc*mean_value(payoff)
        out%sd_vanilla_call=disc*sample_sd(payoff)
        payoff=max(k-terminal,0.0_dp)
        out%vanilla_put=disc*mean_value(payoff)
        out%sd_vanilla_put=disc*sample_sd(payoff)
    end function arithasianmc

    function geomasianmc(s,k,v,r,tt,d,m,numsim,seed) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d
        integer, intent(in) :: m,numsim
        integer, intent(in), optional :: seed
        type(asian_mc_result) :: out
        real(dp), allocatable :: paths(:,:),avg(:),terminal(:),payoff(:)
        real(dp) :: disc
        integer :: j
        type(option_pair) :: exact_price, exact_strike
        call make_paths(s,v,r,tt,d,m,numsim,seed,paths)
        if(.not.allocated(paths))return
        allocate(avg(numsim),terminal(numsim),payoff(numsim))
        avg=1.0_dp
        do j=1,m
            avg=avg*paths(:,j)
        end do
        avg=avg**(1.0_dp/real(m,dp))
        terminal=paths(:,m)
        disc=exp(-r*tt)
        payoff=max(avg-k,0.0_dp)
        out%avg_price_call=disc*mean_value(payoff)
        out%sd_avg_price_call=disc*sample_sd(payoff)
        payoff=max(k-avg,0.0_dp)
        out%avg_price_put=disc*mean_value(payoff)
        out%sd_avg_price_put=disc*sample_sd(payoff)
        payoff=max(terminal-(k/s)*avg,0.0_dp)
        out%avg_strike_call=disc*mean_value(payoff)
        out%sd_avg_strike_call=disc*sample_sd(payoff)
        payoff=max((k/s)*avg-terminal,0.0_dp)
        out%avg_strike_put=disc*mean_value(payoff)
        out%sd_avg_strike_put=disc*sample_sd(payoff)
        payoff=max(terminal-k,0.0_dp)
        out%vanilla_call=disc*mean_value(payoff)
        out%sd_vanilla_call=disc*sample_sd(payoff)
        payoff=max(k-terminal,0.0_dp)
        out%vanilla_put=disc*mean_value(payoff)
        out%sd_vanilla_put=disc*sample_sd(payoff)
        exact_price=geomavgprice(s,k,v,r,tt,d,m)
        exact_strike=geomavgstrike(s,k,v,r,tt,d,m)
        out%exact_avg_price_call=exact_price%call
        out%exact_avg_price_put=exact_price%put
        out%exact_avg_strike_call=exact_strike%call
        out%exact_avg_strike_put=exact_strike%put
    end function geomasianmc

    function arithavgpricecv(s,k,v,r,tt,d,m,numsim,seed) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d
        integer, intent(in) :: m,numsim
        integer, intent(in), optional :: seed
        type(asian_mc_result) :: out
        real(dp), allocatable :: paths(:,:),aavg(:),gavg(:),ap(:),gp(:),corrected(:)
        real(dp) :: truegeom,disc,varpilot
        type(option_pair) :: exact_price
        integer :: nall,j,npilot
        nall=numsim+250
        npilot=min(250,nall)
        call make_paths(s,v,r,tt,d,m,nall,seed,paths)
        if(.not.allocated(paths))return
        allocate(aavg(nall),gavg(nall),ap(nall),gp(nall),corrected(nall))
        aavg=sum(paths,dim=2)/real(m,dp)
        gavg=1.0_dp
        do j=1,m
            gavg=gavg*paths(:,j)
        end do
        gavg=gavg**(1.0_dp/real(m,dp))
        disc=exp(-r*tt)
        ap=disc*max(aavg-k,0.0_dp)
        gp=disc*max(gavg-k,0.0_dp)
        exact_price=geomavgprice(s,k,v,r,tt,d,m)
        truegeom=exact_price%call
        varpilot=variance(gp(1:npilot))
        if(varpilot>0.0_dp)out%beta=covariance(ap(1:npilot),gp(1:npilot))/varpilot
        corrected=ap+out%beta*(truegeom-gp)
        out%avg_price_call=mean_value(corrected)
        out%sd_avg_price_call=sample_sd(corrected)
    end function arithavgpricecv

    subroutine make_paths(s,v,r,tt,d,m,n,seed,paths)
        real(dp), intent(in) :: s,v,r,tt,d
        integer, intent(in) :: m,n
        integer, intent(in), optional :: seed
        real(dp), allocatable, intent(out) :: paths(:,:)
        type(rng_state) :: rng
        real(dp) :: h,zsum
        integer :: i,j,seedval
        if(m<1 .or. n<1 .or. s<=0.0_dp .or. tt<0.0_dp)return
        seedval=13579
        if(present(seed))seedval=seed
        call seed_rng(rng,seedval)
        allocate(paths(n,m))
        h=tt/real(m,dp)
        do i=1,n
            zsum=0.0_dp
            do j=1,m
                zsum=zsum+normal_rng(rng)
                paths(i,j)=s*exp((r-d-0.5_dp*v*v)*h*real(j,dp)+v*sqrt(h)*zsum)
            end do
        end do
    end subroutine make_paths
end module derivmkts_asian_mc
