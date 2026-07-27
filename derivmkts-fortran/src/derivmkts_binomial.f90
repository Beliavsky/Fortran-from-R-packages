! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module derivmkts_binomial
    use derivmkts_kinds, only: dp
    use derivmkts_types, only: binomial_result
    implicit none
    private
    public :: binomopt, binormsdist_discrete
contains

    function binomopt(s,k,v,r,tt,d,nstep,american,putopt,specifyupdn,crr,jarrowrudd,up_in,dn_in,returntrees) result(out)
        real(dp), intent(in) :: s,k,v,r,tt,d
        integer, intent(in), optional :: nstep
        logical, intent(in), optional :: american,putopt,specifyupdn,crr,jarrowrudd,returntrees
        real(dp), intent(in), optional :: up_in,dn_in
        type(binomial_result) :: out
        integer :: n,i,j
        real(dp) :: payoffmult,vnc,eps_move,comb
        logical :: am,isput,spec,usecrr,usejr,trees
        n=10
        if(present(nstep))n=nstep
        am=.true.
        if(present(american))am=american
        isput=.false.
        if(present(putopt))isput=putopt
        spec=.false.
        if(present(specifyupdn))spec=specifyupdn
        usecrr=.false.
        if(present(crr))usecrr=crr
        usejr=.false.
        if(present(jarrowrudd))usejr=jarrowrudd
        trees=.false.
        if(present(returntrees))trees=returntrees
        if(n<1 .or. s<=0.0_dp .or. k<0.0_dp .or. tt<=0.0_dp) return
        out%h=tt/real(n,dp)
        if(spec) then
            if(.not.present(up_in) .or. .not.present(dn_in)) return
            out%up=up_in
            out%down=dn_in
        else if(usecrr) then
            out%up=exp(sqrt(out%h)*v)
            out%down=exp(-sqrt(out%h)*v)
        else if(usejr) then
            out%up=exp((r-d-0.5_dp*v*v)*out%h+sqrt(out%h)*v)
            out%down=exp((r-d-0.5_dp*v*v)*out%h-sqrt(out%h)*v)
        else
            out%up=exp((r-d)*out%h+sqrt(out%h)*v)
            out%down=exp((r-d)*out%h-sqrt(out%h)*v)
        end if
        if(abs(out%up-out%down)<=epsilon(1.0_dp)) return
        out%p=(exp((r-d)*out%h)-out%down)/(out%up-out%down)
        if(out%p<0.0_dp .or. out%p>1.0_dp) return
        allocate(out%stock_tree(n+1,n+1),out%option_tree(n+1,n+1))
        allocate(out%probability_tree(n+1,n+1),out%exercise_tree(n+1,n+1))
        allocate(out%delta_tree(n,n),out%bond_tree(n,n))
        out%stock_tree=0.0_dp
        out%option_tree=0.0_dp
        out%probability_tree=0.0_dp
        out%exercise_tree=.false.
        out%delta_tree=0.0_dp
        out%bond_tree=0.0_dp
        do j=1,n+1
            do i=1,j
                out%stock_tree(i,j)=s*out%up**real(j-i,dp)*out%down**real(i-1,dp)
            end do
        end do
        payoffmult=merge(-1.0_dp,1.0_dp,isput)
        do i=1,n+1
            out%option_tree(i,n+1)=max(payoffmult*(out%stock_tree(i,n+1)-k),0.0_dp)
            out%exercise_tree(i,n+1)=out%option_tree(i,n+1)>0.0_dp
        end do
        do j=n,1,-1
            do i=1,j
                vnc=exp(-r*out%h)*(out%p*out%option_tree(i,j+1)+(1.0_dp-out%p)*out%option_tree(i+1,j+1))
                if(am) then
                    out%option_tree(i,j)=max(vnc,payoffmult*(out%stock_tree(i,j)-k))
                    out%exercise_tree(i,j)=payoffmult*(out%stock_tree(i,j)-k)>vnc .and. &
                        payoffmult*(out%stock_tree(i,j)-k)>0.0_dp
                else
                    out%option_tree(i,j)=vnc
                end if
            end do
        end do
        out%price=out%option_tree(1,1)
        do j=1,n
            do i=1,j
                out%delta_tree(i,j)=exp(-d*out%h)*(out%option_tree(i,j+1)-out%option_tree(i+1,j+1))/ &
                    ((out%up-out%down)*out%stock_tree(i,j))
                out%bond_tree(i,j)=exp(-r*out%h)*(out%up*out%option_tree(i+1,j+1)- &
                    out%down*out%option_tree(i,j+1))/(out%up-out%down)
            end do
        end do
        out%delta=out%delta_tree(1,1)
        if(n>=2) then
            out%gamma=(out%delta_tree(1,2)-out%delta_tree(2,2))/ &
                (out%stock_tree(1,2)-out%stock_tree(2,2))
            eps_move=(out%up*out%down-1.0_dp)*s
            out%theta=(out%option_tree(2,3)-eps_move*out%delta-0.5_dp*eps_move*eps_move*out%gamma- &
                out%option_tree(1,1))/(2.0_dp*out%h)/365.0_dp
        end if
        do j=1,n+1
            do i=1,j
                comb=binomial_coefficient(j-1,i-1)
                out%probability_tree(i,j)=comb*out%p**real(j-i,dp)*(1.0_dp-out%p)**real(i-1,dp)
            end do
        end do
        out%valid=.true.
        if(.not.trees) then
            deallocate(out%stock_tree,out%option_tree,out%probability_tree,out%exercise_tree)
            deallocate(out%delta_tree,out%bond_tree)
        end if
    end function binomopt

    pure real(dp) function binomial_coefficient(n,k) result(c)
        integer,intent(in)::n,k
        integer::i,kk
        if(k<0 .or. k>n) then
            c=0.0_dp
            return
        end if
        kk=min(k,n-k)
        c=1.0_dp
        do i=1,kk
            c=c*real(n-kk+i,dp)/real(i,dp)
        end do
    end function binomial_coefficient

    pure subroutine binormsdist_discrete(n,p,probabilities)
        integer,intent(in)::n
        real(dp),intent(in)::p
        real(dp),intent(out)::probabilities(0:n)
        integer::k
        do k=0,n
            probabilities(k)=binomial_coefficient(n,k)*p**k*(1.0_dp-p)**(n-k)
        end do
    end subroutine binormsdist_discrete
end module derivmkts_binomial
