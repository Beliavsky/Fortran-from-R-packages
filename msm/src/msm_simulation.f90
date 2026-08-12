! SPDX-License-Identifier: GPL-2.0-or-later
module msm_simulation
    use msm_kinds, only : dp
    use msm_stats, only : rand_uniform, rand_exponential
    use msm_emissions, only : emission_model, simulate_emission
    implicit none
    private
    type, public :: ctmc_path
        real(dp), allocatable :: time(:)
        integer, allocatable :: state(:)
    end type ctmc_path
    public :: simulate_ctmc, simulate_hmm_observations
contains
    subroutine simulate_ctmc(q,start_state,end_time,path,max_events)
        real(dp), intent(in) :: q(:,:),end_time
        integer, intent(in) :: start_state
        type(ctmc_path), intent(out) :: path
        integer, intent(in), optional :: max_events
        real(dp), allocatable :: tt(:)
        integer, allocatable :: ss(:)
        real(dp) :: t,rate,u,cum
        integer :: n,st,j,m,cap,ne
        n=size(q,1); cap=10000; if(present(max_events)) cap=max_events
        allocate(tt(cap+1),ss(cap+1)); ne=1; tt(1)=0.0_dp; ss(1)=start_state
        t=0.0_dp; st=start_state
        do while(t<end_time .and. ne<=cap)
            rate=-q(st,st)
            if(rate<=0.0_dp) exit
            t=t+rand_exponential(rate)
            if(t>end_time) exit
            u=rand_uniform(); cum=0.0_dp; j=st
            do m=1,n
                if(m/=st .and. q(st,m)>0.0_dp) then
                    cum=cum+q(st,m)/rate
                    if(u<=cum) then; j=m; exit; end if
                end if
            end do
            st=j; ne=ne+1; tt(ne)=t; ss(ne)=st
        end do
        allocate(path%time(ne),path%state(ne)); path%time=tt(1:ne); path%state=ss(1:ne)
    end subroutine simulate_ctmc

    subroutine simulate_hmm_observations(states,models,obs)
        integer, intent(in) :: states(:)
        type(emission_model), intent(in) :: models(:,:)
        real(dp), allocatable, intent(out) :: obs(:,:)
        integer :: nout,nobs,k,r
        nout=size(models,1); nobs=size(states); allocate(obs(nout,nobs))
        do k=1,nobs
            do r=1,nout
                obs(r,k)=simulate_emission(models(r,states(k)))
            end do
        end do
    end subroutine simulate_hmm_observations
end module msm_simulation
