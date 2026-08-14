! Conditional uniform graph and QAP tests translated from R/sna gtest.R.
! Upstream copyright (C) Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_testing
    use sna_kinds, only : dp
    use sna_types, only : qap_result
    use sna_graph, only : dyad_census
    use sna_random, only : rgraph, rgnm, rguman
    use sna_permutation, only : rmperm
    implicit none
    private

    abstract interface
        function graph_statistic(g) result(v)
            import dp
            real(dp),intent(in)::g(:,:)
            real(dp)::v
        end function graph_statistic
        function graph_stack_statistic(g) result(v)
            import dp
            real(dp),intent(in)::g(:,:,:)
            real(dp)::v
        end function graph_stack_statistic
    end interface

    public :: graph_statistic, graph_stack_statistic, cug_test, qap_test
contains

    function cug_test(dat,fun,mode,cmode,diag,reps) result(out)
        real(dp),intent(in)::dat(:,:)
        procedure(graph_statistic)::fun
        character(len=*),intent(in),optional::mode,cmode
        logical,intent(in),optional::diag
        integer,intent(in),optional::reps
        type(qap_result)::out
        character(len=16)::md,cm
        logical::dg
        integer::nr,r,n,m
        real(dp)::dc(3)
        real(dp),allocatable::gr(:,:),prob(:,:)
        md='digraph'
        if(present(mode))md=trim(mode)
        cm='size'
        if(present(cmode))cm=trim(cmode)
        dg=.false.
        if(present(diag))dg=diag
        nr=1000
        if(present(reps))nr=reps
        n=size(dat,1)
        out%observed=fun(dat)
        allocate(out%simulated(nr))
        select case(trim(cm))
        case('edges')
            if(trim(md)=='graph')then
            m=count_lower_edges(dat,dg)
            else
            m=count_directed_edges(dat,dg)
            end if
        case('dyad.census')
            dc=dyad_census(dat)
        end select
        do r=1,nr
            select case(trim(cm))
            case('edges')
            gr=rgnm(n,m,md,dg)
            case('dyad.census')
            gr=rguman(n,dc(1),dc(2),dc(3),'exact')
            case default
                allocate(prob(n,n))
                prob=0.5_dp
                gr=rgraph(n,prob,md,dg)
                deallocate(prob)
            end select
            out%simulated(r)=fun(gr)
        end do
        call test_pvalues(out)
    end function cug_test

    function qap_test(dat,fun,reps) result(out)
        real(dp),intent(in)::dat(:,:,:)
        procedure(graph_stack_statistic)::fun
        integer,intent(in),optional::reps
        type(qap_result)::out
        real(dp),allocatable::gr(:,:,:)
        integer::nr,r,k
        nr=1000
        if(present(reps))nr=reps
        out%observed=fun(dat)
        allocate(out%simulated(nr),gr(size(dat,1),size(dat,2),size(dat,3)))
        do r=1,nr
            do k=1,size(dat,1)
            gr(k,:,:)=rmperm(dat(k,:,:))
            end do
            out%simulated(r)=fun(gr)
        end do
        call test_pvalues(out)
    end function qap_test

    subroutine test_pvalues(out)
        type(qap_result),intent(inout)::out
        integer::n
        n=size(out%simulated)
        out%p_lower=real(count(out%simulated<=out%observed),dp)/real(n,dp)
        out%p_upper=real(count(out%simulated>=out%observed),dp)/real(n,dp)
        out%p_two_sided=real(count(abs(out%simulated)>=abs(out%observed)),dp)/real(n,dp)
    end subroutine test_pvalues

    integer function count_directed_edges(g,diag) result(m)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::diag
        integer::i,j
        m=0
        do i=1,size(g,1)
        do j=1,size(g,2)
        if(.not.diag.and.i==j)cycle
        if(g(i,j)>0)m=m+1
        end do
        end do
    end function count_directed_edges

    integer function count_lower_edges(g,diag) result(m)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::diag
        integer::i,j
        m=0
        do i=1,size(g,1)
        do j=1,i
        if(.not.diag.and.i==j)cycle
        if(g(i,j)>0.or.g(j,i)>0)m=m+1
        end do
        end do
    end function count_lower_edges
end module sna_testing
