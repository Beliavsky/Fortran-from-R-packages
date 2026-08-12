program test_consensus_pclust
    use clue
    implicit none
    real(dp) :: ens(6,3,3), x(6,2)
    real(dp), allocatable :: m(:,:), cm(:,:)
    integer, allocatable :: ids(:)
    type(pclust_result) :: pc
    ens=0.0_dp
    ens(:,:,1)=membership_from_ids([1,1,2,2,3,3],3)
    ens(:,:,2)=membership_from_ids([2,2,3,3,1,1],3)
    ens(:,:,3)=membership_from_ids([3,3,1,1,2,2],3)
    m=consensus_dwh(ens,k=3)
    ids=class_ids_from_membership(m)
    cm=co_membership(membership_from_ids(ids))
    call check(maxval(abs(cm-co_membership(ens(:,:,1))))<1e-12_dp,'DWH consensus label matching')
    ids=consensus_hard_euclidean(ens,k=3)
    call check(maxval(abs(co_membership(membership_from_ids(ids))-co_membership(ens(:,:,1))))<1e-12_dp,'hard Euclidean consensus')
    m=consensus_soft_manhattan(ens,k=3,maxiter=10)
    call check(maxval(abs(co_membership(m)-co_membership(ens(:,:,1))))<1e-10_dp,'soft Manhattan consensus LP')
    ids=consensus_hard_manhattan(ens,k=3,maxiter=10)
    call check(maxval(abs(co_membership(membership_from_ids(ids))-co_membership(ens(:,:,1))))<1e-12_dp,'hard Manhattan consensus')
    x=reshape([0.0_dp,0.1_dp,-0.1_dp,5.0_dp,5.1_dp,4.9_dp, 0.0_dp,-0.1_dp,0.1_dp,5.0_dp,4.9_dp,5.1_dp],[6,2])
    pc=pclust_euclidean(x,2,m=2.0_dp,start=reshape([0.0_dp,5.0_dp,0.0_dp,5.0_dp],[2,2]))
    call check(pc%objective<0.2_dp,'fuzzy prototype clustering')
    call check(count(pc%class_ids==pc%class_ids(1))==3,'prototype class split')
    print *, 'test_consensus_pclust: PASS'
contains
    subroutine check(ok,msg)
    logical,intent(in)::ok
    character(*),intent(in)::msg
    if(.not.ok)then
    write(*,*)'FAIL: ',trim(msg)
    error stop 1
    end if
    end subroutine
end program
