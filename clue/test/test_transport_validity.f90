program test_transport_validity
    use clue
    implicit none
    real(dp),allocatable::m1(:,:),m2(:,:),s(:)
    real(dp)::v
    integer::stat
    real(dp)::d(4,4)
    m1=membership_from_ids([1,1,2,2])
    m2=membership_from_ids([2,2,1,1])
    v=dissimilarity_mallows(m1,m2,status=stat)
    call check(stat==0 .and. abs(v)<1e-9_dp,'Mallows transport')
    d=reshape([0.0_dp,1.0_dp,5.0_dp,5.0_dp, &
               1.0_dp,0.0_dp,5.0_dp,5.0_dp, &
               5.0_dp,5.0_dp,0.0_dp,1.0_dp, &
               5.0_dp,5.0_dp,1.0_dp,0.0_dp],[4,4])
    v=dissimilarity_accounted_for(m1,d)
    call check(v>0.7_dp,'dissimilarity accounted for')
    s=silhouette_widths([1,1,2,2],d)
    call check(minval(s)>0.7_dp,'silhouette widths')
    print *, 'test_transport_validity: PASS'
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
