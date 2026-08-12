program test_lsap_partition
    use clue
    implicit none
    real(dp) :: c(3,3), v
    integer, allocatable :: p(:), meet(:), join(:)
    integer :: a(6), b(6)
    real(dp), allocatable :: ma(:,:), mb(:,:)
    c = reshape([4.0_dp,2.0_dp,3.0_dp, 1.0_dp,0.0_dp,2.0_dp, 3.0_dp,5.0_dp,2.0_dp],[3,3])
    call solve_lsap(c,p)
    call check(abs(assignment_cost(c,p)-5.0_dp)<1e-12_dp,'LSAP minimum')
    call solve_lsap(c,p,.true.)
    call check(abs(assignment_cost(c,p)-11.0_dp)<1e-12_dp,'LSAP maximum')
    a=[1,1,2,2,3,3]
    b=[2,2,1,1,3,3]
    call check(abs(agreement_rand(a,b)-1.0_dp)<1e-14_dp,'Rand label invariance')
    call check(abs(agreement_adjusted_rand(a,b)-1.0_dp)<1e-14_dp,'adjusted Rand')
    call check(abs(agreement_nmi(a,b)-1.0_dp)<1e-14_dp,'NMI')
    ma=membership_from_ids(a)
    call check(maxval(abs(membership_from_ids([10,10,42,42])-membership_from_ids([1,1,2,2])))<1e-14_dp, &
        'arbitrary class labels canonicalized')
    mb=membership_from_ids(b)
    call check(dissimilarity_euclidean(ma,mb)<1e-14_dp,'matched membership distance')
    call check(dissimilarity_comemberships(ma,mb)<1e-14_dp,'comembership distance')
    meet=partition_meet([1,1,2,2],[1,2,1,2])
    call check(maxval(meet)==4,'partition meet')
    join=partition_join([1,1,2,2],[1,2,1,2])
    call check(maxval(join)==1,'partition join')
    v=fuzziness_pc(ma,.true.)
    call check(abs(v)<1e-14_dp,'hard PC fuzziness')
    print *, 'test_lsap_partition: PASS'
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
