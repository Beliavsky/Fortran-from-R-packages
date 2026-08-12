program partition_demo
    use clue, only: dp, membership_from_ids, agreement_adjusted_rand, agreement_nmi, &
        dissimilarity_euclidean, consensus_dwh, class_ids_from_membership
    implicit none
    integer :: a(8), b(8)
    real(dp) :: ens(8,3,2)
    real(dp), allocatable :: ma(:,:), mb(:,:), consensus(:,:)
    integer, allocatable :: ids(:)

    a = [1,1,1,2,2,3,3,3]
    b = [2,2,2,3,3,1,1,1]
    ma = membership_from_ids(a,3)
    mb = membership_from_ids(b,3)

    print '(a,f8.5)', 'Adjusted Rand: ', agreement_adjusted_rand(a,b)
    print '(a,f8.5)', 'NMI:           ', agreement_nmi(a,b)
    print '(a,es12.4)', 'Matched Euclidean dissimilarity: ', dissimilarity_euclidean(ma,mb)

    ens(:,:,1) = ma
    ens(:,:,2) = mb
    consensus = consensus_dwh(ens,k=3)
    ids = class_ids_from_membership(consensus)
    print '(a,8(i0,1x))', 'Consensus classes: ', ids
end program partition_demo
