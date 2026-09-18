program basic
    use dbscan_api, only : dp, clustering_result, dbscan, ncluster, nnoise
    implicit none

    real(dp) :: x(8, 2)
    type(clustering_result) :: fit
    integer :: status

    x = reshape([ &
        0.0_dp, 0.0_dp, 0.1_dp, 0.0_dp, 0.0_dp, 0.1_dp, &
        3.0_dp, 3.0_dp, 3.1_dp, 3.0_dp, 3.0_dp, 3.1_dp, &
        1.5_dp, 1.5_dp, 6.0_dp, 0.0_dp], shape(x), order = [2, 1])

    call dbscan(x, 0.25_dp, 3, fit, status = status)
    if (status /= 0) error stop 'DBSCAN failed'

    print '(a,*(1x,i0))', 'cluster:', fit%cluster
    print '(a,i0)', 'clusters: ', ncluster(fit)
    print '(a,i0)', 'noise points: ', nnoise(fit)
end program basic
