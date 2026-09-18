program parity_driver
    use dbscan_api
    implicit none

    real(dp) :: x(12, 2)
    type(clustering_result) :: db
    type(hdbscan_result) :: hd
    type(optics_result) :: op
    real(dp), allocatable :: lo(:)
    character(len=256) :: path
    integer :: i
    integer :: status
    integer :: unit

    x = reshape([ &
        0.00_dp, 0.00_dp, 0.10_dp, 0.00_dp, 0.00_dp, 0.12_dp, 0.12_dp, 0.10_dp, 0.20_dp, 0.05_dp, &
        3.00_dp, 3.00_dp, 3.10_dp, 3.00_dp, 3.00_dp, 3.12_dp, 3.12_dp, 3.10_dp, 3.20_dp, 3.05_dp, &
        1.50_dp, 1.50_dp, 6.00_dp, 0.00_dp], shape(x), order = [2, 1])

    path = 'fortran_parity.csv'
    if (command_argument_count() >= 1) call get_command_argument(1, path)
    call dbscan(x, 0.25_dp, 3, db, status = status)
    if (status /= 0) error stop 'DBSCAN failed'
    lo = lof(x, 3, status)
    if (status /= 0) error stop 'LOF failed'
    call hdbscan(x, 3, hd, status)
    if (status /= 0) error stop 'HDBSCAN failed'
    call optics(x, 3, op, eps = 0.5_dp, status = status)
    if (status /= 0) error stop 'OPTICS failed'

    open(newunit = unit, file = trim(path), status = 'replace', action = 'write')
    write(unit, '(a)') 'id,dbscan,lof,hdbscan,hprob,reach,core,predecessor'
    do i = 1, 12
        write(unit, '(i0,a,i0,a,es24.16,a,i0,a,es24.16,a,es24.16,a,es24.16,a,i0)') &
            i, ',', db%cluster(i), ',', lo(i), ',', hd%cluster(i), ',', hd%membership_prob(i), ',', &
            op%reachdist(i), ',', op%coredist(i), ',', op%predecessor(i)
    end do
    close(unit)
end program parity_driver
