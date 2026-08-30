program basic
    use proxy, only: dp, proxy_dist_auto, proxy_pack_dist
    implicit none

    real(dp) :: x(4, 2)
    real(dp), allocatable :: distances(:, :)
    real(dp), allocatable :: packed(:)

    x(1, :) = [0.0_dp, 0.0_dp]
    x(2, :) = [1.0_dp, 0.0_dp]
    x(3, :) = [1.0_dp, 1.0_dp]
    x(4, :) = [2.0_dp, 1.0_dp]

    call proxy_dist_auto(x, 'Euclidean', distances)
    call proxy_pack_dist(distances, packed)

    write (*, '(a)') 'Euclidean distance matrix:'
    write (*, '(4f9.4)') transpose(distances)
    write (*, '(a)') 'R-compatible packed dist vector:'
    write (*, '(*(f9.4,1x))') packed
end program basic
