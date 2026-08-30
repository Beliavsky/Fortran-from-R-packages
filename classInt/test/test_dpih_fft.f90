program test_dpih_fft
    use classint, only: dp, dpih_bandwidth
    use classint_dpih, only: dpih_bandwidth_direct_reference
    implicit none

    real(dp) :: x(80)
    real(dp) :: h_fft
    real(dp) :: h_direct
    real(dp) :: relerr
    integer :: i
    integer :: level

    do i = 1, size(x)
        x(i) = 0.035_dp * real(i, dp) + sin(0.41_dp * real(i, dp)) + &
               0.2_dp * cos(0.17_dp * real(i, dp))
    end do

    do level = 1, 5
        h_fft = dpih_bandwidth(x, "minim", level, 401, [minval(x), maxval(x)], .true.)
        h_direct = dpih_bandwidth_direct_reference(x, "minim", level, 401, [minval(x), maxval(x)], .true.)
        relerr = abs(h_fft - h_direct) / max(1.0_dp, abs(h_direct))
        if (relerr > 2.0e-12_dp) error stop "dpih FFT/direct parity"
    end do

    h_fft = dpih_bandwidth(x, "minim", 2, 8193, [minval(x), maxval(x)], .true.)
    if (.not. (h_fft > 0.0_dp .and. h_fft < maxval(x) - minval(x))) error stop "large-grid FFT bandwidth"

    print *, "test_dpih_fft: ok"
end program test_dpih_fft
