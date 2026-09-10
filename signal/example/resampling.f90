! SPDX-License-Identifier: GPL-2.0-only
program resampling
    use signal
    implicit none

    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    integer :: i

    allocate(x(11))
    do i = 1, size(x)
        x(i) = sin(2.0_dp * signal_pi * real(i - 1, dp) / 5.0_dp)
    end do

    y = resample_signal(x, 1.0_dp, 0.05_dp)
    print '(a,i0)', 'input samples: ', size(x)
    print '(a,i0)', 'resampled samples: ', size(y)
    print '(a,*(1x,f9.5))', 'first samples:', y(1:min(8, size(y)))
end program resampling
