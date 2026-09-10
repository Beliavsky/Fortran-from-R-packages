! SPDX-License-Identifier: GPL-2.0-only
program filter_design
    use signal
    implicit none

    type(arma_filter) :: filt
    type(frequency_response) :: response
    real(dp), allocatable :: impulse(:)
    real(dp), allocatable :: output(:)

    filt = butter_filter(4, [0.2_dp], 'low')
    allocate(impulse(16))
    impulse = 0.0_dp
    impulse(1) = 1.0_dp
    output = filter_signal(filt%b, filt%a, impulse)
    response = digital_frequency_response(filt%b, filt%a, 8)

    print '(a,*(1x,f10.6))', 'b:', filt%b
    print '(a,*(1x,f10.6))', 'a:', filt%a
    print '(a,*(1x,f10.6))', 'first impulse-response samples:', output(1:6)
    print '(a,1x,f10.6)', 'dc magnitude:', abs(response%h(1))
end program filter_design
