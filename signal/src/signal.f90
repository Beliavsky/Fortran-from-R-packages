! SPDX-License-Identifier: GPL-2.0-only
!
! Public API for the modern Fortran translation of CRAN package signal 1.8-1.
module signal
    use signal_kinds, only : dp, signal_pi
    use signal_types, only : arma_filter, zpg_filter, filter_order, fft_filter, median_filter, &
        frequency_response, analog_response, group_delay_response, impulse_response_data, &
        sgolay_filter, spectrogram_data
    use signal_utils, only : convolve, inverse_dft, polyval_real_complex, polyval_complex, &
        polynomial_from_roots, polynomial_roots, sinc_value
    use signal_windows, only : bartlett_window, blackman_window, boxcar_window, &
        chebyshev_window, flattop_window, gaussian_window, hamming_window, hanning_window, &
        kaiser_window, triangular_window
    use signal_filters, only : make_arma, make_ma, make_zpg, arma_from_zpg, zpg_from_arma, &
        make_filter_order, make_fft_filter, make_median_filter, make_spencer_filter, &
        unit_phasor_degrees, filter_signal, fft_filter_signal, zero_phase_filter, &
        median_filter_signal, spencer_smooth, levinson_durbin, unwrap_phase, fractional_difference
    use signal_iir_design, only : bilinear_transform, splane_frequency_transform, butter_filter, &
        butter_order, cheby1_filter, cheby2_filter, cheby1_order, ellip_filter, ellip_order
    use signal_fir_design, only : fir1_filter, fir2_filter, kaiser_order, remez_filter
    use signal_interpolation, only : interp1_signal, pchip_interpolate, interpolate_signal, &
        resample_signal, decimate_signal
    use signal_analysis, only : chirp_signal, analog_frequency_response, digital_frequency_response, &
        digital_frequency_response_at, group_delay, impulse_response, spectrogram
    use signal_sgolay, only : savitzky_golay, savitzky_golay_filter
    implicit none
    public

end module signal
