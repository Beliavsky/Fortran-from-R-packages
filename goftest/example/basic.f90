program basic
    use goftest, only : dp, p_ad, q_ad, p_cvm, q_cvm, ad_test_values, cvm_test_values, gof_result
    implicit none

    real(dp) :: u(8)
    type(gof_result) :: ad, cvm

    u = [0.05_dp, 0.18_dp, 0.29_dp, 0.41_dp, 0.55_dp, 0.68_dp, 0.82_dp, 0.94_dp]
    ad = ad_test_values(u)
    cvm = cvm_test_values(u)

    write(*,'(a,f12.8)') 'AD statistic  = ', ad%statistic
    write(*,'(a,f12.8)') 'AD p-value    = ', ad%p_value
    write(*,'(a,f12.8)') 'CvM statistic = ', cvm%statistic
    write(*,'(a,f12.8)') 'CvM p-value   = ', cvm%p_value
    write(*,'(a,f12.8)') 'qAD(.95)      = ', q_ad(0.95_dp)
    write(*,'(a,f12.8)') 'qCvM(.95)     = ', q_cvm(0.95_dp)
    write(*,'(a,f12.8)') 'pAD(1)        = ', p_ad(1.0_dp)
    write(*,'(a,f12.8)') 'pCvM(.2)      = ', p_cvm(0.2_dp)
end program basic
