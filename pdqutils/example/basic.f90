program basic
    use pdqutils, only : dp, moment2cumulant, dapx_gca, papx_gca, qapx_cf, gca_gamma
    implicit none
    real(dp) :: moments(6),cumulants(6)
    integer :: i,j

    ! First six moments of a chi-square distribution with 30 degrees of freedom.
    do i=1,6
        moments(i)=1.0_dp
        do j=0,i-1
            moments(i)=moments(i)*2.0_dp*(15.0_dp+real(j,dp))
        end do
    end do
    cumulants=moment2cumulant(moments)

    print '(a,es16.8)','GCA density at 30: ', &
        dapx_gca(30.0_dp,moments,basis=gca_gamma,support_lo=0.0_dp)
    print '(a,es16.8)','GCA CDF at 30:     ', &
        papx_gca(30.0_dp,moments,basis=gca_gamma,support_lo=0.0_dp)
    print '(a,es16.8)','CF median:          ', &
        qapx_cf(0.5_dp,cumulants,support_lo=0.0_dp)
end program basic
