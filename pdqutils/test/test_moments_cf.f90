program test_moments_cf
    use pdqutils, only : dp, moment2cumulant, cumulant2moment, as269, as269_orders, qapx_cf
    implicit none
    real(dp), allocatable :: kappa(:), moms(:), ords(:)
    real(dp) :: ec(7), rate
    integer :: r, fails

    fails=0
    kappa=moment2cumulant([0.0_dp,1.0_dp,0.0_dp,3.0_dp,0.0_dp,15.0_dp])
    call check(maxval(abs(kappa-[0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]))<1.0e-13_dp,'normal cumulants')
    moms=cumulant2moment(kappa)
    call check(maxval(abs(moms-[0.0_dp,1.0_dp,0.0_dp,3.0_dp,0.0_dp,15.0_dp]))<1.0e-13_dp,'normal moments')

    rate=0.7_dp
    do r=1,7
        ec(r)=gamma(real(r,dp))/rate**r
    end do
    moms=cumulant2moment(moment2cumulant(ec))
    call check(maxval(abs(moms-ec))<2.0e-11_dp,'moment cumulant roundtrip')

    ords=as269_orders(-2.0_dp,[0.5_dp,1.2_dp,-0.3_dp,2.0_dp])
    call check(maxval(abs(ords-[-1.75_dp,-1.8083333333333333_dp,-1.7469907407407408_dp, &
        -1.6651755401234567_dp]))<2.0e-14_dp,'AS269 all orders')
    call check(abs(as269(0.0_dp,[0.5_dp,1.2_dp,-0.3_dp,2.0_dp])+0.1342746913580247_dp)<2.0e-14_dp,'AS269 final')

    call check(abs(qapx_cf(0.5_dp,ec)-0.9909017104784343_dp)<2.0e-12_dp,'Cornish-Fisher median')
    call check(abs(qapx_cf(0.9_dp,ec)-3.2824886823332515_dp)<2.0e-11_dp,'Cornish-Fisher 90%')
    call check(qapx_cf(0.0_dp,ec,support_lo=0.0_dp)==0.0_dp,'Cornish-Fisher lower endpoint')

    if(fails==0) then
        print '(a)','test_moments_cf: PASS'
    else
        print '(a,i0)','test_moments_cf: FAIL ',fails
        error stop 1
    end if
contains
    subroutine check(ok,name)
        logical,intent(in)::ok
        character(len=*),intent(in)::name
        if(.not.ok) then
            print '(a,a)','FAIL: ',trim(name)
            fails=fails+1
        end if
    end subroutine check
end program test_moments_cf
