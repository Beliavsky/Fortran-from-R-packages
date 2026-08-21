program test_envelope
    use boot_kinds, only : dp
    use boot_envelope
    implicit none
    real(dp)::m(9,3),pl(3),ph(3),ol(3),oh(3),pe,oe
    integer::i,kp,ko
    do i=1,9
        m(i,:)=[real(i,dp),real(i*i,dp),real(10-i,dp)]
    end do
    call confidence_envelope(m,0.8_dp,0.8_dp,pl,ph,ol,oh,pe,oe,kp,ko)
    if(any(pl>ph).or.any(ol>oh))error stop 1
    if(kp<1 .or. ko<1 .or. ko>kp)error stop 2
    if(pe<0.0_dp .or. pe>1.0_dp .or. oe<0.0_dp .or. oe>1.0_dp)error stop 3
    print '(a)', 'test_envelope: PASS'
end program test_envelope
