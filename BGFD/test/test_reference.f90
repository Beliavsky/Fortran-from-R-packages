program test_reference
    use bgfd
    implicit none
    real(dp), parameter :: x = 1.3_dp, prob = 0.37_dp
    integer :: failures
    failures = 0

    call check_close('d_bell_e', d_bell_e(x, 1.2_dp, 0.7_dp), 1.36451103912537064e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_bell_e', p_bell_e(x, 1.2_dp, 0.7_dp), 9.02210964943733540e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_bell_e', q_bell_e(prob, 1.2_dp, 0.7_dp), 1.90881832047162248e-01_dp, 5.0e-12_dp, failures)
    call check_close('d_cbell_e', d_cbell_e(x, 1.2_dp, 0.7_dp), 3.65624386676216662e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_cbell_e', p_cbell_e(x, 1.2_dp, 0.7_dp), 6.22102731299516210e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_cbell_e', q_cbell_e(prob, 1.2_dp, 0.7_dp), 7.22512920194378849e-01_dp, 5.0e-12_dp, failures)

    call check_close('d_bell_ee', d_bell_ee(x, 1.2_dp, 1.4_dp, 0.7_dp), 1.93833588819119385e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_bell_ee', p_bell_ee(x, 1.2_dp, 1.4_dp, 0.7_dp), 8.61561329666649445e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_bell_ee', q_bell_ee(prob, 1.2_dp, 1.4_dp, 0.7_dp), 3.23942392309311489e-01_dp, 5.0e-12_dp, failures)
    call check_close('d_cbell_ee', d_cbell_ee(x, 1.2_dp, 1.4_dp, 0.7_dp), 4.07285099449775978e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_cbell_ee', p_cbell_ee(x, 1.2_dp, 1.4_dp, 0.7_dp), 5.25634096409524521e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_cbell_ee', q_cbell_ee(prob, 1.2_dp, 1.4_dp, 0.7_dp), 9.43050256692747002e-01_dp, 5.0e-12_dp, failures)

    call check_close('d_bell_w', d_bell_w(x, 0.9_dp, 1.6_dp, 0.7_dp), 2.48002147925348193e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_bell_w', p_bell_w(x, 0.9_dp, 1.6_dp, 0.7_dp), 8.77526299085561123e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_bell_w', q_bell_w(prob, 0.9_dp, 1.6_dp, 0.7_dp), 4.25171927547500939e-01_dp, 5.0e-12_dp, failures)
    call check_close('d_cbell_w', d_cbell_w(x, 0.9_dp, 1.6_dp, 0.7_dp), 5.71437909915269882e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_cbell_w', p_cbell_w(x, 0.9_dp, 1.6_dp, 0.7_dp), 5.60763296552227786e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_cbell_w', q_cbell_w(prob, 0.9_dp, 1.6_dp, 0.7_dp), 9.76933993069713069e-01_dp, 5.0e-12_dp, failures)

    call check_close('d_bell_ew', d_bell_ew(x, 0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp), 3.25541915172512131e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_bell_ew', p_bell_ew(x, 0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp), 8.39338697797652022e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_bell_ew', q_bell_ew(prob, 0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp), 5.54013343110647760e-01_dp, 5.0e-12_dp, failures)
    call check_close('d_cbell_ew', d_cbell_ew(x, 0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp), 6.05478104458606836e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_cbell_ew', p_cbell_ew(x, 0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp), 4.81701279943591132e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_cbell_ew', q_cbell_ew(prob, 0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp), &
        1.11518977696466393e+00_dp, 5.0e-12_dp, failures)

    call check_close('d_bell_f', d_bell_f(x, 1.1_dp, 2.2_dp, 0.7_dp), 3.02574088812930964e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_bell_f', p_bell_f(x, 1.1_dp, 2.2_dp, 0.7_dp), 7.76056666625877689e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_bell_f', q_bell_f(prob, 1.1_dp, 2.2_dp, 0.7_dp), 5.93607182028216829e-01_dp, 5.0e-12_dp, failures)
    call check_close('d_cbell_f', d_cbell_f(x, 1.1_dp, 2.2_dp, 0.7_dp), 4.11644755358806402e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_cbell_f', p_cbell_f(x, 1.1_dp, 2.2_dp, 0.7_dp), 3.81024344182112618e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_cbell_f', q_cbell_f(prob, 1.1_dp, 2.2_dp, 0.7_dp), 1.27334226398881301e+00_dp, 5.0e-12_dp, failures)

    call check_close('d_bell_l', d_bell_l(x, 1.5_dp, 2.3_dp, 0.7_dp), 1.10352096841272507e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_bell_l', p_bell_l(x, 1.5_dp, 2.3_dp, 0.7_dp), 8.86819585398187171e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_bell_l', q_bell_l(prob, 1.5_dp, 2.3_dp, 0.7_dp), 1.57077699764316292e-01_dp, 5.0e-12_dp, failures)
    call check_close('d_cbell_l', d_cbell_l(x, 1.5_dp, 2.3_dp, 0.7_dp), 2.68800700008748561e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_cbell_l', p_cbell_l(x, 1.5_dp, 2.3_dp, 0.7_dp), 5.82780109360907184e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_cbell_l', q_cbell_l(prob, 1.5_dp, 2.3_dp, 0.7_dp), 6.86776126352676264e-01_dp, 5.0e-12_dp, failures)

    call check_close('d_bell_b', d_bell_b(x, 1.1_dp, 2.2_dp, 1.4_dp, 0.7_dp), 2.43570752663880369e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_bell_b', p_bell_b(x, 1.1_dp, 2.2_dp, 1.4_dp, 0.7_dp), 8.58586935835261689e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_bell_b', q_bell_b(prob, 1.1_dp, 2.2_dp, 1.4_dp, 0.7_dp), 5.01656470765129425e-01_dp, 5.0e-12_dp, failures)
    call check_close('d_cbell_b', d_cbell_b(x, 1.1_dp, 2.2_dp, 1.4_dp, 0.7_dp), 5.03295685605275378e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_cbell_b', p_cbell_b(x, 1.1_dp, 2.2_dp, 1.4_dp, 0.7_dp), 5.19436368391640979e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_cbell_b', q_cbell_b(prob, 1.1_dp, 2.2_dp, 1.4_dp, 0.7_dp), 1.02582311409353899e+00_dp, 5.0e-12_dp, failures)

    call check_close('d_bell_bx', d_bell_bx(x, 1.3_dp, 0.7_dp), 3.28622019325421677e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_bell_bx', p_bell_bx(x, 1.3_dp, 0.7_dp), 8.89663800344635769e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_bell_bx', q_bell_bx(prob, 1.3_dp, 0.7_dp), 5.91476418397062553e-01_dp, 5.0e-12_dp, failures)
    call check_close('d_cbell_bx', d_cbell_bx(x, 1.3_dp, 0.7_dp), 8.14443250892181836e-01_dp, 5.0e-12_dp, failures)
    call check_close('p_cbell_bx', p_cbell_bx(x, 1.3_dp, 0.7_dp), 5.89768380512192270e-01_dp, 5.0e-12_dp, failures)
    call check_close('q_cbell_bx', q_cbell_bx(prob, 1.3_dp, 0.7_dp), 1.03514283016848063e+00_dp, 5.0e-12_dp, failures)

    if (failures /= 0) error stop 1
    print '(a)', 'test_reference: PASS'
contains

    subroutine check_close(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: got, expected, tol
        integer, intent(inout) :: failures
        if (abs(got-expected) > tol*max(1.0_dp,abs(expected))) then
            print '(a,2es24.15)', trim(name)//' FAIL: ', got, expected
            failures = failures + 1
        end if
    end subroutine check_close

end program test_reference
