program test_models
    use circstats
    implicit none
    real(dp) :: alpha(12),theta(12),wc(80),u,scale
    type(circ_reg_result) :: reg
    type(wrapped_cauchy_fit_result) :: fit
    type(vm_fit_result) :: vm
    integer :: i
    alpha=[0.05_dp,0.4_dp,0.8_dp,1.2_dp,1.7_dp,2.1_dp,2.6_dp,3.0_dp,3.5_dp,4.0_dp,4.6_dp,5.2_dp]
    theta=wrap_2pi(0.7_dp+0.5_dp*sin(alpha)-0.2_dp*cos(alpha)+0.08_dp*sin(2.0_dp*alpha))
    reg=circ_reg(alpha,theta,1)
    call check(reg%rho,0.997383125621789_dp,2.0e-13_dp,"reg rho")
    call check(reg%pvalues(1),0.025310259303224003_dp,2.0e-12_dp,"reg p1")
    call check(reg%pvalues(2),0.011687272086263772_dp,2.0e-12_dp,"reg p2")
    call check(reg%coef(1,1),0.70687827_dp,2.0e-8_dp,"reg coef")
    scale=-log(0.62_dp)
    do i=1,size(wc)
        u=(real(i,dp)-0.5_dp)/real(size(wc),dp)
        wc(i)=wrap_2pi(0.8_dp+scale*tan(pi*(u-0.5_dp)))
    end do
    fit=wrpcauchy_ml(wc,0.75_dp,0.55_dp,1.0e-12_dp)
    if (.not.fit%converged) error stop 1
    if (abs(wrap_pi(fit%mu-0.8_dp))>0.08_dp) error stop 1
    if (abs(fit%rho-0.62_dp)>0.08_dp) error stop 1
    vm=vm_ml(theta)
    if (vm%kappa <= 0.0_dp) error stop 1
    print *, "test_models: PASS"
contains
    subroutine check(got,want,eps,label)
        real(dp),intent(in)::got,want,eps
        character(*),intent(in)::label
        if (abs(got-want)>eps) then
            print *, trim(label),got,want
            error stop 1
        end if
    end subroutine
end program
