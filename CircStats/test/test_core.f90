program test_core
    use circstats
    implicit none
    real(dp) :: x(5),pv,rg
    type(circ_dispersion_result) :: d
    type(circ_summary_result) :: s
    type(trig_moment_result) :: tm
    x=[0.1_dp,0.4_dp,0.9_dp,1.3_dp,5.8_dp]
    call check(circ_mean(x),0.4475203731924086_dp,1.0e-13_dp,"mean")
    call check(est_rho(x),0.818768671012122_dp,1.0e-13_dp,"rho")
    d=circ_disp(x)
    call check(d%rbar,0.818768671012122_dp,1.0e-13_dp,"disp")
    s=circ_summary(x)
    call check(s%mean_dir,0.4475203731924086_dp,1.0e-13_dp,"summary")
    tm=trig_moment(x,2,.true.)
    if (tm%rho < 0.0_dp .or. tm%rho > 1.0_dp) error stop 1
    rg=circ_range(x,pv)
    if (rg <= 0.0_dp .or. rg > twopi) error stop 1
    if (pv < 0.0_dp .or. pv > 1.0_dp) error stop 1
    call check(rad(180.0_dp),pi,1.0e-14_dp,"rad")
    call check(deg(pi),180.0_dp,1.0e-13_dp,"deg")
    if (abs(nck(10,3)-120.0_dp)>1.0e-10_dp) error stop 1
    print *, "test_core: PASS"
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
