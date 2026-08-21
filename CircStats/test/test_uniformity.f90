program test_uniformity
    use circstats
    implicit none
    real(dp) :: x(5),y(6)
    type(test_result) :: r,k,w,rs,w2,v
    x=[0.1_dp,0.4_dp,0.9_dp,1.3_dp,5.8_dp]
    y=[0.2_dp,0.6_dp,1.0_dp,1.6_dp,5.6_dp,6.0_dp]
    r=rayleigh_test(x)
    call check(r%statistic,0.818768671012122_dp,1.0e-13_dp,"ray statistic")
    call check(r%p_value,0.02551030352133594_dp,2.0e-14_dp,"ray p")
    v=v0_test(x,0.2_dp)
    call check(v%statistic,0.7938149860679167_dp,1.0e-13_dp,"v0 statistic")
    call check(v%p_value,0.003973987222075726_dp,2.0e-13_dp,"v0 p")
    k=kuiper_test(x)
    call check(k%statistic,1.7893466501843005_dp,2.0e-13_dp,"kuiper")
    w=watson_uniform_test(x)
    call check(w%statistic,0.1987597537015043_dp,2.0e-13_dp,"watson")
    rs=rao_spacing_test(x,0.05_dp)
    call check(rs%statistic,185.83100780887042_dp,2.0e-11_dp,"rao spacing")
    call check(rs%critical,183.44_dp,1.0e-13_dp,"rao critical")
    if (.not.rs%reject) error stop 1
    w2=watson_two_test(x,y)
    call check(w2%statistic,0.030303030303030307_dp,2.0e-14_dp,"watson two")
    print *, "test_uniformity: PASS"
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
