program test_analysis
    use circstats
    implicit none
    real(dp) :: x(8),data(11)
    integer :: group(11)
    type(change_point_result) :: cp
    type(rao_homogeneity_result) :: rh
    x=[0.1_dp,0.3_dp,0.5_dp,0.7_dp,2.2_dp,2.5_dp,2.8_dp,3.1_dp]
    cp=change_pt(x)
    if (cp%k_r < 1 .or. cp%k_r > size(x)) error stop 1
    if (cp%k_t < 2 .or. cp%k_t > size(x)-2) error stop 1
    if (cp%rmax < 0.0_dp) error stop 1
    data=[0.1_dp,0.3_dp,0.5_dp,0.7_dp,0.9_dp,0.2_dp,0.45_dp,0.65_dp,0.85_dp,1.05_dp,1.2_dp]
    group=[1,1,1,1,1,2,2,2,2,2,2]
    rh=rao_homogeneity_test(data,group)
    call check(rh%polar%statistic,1.0957177596691672_dp,2.0e-12_dp,"rao polar")
    call check(rh%dispersion%statistic,0.42557384362635275_dp,2.0e-12_dp,"rao dispersion")
    call check(rh%polar%p_value,0.2952078014583692_dp,3.0e-12_dp,"rao polar p")
    print *, "test_analysis: PASS"
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
