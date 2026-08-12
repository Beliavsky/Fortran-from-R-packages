program test_pava_medoid
    use clue
    implicit none
    real(dp), allocatable :: y(:)
    real(dp) :: d(4,4)
    type(kmedoids_result) :: km
    y=pava_mean([3.0_dp,2.0_dp,4.0_dp])
    call check(maxval(abs(y-[2.5_dp,2.5_dp,4.0_dp]))<1e-12_dp,'PAVA mean')
    y=pava_median([3.0_dp,2.0_dp,4.0_dp])
    call check(all(y(1:2)>=2.0_dp).and.all(y(1:2)<=3.0_dp),'PAVA median')
    d=reshape([0.0_dp,1.0_dp,10.0_dp,11.0_dp, &
               1.0_dp,0.0_dp,9.0_dp,10.0_dp, &
               10.0_dp,9.0_dp,0.0_dp,1.0_dp, &
               11.0_dp,10.0_dp,1.0_dp,0.0_dp],[4,4])
    km=kmedoids(d,2)
    call check(km%status==0,'exact kmedoids LP status')
    call check(abs(km%criterion-2.0_dp)<1e-8_dp,'exact kmedoids criterion')
    km=pam_kmedoids(d,2)
    call check(abs(km%criterion-2.0_dp)<1e-8_dp,'PAM kmedoids criterion')
    print *, 'test_pava_medoid: PASS'
contains
    subroutine check(ok,msg)
    logical,intent(in)::ok
    character(*),intent(in)::msg
    if(.not.ok)then
    write(*,*)'FAIL: ',trim(msg)
    error stop 1
    end if
    end subroutine
end program
