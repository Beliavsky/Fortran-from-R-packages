program test_simulation
    use msm, only : dp,make_generator,set_random_seed,simulate_ctmc,ctmc_path
    implicit none
    real(dp) :: off(2,2),q(2,2),s
    type(ctmc_path) :: path
    integer :: i,n
    call set_random_seed(12345); off=0.0_dp; off(1,2)=0.5_dp; q=make_generator(off)
    n=4000; s=0.0_dp
    do i=1,n
        call simulate_ctmc(q,1,100.0_dp,path,max_events=4)
        if(size(path%time)>=2) s=s+path%time(2)
        call check(all(path%state>=1.and.path%state<=2),"states in range")
        call check(all(path%time(2:)>=path%time(:size(path%time)-1)),"times ordered")
    end do
    call check(abs(s/real(n,dp)-2.0_dp)<0.08_dp,"mean exponential waiting time")
    print '(a)', 'test_simulation: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; write(*,'(a)') 'FAIL: '//msg; error stop 1; end if
    end subroutine check
end program test_simulation
