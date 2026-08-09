program test_sampling
    use limsolve
    implicit none
    type(sample_result) :: s
    real(dp) :: a0(0,2),b0(0),e(1,2),f(1),g(2,2),h(2)
    integer :: i
    e(1,:)=[1.0_dp,1.0_dp]; f=[1.0_dp]
    g=0.0_dp; g(1,1)=1.0_dp; g(2,2)=1.0_dp; h=0.0_dp
    call run_method('rda',123)
    call run_method('cda',456)
    call run_method('mirror',789)
    print *, 'PASS test_sampling'
contains
    subroutine run_method(method,seed_value)
        character(len=*), intent(in) :: method
        integer, intent(in) :: seed_value
        call xsample(a0,b0,e,f,g,h,s,iter=300,outputlength=60,type=method, &
            seed=seed_value,fulloutput=.true.)
        call check(s%status==LS_SUCCESS,'sampling status '//trim(method))
        do i=1,size(s%x,1)
            call check(abs(sum(s%x(i,:))-1.0_dp)<1.0e-7_dp,'sampling equality')
            call check(minval(s%x(i,:))>=-1.0e-8_dp,'sampling inequality')
        end do
        call check(s%accepted_ratio>0.0_dp,'sampling acceptance')
    end subroutine run_method

    subroutine check(ok,msg)
        logical,intent(in)::ok; character(len=*),intent(in)::msg
        if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
    end subroutine check
end program test_sampling
