program test_sumt
    use clue
    implicit none
    type(sumt_result) :: r
    r=sumt_optimize([-2.0_dp],loss,penalty,grad_L=gl,grad_P=gp,eps=1e-8_dp,max_outer=20)
    call check(r%x(1)>0.99_dp .and. r%x(1)<1.01_dp,'SUMT constrained minimum')
    call check(r%penalty<1e-4_dp,'SUMT penalty')
    print *, 'test_sumt: PASS'
contains
    function loss(x) result(v)
        real(dp),intent(in)::x(:);real(dp)::v;v=x(1)**2
    end function
    function penalty(x) result(v)
        real(dp),intent(in)::x(:);real(dp)::v;v=max(1.0_dp-x(1),0.0_dp)**2
    end function
    subroutine gl(x,g)
        real(dp),intent(in)::x(:);real(dp),intent(out)::g(:);g(1)=2.0_dp*x(1)
    end subroutine
    subroutine gp(x,g)
        real(dp),intent(in)::x(:);real(dp),intent(out)::g(:)
        if(x(1)<1.0_dp)then;g(1)=-2.0_dp*(1.0_dp-x(1));else;g(1)=0.0_dp;end if
    end subroutine
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(*),intent(in)::msg
        if(.not.ok)then;write(*,*)'FAIL: ',trim(msg);error stop 1;end if
    end subroutine
end program
