program test_stats_utils
    use greybox_kinds, only: dp
    use greybox_stats, only: correlation_result, cramer_v, mcor, pcor_matrix, determination
    use greybox_utils, only: xreg_multiplier, xreg_transformer, backshift, temporal_dummy
    implicit none
    integer, parameter :: n=12
    integer :: a(n), b(n), failures, i
    real(dp) :: x(n,3), y(n), pc(3,3), det(3), ex(n,6), tr(n,6), lagged(n), td(n,3)
    type(correlation_result) :: r

    failures=0
    do i=1,n
        a(i)=merge(1,2,i<=6)
        b(i)=a(i)
        x(i,1)=real(i,dp)
        x(i,2)=2.0_dp*x(i,1)+0.1_dp*sin(real(i,dp))
        x(i,3)=sin(2.3_dp*real(i,dp))
        y(i)=1.0_dp+0.8_dp*x(i,1)-0.2_dp*x(i,3)
    end do
    r=cramer_v(a,b)
    call check_true(abs(r%value-1.0_dp)<1.0e-12_dp,'Cramer V perfect')
    r=mcor(x(:,[1,3]),y)
    call check_true(r%value>0.99_dp,'multiple correlation')
    pc=pcor_matrix(x)
    call check_true(maxval(abs(pc-transpose(pc)))<1.0e-12_dp,'partial correlation symmetry')
    det=determination(x)
    call check_true(det(2)>0.99_dp,'determination collinear variable')

    ex=xreg_multiplier(x(:,1:3))
    call check_true(size(ex,2)==6,'xreg multiplier size')
    tr=xreg_transformer(abs(x(:,1:1))+1.0_dp,[1,2,3,4,5])
    call check_true(size(tr,2)==6,'xreg transformer size')
    lagged=backshift(x(:,1),1,0)
    call check_true(abs(lagged(1))<1.0e-14_dp.and.abs(lagged(2)-x(1,1))<1.0e-14_dp,'backshift')
    td=temporal_dummy(n,3)
    call check_true(all(abs(sum(td,dim=2)-1.0_dp)<1.0e-14_dp),'temporal dummies')

    if(failures/=0)then
        write(*,'(a,i0)')'test_stats_utils: FAIL ',failures
        error stop 1
    end if
    write(*,'(a)')'test_stats_utils: PASS'
contains
    subroutine check_true(ok,name)
        logical,intent(in)::ok
        character(len=*),intent(in)::name
        if(.not.ok)then;failures=failures+1;write(*,'(a)')trim(name)//': failed';end if
    end subroutine check_true
end program test_stats_utils
