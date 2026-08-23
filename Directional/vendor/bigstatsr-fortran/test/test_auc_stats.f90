program test_auc_stats
    use bigstatsr
    implicit none
    type(pcor_result) :: pr
    type(auc_boot_result) :: ab
    real(dp) :: x(8),y(8),z(8,1)
    integer :: i
    call check(abs(auc([0.0_dp,0.0_dp],[0,1])-0.5_dp)<1.0e-14_dp,'AUC tie')
    call check(abs(auc([0.2_dp,0.1_dp,1.0_dp],[0,0,1])-1.0_dp)<1.0e-14_dp,'AUC perfect')

    do i=1,8
        z(i,1)=real(i-4,dp)
        x(i)=0.7_dp*z(i,1)+real(mod(3*i,5),dp)
        y(i)=-0.4_dp*z(i,1)+2.0_dp*real(mod(3*i,5),dp)+0.2_dp*real(mod(i,2),dp)
    end do
    pr=pcor(x,y,z)
    call check(pr%info==0,'pcor info')
    call check(pr%r>0.95_dp,'pcor strength')
    call check(pr%lower<pr%r .and. pr%upper>pr%r,'pcor interval')

    ab=auc_bootstrap([0.1_dp,0.2_dp,0.8_dp,0.9_dp],[0,0,1,1],200)
    call check(ab%mean>0.9_dp,'bootstrap auc')
    print *, 'test_auc_stats: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok) then
            print *, 'FAIL: ',trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_auc_stats
