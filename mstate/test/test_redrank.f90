program test_redrank
    use mstate
    implicit none
    integer,parameter::nsub=24,n=2*nsub
    real(dp)::start(n),stop(n),z(n,2),g0(1,2)
    integer::status(n),transid(n),i,info
    type(redrank_result)::rr
    real(dp)::detb

    do i=1,nsub
        start(2*i-1:2*i)=0.0_dp
        stop(2*i-1:2*i)=real(i,dp)
        transid(2*i-1)=1;transid(2*i)=2
        z(2*i-1,1)=real(i-nsub/2,dp)/6.0_dp
        z(2*i-1,2)=sin(real(i,dp))
        z(2*i,:)=z(2*i-1,:)
        status(2*i-1)=merge(1,0,mod(i,3)==0.or.mod(i,5)==0)
        status(2*i)=merge(1,0,mod(i,4)==0.or.mod(i,7)==0)
    end do
    g0=reshape([0.8_dp,-0.6_dp],[1,2])
    call redrank_fit(start,stop,status,transid,z,1,rr,gamma_start=g0,method='breslow', &
                     eps=1.0e-6_dp,maxiter=40,info=info)
    call check(info==0,'redrank info')
    call check(all(shape(rr%beta)==[2,2]),'redrank beta shape')
    detb=rr%beta(1,1)*rr%beta(2,2)-rr%beta(1,2)*rr%beta(2,1)
    call check(abs(detb)<1.0e-8_dp,'rank one beta')
    call check(rr%niter>=1.and.rr%df==3,'redrank metadata')
    call check(all(abs(matmul(rr%alpha,rr%gamma)-rr%beta)<1.0e-12_dp),'factorization')
    print '(a)','test_redrank: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program test_redrank
