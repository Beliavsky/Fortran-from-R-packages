program test_gca
    use pdqutils, only : dp, dapx_gca, papx_gca, gca_normal, gca_gamma, gca_beta, gca_arcsine, gca_wigner
    implicit none
    real(dp) :: rn(5),rg(6),rb(6),rt(6),ra(6),rw(6)
    integer :: i,j,fails
    fails=0
    rn=[2.0_dp,13.0_dp,62.0_dp,475.0_dp,3182.0_dp]
    call near(dapx_gca(1.0_dp,rn,basis=gca_normal),0.12579440923099772_dp,2e-14_dp,'normal pdf')
    call near(papx_gca(1.0_dp,rn,basis=gca_normal),0.36944134018176367_dp,2e-14_dp,'normal cdf')

    do i=1,6
        rg(i)=1.0_dp
        do j=0,i-1
            rg(i)=rg(i)*2.0_dp*(15.0_dp+real(j,dp))
        end do
    end do
    call near(dapx_gca(30.0_dp,rg,basis=gca_gamma,support_lo=0.0_dp), &
        0.05121793333226722_dp,3e-13_dp,'gamma pdf')
    call near(papx_gca(30.0_dp,rg,basis=gca_gamma,support_lo=0.0_dp), &
        0.5343462910559904_dp,3e-13_dp,'gamma cdf')

    do i=1,6
        rb(i)=1.0_dp
        do j=0,i-1
            rb(i)=rb(i)*(20.0_dp+real(j,dp))/(60.0_dp+real(j,dp))
        end do
    end do
    call near(dapx_gca(0.3_dp,rb,basis=gca_beta),5.908419205987833_dp,6e-11_dp,'beta pdf')
    call near(papx_gca(0.3_dp,rb,basis=gca_beta),0.29982167525229364_dp,6e-11_dp,'beta cdf')

    do i=1,6
        rt(i)=1.0_dp
        do j=0,i-1
            rt(i)=rt(i)*(2.0_dp+real(j,dp))/(7.0_dp+real(j,dp))
        end do
    end do
    call near(dapx_gca(0.4_dp,rt,basis=gca_beta,shape1=3.0_dp,shape2=4.0_dp), &
        1.675933286400124_dp,3e-12_dp,'non-parent beta pdf')
    call near(papx_gca(0.4_dp,rt,basis=gca_beta,shape1=3.0_dp,shape2=4.0_dp), &
        0.7781264588799924_dp,3e-12_dp,'non-parent beta cdf')

    ra=[0.0_dp,0.5_dp,0.0_dp,3.0_dp/8.0_dp,0.0_dp,5.0_dp/16.0_dp]
    call near(dapx_gca(0.0_dp,ra,basis=gca_arcsine),0.3183098861837907_dp,3e-14_dp,'arcsine pdf')
    call near(papx_gca(-0.5_dp,ra,basis=gca_arcsine),1.0_dp/3.0_dp,3e-14_dp,'arcsine cdf')

    rw=[0.0_dp,0.25_dp,0.0_dp,0.125_dp,0.0_dp,5.0_dp/64.0_dp]
    call near(dapx_gca(0.0_dp,rw,basis=gca_wigner),0.6366197723675814_dp,3e-14_dp,'wigner pdf')
    call near(papx_gca(0.7_dp,rw,basis=gca_wigner),0.9059397978129061_dp,3e-13_dp,'wigner cdf')

    call check(dapx_gca(-0.1_dp,rb,basis=gca_beta)==0.0_dp,'beta support pdf')
    call check(papx_gca(1.1_dp,rb,basis=gca_beta)==1.0_dp,'beta support cdf')
    call near(papx_gca(0.3_dp,rb,basis=gca_beta)+papx_gca(0.3_dp,rb,basis=gca_beta,lower_tail=.false.), &
        1.0_dp,2e-14_dp,'gca tail complement')

    if(fails==0) then
        print '(a)','test_gca: PASS'
    else
        print '(a,i0)','test_gca: FAIL ',fails
        error stop 1
    end if
contains
    subroutine near(x,y,tol,name)
        real(dp),intent(in)::x,y,tol
        character(len=*),intent(in)::name
        call check(abs(x-y)<=tol,name)
    end subroutine near
    subroutine check(ok,name)
        logical,intent(in)::ok
        character(len=*),intent(in)::name
        if(.not.ok) then
            print '(a,a)','FAIL: ',trim(name); fails=fails+1
        end if
    end subroutine check
end program test_gca
