program test_edgeworth
    use pdqutils, only : dp, dapx_edgeworth, papx_edgeworth
    implicit none
    real(dp) :: kappa(7),p1,p2
    integer :: r,fails
    fails=0
    kappa(1)=10.0_dp
    do r=2,7
        kappa(r)=2.0_dp**(r-1)*gamma(real(r,dp))*10.0_dp
    end do
    call check(abs(dapx_edgeworth(4.0_dp,kappa)-0.04544992358242081_dp)<2.0e-14_dp,'density x=4')
    call check(abs(papx_edgeworth(4.0_dp,kappa)-0.052722527913258495_dp)<2.0e-14_dp,'cdf x=4')
    call check(abs(dapx_edgeworth(10.0_dp,kappa)-0.08773182546165119_dp)<2.0e-14_dp,'density x=10')
    call check(abs(papx_edgeworth(10.0_dp,kappa)-0.5595073831360453_dp)<2.0e-14_dp,'cdf x=10')
    call check(abs(dapx_edgeworth(18.0_dp,kappa)-0.017154576883529258_dp)<2.0e-14_dp,'density x=18')
    call check(abs(papx_edgeworth(18.0_dp,kappa)-0.9454341861582858_dp)<2.0e-14_dp,'cdf x=18')
    p1=papx_edgeworth(12.0_dp,kappa)
    p2=papx_edgeworth(12.0_dp,kappa,lower_tail=.false.)
    call check(abs(p1+p2-1.0_dp)<2.0e-14_dp,'tail complement')
    call check(dapx_edgeworth(-1.0_dp,kappa,support_lo=0.0_dp)==0.0_dp,'density support')
    call check(papx_edgeworth(-1.0_dp,kappa,support_lo=0.0_dp)==0.0_dp,'cdf support')
    if(fails==0) then
        print '(a)','test_edgeworth: PASS'
    else
        print '(a,i0)','test_edgeworth: FAIL ',fails
        error stop 1
    end if
contains
    subroutine check(ok,name)
        logical,intent(in)::ok
        character(len=*),intent(in)::name
        if(.not.ok) then
            print '(a,a)','FAIL: ',trim(name); fails=fails+1
        end if
    end subroutine check
end program test_edgeworth
