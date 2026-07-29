! SPDX-License-Identifier: GPL-2.0-or-later
program test_bivariate
    use evir
    use test_support
    implicit none
    type(gpdbiv_fit_result)::fit,known
    real(dp)::x1(300),x2(300),p,probs(6),u
    integer::i,status,j

    do i=1,300
        p=(real(i,dp)-0.5_dp)/300.0_dp
        x1(i)=qgpd(p,0.15_dp,beta=1.0_dp)
        j=mod(37*(i-1),300)+1
        u=(real(j,dp)-0.5_dp)/300.0_dp
        x2(i)=qgpd(0.65_dp*p+0.35_dp*u,0.10_dp,beta=1.3_dp)
    end do
    fit=gpdbiv(x1,x2,ne1=100,ne2=100,global_fit=.false.)
    call check(fit%alpha>0.0_dp.and.fit%alpha<1.0_dp,'bivariate alpha domain')
    call check(fit%par1(2)>0.0_dp.and.fit%par2(2)>0.0_dp,'bivariate marginal scales')

    known%u1=1.0_dp; known%u2=1.5_dp
    known%lambda1=0.1_dp; known%lambda2=0.15_dp
    known%par1=[0.2_dp,1.0_dp]; known%par2=[0.1_dp,1.2_dp]
    known%alpha=0.6_dp
    call interpret_gpdbiv(known,3.0_dp,3.5_dp,probs,status)
    call check(status==evir_ok,'bivariate interpretation status')
    call check_close(probs(4),probs(1)*probs(2),1.0e-14_dp,'product probability')
    call check_close(probs(5),probs(3)/probs(1),1.0e-14_dp,'conditional probability')
    call check_close(bivariate_survivor(known,3.0_dp,3.5_dp),0.011371663108103669_dp,1.0e-13_dp,'survivor SciPy reference')
    call check_close(bivariate_cdf(known,3.0_dp,3.5_dp),0.9606694725597203_dp,1.0e-13_dp,'CDF SciPy reference')
    call finish_tests()
end program test_bivariate
