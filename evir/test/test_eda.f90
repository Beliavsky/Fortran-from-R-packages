! SPDX-License-Identifier: GPL-2.0-or-later
program test_eda
    use evir
    use test_support
    implicit none
    real(dp)::x(10),times(10),pareto(1000)
    type(xy_result)::xy
    type(records_result)::rr
    type(decluster_result)::dc
    type(band_result)::hb
    type(matrix_result)::ei
    integer::i

    x=[1.0_dp,2.0_dp,4.0_dp,3.0_dp,5.0_dp,2.0_dp,6.0_dp,1.0_dp,7.0_dp,3.0_dp]
    times=[0.0_dp,1.0_dp,2.0_dp,10.0_dp,11.0_dp,20.0_dp,21.0_dp,22.0_dp,40.0_dp,41.0_dp]
    call check_close(findthresh(x,3),4.0_dp,1.0e-14_dp,'find threshold')
    dc=decluster(x,times,3.0_dp)
    call check(dc%status==evir_ok,'decluster status')
    call check(size(dc%values)==4,'decluster count')
    call check_close(dc%values(1),4.0_dp,1.0e-14_dp,'decluster first maximum')
    xy=emplot(x)
    call check(size(xy%x)==10.and.xy%y(1)>xy%y(10),'empirical tail')
    xy=meplot(x,omit=1)
    call check(size(xy%x)>0.and.all(xy%y>=0.0_dp),'mean excess')
    xy=qplot(x,xi=0.1_dp)
    call check(size(xy%x)==10.and.xy%y(10)>xy%y(1),'quantile plot data')
    rr=records(x)
    call check(size(rr%record)==6,'record count')
    call check(rr%trial(6)==9,'last record trial')
    do i=1,size(pareto)
        pareto(i)=(1.0_dp-(real(i,dp)-0.5_dp)/real(size(pareto),dp))**(-0.5_dp)
    end do
    hb=hill(pareto,option='xi',start_k=50,end_k=200)
    call check_close(hb%estimate(1),0.5_dp,3.0e-2_dp,'Hill Pareto xi')
    ei=exindex(x,block_size=2,start_k=1,end_k=5)
    call check(ei%status==evir_ok.and.size(ei%values,2)==5,'extremal index table')
    call finish_tests()
end program test_eda
