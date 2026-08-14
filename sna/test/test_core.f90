program test_core
    use sna
    implicit none

    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp) :: g(3,3), x(1,3,3), y(3,3), dc(3), pt(16), p, lp
    real(dp), allocatable :: tc(:), b(:), cg(:,:), rg(:,:), se(:,:,:), css(:,:,:)
    type(geodist_result) :: gd
    type(network_regression_result) :: fit
    type(brokerage_result) :: brok
    integer :: cl(3), i, j

    g = 0.0_dp
    g(1,2)=1.0_dp; g(2,3)=1.0_dp; g(3,1)=1.0_dp

    gd = geodist(g)
    call assert_close(gd%distance(1,3), 2.0_dp, tol, 'geodist 1->3')
    call assert_close(gd%counts(1,3), 1.0_dp, tol, 'geodesic count')

    dc = dyad_census(g)
    call assert_close(dc(1), 0.0_dp, tol, 'mutual dyads')
    call assert_close(dc(2), 3.0_dp, tol, 'asymmetric dyads')
    call assert_close(dc(3), 0.0_dp, tol, 'null dyads')
    call assert_close(gden(g), 0.5_dp, tol, 'directed density')

    tc = triad_census(g, .true.)
    call assert_close(tc(10), 1.0_dp, tol, '030C triad')
    call assert_close(sum(tc), 1.0_dp, tol, 'triad total')

    b = betweenness(g)
    call assert_close(b(1), 1.0_dp, tol, 'cycle betweenness v1')
    call assert_close(b(2), 1.0_dp, tol, 'cycle betweenness v2')
    call assert_close(b(3), 1.0_dp, tol, 'cycle betweenness v3')

    allocate(se(2,3,3))
    se(1,:,:) = g; se(2,:,:) = g
    cg = centralgraph(se)
    call assert_close(maxval(abs(cg-g)), 0.0_dp, tol, 'central graph')
    call assert_close(structdist(g,g), 0.0_dp, tol, 'structural distance identity')

    rg = rgnm(5, 7, 'digraph', .false.)
    call assert_close(sum(rg), 7.0_dp, tol, 'rgnm exact edge count')

    pt = bn_ptriad(0.2_dp,0.1_dp,0.15_dp,0.25_dp)
    call assert_close(sum(pt), 1.0_dp, 1.0e-9_dp, 'biased-net triad probabilities')
    if (minval(pt) < -tol) error stop 'negative biased-net triad probability'

    lp = bn_lpt_triad(0,0,0,0,0,0,2,1,3,0.23_dp,0.37_dp,0.41_dp,0.11_dp)
    call assert_close(lp,-3.665995484510641_dp,1.0e-14_dp,'biased-net exact triad kernel 000000')
    lp = bn_lpt_triad(1,1,1,1,1,1,2,1,3,0.23_dp,0.37_dp,0.41_dp,0.11_dp)
    call assert_close(lp,-1.3088270601882654_dp,1.0e-14_dp,'biased-net exact triad kernel 111111')

    css=sr2css(g)
    if(stackcount(css)/=3)error stop 'stackcount rank-3 failed'
    if(stackcount(g)/=1)error stop 'stackcount rank-2 failed'

    p = bbnam_probtie([1.0_dp],0.5_dp,[0.1_dp],[0.2_dp])
    call assert_close(p, 0.45_dp/0.55_dp, tol, 'fixed Bayesian tie posterior')

    g=0.0_dp; g(1,2)=1.0_dp; g(2,3)=1.0_dp; cl=[1,1,1]
    brok=brokerage(g,cl)
    call assert_close(brok%raw(2,1),1.0_dp,tol,'brokerage coordinator')
    call assert_close(brok%raw(2,6),1.0_dp,tol,'brokerage total')

    x = 0.0_dp
    x(1,1,2)=1.0_dp; x(1,1,3)=2.0_dp
    x(1,2,1)=3.0_dp; x(1,2,3)=4.0_dp
    x(1,3,1)=5.0_dp; x(1,3,2)=6.0_dp
    y = 0.0_dp
    do i=1,3
        do j=1,3
            if(i/=j)y(i,j)=1.0_dp+2.0_dp*x(1,i,j)
        end do
    end do
    fit=netlm(y,x,intercept=.true.,mode='digraph',diag=.false.,nullhyp='classical')
    if (.not.fit%fit%converged) error stop 'netlm did not converge'
    call assert_close(fit%fit%coef(1),1.0_dp,1.0e-9_dp,'netlm intercept')
    call assert_close(fit%fit%coef(2),2.0_dp,1.0e-9_dp,'netlm slope')

    print '(a)', 'All sna-fortran core tests passed.'

contains
    subroutine assert_close(actual,expected,atol,label)
        real(dp),intent(in)::actual,expected,atol
        character(len=*),intent(in)::label
        if(abs(actual-expected)>atol)then
            write(*,'(a,2(1x,es24.16))') trim(label)//' failed:',actual,expected
            error stop 1
        end if
    end subroutine assert_close
end program test_core
