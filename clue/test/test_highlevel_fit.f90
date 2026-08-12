program test_highlevel_fit
    use clue_kinds, only: dp
    use clue_fit, only: tree_fit_result, fit_ultrametric_sumt, fit_addtree_sumt, &
        fit_l1_ultrametric_sumt, fit_l1_ultrametric_irip, fit_sum_of_ultrametrics
    use clue_trees, only: is_ultrametric, is_additive
    implicit none
    real(dp) :: d(4,4)
    real(dp), allocatable :: terms(:,:,:), total(:,:)
    type(tree_fit_result) :: a,b,c,e

    d = reshape([ &
        0.0_dp,1.0_dp,4.0_dp,5.0_dp, &
        1.0_dp,0.0_dp,3.0_dp,4.0_dp, &
        4.0_dp,3.0_dp,0.0_dp,2.0_dp, &
        5.0_dp,4.0_dp,2.0_dp,0.0_dp ], [4,4])

    a=fit_ultrametric_sumt(d,max_outer=20,max_inner=100)
    call check(is_ultrametric(a%distance,1.0e-7_dp),'SUMT ultrametric')
    call check(all(a%distance>=-1.0e-12_dp),'SUMT ultrametric nonnegative')

    b=fit_addtree_sumt(d,max_outer=20,max_inner=100)
    call check(is_additive(b%distance,1.0e-6_dp),'SUMT addtree')

    c=fit_l1_ultrametric_sumt(d,max_outer=20,max_inner=100)
    call check(is_ultrametric(c%distance,1.0e-7_dp),'L1 SUMT ultrametric')

    e=fit_l1_ultrametric_irip(d,maxiter=3)
    call check(is_ultrametric(e%distance,1.0e-7_dp),'L1 IRIP ultrametric')

    terms=fit_sum_of_ultrametrics(d,2,method='IP',maxiter=3)
    call check(is_ultrametric(terms(:,:,1),1.0e-7_dp),'sum term 1')
    call check(is_ultrametric(terms(:,:,2),1.0e-7_dp),'sum term 2')
    total=sum(terms,dim=3)
    call check(all(total>=-1.0e-12_dp),'sum nonnegative')

    print '(a)', 'test_highlevel_fit: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(*),intent(in)::msg
        if(.not.ok)then
            print '(a)', 'FAIL: '//trim(msg)
            error stop 1
        end if
    end subroutine
end program test_highlevel_fit
