program test_one_shot
    use scip
    implicit none
    real(dp) :: a(2,2), obj(2), b(2)
    real(dp) :: ak(1,6), bk(1), objk(6)
    character(len=2) :: sense(2), sense1(1)
    character(len=1) :: vtype(6)
    type(scip_result) :: res, res_csc, rk
    type(scip_control) :: ctrl
    type(scip_csc_matrix) :: csc

    a = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2,2])
    obj = [-5.0_dp, -4.0_dp]
    b = [6.0_dp, 8.0_dp]
    sense = ['<=', '<=']
    ctrl%verbose = .false.

    res = scip_solve(obj, a, b, sense, control=ctrl)
    call require(trim(res%status) == 'optimal', 'dense LP status')
    call require(abs(res%objval + 22.0_dp) < 1.0e-8_dp, 'dense LP objective')
    call require(maxval(abs(res%x - [10.0_dp/3.0_dp, 4.0_dp/3.0_dp])) < 1.0e-7_dp, &
                 'dense LP solution')

    csc = make_csc_matrix(a)
    res_csc = scip_solve(obj, csc, b, sense, control=ctrl)
    call require(abs(res_csc%objval - res%objval) < 1.0e-12_dp, 'CSC objective')
    call require(maxval(abs(res_csc%x - res%x)) < 1.0e-12_dp, 'CSC solution')

    ak = reshape([7.0_dp, 2.0_dp, 7.0_dp, 5.0_dp, 1.0_dp, 3.0_dp], [1,6])
    objk = -1.0_dp
    bk = [13.0_dp]
    sense1 = ['<=']
    vtype = 'B'
    rk = scip_solve(objk, ak, bk, sense1, vtype=vtype, control=ctrl)
    call require(trim(rk%status) == 'optimal', 'knapsack status')
    call require(abs(rk%objval + 4.0_dp) < 1.0e-8_dp, 'knapsack objective')
    call require(sum(rk%x) > 3.5_dp .and. sum(rk%x) < 4.5_dp, 'knapsack item count')
    call require(dot_product(ak(1,:), rk%x) <= 13.0_dp + 1.0e-8_dp, 'knapsack capacity')

    print *, 'test_one_shot: PASS'
contains
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine
end program
