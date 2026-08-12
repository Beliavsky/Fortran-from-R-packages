program test_trees
    use clue
    implicit none
    real(dp) :: d(4,4), u(4,4), a(4,4)
    integer :: it
    d=0.0_dp
    d(1,2)=1
    d(1,3)=4
    d(1,4)=5
    d(2,3)=2
    d(2,4)=6
    d(3,4)=3
    d=d+transpose(d)
    call check(non_ultrametricity(d,.true.)>0,'non-ultrametric detected')
    u=d
    call fit_ultrametric_ip(u,maxiter=10000,tol=1e-12_dp,iterations=it)
    call check(non_ultrametricity(u,.true.)<1e-10_dp,'ultrametric IP fit')
    u=d
    call fit_ultrametric_ir(u,maxiter=10000,tol=1e-12_dp,iterations=it)
    call check(non_ultrametricity(u,.true.)<1e-8_dp,'ultrametric IR fit')
    a=d
    call fit_addtree_ip(a,maxiter=10000,tol=1e-12_dp,iterations=it)
    call check(non_additivity(a,.true.)<1e-9_dp,'addtree IP fit')
    a=d
    call fit_addtree_ir(a,maxiter=10000,tol=1e-12_dp,iterations=it)
    call check(non_additivity(a,.true.)<1e-8_dp,'addtree IR fit')
    u=ultrametrify(d)
    call check(is_ultrametric(u,1e-12_dp),'single-link ultrametrify')
    print *, 'test_trees: PASS'
contains
    subroutine check(ok,msg)
    logical,intent(in)::ok
    character(*),intent(in)::msg
    if(.not.ok)then
    write(*,*)'FAIL: ',trim(msg)
    error stop 1
    end if
    end subroutine
end program
