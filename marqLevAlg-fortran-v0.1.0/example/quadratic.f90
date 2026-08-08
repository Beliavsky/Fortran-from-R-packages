program quadratic
    use marqlevalg, only : dp, mla_result, marqlev_optimize
    implicit none
    type(mla_result) :: r
    real(dp) :: x0(2)

    x0 = [8.0_dp, 9.0_dp]
    call marqlev_optimize(x0, objective, r)
    print '(a,2f14.8)', 'par = ', r%par
    print '(a,es14.6)', 'f   = ', r%fn_value
contains
    function objective(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = 4.0_dp * (x(1) - 5.0_dp)**2 + (x(2) - 6.0_dp)**2
    end function objective
end program quadratic
