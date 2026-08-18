! SPDX-License-Identifier: GPL-3.0-or-later
program integrate_example
    use cubature, only : dp, cubature_result, hcubature, vegas, vegas_options
    implicit none
    type(cubature_result) :: r
    type(vegas_options) :: opt
    call hcubature(f, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], 1, r, 1.0e-8_dp, 1.0e-12_dp)
    print '(a,es20.12,a,es12.4,a,i0)', 'hcubature = ', r%integral(1), '  error = ', r%error(1), &
        '  neval = ', r%evaluations
    opt%nstart = 4000
    opt%max_eval = 100000
    call vegas(f, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], 1, r, 1.0e-3_dp, 1.0e-10_dp, opt)
    print '(a,es20.12,a,es12.4,a,i0)', 'vegas      = ', r%integral(1), '  error = ', r%error(1), &
        '  neval = ', r%evaluations
contains
    subroutine f(x, value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value(:)
        value(1) = product(cos(x))
    end subroutine f
end program integrate_example
