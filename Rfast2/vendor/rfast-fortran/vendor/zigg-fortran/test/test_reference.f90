program test_reference
    use iso_fortran_env, only : real64
    use zigg, only : dp, zsetseed, zrnorm, zrexp, zrunif
    implicit none

    real(dp), parameter :: tol = 5.0e-15_dp
    real(dp), parameter :: normal_ref(20) = [ &
        1.2977704266575947_dp, -0.69905311826578109_dp, 0.44307945774419810_dp, &
       -0.40529617136667800_dp, -1.7715425529828976_dp, 0.31989625797866905_dp, &
       -1.0571031728518958_dp, 0.22800338862237893_dp, -0.92025508665866851_dp, &
       -0.97783728935397929_dp, 1.0951163803069051_dp, 0.87029331684500943_dp, &
       -0.21738758774612582_dp, -0.20750401602493954_dp, -1.0990403426196333_dp, &
       -1.8673959561522160_dp, -0.58389295227628923_dp, 1.3083155688333912_dp, &
       -0.34079091525222643_dp, 0.93308654446068506_dp ]
    real(dp), parameter :: exp_ref(20) = [ &
        0.40564453632368941_dp, 0.24633474127749741_dp, 0.11836257657620800_dp, &
        1.3193827244380658_dp, 0.68371244257191288_dp, 0.079957301751930582_dp, &
        1.8069859843533540_dp, 0.065858166901015383_dp, 1.2431735506062502_dp, &
        3.1905170608528852_dp, 0.054560495134270963_dp, 0.56004067750450104_dp, &
        0.84087447558381978_dp, 0.18580722702012945_dp, 1.9454944550346160_dp, &
        0.13529754453930168_dp, 1.1048259339220690_dp, 7.8747457100472564_dp, &
        0.29681112238023483_dp, 1.4877357744032622_dp ]
    real(dp), parameter :: unif_ref(20) = [ &
        0.90864950300334801_dp, 0.092843255463961993_dp, 0.71170171781502500_dp, &
        0.22274627798810481_dp, 0.046350645206655217_dp, 0.67500457355606658_dp, &
        0.16126685281570280_dp, 0.59124027973198456_dp, 0.084393795504235425_dp, &
        0.28121379627436038_dp, 0.020424277565642024_dp, 0.89137859289390320_dp, &
        0.97088390910309141_dp, 0.61464155457969194_dp, 0.36242813425668496_dp, &
        0.52611171500127984_dp, 0.49501766707099337_dp, 0.20667617726266280_dp, &
        0.58885906460794879_dp, 0.83724736874554462_dp ]
    real(dp), allocatable :: x(:)

    call zsetseed(12345)
    x = zrnorm(20)
    call check_close(x, normal_ref, 'normal')

    call zsetseed(12345)
    x = zrexp(20)
    call check_close(x, exp_ref, 'exponential')

    call zsetseed(12345)
    x = zrunif(20)
    call check_close(x, unif_ref, 'uniform')

    print '(a)', 'test_reference: PASS'

contains

    subroutine check_close(x, ref, label)
        real(real64), intent(in) :: x(:), ref(:)
        character(len=*), intent(in) :: label
        if (size(x) /= size(ref)) error stop 'size mismatch'
        if (maxval(abs(x - ref)) > tol) then
            print '(a,1x,es24.16)', trim(label)//' max error:', maxval(abs(x - ref))
            error stop 'reference mismatch'
        end if
    end subroutine check_close

end program test_reference
