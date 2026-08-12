program test_global_optimizers
  use iso_fortran_env, only : int64
  use adagio
  implicit none
  type(de_result) :: de
  type(ea_result) :: ea
  type(cmaes_result) :: cm

  de = simple_de(sphere, [-5._dp,-5._dp], [5._dp,5._dp], n_pop=24, nmax=120, seed=1234_int64)
  call check(de%fmin < 1e-6_dp, 'simpleDE convergence')
  call check(de%actual_nfeval == 24 + 24*120, 'simpleDE actual count')
  call check(de%nfeval == 24 + 120*24*24, 'simpleDE source count')

  ea = simple_ea(sphere, [0._dp,0._dp], [4._dp,4._dp], n_pop=40, &
                 eps=1e-4_dp, tol=1e-8_dp, seed=4321_int64)
  call check(ea%val < 1e-2_dp, 'simpleEA convergence')

  cm = pure_cmaes(sphere, [2._dp,-2._dp], [-5._dp,-5._dp], [5._dp,5._dp], &
                  stopfitness=1e-10_dp, stopeval=5000, seed=2024_int64)
  call check(cm%fmin < 1e-6_dp, 'CMA-ES convergence')

  print *, 'test_global_optimizers: PASS'
contains
  function sphere(x) result(f)
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=sum(x*x)
  end function
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine
end program test_global_optimizers
