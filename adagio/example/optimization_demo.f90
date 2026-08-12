program optimization_demo
  use iso_fortran_env, only : int64
  use adagio
  implicit none
  type(opt_result) :: nm
  type(de_result) :: de
  nm = neldermead(fn_rosenbrock, [-1.2_dp,1._dp])
  de = simple_de(fn_rastrigin, [-5.12_dp,-5.12_dp], [5.12_dp,5.12_dp], &
                  n_pop=48, nmax=200, seed=12345_int64)
  print '(a,2f14.7,a,es14.5)', 'Nelder-Mead Rosenbrock x = ',nm%x,'  f = ',nm%f
  print '(a,2f14.7,a,es14.5)', 'simpleDE Rastrigin x     = ',de%xmin,'  f = ',de%fmin
end program optimization_demo
