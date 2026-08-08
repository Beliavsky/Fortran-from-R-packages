program demo
   use neighbours, only: numeric_neighbour_config, init_numeric_neighbour, neighbours_ok
   use neighbours_nmof_adapter, only: nmof_numeric_neighbour
   use nmof_kinds, only: dp, i8
   use nmof_optimization, only: local_search
   use nmof_types, only: optimization_result
   implicit none
   type(numeric_neighbour_config) :: config
   type(optimization_result) :: result
   real(dp) :: x0(5)
   integer :: status

   x0 = [0.50_dp, 0.20_dp, 0.10_dp, 0.10_dp, 0.10_dp]
   call init_numeric_neighbour(config, 5, 0.05_dp, lower=[0.0_dp], upper=[1.0_dp], &
      status=status)
   if (status /= neighbours_ok) error stop 'invalid neighbourhood configuration'

   call local_search(objective, nmof_numeric_neighbour, x0, result, n_steps=20000, &
      seed=12345_i8, context=config)

   print '(a,f12.8)', 'objective = ', result%ofvalue
   print '(a,5f10.6)', 'xbest     = ', result%xbest
contains
   function objective(x, context) result(f)
      real(dp), intent(in) :: x(:)
      class(*), intent(in), optional :: context
      real(dp) :: f

      if (present(context)) continue
      f = sum((x - 0.20_dp)**2)
   end function objective
end program demo
