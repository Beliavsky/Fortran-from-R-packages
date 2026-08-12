program test_velocity_compatibility
   use psoptim, only : dp, ps_control, ps_result, ps_optimize
   implicit none
   type(ps_control) :: legacy, fixed
   type(ps_result) :: r1, r2
   real(dp) :: xmin1(1), xmax1(1), vmax1(1)
   real(dp) :: xmin2(2), xmax2(2), vmax2(2)

   ! Seed 1 makes the second RNG draw (the initial velocity for this 1-D,
   ! one-particle case) strongly negative. With w=3 it crosses -vmax.
   xmin1 = [0.0_dp]
   xmax1 = [10.0_dp]
   vmax1 = [1.0_dp]
   legacy%n = 1
   legacy%max_loop = 1
   legacy%w = 3.0_dp
   legacy%c1 = 0.0_dp
   legacy%c2 = 0.0_dp
   legacy%seed = 1
   call ps_optimize(flat, xmin1, xmax1, vmax1, r1, legacy)
   if (abs(r1%final_velocity(1,1) - 1.0_dp) > 1.0e-14_dp) &
      error stop "legacy negative velocity was not reset to +vmax"

   fixed = legacy
   fixed%legacy_velocity_clip = .false.
   call ps_optimize(flat, xmin1, xmax1, vmax1, r2, fixed)
   if (abs(r2%final_velocity(1,1) + 1.0_dp) > 1.0e-14_dp) &
      error stop "corrected negative velocity was not clipped to -vmax"

   ! R's runif(n*d,min=-vmax,max=vmax) recycles vmax over the linear draw
   ! stream, so with n=3,d=2 a dimension-1 particle can get vmax(2).
   xmin2 = [0.0_dp, 0.0_dp]
   xmax2 = [1.0_dp, 1.0_dp]
   vmax2 = [0.01_dp, 100.0_dp]
   legacy = ps_control()
   legacy%n = 3
   legacy%max_loop = 0
   legacy%seed = 7
   call ps_optimize(flat, xmin2, xmax2, vmax2, r1, legacy)
   if (abs(r1%final_velocity(1,2)) <= 0.01_dp) &
      error stop "legacy velocity initialization did not recycle vmax"

   fixed = legacy
   fixed%legacy_velocity_initialization = .false.
   call ps_optimize(flat, xmin2, xmax2, vmax2, r2, fixed)
   if (maxval(abs(r2%final_velocity(1,:))) > 0.01_dp + 1.0e-14_dp) &
      error stop "corrected velocity initialization ignored dimension limit"

   print *, "test_velocity_compatibility: PASS"
contains
   function flat(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = 0.0_dp*sum(x)
   end function flat
end program test_velocity_compatibility
