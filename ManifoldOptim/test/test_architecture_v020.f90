! SPDX-License-Identifier: GPL-2.0-or-later
program test_architecture_v020
  use manifoldoptim
  implicit none
  type(manifold_domain) :: st, sph, lr, eu
  type(solver_options) :: opt
  type(solver_result) :: res
  real(dp), allocatable :: x(:), eta(:), p(:), y(:), v(:), tv(:)
  real(dp) :: beta, mx, my
  integer :: n, m, r, l1, l2, mode
  logical :: ok
  real(dp) :: target(3), x0(3)
  real(dp), allocatable :: um(:,:), wm(:,:), pum(:,:), pvm(:,:)

  ! Stiefel ParamSet=2: constructed retraction.
  n = 5
  r = 2
  allocate(st%component(1))
  st%component(1) = make_component(MANI_STIEFEL,n,p=r,param_set=2)
  allocate(x(st%length()),eta(st%length()),p(st%length()),y(st%length()))
  call random_manifold_point(st,x)
  call random_number(eta)
  eta = 0.15_dp*(2.0_dp*eta-1.0_dp)
  call project_tangent(st,x,eta,p)
  call retract_point(st,x,p,y,ok)
  if (.not. ok .or. .not. point_is_valid(st,y,2.0e-8_dp)) &
    error stop 'constructed Stiefel retraction'
  deallocate(x,eta,p,y)

  ! Sphere ParamSet=4 QF differentiated retraction / locking beta.
  allocate(sph%component(1))
  sph%component(1) = make_component(MANI_SPHERE,3,param_set=4)
  allocate(x(3),eta(3))
  x = [1.0_dp,0.0_dp,0.0_dp]
  eta = [0.0_dp,0.5_dp,0.0_dp]
  beta = manifold_beta(sph,x,eta)
  if (abs(beta-1.25_dp) > 2.0e-10_dp) then
    write(*,*) 'sphere beta=',beta
    error stop 'sphere locking beta'
  end if
  deallocate(x,eta)

  ! LowRank quotient projection and intrinsic-coordinate transport isometry.
  n = 5
  m = 4
  r = 2
  allocate(lr%component(1))
  lr%component(1) = make_component(MANI_LOWRANK,n,m=m,p=r)
  allocate(x(lr%length()),eta(lr%length()),p(lr%length()),y(lr%length()))
  allocate(v(lr%length()),tv(lr%length()))
  call random_manifold_point(lr,x)
  call random_number(v)
  v = 2.0_dp*v-1.0_dp
  call project_tangent(lr,x,v,p)
  l1 = n*r
  l2 = l1+r*r
  allocate(um(n,r),wm(m,r),pum(n,r),pvm(m,r))
  um = reshape(x(1:l1),[n,r])
  wm = reshape(x(l2+1:),[m,r])
  pum = reshape(p(1:l1),[n,r])
  pvm = reshape(p(l2+1:),[m,r])
  if (maxval(abs(matmul(transpose(um),pum))) > 2.0e-10_dp) &
    error stop 'LowRank U quotient projection'
  if (maxval(abs(matmul(transpose(wm),pvm))) > 2.0e-10_dp) &
    error stop 'LowRank V quotient projection'
  call random_number(eta)
  eta = 0.02_dp*(2.0_dp*eta-1.0_dp)
  call project_tangent(lr,x,eta,p)
  eta = p
  call retract_point(lr,x,eta,y,ok)
  if (.not. ok) error stop 'LowRank retraction in transport test'
  call transport_vector(lr,x,y,p,tv)
  mx = manifold_metric(lr,x,p,p)
  my = manifold_metric(lr,y,tv,tv)
  if (abs(mx-my) > 2.0e-7_dp*max(1.0_dp,mx)) then
    write(*,*) 'lowrank metric before/after=',mx,my
    error stop 'LowRank transport isometry'
  end if
  deallocate(x,eta,p,y,v,tv,um,wm,pum,pvm)

  ! Every built-in line-search mode and the custom callback path.
  allocate(eu%component(1))
  eu%component(1) = make_component(MANI_EUCLIDEAN,3,m=1)
  target = [1.0_dp,-2.0_dp,0.5_dp]
  x0 = [4.0_dp,3.0_dp,-1.0_dp]
  opt = solver_options()
  opt%tolerance = 1.0e-8_dp
  opt%max_iteration = 80
  do mode = LINESEARCH_ARMIJO, LINESEARCH_EXACT
    opt%line_search = mode
    call manifold_optimize(eu,x0,obj,grad,'RSD',res,opt)
    if (maxval(abs(res%xopt-target)) > 2.0e-5_dp) then
      write(*,*) 'line search mode failed',mode,res%xopt
      error stop 'built-in line search'
    end if
  end do
  opt%line_search = LINESEARCH_INPUTFUN
  opt%line_search_proc => custom_line_search
  call manifold_optimize(eu,x0,obj,grad,'RSD',res,opt)
  if (maxval(abs(res%xopt-target)) > 2.0e-10_dp) error stop 'custom line search'

  ! Non-default Broyden member is implemented, not dispatched as BFGS.
  opt = solver_options()
  opt%broyden_phi = 0.35_dp
  opt%tolerance = 1.0e-8_dp
  opt%max_iteration = 100
  call manifold_optimize(eu,x0,obj,grad,'RBroydenFamily',res,opt)
  if (maxval(abs(res%xopt-target)) > 2.0e-5_dp) error stop 'Broyden family phi'

  write(*,*) 'PASS test_architecture_v020'

contains

  subroutine obj(xx,f)
    real(dp), intent(in) :: xx(:)
    real(dp), intent(out) :: f
    f = 0.5_dp*sum((xx-target)**2)
  end subroutine obj

  subroutine grad(xx,g)
    real(dp), intent(in) :: xx(:)
    real(dp), intent(out) :: g(:)
    g = xx-target
  end subroutine grad

  function custom_line_search(xx,dir,initial_step,initial_slope) result(step)
    real(dp), intent(in) :: xx(:), dir(:), initial_step, initial_slope
    real(dp) :: step
    if (size(xx) /= size(dir)) error stop 'custom line search size'
    if (initial_slope >= 0.0_dp) error stop 'custom line search direction'
    step = min(1.0_dp,initial_step)
  end function custom_line_search

end program test_architecture_v020
