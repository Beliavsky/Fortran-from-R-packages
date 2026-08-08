! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
program test_remaining_manifolds
  use manifoldoptim
  implicit none

  call test_grassmann()
  call test_orthgroup()
  call test_spd()
  write(*,*) 'PASS test_remaining_manifolds'

contains

  subroutine test_grassmann()
    type(manifold_domain) :: domain
    type(solver_options) :: opt
    type(solver_result) :: res
    real(dp) :: x0(5), a(5), truth(5)

    allocate(domain%component(1))
    domain%component(1) = make_component(MANI_GRASSMANN, 5, p=1)
    a = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, 1.5_dp]
    truth = a / vecnorm(a)
    x0 = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    opt%tolerance = 1.0e-7_dp
    opt%max_iteration = 400
    call manifold_optimize(domain, x0, grass_obj, grass_grad, 'RCG', res, opt)
    if (.not. point_is_valid(domain, res%xopt, 2.0e-7_dp)) error stop 'Grassmann validity'
    if (abs(dot_product(res%xopt, truth)) < 0.999_dp) error stop 'Grassmann optimum'
  end subroutine test_grassmann

  subroutine grass_obj(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    real(dp) :: a(5)
    a = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, 1.5_dp]
    f = -dot_product(a, x)**2
  end subroutine grass_obj

  subroutine grass_grad(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: a(5)
    a = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, 1.5_dp]
    g = -2.0_dp * dot_product(a, x) * a
  end subroutine grass_grad

  subroutine test_orthgroup()
    integer, parameter :: n = 3
    type(manifold_domain) :: domain
    type(solver_options) :: opt
    type(solver_result) :: res
    real(dp) :: x0(n*n), q(n,n)

    allocate(domain%component(1))
    domain%component(1) = make_component(MANI_ORTHGROUP, n)
    q = 0.0_dp
    q(1,2) = -1.0_dp
    q(2,1) = 1.0_dp
    q(3,3) = 1.0_dp
    x0 = reshape(q, [n*n])
    opt%tolerance = 1.0e-7_dp
    opt%max_iteration = 400
    call manifold_optimize(domain, x0, orth_obj, orth_grad, 'LRBFGS', res, opt)
    if (.not. point_is_valid(domain, res%xopt, 2.0e-7_dp)) error stop 'OrthGroup validity'
    q = reshape(res%xopt, [n,n])
    if (maxval(abs(q - eye_matrix(n))) > 2.0e-3_dp) error stop 'OrthGroup optimum'
  end subroutine test_orthgroup

  subroutine orth_obj(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    real(dp) :: q(3,3)
    q = reshape(x, [3,3])
    f = 0.5_dp * sum((q - eye_matrix(3))**2)
  end subroutine orth_obj

  subroutine orth_grad(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: q(3,3)
    q = reshape(x, [3,3])
    g = reshape(q - eye_matrix(3), [9])
  end subroutine orth_grad

  subroutine test_spd()
    integer, parameter :: n = 3
    type(manifold_domain) :: domain
    type(solver_options) :: opt
    type(solver_result) :: res
    real(dp) :: x0(n*n), q(n,n)

    allocate(domain%component(1))
    domain%component(1) = make_component(MANI_SPD, n)
    q = 0.0_dp
    q(1,1) = 2.0_dp
    q(2,2) = 1.5_dp
    q(3,3) = 0.7_dp
    x0 = reshape(q, [n*n])
    opt%tolerance = 1.0e-7_dp
    opt%max_iteration = 500
    opt%initial_step = 0.5_dp
    call manifold_optimize(domain, x0, spd_obj, spd_grad, 'RSD', res, opt)
    if (.not. point_is_valid(domain, res%xopt, 1.0e-8_dp)) error stop 'SPD validity'
    q = reshape(res%xopt, [n,n])
    if (maxval(abs(q - eye_matrix(n))) > 2.0e-3_dp) then
      write(*,*) res%fval, res%normgf, q
      error stop 'SPD optimum'
    end if
  end subroutine test_spd

  subroutine spd_obj(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    real(dp) :: q(3,3)
    q = reshape(x, [3,3])
    f = 0.5_dp * sum((q - eye_matrix(3))**2)
  end subroutine spd_obj

  subroutine spd_grad(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: q(3,3)
    q = reshape(x, [3,3])
    g = reshape(q - eye_matrix(3), [9])
  end subroutine spd_grad

end program test_remaining_manifolds
