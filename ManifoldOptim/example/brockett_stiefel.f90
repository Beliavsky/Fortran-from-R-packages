! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
program brockett_stiefel
  use manifoldoptim
  implicit none

  integer, parameter :: n = 5, p = 2
  type(manifold_domain) :: domain
  type(solver_options) :: opt
  type(solver_result) :: res
  real(dp) :: x0(n*p), xmat(n,p)
  real(dp) :: amat(n,n), dmat(p,p)
  integer :: i
  logical :: ok

  allocate(domain%component(1))
  domain%component(1) = make_component(MANI_STIEFEL, n, p=p)

  amat = 0.0_dp
  do i = 1, n
    amat(i,i) = real(i, dp)
  end do
  dmat = 0.0_dp
  dmat(1,1) = 1.0_dp
  dmat(2,2) = 2.0_dp

  xmat(:,1) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
  xmat(:,2) = [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp, 0.5_dp]
  call orthonorm(reshape(xmat, [n*p]), n, p, x0, ok)
  if (.not. ok) error stop 'failed to construct Stiefel starting point'

  opt%tolerance = 1.0e-8_dp
  opt%max_iteration = 500
  call manifold_optimize(domain, x0, objective, gradient, 'LRBFGS', res, opt)

  xmat = reshape(res%xopt, [n,p])
  write(*,'(a,1x,es14.6)') 'Brockett objective =', res%fval
  write(*,'(a)') 'X*:'
  do i = 1, n
    write(*,'(*(1x,f10.6))') xmat(i,:)
  end do

contains

  subroutine objective(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    real(dp) :: xm(n,p)
    xm = reshape(x, [n,p])
    f = trace_matrix(matmul(transpose(xm), matmul(amat, matmul(xm, dmat))))
  end subroutine objective

  subroutine gradient(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: xm(n,p), gm(n,p)
    xm = reshape(x, [n,p])
    gm = 2.0_dp * matmul(amat, matmul(xm, dmat))
    g = reshape(gm, [n*p])
  end subroutine gradient

end program brockett_stiefel
