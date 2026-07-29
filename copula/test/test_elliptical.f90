! SPDX-License-Identifier: GPL-3.0-or-later
program test_elliptical
  use copula
  implicit none
  type(copula_model) :: normal_model, student_model, normal3
  real(dp) :: correlation(2,2), correlation3(3,3), u(2), u3(3), tails(2)
  real(dp), allocatable :: sample(:,:)
  logical :: ok

  correlation = reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
  normal_model = normal_copula(correlation)
  u = [0.3_dp,0.7_dp]
  call assert_close(pCopula(u,normal_model),0.2669038488673631_dp,2.0e-9_dp)
  call assert_close(dCopula(u,normal_model),0.8770819376466364_dp,2.0e-12_dp)
  call assert_close(tau(normal_model),1.0_dp/3.0_dp,2.0e-14_dp)
  call assert_close(rho(normal_model),6.0_dp*asin(0.25_dp)/pi,2.0e-14_dp)

  correlation = reshape([1.0_dp,0.4_dp,0.4_dp,1.0_dp],[2,2])
  student_model = t_copula(correlation,5.0_dp)
  call assert_close(pCopula(u,student_model),0.2518449_dp,2.0e-5_dp)
  call assert_close(dCopula(u,student_model),0.893659595147815_dp,2.0e-11_dp)
  tails = lambda(student_model)
  call assert_close(tails(1),tails(2),1.0e-14_dp)
  call assert_true(tails(1) > 0.0_dp .and. tails(1) < 1.0_dp)

  correlation3 = reshape([1.0_dp,0.2_dp,0.1_dp,0.2_dp,1.0_dp,0.3_dp,0.1_dp,0.3_dp,1.0_dp],[3,3])
  normal3 = normal_copula(correlation3)
  u3 = [0.25_dp,0.5_dp,0.8_dp]
  call assert_true(pCopula(u3,normal3) > 0.0_dp)
  call assert_true(pCopula(u3,normal3) < minval(u3))
  call rCopula(2000,normal3,sample,ok,12345_i8)
  call assert_true(ok)
  call assert_true(all(sample > 0.0_dp) .and. all(sample < 1.0_dp))
  call assert_close(sum(sample(:,1))/real(size(sample,1),dp),0.5_dp,0.03_dp)

  print '(a)', 'test_elliptical: PASS'
contains
  subroutine assert_close(actual, reference, tolerance)
    real(dp), intent(in) :: actual, reference, tolerance
    if (abs(actual-reference) > tolerance*(1.0_dp+abs(reference))) then
      print '(a,3es25.16)', 'mismatch: ', actual, reference, abs(actual-reference)
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true
end program test_elliptical
