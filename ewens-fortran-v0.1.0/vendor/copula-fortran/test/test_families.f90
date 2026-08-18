! SPDX-License-Identifier: GPL-3.0-or-later
program test_families
  use copula
  implicit none
  type(copula_model) :: model, rotated
  real(dp) :: u(2), expected, density, tails(2)
  u = [0.3_dp,0.7_dp]

  model = independence_copula(2)
  call assert_close(pCopula(u,model),0.21_dp,1.0e-14_dp)
  call assert_close(dCopula(u,model),1.0_dp,1.0e-14_dp)

  model = clayton_copula(2.0_dp)
  expected = (u(1)**(-2.0_dp)+u(2)**(-2.0_dp)-1.0_dp)**(-0.5_dp)
  call assert_close(pCopula(u,model),expected,1.0e-14_dp)
  density = 3.0_dp*(u(1)*u(2))**(-3.0_dp) * &
    (u(1)**(-2.0_dp)+u(2)**(-2.0_dp)-1.0_dp)**(-2.5_dp)
  call assert_close(dCopula(u,model),density,2.0e-5_dp)
  call assert_close(tau(model),0.5_dp,1.0e-14_dp)
  tails = lambda(model)
  call assert_close(tails(1),sqrt(0.5_dp),1.0e-13_dp)

  model = gumbel_copula(2.0_dp)
  expected = exp(-sqrt((-log(u(1)))**2+(-log(u(2)))**2))
  call assert_close(pCopula(u,model),expected,1.0e-14_dp)
  call assert_close(tau(model),0.5_dp,1.0e-14_dp)

  model = frank_copula(5.0_dp)
  call assert_close(pCopula(u,model), &
    -log(1.0_dp+(exp(-5.0_dp*u(1))-1.0_dp)*(exp(-5.0_dp*u(2))-1.0_dp)/ &
    (exp(-5.0_dp)-1.0_dp))/5.0_dp,1.0e-13_dp)
  call assert_close(iTau(family_frank,tau(model)),5.0_dp,2.0e-5_dp)

  model = fgm_copula(0.6_dp)
  expected = u(1)*u(2)*(1.0_dp+0.6_dp*(1.0_dp-u(1))*(1.0_dp-u(2)))
  call assert_close(pCopula(u,model),expected,1.0e-14_dp)
  density = 1.0_dp+0.6_dp*(1.0_dp-2.0_dp*u(1))*(1.0_dp-2.0_dp*u(2))
  call assert_close(dCopula(u,model),density,2.0e-6_dp)
  call assert_close(tau(model),2.0_dp*0.6_dp/9.0_dp,1.0e-14_dp)
  call assert_close(rho(model),0.2_dp,1.0e-14_dp)

  model = plackett_copula(1.0_dp)
  call assert_close(pCopula(u,model),product(u),1.0e-14_dp)

  model = clayton_copula(2.0_dp)
  rotated = rotated_copula(model,rotation_180)
  expected = sum(u)-1.0_dp+pCopula(1.0_dp-u,model)
  call assert_close(pCopula(u,rotated),expected,1.0e-14_dp)
  call assert_close(tau(rotated),tau(model),1.0e-14_dp)

  model = galambos_copula(2.0_dp)
  call assert_true(pickands_function(0.5_dp,model) >= 0.5_dp)
  call assert_true(pickands_function(0.5_dp,model) <= 1.0_dp)
  call assert_true(pCopula(u,model) >= product(u))

  model = marshall_olkin_copula(0.4_dp,0.7_dp)
  call assert_close(pCopula(u,model),min(u(1)**0.6_dp*u(2),u(1)*u(2)**0.3_dp),1.0e-14_dp)

  print '(a)', 'test_families: PASS'
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
end program test_families
